import sys,time; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
t0=time.time()
N=128; FOV=0.24
def coil_sensitivities(N, FOV, n_coils=8, R_array=0.14, loop_radius=0.05, z_slice=0.0, z_coils=(0.0,), phase0=0.0):
    x=(np.arange(N)-N//2)*FOV/N
    X,Y=np.meshgrid(x,-x)  # y up
    pts=np.stack([X,Y,np.full_like(X,z_slice)],-1)
    C=[]
    for zc in z_coils:
        for j in range(n_coils):
            a=2*np.pi*j/n_coils+phase0
            Rm=mt.Rz(-a)
            p_rot=pts@Rm.T     # rotate points by -a
            B=mt.loop_coil_field(p_rot, loop_radius, center=(R_array,0,zc), normal='x', n_segments=64)
            B=B@mt.Rz(a).T
            C.append(B[...,0]-1j*B[...,1])
    C=np.array(C)
    return C/np.sqrt(np.sum(np.abs(C)**2,0)).max()
C=coil_sensitivities(N,FOV)
print(C.shape, time.time()-t0)
m=mt.shepp_logan(N)
fig,ax=plt.subplots(2,8,figsize=(16,4))
for j in range(8):
    ax[0,j].imshow(np.abs(C[j]));ax[1,j].imshow(np.angle(C[j]),cmap='hsv')
plt.savefig('coils.png')
# SENSE R=2
def sense(alias, C, R, Psi=None):
    Nc,N,_=C.shape
    Psi=np.eye(Nc) if Psi is None else Psi
    Pi=np.linalg.inv(Psi)
    out=np.zeros((N,N),complex); g=np.zeros((N,N))
    Ny=N//R
    for y in range(Ny):
        ys=[y+k*Ny for k in range(R)]
        Cm=C[:,ys,:].transpose(2,0,1)   # (x, Nc, R)
        a=alias[:,y,:].T                 # (x,Nc)
        CH=np.conj(Cm).transpose(0,2,1)
        A=CH@Pi@Cm                       # (x,R,R)
        rhs=CH@Pi@a[...,None]
        A_reg=A+1e-9*np.eye(R)
        sol=np.linalg.solve(A_reg,rhs)[...,0]
        Ainv=np.linalg.inv(A_reg)
        for k,yy in enumerate(ys):
            out[yy,:]=sol[:,k]
            g[yy,:]=np.sqrt(np.abs(Ainv[:,k,k]*A[:,k,k]))
    return out,g
coil_imgs=C*m
K=mt.fft2c(coil_imgs)
for R in [2,3,4]:
    Ku=np.zeros_like(K); Ku[:,::R,:]=K[:,::R,:]  # wait rows are ky axis -2
    Ku=np.zeros_like(K); Ku[:,::R]=K[:,::R]
    al=mt.ifft2c(Ku)*R   # alias full FOV image; take one FOV/R period
    # aliased images: im_alias(y) = sum_k C m(y + k N/R) / R ... with ortho: ifft2c of zero-filled gives (1/R) sum shifted
    alias_small=al[:, :N//R, :]
    # need periodic version: al rows y and y+N/R are identical
    rec,g=sense(al[:, :, :], C, R)
    mask=m>0
    print(R, np.abs(rec-m)[mask].max(), g[mask].mean(), g[mask].max())
print(time.time()-t0)
