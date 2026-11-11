import sys, time
sys.path.insert(0, "/home/dominik/Documents/other/mri")
import numpy as np, mri_tools as mt
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
t0=time.time()
a, Lc = 0.35, 1.2
Nphi, Nz = 96, 160
phi = (np.arange(Nphi)+0.5)*2*np.pi/Nphi
zc = -Lc/2 + (np.arange(Nz)+0.5)*Lc/Nz
PH, ZZ = np.meshgrid(phi, zc, indexing="ij")
dA = a*(2*np.pi/Nphi)*(Lc/Nz)
elem = np.stack([a*np.cos(PH), a*np.sin(PH), ZZ],-1).reshape(-1,3)
phihat = np.stack([-np.sin(PH), np.cos(PH), 0*PH],-1).reshape(-1,3)
zhat = np.array([0,0,1.0])
# target points
rng = np.random.default_rng(1)
pts = rng.normal(size=(400,3)); pts = pts/np.linalg.norm(pts,axis=1,keepdims=True)*0.2*rng.uniform(0.3,1,(400,1))**(1/3)
pts = np.vstack([pts, 0.2*pts[:200]/np.linalg.norm(pts[:200],axis=1,keepdims=True)])
R = pts[:,None,:]-elem[None]
Rn3 = np.linalg.norm(R,axis=-1)**3
def bz_kernel(direction):
    c = np.cross(np.broadcast_to(direction, R.shape), R)
    return mt.MU0/(4*np.pi)*c[...,2]/Rn3*dA
Kphi_ker = bz_kernel(phihat.reshape(1,-1,3))
Kz_ker = bz_kernel(zhat)
M = 12
A = np.zeros((len(pts), M)); Pm = np.zeros((M,M)); Kall=[]
for m in range(1,M+1):
    k = 2*np.pi*m/Lc
    f = np.sin(k*ZZ); fp = k*np.cos(k*ZZ)
    Kp = (-fp*np.cos(PH)).ravel(); Kzv = (-(1/a)*f*np.sin(PH)).ravel()
    A[:,m-1] = Kphi_ker@Kp + Kz_ker@Kzv
    Kall.append(np.concatenate([Kp,Kzv]))
Kall=np.array(Kall); Pm = Kall@Kall.T*dA
target = pts[:,0]*1.0  # G = 1 T/m
for lam in [1e-14, 1e-13,1e-12]:
    c = np.linalg.solve(A.T@A + lam*Pm, A.T@target)
    err = A@c-target
    print(lam, np.abs(err).max()/0.2, "psi max", np.abs(np.sin(2*np.pi*np.arange(1,M+1)[None]*zc[:,None]/Lc)@c).max())
print(time.time()-t0)
lam=1e-13
c = np.linalg.solve(A.T@A + lam*Pm, A.T@target)
fz = np.sin(2*np.pi*np.arange(1,M+1)[None]*zc[:,None]/Lc)@c
psi = fz[None,:]*np.cos(phi)[:,None]
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].plot(zc,fz)
cs=ax[1].contour(np.degrees(phi), zc, psi.T, levels=np.linspace(-abs(psi).max(),abs(psi).max(),22)[1:-1], colors="k")
plt.savefig("fp.png")
print(len(cs.allsegs), [len(s) for s in cs.allsegs])
