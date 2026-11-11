import sys,time; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt, scipy.io as sio
exec(open('h3a.py').read().split("C=coil_sensitivities(N,FOV)")[0])
t0=time.time()
C=coil_sensitivities(N,FOV)
m=mt.shepp_logan(N)
K=mt.fft2c(C*m)
K=mt.add_complex_noise(K,0.002)
# GRAPPA R=2, 24 ACS lines, kernel 2(ky)x3(kx)
R=2; nacs=24
def grappa_r2(Kfull, R=2, nacs=24, lam=1e-3):
    Nc,Ny,Nx=Kfull.shape
    samp=np.zeros(Ny,bool); samp[::R]=True
    acs=np.zeros(Ny,bool); acs[Ny//2-nacs//2:Ny//2+nacs//2]=True
    Ku=Kfull*(samp|acs)[None,:,None]
    # calibration
    A=[];B=[]
    for y in range(Ny//2-nacs//2+1, Ny//2+nacs//2-1):
        for x in range(1,Nx-1):
            A.append(Ku[:,[y-1,y+1]][:,:,x-1:x+2].ravel()); B.append(Ku[:,y,x])
    A=np.array(A);B=np.array(B)
    AH=A.conj().T
    W=np.linalg.solve(AH@A+lam*np.trace(AH@A).real/AH.shape[0]*np.eye(A.shape[1]),AH@B)
    Kr=Ku.copy()
    for y in range(1,Ny-1):
        if samp[y] or acs[y]: continue
        src=np.stack([Kr[:,y-1],Kr[:,y+1]],1)  # Nc,2,Nx
        # build patches for all x
        pad=np.pad(src,((0,0),(0,0),(1,1)))
        P=np.stack([pad[:,:,i:i+Nx] for i in range(3)],-1) # Nc,2,Nx,3
        P=P.transpose(2,0,1,3).reshape(Nx,-1)
        Kr[:,y]=(P@W).T
    return Ku,Kr,samp|acs
Ku,Kr,smask=grappa_r2(K)
rss=lambda X: np.sqrt(np.sum(np.abs(mt.ifft2c(X))**2,0))
ref=rss(K); zf=rss(Ku); gr=rss(Kr)
mask=m>0
print('nrmse zf',np.linalg.norm((zf-ref)[mask])/np.linalg.norm(ref[mask]),'grappa',np.linalg.norm((gr-ref)[mask])/np.linalg.norm(ref[mask]), time.time()-t0)
# Haar
def haar2(x, levels=3):
    c=x.astype(complex).copy(); n=x.shape[0]
    for l in range(levels):
        s=c[:n,:n]
        a=(s[0::2]+s[1::2])/np.sqrt(2); d=(s[0::2]-s[1::2])/np.sqrt(2); s=np.concatenate([a,d],0)
        a=(s[:,0::2]+s[:,1::2])/np.sqrt(2); d=(s[:,0::2]-s[:,1::2])/np.sqrt(2); s=np.concatenate([a,d],1)
        c[:n,:n]=s; n//=2
    return c
def ihaar2(c, levels=3):
    x=c.copy(); n=c.shape[0]>>(levels-1)
    for l in range(levels):
        s=x[:n,:n]; h=n//2
        a,d=s[:,:h],s[:,h:]; t=np.empty_like(s); t[:,0::2]=(a+d)/np.sqrt(2); t[:,1::2]=(a-d)/np.sqrt(2); s=t
        a,d=s[:h],s[h:]; t=np.empty_like(s); t[0::2]=(a+d)/np.sqrt(2); t[1::2]=(a-d)/np.sqrt(2)
        x[:n,:n]=t; n*=2
    return x
x=np.random.randn(64,64); print('haar perfect', np.abs(ihaar2(haar2(x,4),4)-x).max(), np.linalg.norm(haar2(x,4))/np.linalg.norm(x))
# CS ISTA on shepp logan 128 R=4 variable density lines
rng=np.random.default_rng(3)
img=mt.shepp_logan(128)
ky=np.arange(128)-64
pdf=np.minimum(1, 0.08+ 1.0*np.exp(-(ky/18.)**2))
pdf*= (128/4)/pdf.sum()
pdf=np.minimum(pdf,1)
mk=rng.random(128)<pdf; mk[np.abs(ky)<6]=True
print('R=',128/mk.sum())
Mk=np.repeat(mk[:,None],128,1)
y=mt.fft2c(img)*Mk
y=mt.add_complex_noise(y,0.002)*Mk
def soft(z,t): return np.where(np.abs(z)>t, (1-t/np.maximum(np.abs(z),1e-12))*z, 0)
def fista(y,Mk,lam,n_it=100,levels=4):
    x=mt.ifft2c(y); z=x.copy(); tk=1
    hist=[]
    for k in range(n_it):
        g=z+mt.ifft2c(Mk*(y-mt.fft2c(z)*Mk))
        xn=ihaar2(soft(haar2(g,levels),lam),levels)
        tn=(1+np.sqrt(1+4*tk*tk))/2
        z=xn+(tk-1)/tn*(xn-x); x=xn; tk=tn
        hist.append(np.linalg.norm(np.abs(x)-img)/np.linalg.norm(img))
    return x,hist
t1=time.time()
zf=mt.ifft2c(y)
for lam in [0.002,0.01,0.03]:
    xr,h=fista(y,Mk,lam)
    print(lam, 'zf err',np.linalg.norm(np.abs(zf)-img)/np.linalg.norm(img),'cs err',h[-1], time.time()-t1)
