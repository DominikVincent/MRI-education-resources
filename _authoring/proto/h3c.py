exec(open('h3b.py').read().split("# CS ISTA")[0].split("C=coil_sensitivities(N,FOV)\nm=")[0])
import numpy as np, mri_tools as mt, time
exec("def haar2"+open('h3b.py').read().split("def haar2")[1].split("x=np.random")[0])
def soft(z,t): return np.where(np.abs(z)>t, (1-t/np.maximum(np.abs(z),1e-12))*z, 0)
img=mt.shepp_logan(128); rng=np.random.default_rng(3)
k=np.arange(128)-64
def fista(y,Mk,lam,n_it=200,levels=4,keep_coarse=True):
    x=mt.ifft2c(y); z=x.copy(); tk=1; n0=128>>levels
    for it in range(n_it):
        g=z+mt.ifft2c(Mk*(y-mt.fft2c(z)*Mk))
        c=haar2(g,levels); cc=c[:n0,:n0].copy(); c=soft(c,lam)
        if keep_coarse: c[:n0,:n0]=cc
        xn=ihaar2(c,levels)
        tn=(1+np.sqrt(1+4*tk*tk))/2; z=xn+(tk-1)/tn*(xn-x); x=xn; tk=tn
    return x
err=lambda a: np.linalg.norm(np.abs(a)-img)/np.linalg.norm(img)
# 1D lines
pdf=np.exp(-(k/20.)**2)*0.9+0.1; 
for name in ['lines','2d']:
    if name=='lines':
        p=pdf.copy(); p*= (128/4)/p.sum(); p=np.minimum(p,1); mk=rng.random(128)<p; mk[np.abs(k)<5]=True; Mk=np.repeat(mk[:,None],128,1)
    else:
        KX,KY=np.meshgrid(k,k); r=np.sqrt(KX**2+KY**2); p=np.exp(-(r/25.)**2)*0.9+0.05; p*=(128*128/4)/p.sum(); p=np.minimum(p,1)
        Mk=rng.random((128,128))<p; Mk[r<5]=True
    print(name,'R=',Mk.size/Mk.sum())
    y=mt.add_complex_noise(mt.fft2c(img),0.002)*Mk
    print(' zf',err(mt.ifft2c(y)))
    for lam in [0.003,0.01,0.03]:
        for lev in [3,5]:
            t=time.time(); print(' ',lam,lev,err(fista(y,Mk,lam,levels=lev)), round(time.time()-t,2))
