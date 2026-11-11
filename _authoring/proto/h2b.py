import sys,time; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
t0=time.time()
N=96
ph=mt.tissue_phantom(N)
m=ph['M0'].astype(complex)
y,x=np.mgrid[1:-1:N*1j,-1:1:N*1j]
df=150*np.exp(-((x-0.0)**2+(y-0.55)**2)/(2*0.15**2))  # frontal-sinus-like hot spot
fat=ph['labels']==1
df_tot=df+np.where(fat,mt.ppm_to_hz(-3.5,3.0),0)
T2s=np.where(ph['labels']==4,0.2,0.04)
def epi(m,dfmap,t_esp,TE,T2s=None):
    K=np.zeros((N,N),complex)
    for p in range(N):
        t=TE+(p-N//2)*t_esp
        w=m*np.exp(-1j*2*np.pi*dfmap*t)
        if T2s is not None: w=w*np.exp(-t/T2s)
        K[p]=mt.fft2c(w)[p]
    return K
for t_esp in [0.25e-3,0.5e-3]:
    K=epi(m,df_tot,t_esp,30e-3,T2s)
    print(t_esp, 150*N*t_esp, 440*N*t_esp)
img=mt.ifft2c(K)
fig,ax=plt.subplots(1,3,figsize=(12,4)); mt.show_image(m,ax=ax[0]); ax[1].imshow(df_tot); mt.show_image(img,ax=ax[2]); plt.savefig('epi.png')
K0=epi(m,0*df,0.5e-3,30e-3)
for phi in [0.1,0.3]:
    Kg=K0.copy(); Kg[1::2]*=np.exp(1j*phi)
    im=mt.ifft2c(Kg); ghost=np.abs(np.roll(im,N//2,axis=0))
    print(phi, np.abs(np.sin(phi/2))/np.abs(np.cos(phi/2)), np.abs(im[N//2+N//2-5:N//2+N//2-4]).max() )
print(time.time()-t0)
