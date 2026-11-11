import numpy as np, sys, time
sys.path.insert(0,'/home/dominik/Documents/other/mri')
import mri_tools as mt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
GB=mt.GAMMA_BAR
T=10e-3; dt=4e-6; Nt=int(T/dt); n_turns=12; kmax=60.0  # 1/m
t=np.arange(Nt)*dt
# spiral-in, constant angular velocity, ends at k=0 at t=T
u = 1 - (t+dt/2)/T
kx = kmax*u*np.cos(2*np.pi*n_turns*u); ky = kmax*u*np.sin(2*np.pi*n_turns*u)
# gradient: k(t) = -GB int_t^T G => dk/dt = GB G(t)
k = np.stack([kx,ky],1)
G = np.gradient(k, dt, axis=0)/GB
print('Gmax mT/m', np.abs(G).max()*1e3, 'slew', np.abs(np.diff(G,axis=0)/dt).max())
# desired profile: gaussian, FWHM 4 cm -> W(k) = FT
fwhm=0.04; sig=fwhm/2.355
W = np.exp(-2*(np.pi*sig)**2*(kx**2+ky**2))
kdot = np.linalg.norm(np.gradient(k,dt,axis=0),axis=1)
b1 = W*kdot
flip=np.deg2rad(30)
# scale: small tip Mxy(0)=i gamma M0 sum b1 dt ~ sin flip at center -> gamma sum b1 dt = flip
b1 = b1*flip/(mt.GAMMA*np.sum(b1)*dt)
print('peak b1 uT', b1.max()*1e6)
N=64; x=np.linspace(-0.15,0.15,N); X,Y=np.meshgrid(x,x)
r=np.stack([X.ravel(),Y.ravel(),0*X.ravel()],1)
Gf=np.concatenate([G,np.zeros((Nt,1))],1)
t0=time.time()
_,M=mt.bloch_simulate(B1=b1.astype(complex),dt=dt,G=Gf,r=r,return_all=False)
print('time',time.time()-t0)
mxy=(M[:,0]+1j*M[:,1]).reshape(N,N)
plt.subplot(121); plt.imshow(np.abs(mxy),extent=[-15,15,-15,15]); plt.colorbar()
plt.subplot(122); plt.plot(kx,ky,lw=.5); plt.axis('equal')
plt.savefig('j1.png')
print(np.abs(mxy[N//2,::4]).round(3), np.angle(mxy[N//2,N//2]))
