import sys; sys.path.insert(0,'/home/dominik/Documents/other/mri'); sys.path.insert(0,'.')
import numpy as np, mri_tools as mt, slr
from scipy.signal import remez
# HS pulse
T=10e-3; dt=10e-6; n=int(T/dt); tau=(np.arange(n)+0.5)*dt - T/2
kap=800.0; mu=5.0
dw = -mu*kap*np.tanh(kap*tau)            # rad/s RF frequency offset
phi = np.cumsum(dw)*dt
env = 1/np.cosh(kap*tau)
thr = np.sqrt(mu)*kap/mt.GAMMA
print('threshold uT', thr*1e6, 'BW Hz', mu*kap/np.pi)
df=np.linspace(-2000,2000,201)
for A0 in [0.5*thr, thr, 2*thr, 3*thr, 6*thr]:
    b1 = A0*env*np.exp(-1j*phi)
    _,M=mt.bloch_simulate(B1=b1,dt=dt,df=df,return_all=False)
    print('A0 %.1f uT'%(A0*1e6), 'Mz center', M[100,2].round(3), 'Mz at 500Hz', M[np.argmin(abs(df-500)),2].round(3), 'Mz at -500', M[np.argmin(abs(df+500)),2].round(3),'at 1500', M[np.argmin(abs(df-1500)),2].round(3))
# SLR inversion B1 scaling
def dz(n,tbw,p):
    d1=d2=0.01
    if p=='inv': bsf=1; d1=d1/8; d2=np.sqrt(d2/2)
    b=1j*bsf*slr.dzls(n,tbw,d1,d2); return slr.ab2rf(slr.b2a(b),b)
rf=dz(128,6,'inv'); Ts=6e-3; dts=Ts/128
for sc in [0.6,0.8,1,1.2,1.5]:
    _,M=mt.bloch_simulate(B1=sc*rf/(mt.GAMMA*dts),dt=dts,df=[0.0],return_all=False); print('SLR scale',sc, M[0,2].round(3))
# remez even taps
h=remez(60,[0,0.04,0.06,0.5],[1,0],weight=[1,1],fs=1); print('remez ok', h.sum())
