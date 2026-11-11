import sys,time; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
t0=time.time()
def slice_profile(TBW, T_rf=2e-3, dz=5e-3, theta=np.pi/2, z=np.linspace(-15e-3,15e-3,601), dt=10e-6):
    n=int(round(T_rf/dt)); tt=(np.arange(n)-n/2+0.5)*dt
    env=np.sinc(TBW*tt/T_rf)*(0.5+0.5*np.cos(2*np.pi*tt/T_rf))
    b1=env/np.sum(env)*theta/(mt.GAMMA*dt)
    BW=TBW/T_rf; G=BW/(mt.GAMMA_BAR*dz)
    Gw=np.zeros((n,3)); Gw[:,2]=G
    r=np.zeros((len(z),3)); r[:,2]=z
    _,M=mt.bloch_simulate(B1=b1.astype(complex),dt=dt,G=Gw,r=r,return_all=False)
    return np.arccos(np.clip(M[:,2],-1,1))
z=np.linspace(-15e-3,15e-3,601)
for TBW in [2,4,8]:
    th=slice_profile(TBW)
    print(TBW, np.degrees(th[300]), np.degrees(th[np.argmin(abs(z-2.5e-3))]), np.degrees(th[np.argmin(abs(z-5e-3))]))
print(time.time()-t0)
