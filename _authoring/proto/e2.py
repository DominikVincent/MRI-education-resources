import sys; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt, time
def sinc_pulse(TBW, T_rf, dt, theta=np.pi/2):
    n = int(round(T_rf / dt)); t = (np.arange(n) + 0.5) * dt; x = (t - T_rf / 2) / T_rf
    shape = np.sinc(TBW * x) * (0.54 + 0.46 * np.cos(2 * np.pi * x))
    return t, (theta * shape / (mt.GAMMA * shape.sum() * dt)).astype(complex)
dt=10e-6; T=2e-3; TBW=4; dz=5e-3; BW=TBW/T; G=BW/(mt.GAMMA_BAR*dz)
print('G mT/m', G*1e3)
z=np.linspace(-15e-3,15e-3,301); r=np.zeros((301,3)); r[:,2]=z
for th in [30,90]:
    t,b1=sinc_pulse(TBW,T,dt,np.deg2rad(th))
    n=len(b1)
    res=[]
    for frac in np.linspace(0,1,21):
        nrw=int(round(frac*n))
        Gw=np.zeros((n+100,3)); Gw[:n,2]=G
        # rewinder as lobe of amplitude -G*... use 100 steps with area -frac*G*T
        Gw[n:,2]=-frac*G*T/(100*dt)
        _,M=mt.bloch_simulate(B1=np.concatenate([b1,np.zeros(100)]),dt=dt,G=Gw,r=r,return_all=False)
        mxy=M[:,0]+1j*M[:,1]
        res.append(abs(mxy.sum()))
    res=np.array(res); print(th, np.round(res/res.max(),3), 'best', np.linspace(0,1,21)[res.argmax()])
# finer
t,b1=sinc_pulse(TBW,T,dt,np.pi/2); n=len(b1)
fr=np.linspace(0.45,0.56,23); out=[]
for frac in fr:
    Gw=np.zeros((n+100,3)); Gw[:n,2]=G; Gw[n:,2]=-frac*G*T/(100*dt)
    _,M=mt.bloch_simulate(B1=np.concatenate([b1,np.zeros(100)]),dt=dt,G=Gw,r=r,return_all=False)
    out.append(abs((M[:,0]+1j*M[:,1]).sum()))
print('fine best 90', fr[np.argmax(out)])
# chemical shift
for BWx in [1000,4000]:
    Tx=TBW/BWx; Gx=BWx/(mt.GAMMA_BAR*dz); t,b1=sinc_pulse(TBW,Tx,dt,np.pi/2); n=len(b1)
    for dfv in [0,-447]:
        Gw=np.zeros((n+n//2,3)); Gw[:n,2]=Gx; Gw[n:,2]=-Gx
        _,M=mt.bloch_simulate(B1=np.concatenate([b1,np.zeros(n//2)]),dt=dt,G=Gw,r=r,df=dfv,return_all=False)
        p=np.abs(M[:,0]+1j*M[:,1]); c=(p*z).sum()/p.sum()
        print(BWx, dfv, 'centroid mm', c*1e3, 'pred', -dfv/(mt.GAMMA_BAR*Gx)*1e3)
