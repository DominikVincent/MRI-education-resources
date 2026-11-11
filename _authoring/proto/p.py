import sys; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt, time
# M0
N=6.7e28; B0=3; T=310
M0 = N*mt.GAMMA**2*mt.HBAR**2*B0/(4*mt.K_BOLTZMANN*T); print('M0',M0,'A/m', 'mu0M0',mt.MU0*M0)
# magnet design: long solenoid of loops vs optimized pairs
R=0.45; z=np.linspace(-0.25,0.25,101)
zl=np.linspace(-0.8,0.8,161)  # 161 loops
Bsol=sum(mt.loop_field_on_axis(z-zz,R,1.0) for zz in zl)
I=3/Bsol[50]; Bsol*=I
print('solenoid current per loop',I,'ppm pp', (Bsol.max()-Bsol.min())/3*1e6)
# 6 symmetric pairs, least squares
zp=np.array([0.08,0.24,0.40,0.56,0.72,0.85])
A=np.stack([mt.loop_field_on_axis(z-p,R,1)+mt.loop_field_on_axis(z+p,R,1) for p in zp],1)
cur,*_=np.linalg.lstsq(A,3*np.ones_like(z),rcond=None)
B=A@cur; print('currents',cur.round(0),'ppm pp',(B.max()-B.min())/3*1e6)
zt=np.linspace(-0.25,0.25,11)
# eval over wider range
zw=np.linspace(-0.3,0.3,7); Aw=np.stack([mt.loop_field_on_axis(zw-p,R,1)+mt.loop_field_on_axis(zw+p,R,1) for p in zp],1); print(((Aw@cur)/3-1)*1e6)
# loop sensitivities
x=np.linspace(-0.15,0.15,64); X,Y=np.meshgrid(x,x); P=np.stack([X,Y,0*X],-1)
t0=time.time()
for c,nrm in [((0.13,0,0),'x'),((-0.13,0,0),'x'),((0,0.13,0),'y'),((0,-0.13,0),'y')]:
    Bc=mt.loop_coil_field(P,0.05,1.0,center=c,normal=nrm)
print('loops',time.time()-t0)
