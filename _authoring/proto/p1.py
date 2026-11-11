import sys, time; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
# signal EMF estimate
N=6.7e28; T=310; B0=3
M0 = N*(mt.GAMMA*mt.HBAR)**2*B0/(4*mt.K_BOLTZMANN*T)
print('M0 A/m', M0, 'mu0M0 T', mt.MU0*M0)
# Maxwell pair
R=0.3; d=np.sqrt(3)/2*R
z=np.linspace(-0.3,0.3,61); pts=np.stack([0*z,0*z,z],-1)
t=time.time()
B=mt.loop_coil_field(pts,R,1.0,center=(0,0,d))+mt.loop_coil_field(pts,R,-1.0,center=(0,0,-d))
print(time.time()-t)
G=np.gradient(B[:,2],z)[30]; print('G per A',G)
lin = G*z; err=(B[:,2]-lin)/(G*0.15)
print(np.round(err[::5],3))
# finite solenoid homogeneity
def sol(z,R,L): return 0.5*((z+L/2)/np.sqrt(R**2+(z+L/2)**2)-(z-L/2)/np.sqrt(R**2+(z-L/2)**2))
z=np.linspace(-.25,.25,101)
for L in [1.0,1.5,2.0,3.0,5.0]:
    b=sol(z,0.5,L); print(L, (b.max()-b.min())/b[50]*1e6)
