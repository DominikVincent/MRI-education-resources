import sys; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
from scipy.special import ellipk, ellipe
R=0.5
def coefs(zc,a,nm):
    z=np.linspace(-0.05,0.05,201)
    b=mt.loop_field_on_axis(z-zc,a,1)+mt.loop_field_on_axis(z+zc,a,1)
    return np.polyfit(z/0.05,b,2*nm)[::-1][0:2*nm:2]
zs=[0.05,0.16,0.28,0.42,0.58,0.8]; nm=6
A=np.array([coefs(zc,R,nm) for zc in zs]).T; rhs=np.zeros(nm); rhs[0]=1.5; I0=np.linalg.solve(A,rhs)
Rs=0.9; zsh=0.65
A2=np.array([coefs(zc,R,nm) for zc in zs]+[coefs(zsh,Rs,nm)]).T
A2=np.vstack([A2,[2*np.pi*R**2]*6+[2*np.pi*Rs**2]]); rhs=np.zeros(nm+1); rhs[0]=1.5; I=np.linalg.solve(A2,rhs)
zl=zs+[zsh]; Rl=[R]*6+[Rs]
def mutual(R1,R2,d):
    k2=4*R1*R2/((R1+R2)**2+d**2); k=np.sqrt(k2)
    return mt.MU0*np.sqrt(R1*R2)*((2/k-k)*ellipk(k2)-2/k*ellipe(k2))
def energy(zl,Il,Rl,a=0.03):
    zc=np.concatenate([zl,-np.array(zl)]); Ic=np.concatenate([Il,Il]); Rc=np.concatenate([Rl,Rl])
    E=0
    for i in range(len(zc)):
        for j in range(len(zc)):
            M = mt.MU0*Rc[i]*(np.log(8*Rc[i]/a)-1.75) if i==j else mutual(Rc[i],Rc[j],zc[i]-zc[j])
            E+=0.5*Ic[i]*Ic[j]*M
    return E
for a in [0.02,0.03,0.05]:
    print(a,'unshielded E MJ',energy(zs,I0,[R]*6,a)/1e6,' shielded',energy(zl,I,Rl,a)/1e6)
print('ideal estimate', 1.5**2/(2*mt.MU0)*np.pi*0.5**2*1.6/1e6)
