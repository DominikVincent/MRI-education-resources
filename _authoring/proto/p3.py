import sys; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
from scipy.optimize import nnls
R1=0.45; z=np.linspace(-0.25,0.25,101)
zp=np.linspace(0.05,0.6,16)
A=np.stack([mt.loop_field_on_axis(z-p,R1,1)+mt.loop_field_on_axis(z+p,R1,1) for p in zp],1)
cur,_=nnls(A,3*np.ones_like(z))
zs=np.array([-0.75,-0.25,0.25,0.75]); R2=0.85
def Bax(zz,shield):
    B=sum(c*(mt.loop_field_on_axis(zz-p,R1,1)+mt.loop_field_on_axis(zz+p,R1,1)) for c,p in zip(cur,zp))
    if shield:
        Is=-(2*cur.sum()*R1**2)/(len(zs)*R2**2)
        B=B+sum(Is*mt.loop_field_on_axis(zz-p,R2,1) for p in zs)
    return B
zf=np.linspace(0,12,2401)
for sh in [False,True]:
    B=Bax(zf,sh); c=Bax(np.array([0.]),sh)[0]; B=B*3/c
    print(sh,'center',c,'5G at',zf[np.where(np.abs(B)>5e-4)[0].max()])
# Maxwell pair
R=0.35; d=np.sqrt(3)/2*R
zz=np.linspace(-0.3,0.3,7)
f=lambda z: mt.loop_field_on_axis(z-d,R,1)-mt.loop_field_on_axis(z+d,R,1)
h=1e-4; print('grad per A-turn', (f(h)-f(-h))/(2*h))
