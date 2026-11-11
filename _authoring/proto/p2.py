import sys; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
from scipy.optimize import nnls
R=0.45; z=np.linspace(-0.25,0.25,101)
for L in [0.6,0.8]:
  zp=np.linspace(0.05,L,16)
  A=np.stack([mt.loop_field_on_axis(z-p,R,1)+mt.loop_field_on_axis(z+p,R,1) for p in zp],1)
  cur,res=nnls(A,3*np.ones_like(z))
  B=A@cur; print(L,'currents kA-turns',(cur/1e3).round(0),'ppm pp',(B.max()-B.min())/3*1e6, 'total',cur.sum()*2/1e6,'MA-turns')
