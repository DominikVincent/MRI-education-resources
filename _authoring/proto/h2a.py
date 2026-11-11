import sys,time; sys.path.insert(0,'/home/dominik/Documents/other/mri')
import numpy as np, mri_tools as mt
t0=time.time()
def rf_matrix(theta, phase=0.0):
    n=[np.cos(phase),np.sin(phase),0.0]
    return mt.rotation_matrix(n,-theta)
def bssfp_ss(theta, TR, T1, T2, df, TE=None, alt=True):
    TE = TR/2 if TE is None else TE
    df=np.atleast_1d(df)
    beta = 2*np.pi*df*TR + (np.pi if alt else 0)
    E1,E2=np.exp(-TR/T1),np.exp(-TR/T2)
    R=rf_matrix(theta)
    out=[]
    for b in beta:
        A=np.diag([E2,E2,E1])@mt.Rz(-b)
        bb=np.array([0,0,1-E1])
        Mp=np.linalg.solve(np.eye(3)-R@A, R@bb)
        # evolve to TE
        bh=b*TE/TR
        Ah=np.diag([np.exp(-TE/T2)]*2+[np.exp(-TE/T1)])@mt.Rz(-bh)
        M=Ah@Mp+np.array([0,0,1-np.exp(-TE/T1)])
        out.append(M[0]+1j*M[1])
    return np.array(out)
th=np.radians(60); TR=5e-3; T1=0.6;T2=0.08
s=bssfp_ss(th,TR,T1,T2,0.0)
print(abs(s), mt.bssfp_signal(th,TR/2,TR,T1,T2))
df=np.linspace(-1/TR,1/TR,9); print(np.round(abs(bssfp_ss(th,TR,T1,T2,df)),3))
# spoiling sim
def gre_sim(thetas, TR, T1, T2, n_iso=100, n_tr=400, rf_spoil=False, inc=np.radians(117)):
    thetas=np.atleast_1d(thetas); nf=len(thetas)
    M=np.zeros((nf,n_iso,3)); M[...,2]=1
    dphi=2*np.pi*np.arange(n_iso)/n_iso
    E1,E2=np.exp(-TR/T1),np.exp(-TR/T2)
    phase=0.0; Splus=Sminus=None
    for n in range(n_tr):
        if rf_spoil: phase = phase + n*inc
        Rs=np.stack([rf_matrix(t,phase) for t in thetas])  # (nf,3,3)
        M=np.einsum('fij,fsj->fsi',Rs,M)
        Splus=np.abs(np.mean(M[...,0]+1j*M[...,1],axis=1))
        mxy=(M[...,0]+1j*M[...,1])*np.exp(-1j*dphi)*E2
        M=np.stack([mxy.real,mxy.imag,1+(M[...,2]-1)*E1],-1)
        Sminus=np.abs(np.mean(M[...,0]+1j*M[...,1],axis=1))
    return Splus,Sminus
fl=np.radians(np.arange(2,62,4))
for T1,T2 in [(0.6,0.08),(4.0,2.0)]:
    a=gre_sim(fl,10e-3,T1,T2)[0]; b=gre_sim(fl,10e-3,T1,T2,rf_spoil=True)[0]
    c=mt.spoiled_gre_signal(fl,0,10e-3,T1,np.inf)
    print(np.round(a,3));print(np.round(b,3));print(np.round(c,3))
print(time.time()-t0)
