import sys; sys.path.insert(0, "/home/dominik/Documents/other/mri")
import numpy as np, mri_tools as mt
from scipy.linalg import eigh
from scipy.special import jv, hankel2
# mutual inductance of coplanar loops
def loop_pts(a, c, n=200):
    ph = np.linspace(0, 2*np.pi, n+1)
    P = np.stack([c[0]+a*np.cos(ph), c[1]+a*np.sin(ph), 0*ph],-1)
    return 0.5*(P[1:]+P[:-1]), P[1:]-P[:-1]
def mutual(a, dist):
    m1, d1 = loop_pts(a,(0,0)); m2,d2 = loop_pts(a,(dist,0))
    R = np.linalg.norm(m1[:,None]-m2[None],axis=-1)
    return mt.MU0/(4*np.pi)*np.sum((d1@d2.T)/np.maximum(R,1e-4))
a=0.05
ds = np.linspace(0.5*a, 3*a, 51)
Ms = np.array([mutual(a,d) for d in ds])
i = np.where(np.diff(np.sign(Ms)))[0]
print("zero crossing at d/a ~", ds[i]/a, Ms[:3], Ms[-3:])
# birdcage
N=16; Lr=300e-9; Le=40e-9
def circ(d, o):
    M = np.diag(np.full(N, d)); 
    for j in range(N):
        M[j,(j+1)%N]=o; M[j,(j-1)%N]=o
    return M
Lm = circ(2*Le+2*Lr, -Lr)
w0 = 2*np.pi*127.7e6
k=1; s2=np.sin(np.pi*k/N)**2
Cr = 1/(w0**2*(Lr+Le/(2*s2)))
Sm = circ(2/Cr, -1/Cr)
w2, V = eigh(Sm, Lm)
print("low-pass Cr", Cr, np.sqrt(np.abs(w2))/2/np.pi/1e6)
Ce = 1/(w0**2*(Le+2*Lr*s2))
Sm = np.diag(np.full(N,2/Ce)); w2,V=eigh(Sm,Lm); print("high-pass Ce", Ce, np.sqrt(w2)/2/np.pi/1e6)
# tissue k
eps0=8.854e-12
for f, er, sig in [(63.9e6,72,0.69),(127.7e6,63,0.72),(298e6,58,0.77)]:
    w=2*np.pi*f; kk = w*np.sqrt(mt.MU0*eps0*(er - 1j*sig/(w*eps0)))
    print(f/1e6, kk, "lambda", 2*np.pi/kk.real, "atten len", -1/kk.imag, abs(jv(0,kk*0.09)), abs(jv(0,kk*0.17)))
