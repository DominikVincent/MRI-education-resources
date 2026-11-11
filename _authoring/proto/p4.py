import sys, time
sys.path.insert(0, "/home/dominik/Documents/other/mri")
import numpy as np, mri_tools as mt
t0=time.time()
def wire_field(points, path, current=1.0):
    pts = np.asarray(points, float); P = np.asarray(path, float)
    mid = 0.5*(P[1:]+P[:-1]); dl = P[1:]-P[:-1]
    flat = pts.reshape(-1,3); B = np.zeros_like(flat)
    for i in range(0, flat.shape[0], 2048):
        Rv = flat[i:i+2048,None,:]-mid[None]
        Rn = np.maximum(np.linalg.norm(Rv,axis=-1,keepdims=True),1e-4)
        B[i:i+2048] = (mt.MU0*current/(4*np.pi)*np.cross(dl[None],Rv)/Rn**3).sum(1)
    return B.reshape(pts.shape)
# Maxwell
a=0.3
z=np.linspace(-0.3,0.3,61); pts=np.stack([0*z,0*z,z],-1)
for d in [np.sqrt(3)*a, a]:
    B = mt.loop_coil_field(pts,a,1.0,center=(0,0,d/2)) + mt.loop_coil_field(pts,a,-1.0,center=(0,0,-d/2))
    G0 = np.gradient(B[:,2],z)[30]
    print("d/a",d/a,"G0 per A",G0, "theory", 0.6413*mt.MU0/a**2)
# Golay saddle
def saddle(a, z1, z2, phic, span=np.deg2rad(120), n=60):
    ph = np.linspace(phic-span/2, phic+span/2, n)
    arc1 = np.stack([a*np.cos(ph), a*np.sin(ph), np.full(n,z1)],-1)
    arc2 = arc1[::-1].copy(); arc2[:,2]=z2
    return np.vstack([arc1, arc2, arc1[:1]])
x=np.linspace(-0.2,0.2,41); pts=np.stack([x,0*x,0*x],-1)
for signs in [(1,1,-1,-1),(1,-1,-1,1)]:
    B=0
    for s,(phic,zz) in zip(signs,[(0,(0.389,2.57)),(0,(-0.389,-2.57)),(np.pi,(0.389,2.57)),(np.pi,(-0.389,-2.57))]):
        B = B + wire_field(pts, saddle(a, zz[0]*a, zz[1]*a, phic), s)
    print(signs, B[[0,10,20,30,40],2], np.gradient(B[:,2],x)[20])
print(time.time()-t0)
