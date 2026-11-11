import sys; sys.path.insert(0, "/home/dominik/Documents/other/mri")
import numpy as np, mri_tools as mt
a_p, a_s, r_cryo = 0.30, 0.40, 0.45
for nz, zl, zc in [(12,0.6,0.8),(20,0.8,0.8),(24,0.9,0.8),(30,1.0,0.9)]:
    z_shield = np.linspace(-zl, zl, nz)
    z_cryo = np.linspace(-zc, zc, 81)
    cryo_pts = np.stack([np.full_like(z_cryo, r_cryo), 0*z_cryo, z_cryo], -1)
    def primary_B(p):
        return (mt.loop_coil_field(p, a_p, 1.0, center=(0,0,np.sqrt(3)/2*a_p)) + mt.loop_coil_field(p, a_p, -1.0, center=(0,0,-np.sqrt(3)/2*a_p)))
    Bs = np.stack([mt.loop_coil_field(cryo_pts, a_s, 1.0, center=(0,0,zs)) for zs in z_shield], -1)
    for comps in [(0,2),(0,)]:
        A = np.concatenate([Bs[:,c,:] for c in comps]); bb = -np.concatenate([primary_B(cryo_pts)[:,c] for c in comps])
        I = np.linalg.lstsq(A, bb, rcond=None)[0]
        Bt = primary_B(cryo_pts) + np.einsum("pci,i->pc", Bs, I)
        z = np.linspace(-0.05,0.05,11); pz=np.stack([0*z,0*z,z],-1)
        Bz = primary_B(pz)[:,2] + sum(Ii*mt.loop_coil_field(pz,a_s,1.0,center=(0,0,zs))[:,2] for Ii,zs in zip(I,z_shield))
        print(nz, zl, zc, comps, np.abs(primary_B(cryo_pts)).max(), np.abs(Bt).max(), np.abs(Bt[:,0]).max(), np.gradient(Bz,z)[5], np.abs(I).max())
