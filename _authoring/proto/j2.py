import numpy as np, sys, time
sys.path.insert(0,'.')
import slr, mri_tools as mt
from scipy.interpolate import CubicSpline
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

def verse(g, rf):
    g = np.asarray(g, float); rf = np.asarray(rf, complex); m = len(rf)
    k = np.concatenate([[0], np.cumsum(g[:-1])]) + g / 2
    k = m * k / k.max()
    gn = m * g / g.sum()
    x = np.arange(1, m + 1) - 0.5
    return gn * (CubicSpline(x, rf.real)(k) + 1j * CubicSpline(x, rf.imag)(k))

GB = mt.GAMMA_BAR
slewmax = 150.0; gmax = 40e-3
T = 20e-3; TBW = 3; SBW = 6; dz = 5e-3; Nspec = 10; DT = T / Nspec; dt = 10e-6; flip = np.pi / 4
verse_frac = 0.8
Nlobe = int(round(DT / dt)); Nspat = int(round(Nlobe * verse_frac)); Nwait = (Nlobe - Nspat) / 2
rfspec = slr.dzrf(Nspec, TBW)
rfspat = slr.dzrf(Nspat, SBW)
gamp = (SBW / DT) / (GB * dz)
print('gamp mT/m', gamp*1e3)
Nramp = min(int(np.ceil(gmax / slewmax / dt)), Nlobe // 2)
ramp = (np.arange(Nramp) + 0.5) / Nramp
g1 = np.concatenate([ramp, np.ones(Nlobe - 2 * Nramp), ramp[::-1]])
rfspatv = verse(g1, np.concatenate([np.zeros(int(np.floor(Nwait))), rfspat, np.zeros(int(np.ceil(Nwait)))]))
rf = np.concatenate([rfspatv * rfspec[k] for k in range(Nspec)])
g = np.concatenate([(1 if (k+1) % 2 else -1) * gamp * g1 for k in range(Nspec)])
Nrw = int(round((Nlobe - Nramp) / 2 + Nramp))
grw = np.concatenate([ramp, np.ones(Nrw - 2 * Nramp), ramp[::-1]])
g = np.concatenate([g, (1 if (Nspec+2) % 2 else -1) * gamp * grw])
rf = np.concatenate([rf, np.zeros(Nrw)])
rf = flip / rf.sum() * rf
print('len', len(rf), 'peak b1 uT', np.abs(rf).max()/(mt.GAMMA*dt)*1e6, 'Nramp', Nramp)
z = np.linspace(-8e-3, 8e-3, 201); df = np.linspace(-500, 500, 201)
Z, F = np.meshgrid(z, df, indexing='ij')
t0=time.time()
Ns = Z.size
a = np.ones(Ns, complex); b = np.zeros(Ns, complex)
zf, ff = Z.ravel(), F.ravel()
for n in range(len(rf)):
    bz = mt.GAMMA * g[n] * zf * dt + 2*np.pi*ff*dt
    bx, by = rf[n].real, rf[n].imag
    phi = np.sqrt(bx**2 + by**2 + bz**2)
    safe = np.where(phi == 0, 1.0, phi)
    c, s = np.cos(phi / 2), np.sin(phi / 2)/safe
    ar = c + 1j * bz * s; br = 1j * (bx + 1j * by) * s
    a, b = ar * a - np.conj(br) * b, br * a + np.conj(ar) * b
print('time', time.time()-t0)
mxy = (2*np.conj(a)*b).reshape(Z.shape)
plt.imshow(np.abs(mxy), extent=[df[0], df[-1], z[-1]*1e3, z[0]*1e3], aspect='auto'); plt.colorbar()
plt.savefig('j2.png')
print(np.abs(mxy[100, ::20]).round(3))
print(np.abs(mxy[::20, 100]).round(3))
# check bloch on center line
t0=time.time()
_, M = mt.bloch_simulate(B1=rf/(mt.GAMMA*dt), dt=dt, df=df, G=np.stack([0*g,0*g,g],1), r=np.zeros((len(df),3)), return_all=False)
print('bloch line time', time.time()-t0, np.abs(M[:,0]+1j*M[:,1] - mxy[100]).max())
