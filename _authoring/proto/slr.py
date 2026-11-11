import sys, time
sys.path.insert(0, "/home/dominik/Documents/other/mri")
import numpy as np
import mri_tools as mt
from scipy.io import loadmat


def cayley_klein(rf, phase):
    """rf: (Nt,) complex rotation angles [rad] (b1*gamma*dt). phase: (Nt, Ns) precession angles [rad] per step.
    Returns alpha, beta (Ns,). Each step: simultaneous rotation about [Re rf, Im rf, phase] (left-handed)."""
    Ns = phase.shape[1]
    a = np.ones(Ns, complex); b = np.zeros(Ns, complex)
    for n in range(len(rf)):
        bx, by, bz = rf[n].real, rf[n].imag, phase[n]
        phi = np.sqrt(bx**2 + by**2 + bz**2)
        safe = np.where(phi == 0, 1.0, phi)
        nx, ny, nz = bx / safe, by / safe, np.where(phi == 0, 1.0, bz / safe)
        c, s = np.cos(phi / 2), np.sin(phi / 2)
        ar = c + 1j * nz * s
        br = 1j * (nx + 1j * ny) * s
        a, b = ar * a - np.conj(br) * b, br * a + np.conj(ar) * b
    return a, b


def dinf(d1, d2):
    a1, a2, a3, a4, a5, a6 = 5.309e-3, 7.114e-2, -4.761e-1, -2.66e-3, -5.941e-1, -4.278e-1
    l1, l2 = np.log10(d1), np.log10(d2)
    return (a1 * l1**2 + a2 * l1 + a3) * l2 + (a4 * l1**2 + a5 * l1 + a6)


def fir_ls(n, f_pass, f_stop, w_stop=1.0, n_grid=4096):
    """Linear-phase low-pass FIR (n taps) by weighted least squares; f in cycles/sample (0..0.5)."""
    f = np.linspace(0, 0.5, n_grid)
    band_p = f <= f_pass; band_s = f >= f_stop
    ff = np.concatenate([f[band_p], f[band_s]])
    d = np.concatenate([np.ones(band_p.sum()), np.zeros(band_s.sum())])
    w = np.concatenate([np.ones(band_p.sum()), w_stop * np.ones(band_s.sum())])
    m = np.arange(n) - (n - 1) / 2
    A = np.cos(2 * np.pi * np.outer(ff, m))
    h, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], d * np.sqrt(w), rcond=None)
    return h


def dzls(n, tbw, d1=0.01, d2=0.01):
    di = dinf(d1, d2)
    w = di / tbw
    fp = (1 - w) * (tbw / 2) / n   # cycles/sample   ((tbw/2)/(n/2) * 0.5)
    fs = (1 + w) * (tbw / 2) / n
    return fir_ls(n, fp * 1, fs * 1, d1 / d2)


def b2a(b):
    n = len(b)
    M = 16 * 2 ** int(np.ceil(np.log2(n)) + 1)
    B = np.fft.fft(b, M)
    bmx = np.abs(B).max()
    if bmx >= 1:
        B = B / (bmx + 1e-7)
    am = np.sqrt(1 - np.abs(B) ** 2)
    c = np.fft.ifft(np.log(am))
    fold = np.zeros(M); fold[0] = 1; fold[1:M // 2] = 2; fold[M // 2] = 1
    A = np.exp(np.fft.fft(c * fold))
    return np.fft.ifft(A)[:n]


def ab2rf(a, b):
    a = np.array(a, complex); b = np.array(b, complex)
    N = len(a)
    rf = np.zeros(N, complex)
    for j in range(N - 1, -1, -1):
        r = b[0] / a[0]
        C = 1 / np.sqrt(1 + abs(r) ** 2)
        S = C * r
        phi = 2 * np.arctan2(abs(S), C)
        rf[j] = phi * (-1j) * S / abs(S) if abs(S) > 0 else 0
        an = C * a + np.conj(S) * b
        bn = -S * a + C * b
        a = an[:-1]
        b = bn[1:]
    return rf


def dzrf(n, tbw, ptype="st", d1=0.01, d2=0.01):
    if ptype == "st":
        bsf = 1
    elif ptype == "ex":
        bsf = np.sqrt(0.5); d1 = np.sqrt(d1 / 2); d2 = d2 / np.sqrt(2)
    elif ptype == "se":
        bsf = 1; d1 = d1 / 4; d2 = np.sqrt(d2)
    elif ptype == "inv":
        bsf = 1; d1 = d1 / 8; d2 = np.sqrt(d2 / 2)
    elif ptype == "sat":
        bsf = np.sqrt(0.5); d1 = d1 / 2; d2 = np.sqrt(d2)
    b = dzls(n, tbw, d1, d2)
    if ptype == "st":
        return b
    b = bsf * b
    return ab2rf(b2a(b), b)


if __name__ == "__main__":
    ex = loadmat("/home/dominik/Documents/other/mri/data/RFpulse_excitation_tbw6_SLR.mat")["rf_90_SLR"].ravel()
    inv = loadmat("/home/dominik/Documents/other/mri/data/RFpulse_inversion_tbw6_SLR.mat")["rf_inv_SLR"].ravel()
    dt = 1e-4
    rf_ex = dzrf(60, 6, "ex"); rf_inv = dzrf(60, 6, "inv")
    print("sum ex", rf_ex.sum(), "flip", np.abs(rf_ex).sum(), " mat:", mt.GAMMA * dt * np.sum(ex) * 1e-3)
    mat_ex_rad = ex * 1e-3 * mt.GAMMA * dt
    mat_inv_rad = inv * 1e-3 * mt.GAMMA * dt
    print(np.round(rf_ex.real[20:40], 4))
    print(np.round(mat_ex_rad.real[20:40], 4))
    print(np.round(rf_inv.real[20:40], 4))
    print(np.round(mat_inv_rad.real[20:40], 4))
    print("max imag", np.abs(rf_ex.imag).max(), np.abs(rf_inv.imag).max())
    # compare cayley-klein with bloch
    df = np.linspace(-3000, 3000, 101)
    phase = np.tile(2 * np.pi * df * dt, (60, 1))
    a, b = cayley_klein(rf_ex, phase)
    mxy = 2 * np.conj(a) * b
    t0 = time.time()
    _, M = mt.bloch_simulate(B1=rf_ex / (mt.GAMMA * dt), dt=dt, df=df, return_all=False)
    print("bloch time", time.time() - t0)
    print("max diff", np.abs(mxy - (M[:, 0] + 1j * M[:, 1])).max(), np.abs((1 - 2 * abs(b) ** 2) - M[:, 2]).max())
    print(np.round(np.abs(mxy[::5]), 3))
    # refocusing check: start M=(1,0,0) and (0,1,0)
    for m0 in [(1, 0, 0), (0, 1, 0)]:
        _, M = mt.bloch_simulate(B1=rf_inv / (mt.GAMMA * dt), dt=dt, df=df, M_init=m0, return_all=False)
        a, b = cayley_klein(rf_inv, phase)
        m = m0[0] + 1j * m0[1]
        pred = np.conj(a) ** 2 * m - b ** 2 * np.conj(m)
        print("refoc diff", np.abs(pred - (M[:, 0] + 1j * M[:, 1])).max())
