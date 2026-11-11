"""
mri_tools.py -- shared helper functions for the "MRI: from intuition to engineering" notebooks.

Every notebook in this folder imports from here, so that all chapters use the
SAME conventions and units. Read this docstring once; it is the contract.

Units (SI everywhere unless a function says otherwise)
------------------------------------------------------
    time            s        (seconds; 1 ms = 1e-3 s)
    magnetic field  T        (tesla;   1 mT = 1e-3 T, 1 uT = 1e-6 T)
    gradient        T/m      (scanner specs are quoted in mT/m = 1e-3 T/m)
    frequency       Hz       (cycles per second)
    angular freq.   rad/s    (omega = 2*pi*f)
    position        m
    k-space         1/m      (cycles per metre)
    magnetization   normalized, i.e. M0 = 1 unless stated otherwise

Gyromagnetic ratio: two versions are used, and confusing them is the #1 source
of factor-2*pi errors in MRI:
    GAMMA_BAR  [Hz/T]     = 42.577e6 for 1H   ("gamma-bar")   f     = GAMMA_BAR * B
    GAMMA      [rad/(s T)] = 2*pi*GAMMA_BAR                    omega = GAMMA * B

Rotation convention
-------------------
The Bloch equation dM/dt = GAMMA * (M x B) makes M precess CLOCKWISE (seen
from the tip of B looking back), i.e. a LEFT-handed rotation about B.
For B along +z this means  M_xy(t) = M_xy(0) * exp(-i * omega * t),
where M_xy = M_x + i*M_y. All Bloch-simulation helpers below follow this.
`rotation_matrix(axis, angle)` itself is the ordinary RIGHT-handed rotation
from linear algebra; `bloch_rotate` simply uses angle -> -angle.

Fourier convention (matches MRI literature, e.g. Nishimura)
-----------------------------------------------------------
    F(k) = integral f(x) exp(-i 2 pi k x) dx
`fftc`/`ifftc` are "centered" FFTs: the zero position / zero frequency sits at
index N//2 (the middle of the array) in BOTH domains, which is how images and
k-space are displayed.

Content adapted and extended from P. Larson, "Principles of MRI",
https://github.com/LarsonLab/MRI-education-resources (MIT License).
"""

from __future__ import annotations

import numpy as np

# ----------------------------------------------------------------------------
# Physical constants
# ----------------------------------------------------------------------------
H_PLANCK = 6.62607015e-34       # Planck constant h [J s]  (J/Hz)
HBAR = H_PLANCK / (2 * np.pi)   # reduced Planck constant hbar [J s]
K_BOLTZMANN = 1.380649e-23      # Boltzmann constant k_B [J/K]
MU0 = 4e-7 * np.pi              # vacuum permeability mu_0 [T m / A] (~1.2566e-6)
AVOGADRO = 6.02214076e23        # [1/mol]

#: gyromagnetic ratios gamma-bar [Hz/T] of common MRI nuclei
GAMMA_BAR_NUCLEI = {
    "1H": 42.577e6,
    "2H": 6.536e6,
    "13C": 10.708e6,
    "19F": 40.078e6,
    "23Na": 11.262e6,
    "31P": 17.235e6,
    "129Xe": -11.777e6,
}
GAMMA_BAR = GAMMA_BAR_NUCLEI["1H"]   # [Hz/T]
GAMMA = 2 * np.pi * GAMMA_BAR        # [rad/(s T)]  ~ 2.675e8

FAT_SHIFT_PPM = -3.5   # main fat peak relative to water [ppm] (fat precesses slower)


def larmor_frequency(B0, gamma_bar=GAMMA_BAR):
    """Larmor (resonance) frequency f0 = gamma_bar * B0, in Hz. B0 in tesla."""
    return gamma_bar * np.asarray(B0)


def ppm_to_hz(ppm, B0, gamma_bar=GAMMA_BAR):
    """Convert a chemical shift in parts-per-million to Hz at field B0 [T]."""
    return np.asarray(ppm) * 1e-6 * gamma_bar * B0


# ----------------------------------------------------------------------------
# Rotations
# ----------------------------------------------------------------------------
def Rx(angle):
    """Right-handed 3x3 rotation matrix about the x-axis by `angle` [rad]."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)


def Ry(angle):
    """Right-handed 3x3 rotation matrix about the y-axis by `angle` [rad]."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=float)


def Rz(angle):
    """Right-handed 3x3 rotation matrix about the z-axis by `angle` [rad]."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)


def rotation_matrix(axis, angle):
    """Right-handed rotation by `angle` [rad] about the (not necessarily unit) `axis`.

    Uses Rodrigues' formula R = I cos(a) + sin(a) [n]_x + (1 - cos(a)) n n^T.
    """
    n = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(n)
    if norm == 0:
        return np.eye(3)
    n = n / norm
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) * np.cos(angle) + np.sin(angle) * K + (1 - np.cos(angle)) * np.outer(n, n)


def _rodrigues_many(M, n, phi):
    """Rotate many vectors M (..., 3) about unit axes n (..., 3) by angles phi (...) (right-handed)."""
    c = np.cos(phi)[..., None]
    s = np.sin(phi)[..., None]
    ndotM = np.sum(n * M, axis=-1, keepdims=True)
    return M * c + np.cross(n, M) * s + n * ndotM * (1 - c)


def bloch_rotate(M, dt, B, gamma=GAMMA):
    """Precess magnetization M around field B for a time dt (no relaxation).

    Parameters
    ----------
    M : array (..., 3)   magnetization vector(s) [Mx, My, Mz]
    dt : float           duration [s]
    B : array (..., 3)   magnetic field vector(s) [Bx, By, Bz] in tesla,
                         broadcastable against M (one field per spin is allowed)
    gamma : float        gyromagnetic ratio [rad/(s T)]

    Returns
    -------
    M_new : array (..., 3)

    The rotation is LEFT-handed (clockwise) by the angle gamma*|B|*dt, which is
    exactly what dM/dt = gamma M x B prescribes.
    """
    M = np.asarray(M, dtype=float)
    B = np.asarray(B, dtype=float)
    B, M = np.broadcast_arrays(B, M)
    Bnorm = np.linalg.norm(B, axis=-1)
    n = B / np.where(Bnorm == 0, 1.0, Bnorm)[..., None]
    phi = gamma * Bnorm * dt
    return _rodrigues_many(M, n, -phi)


def relax(M, dt, T1, T2, M0=1.0):
    """Apply T1 and T2 relaxation for a time dt [s] (no precession).

        Mx, My  ->  Mx, My * exp(-dt/T2)
        Mz      ->  M0 + (Mz - M0) * exp(-dt/T1)

    T1, T2, M0 may be scalars or arrays broadcastable to M[..., 0]. Use np.inf
    to switch a relaxation process off.
    """
    M = np.array(M, dtype=float, copy=True)
    E1 = np.exp(-dt / np.asarray(T1, dtype=float))
    E2 = np.exp(-dt / np.asarray(T2, dtype=float))
    M[..., 0] = M[..., 0] * E2
    M[..., 1] = M[..., 1] * E2
    M[..., 2] = M0 + (M[..., 2] - M0) * E1
    return M


def bloch_simulate(B1=None, dt=1e-6, n_steps=None, M_init=(0.0, 0.0, 1.0), T1=np.inf, T2=np.inf,
                   M0=1.0, df=0.0, G=None, r=None, B0=0.0, return_all=True, gamma=GAMMA):
    """Step-by-step Bloch simulation for one or many spins ("isochromats").

    The field seen by spin j at time step n is
        B = [Re b1[n], Im b1[n], B0 + df_j/gamma_bar + G[n] . r_j]
    so in the ROTATING FRAME use B0 = 0 (the default) and describe the RF by
    its complex envelope b1(t) = B1x(t) + i*B1y(t).

    Parameters
    ----------
    B1 : complex array (Nt,) or None   RF envelope b1(t) [T]; None = no RF
    dt : float                         time step [s]
    n_steps : int                      number of steps if B1 and G are None
    M_init : array (3,) or (Ns, 3)     starting magnetization
    T1, T2 : float or array (Ns,)      relaxation times [s] (np.inf = none)
    M0 : float or array (Ns,)          equilibrium magnetization
    df : float or array (Ns,)          off-resonance of each spin [Hz]
    G : array (Nt, 3) or None          gradient waveform [T/m]
    r : array (Ns, 3) or None          spin positions [m] (needed if G given)
    B0 : float                         static field along z [T] (0 in rotating frame)
    return_all : bool                  return every time step or only the end

    Returns
    -------
    t : array (Nt+1,)   times [s], starting at 0
    M : array (Nt+1, Ns, 3) if return_all else (Ns, 3)
    """
    lengths = [len(x) for x in (B1, G) if x is not None]
    Nt = lengths[0] if lengths else int(n_steps)
    b1 = np.zeros(Nt, complex) if B1 is None else np.asarray(B1, complex)

    df = np.atleast_1d(np.asarray(df, dtype=float))
    M = np.atleast_2d(np.asarray(M_init, dtype=float))
    Ns = max(M.shape[0], df.shape[0], 1 if r is None else np.atleast_2d(r).shape[0])
    M = np.broadcast_to(M, (Ns, 3)).copy()
    df = np.broadcast_to(df, (Ns,))
    T1 = np.broadcast_to(np.asarray(T1, float), (Ns,))
    T2 = np.broadcast_to(np.asarray(T2, float), (Ns,))
    M0 = np.broadcast_to(np.asarray(M0, float), (Ns,))
    if G is not None:
        G = np.asarray(G, dtype=float).reshape(Nt, 3)
        r = np.broadcast_to(np.atleast_2d(np.asarray(r, dtype=float)), (Ns, 3))

    out = np.empty((Nt + 1, Ns, 3)) if return_all else None
    if return_all:
        out[0] = M
    Bz_static = B0 + df / (gamma / (2 * np.pi))
    for n in range(Nt):
        Bz = Bz_static + (r @ G[n] if G is not None else 0.0)
        B = np.stack([np.full(Ns, b1[n].real), np.full(Ns, b1[n].imag), Bz], axis=-1)
        M = bloch_rotate(M, dt, B, gamma=gamma)
        M = relax(M, dt, T1, T2, M0)
        if return_all:
            out[n + 1] = M
    t = np.arange(Nt + 1) * dt
    return t, (out if return_all else M)


def flip_angle(b1, dt, gamma=GAMMA):
    """Flip angle [rad] of an on-resonance RF pulse with fixed phase axis: theta = gamma * |integral b1(t) dt|.

    The integral is signed: negative lobes (e.g. of a sinc) rotate back and reduce the flip angle.
    """
    return gamma * np.abs(np.sum(np.asarray(b1))) * dt


def b1_for_flip(theta, T_rf, gamma=GAMMA):
    """Amplitude [T] of a rectangular ("hard") pulse of length T_rf [s] giving flip angle theta [rad]."""
    return theta / (gamma * T_rf)


# ----------------------------------------------------------------------------
# Signal equations for common sequences (all times in seconds, flip in radians)
# ----------------------------------------------------------------------------
def spin_echo_signal(TE, TR, T1, T2, M0=1.0):
    """90-180 spin echo (assuming TE << TR):  S = M0 (1 - exp(-TR/T1)) exp(-TE/T2)."""
    return M0 * (1 - np.exp(-TR / T1)) * np.exp(-TE / T2)


def spoiled_gre_signal(flip, TE, TR, T1, T2star, M0=1.0):
    """Spoiled gradient echo (FLASH/SPGR) steady state:
        S = M0 sin(theta) (1 - E1) / (1 - cos(theta) E1) * exp(-TE/T2*),   E1 = exp(-TR/T1)."""
    E1 = np.exp(-TR / T1)
    return M0 * np.sin(flip) * (1 - E1) / (1 - np.cos(flip) * E1) * np.exp(-TE / T2star)


def ernst_angle(TR, T1):
    """Flip angle [rad] maximizing spoiled-GRE signal: theta_E = arccos(exp(-TR/T1))."""
    return np.arccos(np.exp(-TR / T1))


def inversion_recovery_signal(TI, TR, TE, T1, T2, M0=1.0):
    """Inversion recovery (180 - TI - 90 readout), assuming TE << TR:
        S = M0 (1 - 2 exp(-TI/T1) + exp(-TR/T1)) exp(-TE/T2).  Can be negative."""
    return M0 * (1 - 2 * np.exp(-TI / T1) + np.exp(-TR / T1)) * np.exp(-TE / T2)


def bssfp_signal(flip, TE, TR, T1, T2, M0=1.0):
    """Balanced SSFP on resonance (Haacke et al. 1999), typically TE = TR/2:
        S = M0 sin(theta)(1 - E1) / (1 - (E1 - E2) cos(theta) - E1 E2) * exp(-TE/T2)."""
    E1, E2 = np.exp(-TR / T1), np.exp(-TR / T2)
    return M0 * np.sin(flip) * (1 - E1) / (1 - (E1 - E2) * np.cos(flip) - E1 * E2) * np.exp(-TE / T2)


# ----------------------------------------------------------------------------
# Fourier helpers
# ----------------------------------------------------------------------------
def fftc(x, axis=-1, norm="ortho"):
    """Centered 1D FFT: zero position and zero frequency both at index N//2."""
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x, axes=axis), axis=axis, norm=norm), axes=axis)


def ifftc(X, axis=-1, norm="ortho"):
    """Centered 1D inverse FFT (inverse of `fftc`)."""
    return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(X, axes=axis), axis=axis, norm=norm), axes=axis)


def fft2c(x, norm="ortho"):
    """Centered 2D FFT over the last two axes (image -> k-space)."""
    ax = (-2, -1)
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x, axes=ax), axes=ax, norm=norm), axes=ax)


def ifft2c(X, norm="ortho"):
    """Centered 2D inverse FFT over the last two axes (k-space -> image)."""
    ax = (-2, -1)
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(X, axes=ax), axes=ax, norm=norm), axes=ax)


def sinc(x):
    """Normalized sinc, sin(pi x)/(pi x) (same as numpy.sinc). FT of a unit-width rect."""
    return np.sinc(x)


def rect(x):
    """Unit rectangle: 1 for |x| < 1/2, 1/2 at |x| = 1/2, 0 otherwise."""
    x = np.abs(np.asarray(x, dtype=float))
    return np.where(x < 0.5, 1.0, np.where(x == 0.5, 0.5, 0.0))


def jinc(x):
    """jinc(x) = J1(pi x) / (pi x)  (as in Larson's code; equals 1/2 at x = 0).

    The 2D Fourier transform of a disk is a jinc, the circular analogue of sinc.
    """
    from scipy.special import j1
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-12
    out[small] = 0.5
    px = np.pi * x[~small]
    out[~small] = j1(px) / px
    return out


# ----------------------------------------------------------------------------
# Phantoms
# ----------------------------------------------------------------------------
def disk_phantom_kspace(N=256, centers=((-0.3, -0.25), (0.0, 0.25), (0.3, -0.25)),
                        radius=0.125, values=(1.0, 1.0, 1.0)):
    """Analytic k-space of a few disks (port of Larson's contrast_phantom.m).

    centers and radius are in units of the field of view (FOV = 1, image spans
    -0.5..0.5). values = signal of each disk (e.g. from a signal equation).
    Returns the N x N k-space matrix M(kx, ky); image = ifft2c(M).
    Computing k-space analytically (instead of FFT of a pixel image) gives a
    realistic, non-periodic object: useful for truncation/ringing demos.
    """
    k = (np.arange(N) - N // 2) / N            # cycles per pixel
    kx, ky = np.meshgrid(k, k)
    M = np.zeros((N, N), complex)
    r_pix = radius * N
    for (xc, yc), v in zip(centers, values):
        # FT of a disk of radius R:  pi R^2 * 2 J1(2 pi R k)/(2 pi R k) = pi R^2 * 2 jinc(2 R k)
        kr = np.sqrt(kx**2 + ky**2)
        disk_ft = np.pi * r_pix**2 * 2 * jinc(2 * r_pix * kr)
        M += v * disk_ft * np.exp(-1j * 2 * np.pi * (kx * xc * N + ky * yc * N))
    return M / N   # scale so that ifft2c(..., norm="ortho") gives ~unit amplitude


def shepp_logan(N=256):
    """Modified Shepp-Logan head phantom (Toft), N x N, values in [0, 1]."""
    ellipses = [  # A, a, b, x0, y0, phi(deg)
        (1.0, .69, .92, 0, 0, 0),
        (-.8, .6624, .8740, 0, -.0184, 0),
        (-.2, .1100, .3100, .22, 0, -18),
        (-.2, .1600, .4100, -.22, 0, 18),
        (.1, .2100, .2500, 0, .35, 0),
        (.1, .0460, .0460, 0, .1, 0),
        (.1, .0460, .0460, 0, -.1, 0),
        (.1, .0460, .0230, -.08, -.605, 0),
        (.1, .0230, .0230, 0, -.606, 0),
        (.1, .0230, .0460, .06, -.605, 0),
    ]
    y, x = np.mgrid[1:-1:N * 1j, -1:1:N * 1j]
    img = np.zeros((N, N))
    for A, a, b, x0, y0, phi in ellipses:
        p = np.deg2rad(phi)
        xr = (x - x0) * np.cos(p) + (y - y0) * np.sin(p)
        yr = -(x - x0) * np.sin(p) + (y - y0) * np.cos(p)
        img[(xr / a) ** 2 + (yr / b) ** 2 <= 1] += A
    return img


def tissue_phantom(N=256):
    """Shepp-Logan-shaped phantom with tissue-like (M0, T1, T2) maps at ~1.5 T.

    Returns dict with 'M0', 'T1' [s], 'T2' [s] (each N x N) and 'labels' (int map:
    0 background, 1 skull/fat-like, 2 white matter, 3 gray matter, 4 CSF, 5 lesion).
    Values are textbook-typical, meant for contrast demos, not quantitative work.
    """
    sl = np.round(shepp_logan(N), 3)
    y, x = np.mgrid[1:-1:N * 1j, -1:1:N * 1j]
    head = (x / .69) ** 2 + (y / .92) ** 2 <= 1
    labels = np.zeros((N, N), int)
    labels[head & (sl == 1.0)] = 1                  # outer ring (fat/scalp-like)
    labels[head & (sl == 0.2)] = 2                  # white-matter-like
    labels[head & (sl == 0.3)] = 3                  # gray-matter-like ellipses
    labels[head & (sl == 0.0)] = 4                  # dark ellipses = CSF-filled ventricles
    labels[head & ((sl == 0.1) | (sl == 0.4))] = 5  # small spots = "lesions"
    #           bg     fat   WM    GM    CSF   lesion
    M0_tab = [0.0, 0.9, 0.7, 0.8, 1.0, 0.85]
    T1_tab = [1.0, 0.26, 0.60, 0.95, 4.0, 1.3]
    T2_tab = [1e-3, 0.08, 0.08, 0.10, 2.0, 0.15]
    M0 = np.take(M0_tab, labels)
    T1 = np.take(T1_tab, labels)
    T2 = np.take(T2_tab, labels)
    return {"M0": M0, "T1": T1, "T2": T2, "labels": labels}


# ----------------------------------------------------------------------------
# Magnetic fields of coils
# ----------------------------------------------------------------------------
def loop_field_on_axis(z, radius, current):
    """|B| [T] on the axis of a circular loop: mu0 I R^2 / (2 (R^2 + z^2)^(3/2))."""
    z = np.asarray(z, dtype=float)
    return MU0 * current * radius**2 / (2 * (radius**2 + z**2) ** 1.5)


def loop_coil_field(points, radius, current=1.0, center=(0, 0, 0), normal="z", n_segments=200):
    """Biot-Savart field [T] of a circular current loop at arbitrary points.

    Parameters
    ----------
    points : array (..., 3)   evaluation positions [m]
    radius : float            loop radius [m]
    current : float           current [A]
    center : (3,)             loop center [m]
    normal : 'x' | 'y' | 'z'  loop axis direction
    n_segments : int          number of straight segments approximating the circle

    Returns B with the same shape as points.  dB = mu0 I / (4 pi) * dl x R / |R|^3.
    """
    pts = np.asarray(points, dtype=float)
    phi = np.linspace(0, 2 * np.pi, n_segments + 1)
    ring = np.stack([radius * np.cos(phi), radius * np.sin(phi), np.zeros_like(phi)], axis=-1)
    perm = {"z": [0, 1, 2], "x": [2, 0, 1], "y": [1, 2, 0]}[normal]
    ring = ring[:, perm] + np.asarray(center, dtype=float)
    mid = 0.5 * (ring[1:] + ring[:-1])            # segment midpoints (S, 3)
    dl = ring[1:] - ring[:-1]                      # segment vectors   (S, 3)
    flat = pts.reshape(-1, 3)
    B = np.zeros_like(flat)
    for i in range(0, flat.shape[0], 4096):       # chunk to bound memory
        Rv = flat[i:i + 4096, None, :] - mid[None]                  # (P, S, 3)
        Rn = np.linalg.norm(Rv, axis=-1, keepdims=True)
        Rn = np.maximum(Rn, 1e-3 * radius)                          # avoid singularity on the wire
        B[i:i + 4096] = (MU0 * current / (4 * np.pi) * np.cross(dl[None], Rv) / Rn**3).sum(axis=1)
    return B.reshape(pts.shape)


# ----------------------------------------------------------------------------
# Noise / SNR helpers
# ----------------------------------------------------------------------------
def add_complex_noise(x, sigma, rng=None):
    """Add independent Gaussian noise of std `sigma` to real and imaginary parts."""
    rng = np.random.default_rng(0) if rng is None else rng
    x = np.asarray(x)
    return x + sigma * (rng.standard_normal(x.shape) + 1j * rng.standard_normal(x.shape))


# ----------------------------------------------------------------------------
# Plotting helpers (matplotlib is imported lazily so the module works headless)
# ----------------------------------------------------------------------------
def setup_plots():
    """Consistent figure style for all notebooks."""
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.figsize": (8, 4),
        "figure.dpi": 100,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 11,
        "legend.frameon": False,
    })


def plot_magnetization(t, M, ax=None, title=None, time_unit="ms", components="xyz"):
    """Plot Mx, My, Mz (and |Mxy|) versus time. M has shape (Nt, 3) or (Nt, 1, 3)."""
    import matplotlib.pyplot as plt
    M = np.asarray(M).reshape(len(t), -1, 3)[:, 0, :]
    scale = {"s": 1, "ms": 1e3, "us": 1e6}[time_unit]
    ax = plt.gca() if ax is None else ax
    names = {"x": r"$M_x$", "y": r"$M_y$", "z": r"$M_z$"}
    for i, c in enumerate("xyz"):
        if c in components:
            ax.plot(t * scale, M[:, i], label=names[c])
    if "a" in components:
        ax.plot(t * scale, np.abs(M[:, 0] + 1j * M[:, 1]), "k--", label=r"$|M_{xy}|$")
    ax.set_xlabel(f"time [{time_unit}]")
    ax.set_ylabel(r"magnetization [$M_0$]")
    ax.legend(loc="best")
    if title:
        ax.set_title(title)
    return ax


def plot_magnetization_3d(M, ax=None, title=None, color="C0"):
    """3D trajectory of the tip of M (shape (Nt, 3)), with the final vector drawn as an arrow."""
    import matplotlib.pyplot as plt
    M = np.asarray(M).reshape(-1, 3)
    if ax is None:
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(projection="3d")
    ax.plot(M[:, 0], M[:, 1], M[:, 2], color=color, lw=1)
    ax.quiver(0, 0, 0, *M[-1], color=color, lw=2, arrow_length_ratio=0.1)
    for v, lab in [((1, 0, 0), "x"), ((0, 1, 0), "y"), ((0, 0, 1), "z")]:
        ax.quiver(0, 0, 0, *v, color="gray", lw=0.8, arrow_length_ratio=0.05)
        ax.text(*(1.1 * np.array(v)), lab, color="gray")
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_xlabel(r"$M_x$"); ax.set_ylabel(r"$M_y$"); ax.set_zlabel(r"$M_z$")
    if title:
        ax.set_title(title)
    return ax


def show_image(img, ax=None, title=None, cmap="gray", vmin=None, vmax=None, colorbar=False, extent=None):
    """imshow with MRI-friendly defaults (magnitude of complex input, no ticks)."""
    import matplotlib.pyplot as plt
    img = np.asarray(img)
    if np.iscomplexobj(img):
        img = np.abs(img)
    ax = plt.gca() if ax is None else ax
    im = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent, origin="upper")
    if extent is None:
        ax.set_xticks([]); ax.set_yticks([])
    ax.grid(False)
    if title:
        ax.set_title(title)
    if colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return im


def show_kspace(K, ax=None, title=None, log=True):
    """Display log-magnitude of k-space so that the outer (weak) samples remain visible."""
    K = np.abs(np.asarray(K))
    return show_image(np.log10(K + 1e-6 * K.max()) if log else K, ax=ax, title=title, cmap="gray")
