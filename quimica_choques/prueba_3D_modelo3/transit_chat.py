import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

# ==========================================================
# PARÁMETROS DEL SISTEMA
# ==========================================================

Rjup = 7.1492e9
Rsun = 6.957e10
AU = 1.495978707e13

Rp = 0.4466 * Rjup
Rstar = 0.76 * Rsun
a = 0.0532 * AU

# parámetro de impacto
b = 0.0

# ==========================================================
# CAJA DE LA SIMULACIÓN
# ==========================================================

nx = 256
ny = 256
nv = 250
nsteps = 81

Lbox = 2.0 * 5.09421794574017 * Rp

dx = Lbox / nx
dy = Lbox / ny
dA = dx * dy

# ==========================================================
# COORDENADAS DE LA CAJA
# ==========================================================

x = (np.arange(nx) + 0.5 - nx/2) * dx
y = (np.arange(ny) + 0.5 - ny/2) * dy

Xbox, Ybox = np.meshgrid(x, y, indexing='ij')

# ==========================================================
# FASES
# ==========================================================

phase_deg = np.linspace(-5.0, 5.0, nsteps)
phase = np.deg2rad(phase_deg)

# ==========================================================
# ÁREA ESTELAR
# ==========================================================

Astar = np.pi * Rstar**2

# ==========================================================
# SALIDAS
# ==========================================================

spectra = np.ones((nsteps, nv))
planet_lc = np.ones(nsteps)

# ==========================================================
# BUCLE SOBRE LAS FASES
# ==========================================================

for i, th in enumerate(phase):

    print(f'{i+1}/{nsteps}')

    # ------------------------------------------------------
    # posición del centro del planeta
    # ------------------------------------------------------

    xp = a * np.sin(th)
    yp = b * Rstar

    # ======================================================
    # TRÁNSITO GEOMÉTRICO DEL PLANETA
    # ======================================================

    d = np.sqrt(xp**2 + yp**2)

    if d >= Rstar + Rp:

        planet_lc[i] = 1.0

    elif d <= Rstar - Rp:

        planet_lc[i] = 1.0 - (Rp/Rstar)**2

    else:

        alpha = np.arccos(
            (d**2 + Rp**2 - Rstar**2)
            /(2*d*Rp)
        )

        beta = np.arccos(
            (d**2 + Rstar**2 - Rp**2)
            /(2*d*Rstar)
        )

        overlap = (
            Rp**2 * alpha
            + Rstar**2 * beta
            - 0.5*np.sqrt(
                (-d+Rp+Rstar)
                *(d+Rp-Rstar)
                *(d-Rp+Rstar)
                *(d+Rp+Rstar)
            )
        )

        planet_lc[i] = 1.0 - overlap/Astar

    # ======================================================
    # SI LA CAJA NO TOCA LA ESTRELLA
    # ======================================================

    if d > Rstar + np.sqrt(2)*Lbox/2:
        continue

    # ======================================================
    # LEER MAPA TAU
    # ======================================================

    fname = f'helines/HE_tau_3pj2_{i:03d}.bin'

    f = FortranFile(fname, 'r')
    tau = f.read_reals(np.float64)
    f.close()

    tau = tau.reshape((nx, ny, nv), order='F')

    # ======================================================
    # POSICIÓN ABSOLUTA DE CADA PÍXEL
    # ======================================================

    X = Xbox + xp
    Y = Ybox + yp

    stellar = (X**2 + Y**2) <= Rstar**2

    if stellar.sum() == 0:
        continue

    # ======================================================
    # ESPECTRO TRANSMITIDO
    # ======================================================

    for k in range(nv):

        absorption = 1.0 - np.exp(-tau[:, :, k])

        absorbed_area = (
            absorption[stellar].sum()
            * dA
        )

        spectra[i, k] = (
            1.0
            - absorbed_area/Astar
        )

# ==========================================================
# CURVA DE LUZ DEL GAS
# ==========================================================

vel = np.linspace(-150.0, 150.0, nv)

mask = (vel >= -20.0) & (vel <= 20.0)

gas_lc = np.trapezoid(
    1.0 - spectra[:, mask],
    vel[mask],
    axis=1
)

gas_lc = 1.0 - gas_lc/(40.0)

# ==========================================================
# GRAFICAR
# ==========================================================

plt.figure(figsize=(8,6))

plt.plot(
    phase_deg,
    planet_lc,
    lw=3,
    label='Planetary disk'
)

plt.plot(
    phase_deg,
    gas_lc,
    lw=3,
    label='He 10830'
)

plt.xlabel('Orbital phase [deg]')
plt.ylabel('Normalized flux')

plt.legend()
plt.tight_layout()
plt.savefig('transit_He_chat.png')
plt.show()