import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from matplotlib.colors import LogNorm

# ==========================================================
# PARAMETROS FISICOS
# ==========================================================

Rjup = 7.1492e9
Rsun = 6.957e10

Rp = 0.4466 * Rjup
Rstar = 0.76 * Rsun

Rstar_Rp = Rstar / Rp

Lbox = 2.0 * 5.09421794574017 * Rp
Lbox_Rp = Lbox / Rp

print(f"Rstar/Rp = {Rstar_Rp:.3f}")
print(f"Lbox/Rp  = {Lbox_Rp:.3f}")

# ==========================================================
# DIMENSIONES
# ==========================================================

nx = 256
ny = 256
nv = 250

# ==========================================================
# LEER MAPA DE TAU
# ==========================================================

fname = 'helines/HE_tau_3pj2_040.bin'

f = FortranFile(fname, 'r')

tau = f.read_reals(np.float64)
tau = tau.reshape((nx, ny, nv), order='F')

print('Tau leida correctamente.')

# ==========================================================
# GRAFICAR TAU CENTRAL
# ==========================================================

iv = nv // 2

plt.figure(figsize=(6, 6))
plt.imshow(
    tau[:, :, iv].T,
    origin='lower',
    cmap='inferno',
    norm=LogNorm()
    
)
plt.colorbar(label=r'$\tau$')
plt.title(f'HE_tau_3pj2_040, canal {iv}')
plt.tight_layout()
plt.savefig('tau_map_040.png', dpi=300)
plt.close()

# ==========================================================
# MALLA DEL MAPA DE TAU
# ==========================================================

xtau = np.linspace(
    -Lbox_Rp/2,
    +Lbox_Rp/2,
    nx
)

ytau = np.linspace(
    -Lbox_Rp/2,
    +Lbox_Rp/2,
    ny
)

dx = xtau[1] - xtau[0]
dA = dx**2

Xtau, Ytau = np.meshgrid(
    xtau,
    ytau,
    indexing='ij'
)

# ==========================================================
# DISCO OPACO DEL PLANETA
# ==========================================================

planet_mask = (
    Xtau**2 + Ytau**2
    <= 1.0
)

# ==========================================================
# AREA DE LA ESTRELLA
# ==========================================================

Astar = np.pi * Rstar_Rp**2

# ==========================================================
# TRANSITO
# ==========================================================

nphase = 600

Ltransit = 5.0 * Rstar_Rp

xplanet = np.linspace(
    -Ltransit/2,
    +Ltransit/2,
    nphase
)

lightcurve = np.zeros((nphase, nv))

# ==========================================================
# LOOP DE TRANSITO
# ==========================================================

for i in range(nphase):

    if i % 20 == 0:
        print(f'phase {i}/{nphase}')

    xp = xplanet[i]

    Xs = Xtau + xp

    star_mask = (
        Xs**2 + Ytau**2
        <= Rstar_Rp**2
    )

    Acovered = star_mask.sum() * dA

    for iv in range(nv):

        transmission = np.exp(
            -tau[:, :, iv]
        )

        transmission[planet_mask] = 0.0

        absorbed_flux = (
            transmission[star_mask].sum()
            * dA
        )

        lightcurve[i, iv] = (
            absorbed_flux
            + (Astar - Acovered)
        ) / Astar

# ==========================================================
# CURVA INTEGRADA
# ==========================================================

lc = lightcurve.mean(axis=1)

# ==========================================================
# GUARDAR DATOS
# ==========================================================

np.savetxt(
    'lightcurve_040.dat',
    np.c_[xplanet, lc],
    header='xplanet(Rp) Flux'
)

# ==========================================================
# GRAFICA DEL TRANSITO
# ==========================================================

plt.figure(figsize=(8, 5))

plt.plot(
    xplanet,
    lc,
    lw=2
)

plt.xlabel(r'$x_p\ (R_p)$')
plt.ylabel('Normalized flux')
plt.title('Synthetic transit using HE_tau_3pj2_040')
plt.grid()

plt.tight_layout()
plt.savefig(
    'transit_040.png',
    dpi=300
)

plt.show()