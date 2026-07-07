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

Lbox = 15.6822 * Rp
Lbox_Rp = Lbox / Rp

print(f'Rstar/Rp = {Rstar_Rp:.3f}')
print(f'Lbox/Rp  = {Lbox_Rp:.3f}')

# ==========================================================
# DIMENSIONES
# ==========================================================

nx = 256
ny = 256
nv = 250

# ==========================================================
# LEER LAS TRES LINEAS
# ==========================================================

print('Leyendo mapas de tau...')

f0 = FortranFile('helines/HE_tau_3pj0_040.bin', 'r')
tau0 = f0.read_reals(np.float64)
tau0 = tau0.reshape((nx, ny, nv), order='F')

f1 = FortranFile('helines/HE_tau_3pj1_040.bin', 'r')
tau1 = f1.read_reals(np.float64)
tau1 = tau1.reshape((nx, ny, nv), order='F')

f2 = FortranFile('helines/HE_tau_3pj2_040.bin', 'r')
tau2 = f2.read_reals(np.float64)
tau2 = tau2.reshape((nx, ny, nv), order='F')

print('Mapas leidos correctamente.')

# ==========================================================
# GRAFICAR MAPA CENTRAL
# ==========================================================

iv = nv // 2

fig, ax = plt.subplots(1, 3, figsize=(15, 5))

ax[0].imshow(
    tau0[:, :, iv].T,
    origin='lower',
    cmap='inferno',
    norm=LogNorm()
)
ax[0].set_title('3P0')

ax[1].imshow(
    tau1[:, :, iv].T,
    origin='lower',
    cmap='inferno',
    norm=LogNorm()
)
ax[1].set_title('3P1')

im = ax[2].imshow(
    tau2[:, :, iv].T,
    origin='lower',
    cmap='inferno',
    norm=LogNorm()
)
ax[2].set_title('3P2')


cbar = fig.colorbar(
    im,
    ax=ax,
    orientation='horizontal',
    pad=0.1,      # separación respecto a los paneles
    fraction=0.09 # grosor de la barra
)

#plt.tight_layout()
plt.savefig('tau_maps_040.png', dpi=300,bbox_inches='tight')
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
# DISCO DEL PLANETA
# ==========================================================

planet_mask = (
    Xtau**2 + Ytau**2
    <= 1.0
)

# ==========================================================
# AREA ESTELAR
# ==========================================================

Astar = np.pi * Rstar_Rp**2

# ==========================================================
# RECORRIDO DEL PLANETA
# ==========================================================

nphase = 600

Ltransit = 3.5 * Rstar_Rp

xplanet = np.linspace(
    -Ltransit/2,
    +Ltransit/2,
    nphase
)

# ==========================================================
# ARRAYS DE SALIDA
# ==========================================================

lightcurve0 = np.zeros((nphase, nv))
lightcurve1 = np.zeros((nphase, nv))
lightcurve2 = np.zeros((nphase, nv))

lc_planet = np.zeros(nphase)

# ==========================================================
# LOOP DE TRANSITO
# ==========================================================

for i in range(nphase):

    if i % 50 == 0:
        print(f'{i}/{nphase}')

    xp = xplanet[i]

    Xs = Xtau + xp

    star_mask = (
        Xs**2 + Ytau**2
        <= Rstar_Rp**2
    )

    Acovered = star_mask.sum() * dA

    # ------------------------------------------------------
    # transito geometrico
    # ------------------------------------------------------

    occulted = (
        star_mask &
        planet_mask
    )

    Aplanet = occulted.sum() * dA

    lc_planet[i] = (
        Astar - Aplanet
    ) / Astar

    # ------------------------------------------------------
    # lineas espectrales
    # ------------------------------------------------------

    for iv in range(nv):

        trans0 = np.exp(-tau0[:, :, iv])
        trans1 = np.exp(-tau1[:, :, iv])
        trans2 = np.exp(-tau2[:, :, iv])

        trans0[planet_mask] = 0.0
        trans1[planet_mask] = 0.0
        trans2[planet_mask] = 0.0

        flux0 = (
            trans0[star_mask].sum() * dA
            + (Astar - Acovered)
        ) / Astar

        flux1 = (
            trans1[star_mask].sum() * dA
            + (Astar - Acovered)
        ) / Astar

        flux2 = (
            trans2[star_mask].sum() * dA
            + (Astar - Acovered)
        ) / Astar

        lightcurve0[i, iv] = flux0
        lightcurve1[i, iv] = flux1
        lightcurve2[i, iv] = flux2

# ==========================================================
# INTEGRAR EN VELOCIDAD
# ==========================================================

lc0 = lightcurve0.mean(axis=1)
lc1 = lightcurve1.mean(axis=1)
lc2 = lightcurve2.mean(axis=1)

# ==========================================================
# SUMAR ABSORCIONES
# ==========================================================

abs0 = 1.0 - lc0
abs1 = 1.0 - lc1
abs2 = 1.0 - lc2

lc_triplet = 1.0 - (
      abs0
    + abs1
    + abs2
)

# ==========================================================
# GUARDAR DATOS
# ==========================================================

np.savetxt(
    'lightcurve_triplet.dat',
    np.c_[
        xplanet,
        lc_planet,
        lc0,
        lc1,
        lc2,
        lc_triplet
    ],
    header=
    'xplanet(Rp) '
    'planet '
    '3pj0 '
    '3pj1 '
    '3pj2 '
    'triplet'
)

# ==========================================================
# GRAFICA
# ==========================================================

plt.figure(figsize=(9, 6))

plt.plot(
    xplanet / Rstar_Rp,
    lc_planet,
    'k--',
    lw=3,
    label='Planet'
)

plt.plot(
    xplanet / Rstar_Rp,
    lc0,
    lw=2,
    label='3P0'
)

plt.plot(
    xplanet / Rstar_Rp,
    lc1,
    lw=2,
    label='3P1'
)

plt.plot(
    xplanet / Rstar_Rp,
    lc2,
    lw=2,
    label='3P2'
)

plt.plot(
    xplanet / Rstar_Rp,
    lc_triplet,
    'k',
    lw=3,
    label='Triplet'
)

plt.xlabel(r'$x_p/R_\star$')
plt.ylabel('Normalized flux')
plt.legend()
plt.grid()

plt.tight_layout()

plt.savefig(
    'transit_triplet_040.png',
    dpi=300
)

plt.show()
