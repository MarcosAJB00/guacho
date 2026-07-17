import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LogNorm
from scipy.io import FortranFile

# ==========================================================
# PARAMETROS FISICOS
# ==========================================================

Rjup = 7.1492e9
Rsun = 6.957e10

Rp = 0.4466 * Rjup
Rstar = 0.76 * Rsun

Rstar_Rp = Rstar / Rp

Lbox = 2.0*4.77 * Rp
Lbox_Rp = Lbox / Rp

nx = 256
ny = 256
nv = 250

# ==========================================================
# LEER CURVA DE LUZ
# ==========================================================

data = np.loadtxt('lightcurve_triplet_sinchoque2.dat')

xplanet = data[:, 0]
flux = data[:, -1]

nframes = len(xplanet)

# ==========================================================
# LEER TAUS
# ==========================================================

print('Leyendo mapas de tau...')

f0 = FortranFile(
    'helines/HE_tau_3pj0_002.bin',
    'r'
)
tau0 = f0.read_reals(np.float64)
tau0 = tau0.reshape((nx, ny, nv), order='F')

f1 = FortranFile(
    'helines/HE_tau_3pj1_002.bin',
    'r'
)
tau1 = f1.read_reals(np.float64)
tau1 = tau1.reshape((nx, ny, nv), order='F')

f2 = FortranFile(
    'helines/HE_tau_3pj2_002.bin',
    'r'
)
tau2 = f2.read_reals(np.float64)
tau2 = tau2.reshape((nx, ny, nv), order='F')

print('Tau leidas.')

# ==========================================================
# TAU INTEGRADA EN VELOCIDAD
# ==========================================================

tau_int = (
      tau0.sum(axis=2)
    + tau1.sum(axis=2)
    + tau2.sum(axis=2)
)

tau_masked = np.ma.masked_where(
    tau_int <= 0,
    tau_int
)

cmap = plt.cm.magma.copy()
cmap.set_bad(alpha=0)

vmin = np.min(tau_int[tau_int > 0])
vmax = tau_int.max()

# ==========================================================
# MALLA
# ==========================================================

xtau = np.linspace(
    -Lbox_Rp/2,
     Lbox_Rp/2,
     nx
)

ytau = np.linspace(
    -Lbox_Rp/2,
     Lbox_Rp/2,
     ny
)

extent0 = [
    xtau.min(),
    xtau.max(),
    ytau.min(),
    ytau.max()
]

# ==========================================================
# FIGURA
# ==========================================================

fig = plt.figure(
    figsize=(10, 10),
    facecolor='black'
)

gs = fig.add_gridspec(
    2,
    1,
    height_ratios=[3, 1]
)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

ax1.set_facecolor('black')
ax2.set_facecolor('black')

# ==========================================================
# ESTRELLA
# ==========================================================

star = patches.Circle(
    (0, 0),
    Rstar_Rp,
    facecolor='orange',
    edgecolor='darkorange',
    alpha=0.5,
    zorder=1
)

ax1.add_patch(star)

# ==========================================================
# MAPA DE TAU
# ==========================================================

im = ax1.imshow(
    tau_masked.T,
    origin='lower',
    extent=extent0,
    cmap=cmap,
    norm=LogNorm(
        vmin=vmin,
        vmax=vmax
    ),
    zorder=3
)

# ==========================================================
# PLANETA
# ==========================================================

planet = patches.Circle(
    (0, 0),
    1.0,
    facecolor='deepskyblue',
    edgecolor='white',
    lw=1.5,
    zorder=4
)

ax1.add_patch(planet)

# ==========================================================
# PANEL SUPERIOR
# ==========================================================

ax1.set_xlim(
    -1.2*Rstar_Rp,
     1.2*Rstar_Rp
)

ax1.set_ylim(
    -1.2*Rstar_Rp,
     1.2*Rstar_Rp
)

ax1.set_aspect('equal')

ax1.set_xlabel(
    r'x [$R_p$]',
    color='white'
)

ax1.set_ylabel(
    r'y [$R_p$]',
    color='white'
)

ax1.tick_params(colors='white')

title = ax1.set_title(
    '',
    color='white',
    fontsize=15
)

cbar = fig.colorbar(
    im,
    ax=ax1,
    pad=0.02
)

cbar.set_label(
    r'$\tau$ integrada',
    color='white'
)

cbar.ax.tick_params(
    colors='white'
)

# ==========================================================
# PANEL INFERIOR
# ==========================================================

ax2.plot(
    xplanet/Rstar_Rp,
    flux,
    color='0.3',
    lw=3
)

trail, = ax2.plot(
    [],
    [],
    color='tomato',
    lw=3
)

point, = ax2.plot(
    [],
    [],
    'o',
    color='tomato',
    ms=8
)

vline = ax2.axvline(
    0,
    color='white',
    ls='--',
    alpha=0.5
)

ax2.set_xlim(
    (xplanet/Rstar_Rp).min(),
    (xplanet/Rstar_Rp).max()
)

ax2.set_ylim(
    flux.min()-0.01,
    1.01
)

ax2.set_xlabel(
    r'$x_p/R_\star$',
    color='white'
)

ax2.set_ylabel(
    'Flujo normalizado',
    color='white'
)

ax2.tick_params(colors='white')
ax2.grid(alpha=0.2)

# ==========================================================
# ANIMACION
# ==========================================================

def update(i):

    xp = xplanet[i]

    im.set_extent([
        xtau.min() + xp,
        xtau.max() + xp,
        ytau.min(),
        ytau.max()
    ])

    planet.center = (xp, 0)

    trail.set_data(
        xplanet[:i+1]/Rstar_Rp,
        flux[:i+1]
    )

    point.set_data(
        [xplanet[i]/Rstar_Rp],
        [flux[i]]
    )

    xv = xplanet[i]/Rstar_Rp

    vline.set_xdata(
        [xv, xv]
    )

    title.set_text(
        f'x = {xp:.2f} Rp'
    )

    return (
        im,
        planet,
        trail,
        point,
        vline,
        title
    )

# ==========================================================
# CREAR GIF
# ==========================================================

print('Generando GIF...')

ani = FuncAnimation(
    fig,
    update,
    frames=nframes,
    interval=10,
    blit=False
)

ani.save(
    'transit_live.gif',
    writer=PillowWriter(fps=25)
)

print('GIF guardado:')
print('transit_live.gif')