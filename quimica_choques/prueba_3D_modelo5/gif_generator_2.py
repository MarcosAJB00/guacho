import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as pe
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from scipy.io import FortranFile

# ==========================================================
# PARAMETROS FISICOS  (sin cambios)
# ==========================================================

Rjup = 7.1492e9
Rsun = 6.957e10

Rp    = 0.4466 * Rjup
Rstar = 0.76   * Rsun

Rstar_Rp = Rstar / Rp

Lbox    = 14.73812 * Rp
Lbox_Rp = Lbox / Rp

nx = 256
ny = 256
nv = 250

# ==========================================================
# LEER CURVA DE LUZ  (sin cambios)
# ==========================================================

data = np.loadtxt('lightcurve_triplet.dat')

xplanet = data[:, 0]
flux    = data[:, -1]

nframes = len(xplanet)

# ==========================================================
# LEER TAUS  (sin cambios)
# ==========================================================

print('Leyendo mapas de tau...')

f0  = FortranFile('helines/HE_tau_3pj0_040.bin', 'r')
tau0 = f0.read_reals(np.float64).reshape((nx, ny, nv), order='F')

f1  = FortranFile('helines/HE_tau_3pj1_040.bin', 'r')
tau1 = f1.read_reals(np.float64).reshape((nx, ny, nv), order='F')

f2  = FortranFile('helines/HE_tau_3pj2_040.bin', 'r')
tau2 = f2.read_reals(np.float64).reshape((nx, ny, nv), order='F')

print('Tau leidas.')

# ==========================================================
# TAU INTEGRADA EN VELOCIDAD  (sin cambios)
# ==========================================================

tau_int = (
      tau0.sum(axis=2)
    + tau1.sum(axis=2)
    + tau2.sum(axis=2)
)

tau_masked = np.ma.masked_where(tau_int <= 0, tau_int)

vmin = np.min(tau_int[tau_int > 0])
vmax = tau_int.max()

# ==========================================================
# MALLA  (sin cambios)
# ==========================================================

xtau = np.linspace(-Lbox_Rp / 2,  Lbox_Rp / 2, nx)
ytau = np.linspace(-Lbox_Rp / 2,  Lbox_Rp / 2, ny)

extent0 = [xtau.min(), xtau.max(), ytau.min(), ytau.max()]

# ==========================================================
# PALETA Y ESTILO VISUAL
# ==========================================================

# Fondo negro profundo tipo espacio
BG      = '#04060e'
STAR_C  = '#ffb347'   # naranja para la estrella
STAR_E  = '#e07000'
CURVE_C = '#ff4d6d'
TRAIL_C = '#ff4d6d'
GRID_C  = '#1a2340'
TEXT_C  = '#ccd6f6'
TICK_C  = '#7a8ab0'
TITLE_C = '#e8f0ff'
LABEL_C = '#8892b0'

# Colormap inferno para el planeta y su material
cmap_he = plt.cm.inferno.copy()
cmap_he.set_bad(alpha=0)

# ==========================================================
# FIGURA
# ==========================================================

plt.rcParams.update({
    'font.family'      : 'DejaVu Sans',
    'font.size'        : 13,
    'axes.labelsize'   : 14,
    'axes.titlesize'   : 15,
    'xtick.labelsize'  : 11,
    'ytick.labelsize'  : 11,
    'axes.linewidth'   : 0.6,
    'xtick.color'      : TICK_C,
    'ytick.color'      : TICK_C,
})

fig = plt.figure(figsize=(9, 11), facecolor=BG, dpi=110)

gs = GridSpec(
    2, 1,
    height_ratios=[3, 1],
    hspace=0.06,
    left=0.10, right=0.92,
    top=0.93,  bottom=0.08,
)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

for ax in (ax1, ax2):
    ax.set_facecolor(BG)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID_C)
        spine.set_linewidth(0.8)

# ==========================================================
# PANEL SUPERIOR — estrella + tau + planeta
# ==========================================================

# Estrella: un solo círculo limpio
star = patches.Circle(
    (0, 0), Rstar_Rp,
    facecolor=STAR_C, edgecolor=STAR_E,
    alpha=0.85, linewidth=1.2, zorder=1,
)
ax1.add_patch(star)

# Mapa de tau
im = ax1.imshow(
    tau_masked.T,
    origin='lower',
    extent=extent0,
    cmap=cmap_he,
    norm=LogNorm(vmin=vmin, vmax=vmax),
    zorder=3,
    interpolation='bilinear',
)

# Planeta: negro/oscuro sobre el material inferno
planet = patches.Circle(
    (0, 0), 1.0,
    facecolor='#0a0a0a', edgecolor='#555555',
    linewidth=1.2, zorder=4,
)
ax1.add_patch(planet)

# Límites y aspecto  (sin cambios)
ax1.set_xlim(-1.2*Rstar_Rp,  1.2*Rstar_Rp)
ax1.set_ylim(-1.2*Rstar_Rp,  1.2*Rstar_Rp)
ax1.set_aspect('equal')

ax1.set_xlabel(r'$x\ [R_\mathrm{p}]$',   color=TEXT_C, labelpad=6)
ax1.set_ylabel(r'$y\ [R_\mathrm{p}]$',   color=TEXT_C, labelpad=6)
ax1.tick_params(colors=TICK_C, direction='in', length=4, width=0.7)

# Colorbar
cbar = fig.colorbar(im, ax=ax1, pad=0.02, fraction=0.035)
cbar.set_label(
    r'$\tau_\mathrm{He\,10830}$  (integrado en $v$)',
    color=TEXT_C, fontsize=12,
)
cbar.ax.tick_params(colors=TICK_C, labelsize=10)
cbar.outline.set_edgecolor(GRID_C)

# Título con path effects para legibilidad
title = ax1.set_title(
    '',
    color=TITLE_C, fontsize=14, fontweight='bold', pad=10,
    path_effects=[pe.withStroke(linewidth=3, foreground=BG)],
)

# Anotaciones fijas
ax1.text(
    0.02, 0.97,
    r'HD 189733 b — He\,{\sc i} 10830\,Å',
    transform=ax1.transAxes,
    color=TEXT_C, fontsize=11, va='top',
    path_effects=[pe.withStroke(linewidth=2, foreground=BG)],
)
ax1.text(
    0.02, 0.90,
    r'Triplete: $^3P_0$, $^3P_1$, $^3P_2$',
    transform=ax1.transAxes,
    color=LABEL_C, fontsize=10, va='top',
    path_effects=[pe.withStroke(linewidth=2, foreground=BG)],
)

# ==========================================================
# PANEL INFERIOR — curva de luz  (sin cambios en los datos)
# ==========================================================

ax2.plot(
    xplanet / Rstar_Rp, flux,
    color=GRID_C, lw=2.5, zorder=1,
)

# Relleno bajo la curva completa (muy sutil)
ax2.fill_between(
    xplanet / Rstar_Rp, flux, flux.min() - 0.01,
    color=CURVE_C, alpha=0.04, zorder=0,
)

# Líneas de contacto (I, II, III, IV) — opcionales, si querés
# ax2.axvline(x_contacto, color=TICK_C, ls=':', lw=0.8)

trail, = ax2.plot([], [], color=TRAIL_C, lw=2.5, zorder=3)
point, = ax2.plot([], [], 'o', color=TRAIL_C, ms=7, zorder=4,
                  markeredgecolor='white', markeredgewidth=0.8)
vline  = ax2.axvline(0, color='white', ls='--', lw=0.8, alpha=0.40)

ax2.set_xlim((xplanet / Rstar_Rp).min(), (xplanet / Rstar_Rp).max())
ax2.set_ylim(flux.min() - 0.01, 1.01)

ax2.set_xlabel(r'$x_\mathrm{p}\ /\ R_\star$', color=TEXT_C, labelpad=6)
ax2.set_ylabel('Flujo norm.',                  color=TEXT_C, labelpad=6)

ax2.tick_params(colors=TICK_C, direction='in', length=4, width=0.7)

ax2.grid(color=GRID_C, linewidth=0.5, alpha=0.6)

ax2.text(
    0.98, 0.95,
    'Curva de luz He 10830',
    transform=ax2.transAxes,
    color=LABEL_C, fontsize=10, ha='right', va='top',
)

# ==========================================================
# ANIMACION  (sin cambios en la lógica)
# ==========================================================

def update(i):

    xp = xplanet[i]

    # --- mapa de tau (sin cambios) ---
    im.set_extent([
        xtau.min() + xp,
        xtau.max() + xp,
        ytau.min(),
        ytau.max(),
    ])

    # --- planeta ---
    planet.center = (xp, 0)

    # --- curva de luz (sin cambios) ---
    trail.set_data(xplanet[:i+1] / Rstar_Rp, flux[:i+1])
    point.set_data([xplanet[i]   / Rstar_Rp], [flux[i]])

    xv = xplanet[i] / Rstar_Rp
    vline.set_xdata([xv, xv])

    # --- título con absorción actual ---
    abs_pct = (1.0 - flux[i]) * 100.0
    title.set_text(
        rf'$x_\mathrm{{p}}$ = {xp:.1f} $R_\mathrm{{p}}$'
        rf'     absorción = {abs_pct:.2f}%'
    )

    return (im, planet, trail, point, vline, title)

# ==========================================================
# CREAR GIF  — fps=30 → más rápido que el original (25 ms → ~33 fps)
# ==========================================================

print('Generando GIF...')

ani = FuncAnimation(
    fig, update,
    frames=nframes,
    interval=25,     # sin cambios
    blit=False,
)

ani.save(
    'transit_live.gif',
    writer=PillowWriter(fps=60),
)

print('GIF guardado: transit_live.gif')