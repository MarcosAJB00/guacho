"""
planetary_wind_transit.py
=========================
Visualizacion 3D de viento planetario con region de choque (bow shock)
y proyeccion en el plano del cielo para calcular transito.

ESTRUCTURA DEL CODIGO:
  1. PARAMETROS  <-- modifica esta seccion
  2. Carga de perfiles (archivo o analitico)
  3. Construccion de la grilla 3D
  4. Visualizacion 3D de la densidad
  5. Proyeccion en el plano del cielo (column density)
  6. Curva de luz del transito

FORMATOS DE ARCHIVO ACEPTADOS:
  - Dos columnas separadas por espacios/tabs/comas: r   densidad
  - Lineas que empiezan con # se ignoran (comentarios/header)
  - Ejemplo:
      # r[Rjup]   n[cm-3]
      0.12   1.0e12
      0.20   4.5e11
      ...
      1.00   2.1e9
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d, RegularGridInterpolator
import os

# ============================================================
#  1. PARAMETROS  (modifica esta seccion)
# ============================================================

# --- Geometria del sistema ---
R_PLANET   = 0.12    # Radio del planeta en unidades arbitrarias (ej. R_Jupiter)
R_HILL     = 1.0     # Radio de Hill (define la esfera, unidad de referencia)
R_STAR     = 1.0     # Radio de la estrella EN LAS MISMAS UNIDADES (para la curva de luz)

# --- Perfil de viento planetario ---
# Opcion A: archivo con columnas (r, densidad)
WIND_FILE  = None    # ej. "perfil_viento.dat"  -- pone None para usar perfil analitico

# Opcion B: perfil analitico  rho(r) = RHO0_WIND * (r / R_HILL)^ALPHA_WIND
RHO0_WIND  = 1.0     # densidad en R_HILL (unidades normalizadas)
ALPHA_WIND = -2.0    # indice de ley de potencia (negativo = decrece hacia afuera)

# --- Region de choque (bow shock) ---
# Opcion A: archivo con columnas (r, densidad)  r va desde R_HILL hasta R_EXT_SHOCK
SHOCK_FILE = None    # ej. "perfil_choque.dat"  -- pone None para usar perfil analitico

# Opcion B: perfil analitico  rho(r) = RHO0_SHOCK * (r / R_HILL)^BETA_SHOCK
RHO0_SHOCK = 5.0     # densidad en R_HILL (relativa al viento)
BETA_SHOCK = -3.0    # indice del choque (mas empinado que el viento)
R_EXT_SHOCK = 1.5    # radio externo del choque en unidades de R_HILL

# --- Geometria del parche de choque ---
DELTA_PHI   = 40.0   # semi-angulo del parche en grados (el parche subtiende 2*DELTA_PHI)
THETA_PATCH = 90.0   # colatitud del centro del parche en grados (90 = ecuador)
PHI0_PATCH  = 0.0    # azimut del centro del parche en grados

# --- Grilla de calculo ---
N_R     = 60         # puntos radiales
N_THETA = 80         # puntos en colatitud
N_PHI   = 160        # puntos en azimut

# --- Proyeccion en el plano del cielo ---
# La linea de vision esta a lo largo del eje Z del sistema 3D.
# El planeta se mueve en la direccion X del plano del cielo durante el transito.
N_IMPACT = 100       # resolucion de la grilla del plano del cielo (NxN)
N_LOS    = 300       # puntos de integracion a lo largo de la linea de vision

# --- Curva de luz ---
N_TRANSIT_STEPS = 200   # pasos temporales del transito
B_IMPACT = 0.0          # parametro de impacto del transito en Y_sky (unidades de R_Hill)
                         # 0 = transito central, +/- = transito no central

# --- Salida ---
SAVE_FIGURES = True      # guarda los plots como PNG
OUTPUT_DIR   = "."       # directorio de salida


# ============================================================
#  2. CARGA DE PERFILES
# ============================================================

def load_profile(filename):
    """Carga un archivo de dos columnas (r, densidad). Ignora lineas con #."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                parts = line.replace(',', ' ').split()
                data.append((float(parts[0]), float(parts[1])))
            except (ValueError, IndexError):
                continue
    data = np.array(data)
    return data[:, 0], data[:, 1]


def make_wind_profile():
    if WIND_FILE and os.path.exists(WIND_FILE):
        r_data, rho_data = load_profile(WIND_FILE)
        print(f"  Viento: cargado desde '{WIND_FILE}' ({len(r_data)} puntos)")
        log_interp = interp1d(np.log(r_data), np.log(rho_data),
                              kind='linear', fill_value='extrapolate')
        return lambda r: np.exp(log_interp(np.log(np.clip(r, r_data.min(), r_data.max()))))
    else:
        if WIND_FILE:
            print(f"  AVISO: '{WIND_FILE}' no encontrado, usando perfil analitico.")
        else:
            print(f"  Viento: perfil analitico rho ~ r^{ALPHA_WIND:.1f}")
        return lambda r: RHO0_WIND * (np.asarray(r) / R_HILL) ** ALPHA_WIND


def make_shock_profile():
    if SHOCK_FILE and os.path.exists(SHOCK_FILE):
        r_data, rho_data = load_profile(SHOCK_FILE)
        print(f"  Choque: cargado desde '{SHOCK_FILE}' ({len(r_data)} puntos)")
        log_interp = interp1d(np.log(r_data), np.log(rho_data),
                              kind='linear', fill_value='extrapolate')
        return lambda r: np.exp(log_interp(np.log(np.clip(r, r_data.min(), r_data.max()))))
    else:
        if SHOCK_FILE:
            print(f"  AVISO: '{SHOCK_FILE}' no encontrado, usando perfil analitico.")
        else:
            print(f"  Choque: perfil analitico rho ~ r^{BETA_SHOCK:.1f}")
        return lambda r: RHO0_SHOCK * (np.asarray(r) / R_HILL) ** BETA_SHOCK


# ============================================================
#  3. FUNCION DE DENSIDAD VECTORIZADA
# ============================================================

def is_in_patch(vx, vy, vz):
    """Devuelve mascara booleana: True donde la direccion cae en el parche del choque."""
    dphi_rad  = np.radians(DELTA_PHI)
    theta_rad = np.radians(THETA_PATCH)
    phi0_rad  = np.radians(PHI0_PATCH)
    cx = np.sin(theta_rad) * np.cos(phi0_rad)
    cy = np.cos(theta_rad)
    cz = np.sin(theta_rad) * np.sin(phi0_rad)
    dot = vx * cx + vy * cy + vz * cz
    return dot >= np.cos(dphi_rad)


def density_at(r, vx, vy, vz, wind_func, shock_func):
    """Densidad en un punto dado su radio y vector unitario de direccion."""
    r = np.asarray(r)
    vx, vy, vz = np.asarray(vx), np.asarray(vy), np.asarray(vz)
    rho = np.zeros_like(r, dtype=float)

    wind_mask  = r <= R_HILL
    shock_mask = (r > R_HILL) & (r <= R_EXT_SHOCK * R_HILL) & is_in_patch(vx, vy, vz)

    if np.any(wind_mask):
        rho[wind_mask] = wind_func(r[wind_mask])
    if np.any(shock_mask):
        rho[shock_mask] = shock_func(r[shock_mask])
    return rho


# ============================================================
#  4. VISUALIZACION 3D (cortes y perfiles)
# ============================================================

def plot_3d_density(wind_func, shock_func):
    print("\nGenerando visualizacion de cortes 3D...")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    r_arr    = np.linspace(R_PLANET * 0.5, R_EXT_SHOCK * R_HILL * 1.05, 300)
    th_arr   = np.linspace(0, np.pi, 300)
    R2d, TH2d = np.meshgrid(r_arr, th_arr, indexing='ij')

    def make_slice(phi_val):
        VX = np.sin(TH2d) * np.cos(phi_val)
        VY = np.cos(TH2d)
        VZ = np.sin(TH2d) * np.sin(phi_val)
        return density_at(R2d.ravel(), VX.ravel(), VY.ravel(), VZ.ravel(),
                          wind_func, shock_func).reshape(R2d.shape)

    for ax, phi_val, xlabel, ylabel, title in [
        (axes[0], 0.0,       'X [R$_{Hill}$]', 'Y [R$_{Hill}$]', 'Corte φ=0° (plano XY)'),
        (axes[1], np.pi/2,   'X [R$_{Hill}$]', 'Z [R$_{Hill}$]', 'Corte φ=90° (plano XZ)'),
    ]:
        rho2d = make_slice(phi_val)
        Xc = R2d * np.sin(TH2d)
        Yc = R2d * np.cos(TH2d)
        rho_masked = np.ma.masked_where(rho2d <= 0, rho2d)
        vm = rho_masked.max()
        im = ax.pcolormesh(Xc, Yc, rho_masked,
                           norm=LogNorm(vmin=vm * 1e-4, vmax=vm),
                           cmap='plasma', shading='auto')
        for r_c, col, ls in [(R_HILL, 'white', '--'), (R_EXT_SHOCK*R_HILL, 'yellow', ':')]:
            ax.add_patch(plt.Circle((0, 0), r_c, fill=False, color=col, lw=1.2, ls=ls, alpha=0.8))
        ax.add_patch(plt.Circle((0, 0), R_PLANET, color='cyan', zorder=5))
        lim = R_EXT_SHOCK * R_HILL * 1.1
        ax.set_xlim(0, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect('equal')
        ax.set_xlabel(xlabel, fontsize=10); ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=10)
        plt.colorbar(im, ax=ax, label='Densidad', fraction=0.046, pad=0.04)

    # Perfiles 1D
    ax = axes[2]
    r_p = np.linspace(R_PLANET, R_EXT_SHOCK * R_HILL, 500)
    ax.semilogy(r_p[r_p <= R_HILL], wind_func(r_p[r_p <= R_HILL]),
                'b-', lw=2.5, label='Viento planetario')
    ax.semilogy(r_p[r_p >= R_HILL], shock_func(r_p[r_p >= R_HILL]),
                'r--', lw=2.5, label='Choque (bow shock)')
    ax.axvline(R_HILL, color='gray', ls=':', lw=1.2, alpha=0.7, label=f'R$_{{Hill}}$={R_HILL}')
    ax.axvline(R_PLANET, color='cyan', ls=':', lw=1.2, alpha=0.8, label=f'R$_p$={R_PLANET}')
    ax.set_xlabel('r [R$_{Hill}$]', fontsize=10)
    ax.set_ylabel('Densidad', fontsize=10)
    ax.set_title('Perfiles 1D', fontsize=10)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    plt.suptitle(
        f'Densidad 3D  |  2Δφ={2*DELTA_PHI:.0f}°  θ_patch={THETA_PATCH:.0f}°  '
        f'φ₀={PHI0_PATCH:.0f}°  R_ext={R_EXT_SHOCK}R_H',
        fontsize=11)
    plt.tight_layout()
    if SAVE_FIGURES:
        fname = os.path.join(OUTPUT_DIR, 'fig1_densidad_3d.png')
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"  Guardado: {fname}")
    plt.show()


# ============================================================
#  5. PROYECCION EN EL PLANO DEL CIELO
# ============================================================

def compute_column_density(wind_func, shock_func):
    """
    Integra la densidad a lo largo de la linea de vision (eje Z).
    Devuelve grilla 2D N(x_sky, y_sky).
    """
    print("\nCalculando densidad columnar...")
    xy_range = R_EXT_SHOCK * R_HILL * 1.15
    x_sky = np.linspace(-xy_range, xy_range, N_IMPACT)
    y_sky = np.linspace(-xy_range, xy_range, N_IMPACT)
    XX, YY = np.meshgrid(x_sky, y_sky, indexing='ij')
    N_col  = np.zeros_like(XX)

    z_los = np.linspace(-xy_range, xy_range, N_LOS)
    dz = z_los[1] - z_los[0]

    for ix in range(N_IMPACT):
        x0 = x_sky[ix]
        for iy in range(N_IMPACT):
            y0 = y_sky[iy]
            r_arr = np.sqrt(x0**2 + y0**2 + z_los**2)
            valid = (r_arr >= R_PLANET) & (r_arr <= R_EXT_SHOCK * R_HILL * 1.02)
            if not np.any(valid):
                continue
            rv = r_arr[valid]
            zv = z_los[valid]
            vx = x0 / rv; vy = y0 / rv; vz_v = zv / rv
            rho_v = density_at(rv, vx, vy, vz_v, wind_func, shock_func)
            N_col[ix, iy] = np.sum(rho_v) * dz
        if (ix + 1) % 20 == 0:
            print(f"  {ix+1}/{N_IMPACT} columnas procesadas")

    print(f"  N_col rango: [{N_col[N_col>0].min():.3e}, {N_col.max():.3e}]")
    return x_sky, y_sky, N_col


def plot_column_density(x_sky, y_sky, N_col):
    print("\nGenerando mapa de densidad columnar...")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    nc_pos = np.where(N_col > 0, N_col, np.nan)
    vm = np.nanmax(nc_pos)
    im = ax.pcolormesh(x_sky, y_sky, nc_pos.T,
                       norm=LogNorm(vmin=vm * 1e-4, vmax=vm),
                       cmap='inferno', shading='auto')
    ax.add_patch(plt.Circle((0, 0), R_STAR,
                 fill=False, color='yellow', lw=1.5, label=f'Estrella R={R_STAR}'))
    ax.add_patch(plt.Circle((0, 0), R_HILL,
                 fill=False, color='white', lw=1, ls='--', label=f'R$_{{Hill}}$={R_HILL}'))
    ax.set_aspect('equal')
    ax.set_xlabel('X$_{sky}$ [R$_{Hill}$]', fontsize=11)
    ax.set_ylabel('Y$_{sky}$ [R$_{Hill}$]', fontsize=11)
    ax.set_title('Densidad columnar N(x,y)  [log]', fontsize=11)
    plt.colorbar(im, ax=ax, label='N  [densidad × longitud]', fraction=0.046, pad=0.04)
    ax.legend(fontsize=8, loc='upper right')

    ax2 = axes[1]
    iy_mid = np.argmin(np.abs(y_sky - B_IMPACT))
    ax2.semilogy(x_sky, N_col[:, iy_mid], 'b-', lw=2,
                 label=f'N(x, y={y_sky[iy_mid]:.2f})')
    ax2.axvline(-R_HILL, color='gray', ls='--', alpha=0.6)
    ax2.axvline(R_HILL,  color='gray', ls='--', alpha=0.6, label='±R$_{Hill}$')
    ax2.axvspan(-R_STAR, R_STAR, alpha=0.1, color='yellow', label='Disco estelar')
    ax2.set_xlabel('X$_{sky}$ [R$_{Hill}$]', fontsize=11)
    ax2.set_ylabel('N(x, y=B)', fontsize=11)
    ax2.set_title('Perfil columnar en y = B_impact', fontsize=11)
    ax2.legend(fontsize=9); ax2.grid(True, alpha=0.3)

    plt.suptitle('Proyeccion en el plano del cielo', fontsize=12)
    plt.tight_layout()
    if SAVE_FIGURES:
        fname = os.path.join(OUTPUT_DIR, 'fig2_columna_densidad.png')
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"  Guardado: {fname}")
    plt.show()


# ============================================================
#  6. CURVA DE LUZ DEL TRANSITO
# ============================================================

def compute_transit_lightcurve(x_sky, y_sky, N_col):
    """
    Mueve el planeta frente a la estrella (en X_sky) y calcula
    el flujo bloqueado en cada paso.
    La estrella tiene disco uniforme (Lambertian flat).
    """
    print("\nCalculando curva de luz del transito...")

    # Kappa: normalizado para que tau_max ~ 2
    tau_max_target = 2.0
    N_max = N_col.max()
    kappa = tau_max_target / N_max if N_max > 0 else 1.0
    print(f"  kappa = {kappa:.4e}  (tau_max ~ {tau_max_target:.1f})")

    # Interpolador del mapa de columna
    N_interp = RegularGridInterpolator(
        (x_sky, y_sky), N_col,
        method='linear', bounds_error=False, fill_value=0.0)

    # Grilla de la estrella
    N_STAR = 100
    xs = np.linspace(-R_STAR, R_STAR, N_STAR)
    ys = np.linspace(-R_STAR, R_STAR, N_STAR)
    XS, YS = np.meshgrid(xs, ys)
    on_star = (XS**2 + YS**2) <= R_STAR**2
    total_star = on_star.sum()

    x_transit = np.linspace(-R_STAR * 2.5, R_STAR * 2.5, N_TRANSIT_STEPS)
    flux = np.zeros(N_TRANSIT_STEPS)

    for i, xp in enumerate(x_transit):
        x_rel = (XS[on_star] - xp).ravel()
        y_rel = (YS[on_star] - B_IMPACT).ravel()
        pts   = np.column_stack([x_rel, y_rel])
        N_los = N_interp(pts)
        tau   = kappa * N_los
        flux[i] = np.sum(np.exp(-tau)) / total_star
        if (i + 1) % 50 == 0:
            print(f"  Paso {i+1}/{N_TRANSIT_STEPS}, flujo={flux[i]:.5f}")

    return x_transit, flux


def plot_lightcurve(x_transit, flux):
    print("\nGenerando curva de luz...")
    fig, ax = plt.subplots(figsize=(10, 5))
    t_norm = x_transit / R_STAR
    ax.plot(t_norm, flux * 100, 'steelblue', lw=2)
    ax.axhline(100, color='gray', ls='--', lw=1, alpha=0.6, label='Sin transito')
    depth = (1 - flux.min()) * 100
    i_min = np.argmin(flux)
    ax.annotate(f'Profundidad: {depth:.3f}%',
                xy=(t_norm[i_min], flux[i_min]*100),
                xytext=(0.65, 0.2), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
                fontsize=11, color='red')
    ax.set_xlabel('X$_{sky}$ / R$_*$', fontsize=12)
    ax.set_ylabel('Flujo normalizado [%]', fontsize=12)
    ax.set_title(
        f'Curva de luz del transito\n'
        f'2Δφ={2*DELTA_PHI:.0f}°  B_impact={B_IMPACT:.2f} R_H  '
        f'R_ext={R_EXT_SHOCK} R_H', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    plt.tight_layout()
    if SAVE_FIGURES:
        fname = os.path.join(OUTPUT_DIR, 'fig3_curva_luz.png')
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"  Guardado: {fname}")
    plt.show()


# ============================================================
#  MAIN
# ============================================================

if __name__ == '__main__':
    print("=" * 55)
    print("  Viento Planetario 3D + Bow Shock + Transito")
    print("=" * 55)
    print(f"\nParametros:")
    print(f"  R_planeta={R_PLANET}  R_Hill={R_HILL}  R_estrella={R_STAR}")
    print(f"  Parche: 2Δφ={2*DELTA_PHI:.0f}°  θ={THETA_PATCH:.0f}°  φ₀={PHI0_PATCH:.0f}°")
    print(f"  R_ext_choque={R_EXT_SHOCK} R_Hill  B_impact={B_IMPACT}")

    print("\nCargando perfiles...")
    wind_func  = make_wind_profile()
    shock_func = make_shock_profile()

    # Fig 1: cortes 3D
    plot_3d_density(wind_func, shock_func)

    # Fig 2: proyeccion plano del cielo
    x_sky, y_sky, N_col = compute_column_density(wind_func, shock_func)
    plot_column_density(x_sky, y_sky, N_col)

    # Fig 3: curva de luz
    x_transit, flux = compute_transit_lightcurve(x_sky, y_sky, N_col)
    plot_lightcurve(x_transit, flux)

    print("\nListo!")