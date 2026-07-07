#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
transit_lightcurve.py
======================

Calcula la curva de luz (banda ancha y resuelta en velocidad) del transito
de un exoplaneta en la linea del triplete de Helio 10830 A (componente 3P-2,
"3pj2"), a partir de los mapas de profundidad optica (tau) generados por el
codigo Fortran `helines_tau_map` (modulo tau_utilities, version "primitivas").

------------------------------------------------------------------------------
QUE HACE
------------------------------------------------------------------------------
1. Lee los 81 archivos binarios HE_tau_3pj2_000.bin ... HE_tau_3pj2_080.bin
   (uno por paso de rotacion theta_y, formato Fortran "unformatted" con
   record markers, real(kind=8), shape (nxmap, nymap, nvmap) en orden Fortran).

2. Convierte el indice de fase i -> angulo orbital theta_y -> tiempo real,
   usando el periodo P y el semieje a (orbita circular, transito central).

3. Construye una grilla 2D que SI contiene a la estrella (con limb darkening
   opcional) usando el mismo tamano de pixel que la caja del planeta, para
   poder "pegar" el mapa de tau por indices, sin interpolar.

4. Para cada fase: ubica el mapa de tau del planeta en la posicion orbital
   correspondiente sobre la grilla estelar, agrega el disco solido opaco
   del planeta (radio Rp), convierte tau -> transmitancia exp(-tau), integra
   el flujo total y lo normaliza por el flujo fuera de transito.

5. Genera:
   - Curva de luz integrada en banda ancha de la linea (flujo vs tiempo)
   - Espectro de transmision resuelto en velocidad (mapa tiempo x velocidad)
   - Archivos .csv con los resultados numericos
   - Figuras .png

------------------------------------------------------------------------------
COMO USAR
------------------------------------------------------------------------------
1. Editar el bloque "PARAMETROS DEL USUARIO" mas abajo.
2. Poner los archivos HE_tau_3pj2_XXX.bin en `BIN_DIR`.
3. Correr: python3 transit_lightcurve.py

Si los archivos NO se leen bien (valores absurdos / NaN / shape mismatch),
lo mas probable es que el record-marker de Fortran no sea de 4 bytes (ver
funcion `read_fortran_unformatted_3d`, parametro `marker_bytes`).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# =============================================================================
# PARAMETROS DEL USUARIO
# =============================================================================

# ---- Rutas -------------------------------------------------------------
BIN_DIR   = "./helines"          # carpeta donde estan HE_tau_3pj2_XXX.bin
OUT_DIR   = "./resultados_transito"
FILE_TEMPLATE = "HE_tau_3pj2_{:03d}.bin"

# ---- Grilla de la caja del planeta (simulacion) -------------------------
NX_BOX, NY_BOX, NZ_BOX = 256, 256, 256     # nxtot, nytot, nztot
NVMAP = 250                                  # canales de velocidad
VMIN, VMAX = -150.0e5, 150.0e5               # cm/s (igual que en el Fortran)

# Diametro total de la caja simulada = 10 Rp (ASUNCION, ver mensaje).
# Si en realidad xmax es el SEMI-lado (caja de 20 Rp de diametro), cambiar
# BOX_DIAMETER_IN_RP a 20.0.
BOX_DIAMETER_IN_RP = 2.0*5.09422

# ---- Parametros fisicos del planeta y la estrella ------------------------
RJUP_CM   = 6.9911e9          # radio de Jupiter en cm
RSUN_CM   = 6.957e10          # radio solar en cm
AU_CM     = 1.495978707e13    # unidad astronomica en cm

RP_RJUP   = 0.4466            # radio del planeta [Rjup]
RSTAR_RSUN = 0.76             # radio de la estrella [Rsun]

PERIOD_DAYS = 4.887802443     # periodo orbital [dias]
A_AU        = 0.0532          # semieje orbital [AU]

# Geometria orbital: se asume orbita circular y transito central (b=0),
# consistente con que el Fortran solo rota la caja en theta_y (no la
# traslada en Y). Si tu planeta tiene un parametro de impacto b != 0,
# avisame y lo agrego (offset constante en Y).
IMPACT_PARAMETER_RSTAR = 0.0

# Rango de fases que escribio el Fortran: theta_y en [THETA_MIN, THETA_MAX]
# grados, en NSTEPS+1 pasos (i = 0..NSTEPS).
THETA_MIN_DEG, THETA_MAX_DEG = -5.0, 5.0
NSTEPS = 80

# ---- Limb darkening de la estrella ----------------------------------------
# Ley lineal:  I(mu)/I(1) = 1 - u*(1-mu),   mu = cos(angulo heliocentrico)
# u=0 -> disco uniforme. u~0.5-0.7 es tipico para una enana K en el IR
# cercano (10830 A esta en el infrarrojo cercano, el LD ahi es mas chico
# que en el visible). Cambiar aqui si tenes un valor mejor para tu estrella.
LIMB_DARKENING_U = 0.3

# ---- Resolucion / region de la grilla estelar -----------------------------
# La grilla estelar usa el MISMO tamano de pixel que la caja del planeta
# (dx_box = diametro_caja / NX_BOX), para poder pegar el mapa de tau por
# indices sin interpolar. El ancho total de la grilla se calcula
# automaticamente para cubrir la estrella + el recorrido orbital + la caja,
# con un margen de seguridad MARGIN_FACTOR.
MARGIN_FACTOR = 1.3

# Si la grilla resultante es muy grande y el calculo es lento, se puede
# fijar un limite superior de pixeles por lado (recorta el margen, no la
# fisica) cambiando MAX_GRID_SIDE. None = sin limite.
MAX_GRID_SIDE = 2200

# ---- Record marker de Fortran (bytes) -------------------------------------
# gfortran / ifort por defecto usan 4 bytes. Cambiar a 8 si tu compilador
# fue configurado con record markers de 8 bytes (poco comun).
FORTRAN_MARKER_BYTES = 4

# =============================================================================
# FUNCIONES
# =============================================================================

def read_fortran_unformatted_3d(filepath, shape, dtype=np.float64,
                                 marker_bytes=4, order="F"):
    """
    Lee un archivo escrito por Fortran con `form='unformatted'` (acceso
    secuencial clasico, NO 'stream'), conteniendo un unico WRITE de un
    array 3D. Estos archivos tienen un "record marker" (entero de
    `marker_bytes` bytes) antes y despues del bloque de datos, indicando
    el tamano del record en bytes.

    Parameters
    ----------
    filepath : str
    shape : tuple (nx, ny, nz)
    dtype : numpy dtype (float64 para real(kind=8))
    marker_bytes : 4 u 8
    order : 'F' porque Fortran es column-major

    Returns
    -------
    np.ndarray con la forma `shape`
    """
    nbytes_expected = int(np.prod(shape)) * np.dtype(dtype).itemsize
    marker_dtype = np.dtype(f"<u{marker_bytes}")

    with open(filepath, "rb") as f:
        raw = f.read()

    if len(raw) < 2 * marker_bytes:
        raise ValueError(f"Archivo demasiado chico: {filepath}")

    n1 = np.frombuffer(raw[:marker_bytes], dtype=marker_dtype)[0]
    n2 = np.frombuffer(raw[-marker_bytes:], dtype=marker_dtype)[0]

    if n1 != nbytes_expected or n2 != nbytes_expected:
        raise ValueError(
            f"Record marker no coincide en {filepath}: "
            f"marker_inicio={n1}, marker_fin={n2}, esperado={nbytes_expected} bytes. "
            f"Proba cambiar FORTRAN_MARKER_BYTES (actual={marker_bytes}) "
            f"a 8, o revisa NX_BOX/NY_BOX/NVMAP."
        )

    data = np.frombuffer(raw[marker_bytes:marker_bytes + nbytes_expected],
                          dtype=dtype)
    return data.reshape(shape, order=order)


def build_phase_time_array(nsteps, theta_min_deg, theta_max_deg, period_days):
    """
    Convierte el indice de fase i=0..nsteps en angulo orbital theta_y (rad)
    y en tiempo real respecto del centro del transito (theta_y=0), para una
    orbita circular: t = theta_y / (2*pi/P) = theta_y * P / (2*pi).
    """
    i_arr = np.arange(nsteps + 1)
    theta_deg = (theta_max_deg - theta_min_deg) * i_arr / nsteps + theta_min_deg
    theta_rad = np.deg2rad(theta_deg)
    t_days = theta_rad * period_days / (2.0 * np.pi)
    return i_arr, theta_rad, t_days


def build_star_grid(rstar_cm, dx_pix_cm, planet_x_max_cm, box_halfwidth_cm,
                     margin_factor=1.3, max_side=None):
    """
    Construye una grilla cuadrada centrada en la estrella, con tamano de
    pixel dx_pix_cm (igual al de la caja del planeta), suficientemente
    grande para contener la estrella completa y el recorrido orbital del
    planeta (incluyendo su propia caja).

    Returns
    -------
    xc, yc : arrays 1D con las coordenadas (cm) del CENTRO de cada pixel
    n_side : int, numero de pixeles por lado
    """
    halfwidth_needed = margin_factor * max(rstar_cm,
                                            planet_x_max_cm + box_halfwidth_cm)
    n_side = int(np.ceil(2.0 * halfwidth_needed / dx_pix_cm))
    if n_side % 2 == 0:
        n_side += 1  # impar para tener un pixel central exacto en 0

    if max_side is not None and n_side > max_side:
        n_side = max_side if max_side % 2 == 1 else max_side - 1

    half = n_side // 2
    idx = np.arange(-half, half + 1)
    xc = idx * dx_pix_cm
    yc = idx * dx_pix_cm
    return xc, yc, n_side


def stellar_brightness_map(xc, yc, rstar_cm, u_ld):
    """
    Mapa de brillo estelar I0(x,y) con limb darkening lineal:
        I(mu) = I0 * (1 - u*(1-mu)),  mu = sqrt(1-(r/Rstar)^2)
    Fuera del disco estelar, I0 = 0.
    Normalizado de forma arbitraria (solo importan razones In/Out).
    """
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    r2 = X**2 + Y**2
    on_disk = r2 <= rstar_cm**2
    mu = np.zeros_like(r2)
    mu[on_disk] = np.sqrt(1.0 - r2[on_disk] / rstar_cm**2)
    I0 = np.zeros_like(r2)
    I0[on_disk] = 1.0 - u_ld * (1.0 - mu[on_disk])
    return I0, on_disk


def add_opaque_planet_disk(transmittance_grid, xc, yc, x0_cm, y0_cm, rp_cm):
    """
    Pone transmitancia = 0 dentro del disco solido del planeta (radio rp_cm),
    centrado en (x0_cm, y0_cm). Se aplica DESPUES de pegar el mapa de tau de
    la atmosfera extendida, asi el cuerpo opaco gana donde se superponen.
    """
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    mask = (X - x0_cm)**2 + (Y - y0_cm)**2 <= rp_cm**2
    transmittance_grid[mask] = 0.0


def overlap_indices(x0_cm, y0_cm, dx_pix_cm, n_side_grid, nx_box, ny_box):
    """
    Calcula los indices de recorte (grilla y caja) para la superposicion
    de la caja del planeta (nx_box, ny_box) centrada en (x0_cm, y0_cm) sobre
    la grilla estelar (n_side_grid, n_side_grid), igual que hace
    `place_planet_map_on_grid`, pero sin tocar ningun array: solo devuelve
    los rangos de indices. Permite vectorizar operaciones por canal sin
    repetir el trabajo de recorte 250 veces.

    Returns
    -------
    (gi0, gi1, gj0, gj1) : rango de indices en la grilla estelar
    (bi0, bi1, bj0, bj1) : rango de indices correspondiente en la caja
    """
    half_grid = n_side_grid // 2
    i0 = int(round(x0_cm / dx_pix_cm)) + half_grid
    j0 = int(round(y0_cm / dx_pix_cm)) + half_grid

    half_box_x = nx_box // 2
    half_box_y = ny_box // 2

    gi0, gi1 = i0 - half_box_x, i0 - half_box_x + nx_box
    gj0, gj1 = j0 - half_box_y, j0 - half_box_y + ny_box
    bi0, bi1 = 0, nx_box
    bj0, bj1 = 0, ny_box

    if gi0 < 0:
        bi0 -= gi0
        gi0 = 0
    if gj0 < 0:
        bj0 -= gj0
        gj0 = 0
    if gi1 > n_side_grid:
        bi1 -= (gi1 - n_side_grid)
        gi1 = n_side_grid
    if gj1 > n_side_grid:
        bj1 -= (gj1 - n_side_grid)
        gj1 = n_side_grid

    return (gi0, gi1, gj0, gj1), (bi0, bi1, bj0, bj1)


# =============================================================================
# PROGRAMA PRINCIPAL
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- Constantes fisicas derivadas -------------------------------------
    Rp_cm = RP_RJUP * RJUP_CM
    Rstar_cm = RSTAR_RSUN * RSUN_CM
    a_cm = A_AU * AU_CM

    box_diam_cm = BOX_DIAMETER_IN_RP * Rp_cm
    dx_box_cm = box_diam_cm / NX_BOX           # = dxT del Fortran
    box_halfwidth_cm = box_diam_cm / 2.0

    print("="*70)
    print("Parametros derivados")
    print("="*70)
    print(f"  Rp          = {Rp_cm:.4e} cm  ({RP_RJUP} Rjup)")
    print(f"  Rstar       = {Rstar_cm:.4e} cm  ({RSTAR_RSUN} Rsun)")
    print(f"  a (semieje) = {a_cm:.4e} cm  ({A_AU} AU)  = {a_cm/Rstar_cm:.3f} Rstar")
    print(f"  Diametro caja planeta = {box_diam_cm:.4e} cm = {box_diam_cm/Rp_cm:.1f} Rp"
          f" = {box_diam_cm/Rstar_cm:.4f} Rstar")
    print(f"  dx_box (pixel)        = {dx_box_cm:.4e} cm = {dx_box_cm/Rstar_cm:.5f} Rstar")

    # ---- Fase orbital -> tiempo ---------------------------------------------
    i_arr, theta_rad, t_days = build_phase_time_array(
        NSTEPS, THETA_MIN_DEG, THETA_MAX_DEG, PERIOD_DAYS)

    # Posicion orbital proyectada del planeta en el plano del cielo.
    # x: direccion del movimiento orbital (a lo largo del transito)
    # y: parametro de impacto (constante, asumido 0 salvo que se indique)
    x_planet_cm = a_cm * np.sin(theta_rad)
    y_planet_cm = np.full_like(x_planet_cm, IMPACT_PARAMETER_RSTAR * Rstar_cm)

    x_planet_max_cm = np.max(np.abs(x_planet_cm))
    print(f"  Desplazamiento maximo del planeta en el cielo = "
          f"{x_planet_max_cm:.4e} cm = {x_planet_max_cm/Rstar_cm:.3f} Rstar")
    print(f"  Duracion total simulada = {t_days[-1]-t_days[0]:.4f} dias "
          f"({(t_days[-1]-t_days[0])*24:.2f} h)")

    # ---- Grilla estelar -----------------------------------------------------
    xc, yc, n_side = build_star_grid(
        Rstar_cm, dx_box_cm, x_planet_max_cm, box_halfwidth_cm,
        margin_factor=MARGIN_FACTOR, max_side=MAX_GRID_SIDE)
    print(f"  Grilla estelar: {n_side} x {n_side} pixeles "
          f"({n_side*dx_box_cm/Rstar_cm:.2f} Rstar de lado)")

    I0_grid, on_disk = stellar_brightness_map(xc, yc, Rstar_cm, LIMB_DARKENING_U)
    F_out = I0_grid.sum()  # flujo total fuera de transito (referencia)
    if F_out <= 0:
        raise RuntimeError("El flujo estelar de referencia es cero: revisar "
                            "Rstar / tamano de grilla.")

    # ---- Arrays de salida ----------------------------------------------------
    n_phases = NSTEPS + 1
    flux_broadband = np.ones(n_phases)             # curva de luz integrada
    flux_per_channel = np.ones((n_phases, NVMAP))  # espectro resuelto en vel.
    velocity_axis = (np.arange(NVMAP) + 0.5) * (VMAX - VMIN) / NVMAP + VMIN

    print("="*70)
    print("Procesando fases orbitales...")
    print("="*70)

    for k, i in enumerate(i_arr):
        fname = os.path.join(BIN_DIR, FILE_TEMPLATE.format(i))
        if not os.path.isfile(fname):
            print(f"  [AVISO] no encontrado: {fname} -> se omite esta fase "
                  f"(flujo queda en 1.0, fuera de transito)")
            continue

        tau_map = read_fortran_unformatted_3d(
            fname, shape=(NX_BOX, NY_BOX, NVMAP),
            dtype=np.float64, marker_bytes=FORTRAN_MARKER_BYTES, order="F")

        # --- Indices de superposicion caja-grilla (iguales para banda ----
        # ancha y para cada canal, ya que la caja no cambia de tamano) ---
        (gi0, gi1, gj0, gj1), (bi0, bi1, bj0, bj1) = overlap_indices(
            x_planet_cm[k], y_planet_cm[k], dx_box_cm, n_side, NX_BOX, NY_BOX)

        I0_sub = I0_grid[gi0:gi1, gj0:gj1]            # brillo estelar bajo la caja
        F_fixed_outside = F_out - I0_sub.sum()        # resto de la estrella, constante

        # mascara del disco opaco del planeta, recortada a la zona de la caja
        Xc_sub, Yc_sub = np.meshgrid(xc[gi0:gi1], yc[gj0:gj1], indexing="ij")
        opaque_mask_sub = ((Xc_sub - x_planet_cm[k])**2 +
                            (Yc_sub - y_planet_cm[k])**2) <= Rp_cm**2

        # --- Banda ancha: tau total = suma sobre canales de velocidad ---
        tau_broadband = tau_map.sum(axis=2)           # shape (NX_BOX, NY_BOX)
        transmittance_bb_sub = np.exp(-tau_broadband[bi0:bi1, bj0:bj1])
        transmittance_bb_sub[opaque_mask_sub] = 0.0
        flux_broadband[k] = (F_fixed_outside +
                              (I0_sub * transmittance_bb_sub).sum()) / F_out

        # --- Resuelto en velocidad: vectorizado sobre los NVMAP canales ---
        tau_sub_allv = tau_map[bi0:bi1, bj0:bj1, :]            # (nbx, nby, NVMAP)
        transmittance_sub_allv = np.exp(-tau_sub_allv)
        transmittance_sub_allv[opaque_mask_sub, :] = 0.0
        # suma ponderada por brillo estelar, para cada canal a la vez
        contrib_v = np.einsum('ij,ijk->k', I0_sub, transmittance_sub_allv)
        flux_per_channel[k, :] = (F_fixed_outside + contrib_v) / F_out

        print(f"  fase i={i:3d}  t={t_days[k]:+.4f} d  "
              f"theta_y={np.rad2deg(theta_rad[k]):+.2f} deg  "
              f"flujo_banda_ancha={flux_broadband[k]:.6f}")

    # ---- Guardar resultados ---------------------------------------------------
    out_csv_bb = os.path.join(OUT_DIR, "lightcurve_broadband.csv")
    header_bb = "i,theta_y_deg,t_days,t_hours,flux_relative"
    np.savetxt(out_csv_bb,
               np.column_stack([i_arr, np.rad2deg(theta_rad), t_days,
                                 t_days*24.0, flux_broadband]),
               header=header_bb, delimiter=",", comments="")

    out_csv_vel = os.path.join(OUT_DIR, "lightcurve_velocity_resolved.csv")
    vel_header = "t_days," + ",".join(f"v{vv/1e5:.2f}kms" for vv in velocity_axis)
    np.savetxt(out_csv_vel,
               np.column_stack([t_days, flux_per_channel]),
               header=vel_header, delimiter=",", comments="")

    print(f"\nGuardado: {out_csv_bb}")
    print(f"Guardado: {out_csv_vel}")

    # ---- Figura 1: curva de luz en banda ancha --------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(t_days*24.0, flux_broadband, "-o", ms=3, color="firebrick")
    ax.set_xlabel("Tiempo desde medio-transito [horas]")
    ax.set_ylabel("Flujo relativo (linea He 10830, banda ancha)")
    ax.set_title("Curva de luz - triplete de Helio 10830 A (3p-j2)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "lightcurve_broadband.png"), dpi=150)
    plt.close(fig)

    # ---- Figura 2: espectro de transmision resuelto en velocidad -------------
    fig, ax = plt.subplots(figsize=(8, 5))
    extent = [velocity_axis[0]/1e5, velocity_axis[-1]/1e5,
              t_days[0]*24.0, t_days[-1]*24.0]
    im = ax.imshow(flux_per_channel, aspect="auto", origin="lower",
                    extent=extent, cmap="inferno",
                    norm=Normalize(vmin=flux_per_channel.min(), vmax=1.0))
    ax.set_xlabel("Velocidad [km/s]")
    ax.set_ylabel("Tiempo desde medio-transito [horas]")
    ax.set_title("Espectro de transmision resuelto en velocidad (He 10830, 3p-j2)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Flujo relativo")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "transmission_spectrum_time_velocity.png"),
                dpi=150)
    plt.close(fig)

    print(f"Guardado: {os.path.join(OUT_DIR, 'lightcurve_broadband.png')}")
    print(f"Guardado: "
          f"{os.path.join(OUT_DIR, 'transmission_spectrum_time_velocity.png')}")
    print("\nListo.")


if __name__ == "__main__":
    main()
