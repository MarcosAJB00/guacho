#!/usr/bin/env python3
"""
Curva de luz del triplete He 10830 con geometría del tránsito correcta.

Dos grillas independientes:
  - Grilla simulación (dx_sim): la caja del modelo 3D, 256x256 px, ~10 Rp de lado
  - Grilla tránsito  (dx_T)  : mapa del cielo, lo suficientemente grande para
                                contener la estrella completa

Para cada fase orbital:
  - Si el píxel del mapa tránsito está dentro del planeta → emisión = 0
  - Si está fuera del planeta pero dentro de la caja sim  → usa emtaunu del .bin
  - Si está fuera de la caja sim (material estelar sin modelo) → emisión = 1
"""

from scipy.io import FortranFile
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Parámetros físicos
# ---------------------------------------------------------------------------
Rjup  = 7.1492e9        # cm
Rsun  = 6.955e10        # cm
Rstar = 0.76  * Rsun   # cm
Rp    = 0.4466 * Rjup  # cm
rorb  = 0.0532 * 1.49e13  # cm  (semieje orbital)

# ---------------------------------------------------------------------------
# Grilla de la simulación (la caja del modelo 3D)
# ---------------------------------------------------------------------------
nxmap_sim = 256
nymap_sim = 256
nvmap     = 250
dx_sim    = 2.0 * 5.09421794574017 * Rp / float(nxmap_sim)

# Verificación
print("="*55)
print("Grilla simulación")
print(f"  dx_sim     = {dx_sim:.4e} cm")
print(f"  caja       = {nxmap_sim*dx_sim/Rp:.2f} Rp de lado")
print(f"  rstar_sim  = {Rstar/dx_sim:.1f} px  (MAYOR que {nxmap_sim//2} → no cabe)")
print(f"  rplan_sim  = {Rp/dx_sim:.2f} px")

# ---------------------------------------------------------------------------
# Grilla del tránsito (mapa del cielo, independiente de la simulación)
# Elegimos que quepa 2.5 diámetros estelares en cada dirección
# ---------------------------------------------------------------------------
nxmap_T = 1032
nymap_T = 1032 
dx_T    = 2.5 * 2.0 * Rstar / float(nxmap_T)

rstar_T = Rstar / dx_T
rplan_T = Rp    / dx_T

print()
print("Grilla tránsito")
print(f"  dx_T       = {dx_T:.4e} cm")
print(f"  caja       = {nxmap_T*dx_T/Rstar:.2f} R★ de lado")
print(f"  rstar_T    = {rstar_T:.1f} px")
print(f"  rplan_T    = {rplan_T:.2f} px")
print(f"  DF_geom    = (Rp/Rstar)^2 = {(Rp/Rstar)**2:.6f}")
print("="*55)

# ---------------------------------------------------------------------------
# Parámetros espectrales
# ---------------------------------------------------------------------------
velinit  = -150e5   # cm/s
velfinal =  150e5   # cm/s
vel_axis = np.linspace(velinit, velfinal, nvmap)

typeouts = ['3pj0', '3pj1', '3pj2']

# ---------------------------------------------------------------------------
# Directorios
# ---------------------------------------------------------------------------
path     = Path('./helines/')
path_dat = path / 'dat'
path_dat.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Loop principal
# ---------------------------------------------------------------------------
run = 0

for typeout in typeouts:
    print(f"\nProcesando línea: {typeout}")

    fout1 = open(path_dat / f'HE_tau_{run:03d}_{typeout}.dat', 'w')

    for nout in range(0, 81):

        # --- Leer cubo de tau de la simulación ---
        filein = path / f'HE_tau_{typeout}_{nout:03d}.bin'
        f = FortranFile(str(filein), 'r')
        data   = f.read_reals('d').reshape(nxmap_sim, nymap_sim, nvmap, order='F').T
        f.close()
        emtaunu = np.exp(-data)   # shape: (nvmap, nymap_sim, nxmap_sim)

        # --- Posición del planeta en el mapa del tránsito ---
        theta_plan = (nout - 40) * np.pi * 0.125 / 180.0
        xplan_T    = rorb / dx_T * np.tan(theta_plan)   # en px del mapa tránsito
        yplan_T    = 0.0                                 # b = 0, tránsito central

        EmissionVel   = np.zeros(nvmap)
        TotalEmission = 0.0
        Ref           = 0.0
        n_planet      = 0
        n_sim         = 0
        n_outside     = 0

        # --- Loop sobre la grilla del tránsito ---
        for i in range(nxmap_T):
            for j in range(nymap_T):

                # Coordenadas centradas en el mapa del tránsito (en px)
                x_T = (float(i - nxmap_T / 2) + 0.5)
                y_T = (float(j - nymap_T / 2) + 0.5)

                # Distancia al centro de la estrella
                r1 = np.sqrt(x_T**2 + y_T**2)

                if r1 > rstar_T:
                    continue   # fuera del disco estelar, no contribuye

                Ref += 1.0

                # Distancia al centro del planeta
                r2 = np.sqrt((x_T - xplan_T)**2 + (y_T - yplan_T)**2)

                if r2 < rplan_T:
                    # Píxel tapado por el planeta opaco → emisión = 0
                    n_planet += 1
                    # (no sumamos nada a EmissionVel ni TotalEmission)

                else:
                    # Convertir posición del mapa tránsito a posición física (cm)
                    x_cm = x_T * dx_T
                    y_cm = y_T * dx_T

                    # Índice correspondiente en la grilla de simulación
                    i_sim = int(x_cm / dx_sim + nxmap_sim / 2)
                    j_sim = int(y_cm / dx_sim + nymap_sim / 2)

                    if 0 <= i_sim < nxmap_sim and 0 <= j_sim < nymap_sim:
                        # Dentro de la caja del modelo: usamos emtaunu real
                        EmissionVel   += emtaunu[:, j_sim, i_sim]
                        TotalEmission += np.sum(emtaunu[:, j_sim, i_sim])
                        n_sim += 1
                    else:
                        # Fuera de la caja del modelo: superficie estelar sin
                        # absorción del viento planetario → emtaunu = 1
                        EmissionVel   += 1.0
                        TotalEmission += nvmap
                        n_outside += 1

        # --- Normalización ---
        EmissionVel   = EmissionVel   / Ref
        TotalEmission = TotalEmission / Ref / float(nvmap)

        absorcion_pct = (1.0 - TotalEmission) * 100.0

        print(f"  fase {nout:2d}  absorción = {absorcion_pct:.3f}%"
              f"  [planeta={n_planet}px  sim={n_sim}px  fuera={n_outside}px"
              f"  ref={int(Ref)}px]")

        # --- Guardar perfil de velocidades ---
        fout2 = open(path_dat / f'HE_tau-{run:03d}-abs-{nout:03d}_{typeout}.dat', 'w')
        for l in range(nvmap):
            fout2.write(f'{vel_axis[l]:.6e}  {EmissionVel[l]:.8f}\n')
        fout2.close()

        # --- Guardar resumen de fase ---
        fout1.write(f'{nout}  {np.mean(EmissionVel):.8f}\n')

        # --- Mapa de tau integrado en velocidad ---
        fig, ax = plt.subplots(figsize=(5, 5))
        im = ax.imshow(
            np.sum(emtaunu, axis=0) / nvmap,
            origin='lower',
            norm=plt.matplotlib.colors.LogNorm(vmin=9e-1, vmax=1.0),
            extent=[-nxmap_sim/2, nxmap_sim/2, -nymap_sim/2, nymap_sim/2]
        )
        # Dibujar estrella y planeta en coordenadas de la caja sim
        star_in_sim  = plt.Circle((0, 0), Rstar/dx_sim,
                                   color='yellow', fill=False, lw=0.8, ls='--')
        plan_in_sim  = plt.Circle((xplan_T * dx_T / dx_sim, 0), Rp/dx_sim,
                                   color='cyan',   fill=False, lw=0.8)
        ax.add_patch(star_in_sim)
        ax.add_patch(plan_in_sim)
        ax.set_title(f'{typeout}  fase {nout:02d}  abs={absorcion_pct:.3f}%', fontsize=9)
        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        plt.savefig(path_dat / f'{typeout}_tau_map-{run:03d}-abs-{nout:03d}.png', dpi=90)
        plt.close()

    fout1.close()

print("\nListo.")