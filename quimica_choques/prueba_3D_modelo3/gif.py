#!/usr/bin/python
import os
import numpy as np
from scipy.io import FortranFile
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm

# --- Configuración de parámetros originales ---
endian = '<'      # little endian
kind = 'd'        # double precision

nxmap = 256
nymap = 256
nvmap = 250

rorb = 0.0532 * 1.49e13
torb = 4.887802443 * 86400.
Rjup = 7.1492E9
Rsun = 6.955e10
Rstar = 0.76 * Rsun
Rp = 0.4466 * Rjup
dx = 2.0 * 5.09421794574017 * Rp / float(nxmap)

rstar = Rstar / dx     # Radio estelar en píxeles
rplan = Rp / dx       # Radio del planeta en píxeles
miss_pix = 144044

typeouts = ['3pj0', '3pj1', '3pj2']
path = './helines/'
num_phases = 81

# --- Listas para almacenar datos precalculados ---
phases = np.arange(num_phases)
xplan_all = np.zeros(num_phases)
absorption_total = np.zeros(num_phases)
maps_visual = np.zeros((num_phases, nymap, nxmap))

print("Cargando archivos y precalculando curvas de luz...")

# --- Etapa 1: Procesamiento de datos en memoria ---
Emission_missing_pixel = np.ones(nvmap)

for nout in range(num_phases):
    theta_plan = (nout - 40) * np.pi * 0.125 / 180
    xplan = rorb / dx * np.tan(theta_plan)
    xplan_all[nout] = xplan
    
    abs_lines = []
    map_lines = []
    
    for typeout in typeouts:
        filein = f"{path}HE_tau_{typeout}_{str(nout).zfill(3)}.bin"
        
        if os.path.exists(filein):
            f = FortranFile(filein, 'r')
            data = f.read_reals(kind).reshape(nxmap, nymap, nvmap, order='F').T
            emtaunu = np.exp(-data)
        else:
            # Si falta algún archivo, usamos matriz vacía (transmisión 1)
            emtaunu = np.ones((nvmap, nymap, nxmap))
            
        EmissionVel = np.zeros(nvmap)
        TotalEmission = 0.
        Ref = 0.
        
        # Máscara del núcleo opaco del planeta
        for i in range(nxmap):
            for j in range(nymap):
                x = (float(i - nxmap / 2) + 0.5)
                y = (float(j - nymap / 2) + 0.5)
                r1 = np.sqrt(x**2 + y**2)
                
                # Nota: r1 <= rstar siempre se cumple en la caja central de 256x256
                r2 = np.sqrt((x - xplan)**2 + y**2)
                if r2 < rplan:
                    emtaunu[:, j, i] = 0.
                else:
                    EmissionVel += emtaunu[:, j, i]
                    TotalEmission += np.sum(emtaunu[:, j, i])
                Ref += 1.
                
        Ref += miss_pix
        TotalEmission = (TotalEmission + miss_pix * 250) / Ref / float(nvmap)
        
        # Guardamos la absorción de esta línea (%) -> (1 - Transmisión) * 100
        abs_lines.append((1.0 - TotalEmission) * 100)
        
        # Guardamos el mapa promediado en velocidades para visualización
        map_lines.append(np.sum(emtaunu, axis=0) / nvmap)
        
    # Sumamos las absorciones de las 3 líneas para la curva de luz final
    absorption_total[nout] = sum(abs_lines)
    # Promedio visual de transmisión de las 3 líneas para el mapa central
    maps_visual[nout] = np.mean(map_lines, axis=0)

print("¡Datos cargados con éxito! Generando animación...")

# --- Etapa 2: Configuración de los Gráficos para la Animación ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 11))
plt.subplots_adjust(hspace=0.3)

# --- PANEL SUPERIOR: Geometría del Tránsito ---
ax1.set_title("Tránsito Planetario y Mapa de Transmisión de He I", fontsize=14, pad=10)
# Dibujamos el disco completo de la estrella
star_disk = plt.Circle((0, 0), rstar, color='#FFA500', alpha=0.4, label='Superficie Estelar')
ax1.add_patch(star_disk)

# Dibujamos el límite de tu caja de simulación de 10 Rp en el centro
box_rect = plt.Rectangle((-nxmap/2, -nymap/2), nxmap, nymap, fill=False, 
                         edgecolor='black', linestyle='--', alpha=0.6, label='Caja de Simulación (~10 $R_p$)')
ax1.add_patch(box_rect)

# Inicializamos el mapa de calor (imshow) dentro de la caja central
im_extent = [-nxmap/2, nxmap/2, -nymap/2, nymap/2]
im = ax1.imshow(maps_visual[0], origin='lower', extent=im_extent, 
                norm=LogNorm(vmin=0.9, vmax=1.0), cmap='inferno', zorder=2)

# Inicializamos el planeta en su posición física (puede estar fuera de la caja)
planet_circle = plt.Circle((xplan_all[0], 0), rplan, color='black', zorder=3, label='Núcleo del Planeta')
ax1.add_patch(planet_circle)

# Ajustamos límites para ver TODA la estrella y el trayecto del planeta
ax1.set_xlim(-rstar * 1.4, rstar * 1.4)
ax1.set_ylim(-rstar * 1.4, rstar * 1.4)
ax1.set_aspect('equal')
ax1.set_xlabel("Posición X [píxeles]", fontsize=11)
ax1.set_ylabel("Posición Y [píxeles]", fontsize=11)
ax1.legend(loc='upper right', framealpha=0.9)

# Barra de color para el mapa de transmisión
cbar = fig.colorbar(im, ax=ax1, orientation='vertical', pad=0.02, shrink=0.8)
cbar.set_label('Transmisión Promedio Espectral', fontsize=11)

# --- PANEL INFERIOR: Curva de Luz 'En Vivo' ---
ax2.set_title("Curva de Absorción Total Sumada (3 Líneas)", fontsize=14, pad=10)
# Línea tenue de fondo para mostrar la trayectoria completa que se va a rellenar
ax2.plot(phases, absorption_total, color='gray', linestyle=':', alpha=0.5)

# Líneas dinámicas que se actualizarán cuadro por cuadro
line_lc, = ax2.plot([], [], color='crimson', lw=2.5, label='Absorción Combinada')
point_lc, = ax2.plot([], [], color='crimson', marker='o', markersize=8)

ax2.set_xlim(0, num_phases - 1)
ax2.set_ylim(-0.1, np.max(absorption_total) * 1.2)
ax2.set_xlabel("Fase Orbital (Índice de Salida)", fontsize=11)
ax2.set_ylabel("Absorción Total Sumada [%]", fontsize=11)
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend(loc='upper left')

# --- Función de actualización para cada fotograma (frame) ---
def update(frame):
    # 1. Actualizar panel superior
    im.set_array(maps_visual[frame])
    planet_circle.center = (xplan_all[frame], 0)
    
    # 2. Actualizar panel inferior (Curva de luz acumulativa)
    line_lc.set_data(phases[:frame+1], absorption_total[:frame+1])
    point_lc.set_data([phases[frame]], [absorption_total[frame]])
    
    return im, planet_circle, line_lc, point_lc

# --- Creación de la animación ---
ani = animation.FuncAnimation(fig, update, frames=num_phases, interval=100, blit=True)

# --- Guardar el archivo ---
# Opción A: Guardar como GIF (Ideal para Powerpoint/Keynote y no requiere software externo)
output_gif = "transito_helio_completo.gif"
print(f"Guardando animación en {output_gif}...")
ani.save(output_gif, writer='pillow', fps=10)

# Opción B: Guardar como MP4 de alta calidad (Descomenta las líneas de abajo si tienes ffmpeg instalado)
# output_mp4 = "transito_helio_completo.mp4"
# print(f"Guardando animación en {output_mp4}...")
# ani.save(output_mp4, writer='ffmpeg', fps=10, dpi=200)

plt.close()
print("¡Proceso finalizado con éxito!")