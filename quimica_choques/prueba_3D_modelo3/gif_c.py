#!/usr/bin/env python3
"""
Animación de tránsito planetario con curva de luz del triplete He 10830
Lee archivos .bin (FortranFile) o .dat pre-procesados generados por el código original.

Uso:
    python transit_animation.py               # Lee .dat por defecto
    python transit_animation.py --from-bin    # Lee .bin y computa la absorción
    python transit_animation.py --save        # Guarda un GIF en lugar de mostrar
    python transit_animation.py --save --output mi_video.gif
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec
from pathlib import Path

# ---------------------------------------------------------------------------
# Parámetros físicos (del código original)
# ---------------------------------------------------------------------------
Rsun    = 6.955e10       # cm
Rjup    = 7.1492e9       # cm
Rstar   = 0.76 * Rsun   # cm
Rp      = 0.4466 * Rjup # cm
rorb    = 0.0532 * 1.49e13  # cm

nphases = 81
nxmap   = 256
nymap   = 256
nvmap   = 250
velinit = -150e5
velfinal=  150e5

typeouts = ['3pj0', '3pj1', '3pj2']

# Tamaño del pixel de la caja de simulación
dx = 2*5.09421794574017 * Rp / float(nxmap)  # cm/pixel
rstar_px = Rstar / dx    # radio estelar en píxeles de la caja
rplan_px = Rp / dx       # radio planetario en píxeles de la caja
rorb_px  = rorb / dx     # radio orbital en píxeles

miss_pix = 144044        # píxeles faltantes (bordes de la estrella)

DF = (Rp / Rstar) ** 2  # profundidad de tránsito geométrica

print(f"Rstar  = {Rstar/Rsun:.3f} Rsun")
print(f"Rp     = {Rp/Rjup:.4f} Rjup")
print(f"DF     = (Rp/Rstar)^2 = {DF*100:.3f}%")
print(f"rstar_px = {rstar_px:.1f} px  |  rplan_px = {rplan_px:.1f} px")


# ---------------------------------------------------------------------------
# Funciones de carga de datos
# ---------------------------------------------------------------------------

def load_from_bin(path='./helines/', run=0):
    """
    Lee los archivos .bin originales (FortranFile) y computa la absorción
    integrada en velocidad para cada fase y cada línea del triplete.
    Devuelve flux[nphases] = flujo normalizado sumado sobre las 3 líneas.
    """
    try:
        from scipy.io import FortranFile
    except ImportError:
        raise ImportError("scipy no encontrado. Instalá con: pip install scipy")

    Emission_missing_pixel = np.ones(nvmap)
    flux_per_line = np.zeros((len(typeouts), nphases))

    for li, typeout in enumerate(typeouts):
        print(f"\nLeyendo línea {typeout}...")
        for nout in range(nphases):
            filein = Path(path) / f'HE_tau_{typeout}_{str(nout).zfill(3)}.bin'
            if not filein.exists():
                # Intentar formato alternativo con run
                filein = Path(path) / f'HE_tau_{typeout}_{str(run).zfill(3)}_{str(nout).zfill(3)}.bin'
            if not filein.exists():
                print(f"  Archivo no encontrado: {filein}, usando flujo=1")
                flux_per_line[li, nout] = 1.0
                continue

            f = FortranFile(str(filein), 'r')
            data = f.read_reals('d').reshape(nxmap, nymap, nvmap, order='F').T
            f.close()
            emtaunu = np.exp(-data)

            EmissionVel  = np.zeros(nvmap)
            TotalEmission = 0.0
            Ref = 0.0

            theta = (nout - 40) * np.pi * 0.125 / 180.0
            xplan = rorb_px * np.tan(theta)
            yplan = 0.0

            for i in range(nxmap):
                for j in range(nymap):
                    x = float(i - nxmap/2) + 0.5
                    y = float(j - nymap/2) + 0.5
                    r1 = np.sqrt(x**2 + y**2)

                    if r1 <= rstar_px:
                        r2 = np.sqrt((x - xplan)**2 + (y - yplan)**2)
                        if r2 < rplan_px:
                            emtaunu[:, j, i] = 0.0
                        else:
                            EmissionVel    += emtaunu[:, j, i]
                            TotalEmission  += np.sum(emtaunu[:, j, i])
                        Ref += 1.0

            Ref += miss_pix
            EmissionVel   = (EmissionVel + miss_pix * Emission_missing_pixel) / Ref
            TotalEmission = (TotalEmission + miss_pix * nvmap) / Ref / float(nvmap)

            flux_per_line[li, nout] = TotalEmission
            print(f"  fase {nout:3d}  absorción = {(1-TotalEmission)*100:.3f}%")

    # Flujo total = promedio (o suma normalizada) de las 3 líneas
    flux_total = np.mean(flux_per_line, axis=0)
    return flux_total, flux_per_line


def load_from_dat(path='./helines/dat/', run=0):
    """
    Lee los archivos .dat pre-procesados por el código original.
    Formato esperado: HE_tau-RRR-abs-NNN_TIPO.dat  con columnas (velocidad, flujo)
    Devuelve flux[nphases] = flujo integrado sumado sobre las 3 líneas.
    """
    flux_per_line = np.zeros((len(typeouts), nphases))

    for li, typeout in enumerate(typeouts):
        print(f"Leyendo .dat para línea {typeout}...")
        for nout in range(nphases):
            fname = Path(path) / f'HE_tau-{str(run).zfill(3)}-abs-{str(nout).zfill(3)}_{typeout}.dat'
            if not fname.exists():
                print(f"  No encontrado: {fname}")
                flux_per_line[li, nout] = 1.0
                continue

            d = np.loadtxt(fname)
            vel, emv = d[:, 0], d[:, 1]
            # Flujo integrado = promedio del perfil de emisión normalizado
            flux_per_line[li, nout] = np.mean(emv)

    flux_total = np.mean(flux_per_line, axis=0)
    return flux_total, flux_per_line


def load_summary_dat(path='./helines/dat/', run=0):
    """
    Alternativa: lee los archivos resumen HE_tau_RRR_TIPO.dat
    Formato: columnas (nout, mean_EmissionVel)
    """
    flux_per_line = np.zeros((len(typeouts), nphases))

    for li, typeout in enumerate(typeouts):
        fname = Path(path) / f'HE_tau_{str(run).zfill(3)}_{typeout}.dat'
        if not fname.exists():
            print(f"  Resumen no encontrado: {fname}")
            flux_per_line[li, :] = 1.0
            continue
        d = np.loadtxt(fname)
        phases_in = d[:, 0].astype(int)
        vals      = d[:, 1]
        for idx, ph in enumerate(phases_in):
            if 0 <= ph < nphases:
                flux_per_line[li, ph] = vals[idx]

    flux_total = np.mean(flux_per_line, axis=0)
    return flux_total, flux_per_line


# ---------------------------------------------------------------------------
# Geometría del tránsito (para la visualización)
# ---------------------------------------------------------------------------

def planet_screen_pos(nout, scale_px):
    """
    Posición en pantalla (x, y) del centro del planeta para la fase `nout`.
    scale_px: factor de escala de la caja de simulación a píxeles de pantalla.
    Parámetro de impacto b = 0 (tránsito central).
    """
    theta = (nout - 40) * np.pi * 0.125 / 180.0
    x_sim = rorb_px * np.tan(theta)   # posición en píxeles de la caja
    y_sim = 0.0
    return x_sim * scale_px, y_sim * scale_px


# ---------------------------------------------------------------------------
# Animación
# ---------------------------------------------------------------------------

def make_animation(flux_total, flux_per_line, save_path=None, fps=12):

    phases = np.arange(nphases)

    # Escala visual: el radio estelar ocupa `star_display_r` puntos en pantalla
    star_display_r = 1.0   # unidades normalizadas (Rstar = 1)
    planet_display_r = Rp / Rstar
    orbit_display_r  = rorb / Rstar

    # -----------------------------------------------------------------------
    # Figura
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=(11, 9), facecolor='#0d1117')
    gs  = GridSpec(2, 1, height_ratios=[1.5, 1], hspace=0.08,
                   top=0.95, bottom=0.09, left=0.10, right=0.97)

    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1])

    for ax in [ax_top, ax_bot]:
        ax.set_facecolor('#0d1117')
        for spine in ax.spines.values():
            spine.set_edgecolor('#30363d')

    # ---- Panel superior: escena del tránsito ----
    ax_top.set_xlim(-3.5, 3.5)
    ax_top.set_ylim(-2.2, 2.2)
    ax_top.set_aspect('equal')
    ax_top.tick_params(colors='#8b949e', labelsize=8)
    ax_top.set_xlabel('x / R★', color='#8b949e', fontsize=9)
    ax_top.set_ylabel('y / R★', color='#8b949e', fontsize=9)
    ax_top.set_title('Tránsito He 10830 — HD 189733b', color='#e6edf3',
                     fontsize=11, fontweight='bold', pad=8)

    # Estrellas de fondo (decorativas)
    rng = np.random.default_rng(42)
    bx = rng.uniform(-3.4, 3.4, 120)
    by = rng.uniform(-2.1, 2.1, 120)
    bs = rng.uniform(0.3, 1.2, 120)
    ax_top.scatter(bx, by, s=bs, color='white', alpha=0.4, zorder=0)

    # Órbita (arco punteado)
    theta_arr = np.linspace(-40*0.125*np.pi/180, 40*0.125*np.pi/180, 200)
    orb_x = orbit_display_r * np.sin(theta_arr)
    orb_y = np.zeros_like(orb_x)
    ax_top.plot(orb_x, orb_y, '--', color='#30363d', lw=0.8, zorder=1)

    # Estrella (manchas solares opcionales + degradado radial via scatter)
    star_base = plt.matplotlib.patches.Circle(
        (0, 0), star_display_r, color='#FDB813', zorder=2, linewidth=0)
    ax_top.add_patch(star_base)
    # Oscurecimiento de borde (limb darkening) simulado con círculo más oscuro
    for r_frac, alpha in [(1.0, 0.30), (0.85, 0.12), (0.60, 0.05)]:
        ld = plt.matplotlib.patches.Circle(
            (0, 0), star_display_r * r_frac,
            color='#7a3800', alpha=alpha, zorder=3, linewidth=0)
        ax_top.add_patch(ld)

    # Label estrella
    ax_top.text(0, -star_display_r - 0.18, 'HD 189733 (0.76 R☉)',
                ha='center', va='top', color='#FDB813', fontsize=7.5, zorder=10)

    # Exosferas (3 líneas, radios distintos)
    exo_factors = [2.2, 2.5, 2.8]
    exo_colors  = ['#4a9eff', '#7ab8ff', '#aad0ff']
    exo_patches = []
    for fi, (fc, ec) in enumerate(zip(exo_factors, exo_colors)):
        ep = plt.matplotlib.patches.Circle(
            (orbit_display_r, 0), planet_display_r * fc,
            color=ec, alpha=0.10 + 0.03*(2-fi), zorder=4, linewidth=0)
        ax_top.add_patch(ep)
        exo_patches.append(ep)

    # Planeta
    planet_patch = plt.matplotlib.patches.Circle(
        (orbit_display_r, 0), planet_display_r,
        color='#111827', zorder=5, linewidth=0.8,
        edgecolor='#4a9eff')
    ax_top.add_patch(planet_patch)

    # Label planeta
    planet_label = ax_top.text(
        orbit_display_r, planet_display_r + 0.12, 'HD 189733b',
        ha='center', va='bottom', color='#4a9eff', fontsize=7, zorder=11)

    # Texto de fase y ángulo
    phase_text = ax_top.text(
        -3.35, 1.95, '', ha='left', va='top', color='#e6edf3',
        fontsize=9, fontfamily='monospace', zorder=12)
    abs_text = ax_top.text(
        3.35, 1.95, '', ha='right', va='top', color='#f47067',
        fontsize=9, fontfamily='monospace', zorder=12)

    # Leyenda de exosferas
    for fi, (typeout, ec) in enumerate(zip(typeouts, exo_colors)):
        ax_top.plot([], [], 'o', color=ec, alpha=0.6,
                    markersize=5, label=f'Exosfera {typeout}')
    ax_top.legend(loc='lower right', fontsize=7, facecolor='#161b22',
                  edgecolor='#30363d', labelcolor='#8b949e', markerscale=0.9)

    # ---- Panel inferior: curva de luz ----
    ax_bot.set_xlim(-0.5, nphases - 0.5)
    ypad = 0.002
    ymin = np.min(flux_total) - ypad
    ymax = 1.0 + ypad
    ax_bot.set_ylim(ymin, ymax)
    ax_bot.set_xlabel('Fase orbital', color='#8b949e', fontsize=9)
    ax_bot.set_ylabel('Flujo normalizado', color='#8b949e', fontsize=9)
    ax_bot.tick_params(colors='#8b949e', labelsize=8)
    ax_bot.axhline(1.0, color='#30363d', lw=0.8, ls='--', zorder=1)

    # Curvas individuales por línea (finas, en fondo)
    line_colors_lc = ['#58a6ff', '#3fb950', '#f78166']
    for li, (typeout, lc) in enumerate(zip(typeouts, line_colors_lc)):
        ax_bot.plot(phases, flux_per_line[li], color=lc, lw=0.7,
                    alpha=0.35, label=typeout, zorder=2)

    # Curva total (suma de 3 líneas) — se va dibujando
    lc_future, = ax_bot.plot(phases, flux_total, color='#e6edf3',
                              lw=0.6, alpha=0.15, ls='--', zorder=3)
    lc_drawn,  = ax_bot.plot([], [], color='#f47067', lw=2.0, zorder=5)
    lc_dot,    = ax_bot.plot([], [], 'o', color='#f47067',
                              markersize=5, zorder=6)
    vline = ax_bot.axvline(0, color='#f47067', lw=0.8, alpha=0.4, zorder=4)

    ax_bot.legend(loc='lower center', fontsize=7, ncol=4,
                  facecolor='#161b22', edgecolor='#30363d', labelcolor='#8b949e')

    ax_bot.text(0.01, 0.97, 'Suma 3pj0 + 3pj1 + 3pj2',
                transform=ax_bot.transAxes, ha='left', va='top',
                color='#f47067', fontsize=7.5)

    # -----------------------------------------------------------------------
    # Función de actualización
    # -----------------------------------------------------------------------
    def update(nout):
        # Posición del planeta
        theta = (nout - 40) * np.pi * 0.125 / 180.0
        xp = orbit_display_r * np.sin(theta)
        yp = 0.0

        planet_patch.set_center((xp, yp))
        for ep, fc in zip(exo_patches, exo_factors):
            ep.set_center((xp, yp))

        planet_label.set_position((xp, yp + planet_display_r + 0.12))

        # Color del planeta: negro si está sobre la estrella, azul si no
        sep = np.sqrt(xp**2 + yp**2)
        if sep < star_display_r + planet_display_r:
            planet_patch.set_facecolor('#050a10')
        else:
            planet_patch.set_facecolor('#1a3a5c')

        # Textos
        theta_deg = (nout - 40) * 0.125
        phase_text.set_text(f'fase: {nout:02d}/80   θ={theta_deg:+.3f}°')
        absorcion = (1.0 - flux_total[nout]) * 100
        abs_text.set_text(f'absorción: {absorcion:.3f}%')

        # Curva de luz
        lc_drawn.set_data(phases[:nout+1], flux_total[:nout+1])
        lc_dot.set_data([phases[nout]], [flux_total[nout]])
        vline.set_xdata([nout])

        return (planet_patch, *exo_patches, planet_label,
                phase_text, abs_text, lc_drawn, lc_dot, vline)

    # -----------------------------------------------------------------------
    # Crear animación
    # -----------------------------------------------------------------------
    ani = animation.FuncAnimation(
        fig, update, frames=nphases,
        interval=1000/fps, blit=True)

    if save_path:
        ext = Path(save_path).suffix.lower()
        print(f"\nGuardando animación en: {save_path}")
        if ext == '.gif':
            writer = animation.PillowWriter(fps=fps)
            ani.save(save_path, writer=writer, dpi=130)
        elif ext == '.mp4':
            writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
            ani.save(save_path, writer=writer, dpi=130)
        else:
            ani.save(save_path, dpi=130)
        print("¡Listo!")
    else:
        plt.show()

    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Animación del tránsito He 10830 de HD 189733b')
    parser.add_argument('--from-bin', action='store_true',
                        help='Leer archivos .bin crudos (más lento)')
    parser.add_argument('--from-summary', action='store_true',
                        help='Leer archivos resumen HE_tau_RRR_TIPO.dat')
    parser.add_argument('--path', default='./helines/',
                        help='Directorio raíz de los datos (default: ./helines/)')
    parser.add_argument('--dat-path', default=None,
                        help='Directorio de los .dat (default: PATH/dat/)')
    parser.add_argument('--run', type=int, default=0,
                        help='Número de run de GUACHO (default: 0)')
    parser.add_argument('--save', action='store_true',
                        help='Guardar como archivo en lugar de mostrar')
    parser.add_argument('--output', default='transit_he10830.gif',
                        help='Archivo de salida (default: transit_he10830.gif)')
    parser.add_argument('--fps', type=int, default=12,
                        help='Fotogramas por segundo (default: 12)')
    args = parser.parse_args()

    dat_path = args.dat_path or str(Path(args.path) / 'dat')

    print("=" * 55)
    print("  Tránsito He 10830 — HD 189733b")
    print("=" * 55)

    if args.from_bin:
        print(f"\nModo: leyendo .bin desde {args.path}")
        flux_total, flux_per_line = load_from_bin(args.path, args.run)
    elif args.from_summary:
        print(f"\nModo: leyendo resumen .dat desde {dat_path}")
        flux_total, flux_per_line = load_summary_dat(dat_path, args.run)
    else:
        print(f"\nModo: leyendo .dat por fase desde {dat_path}")
        flux_total, flux_per_line = load_from_dat(dat_path, args.run)

    print(f"\nFlujo mínimo: {np.min(flux_total):.5f}")
    print(f"Absorción máxima: {(1-np.min(flux_total))*100:.3f}%")
    print(f"Fase de mínimo: {np.argmin(flux_total)}")

    save_path = args.output if args.save else None
    make_animation(flux_total, flux_per_line, save_path=save_path, fps=args.fps)


if __name__ == '__main__':
    main()