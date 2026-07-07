import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# --- Datos ---
ATES    = np.loadtxt('lightcurve_triplet_sinchoque2.dat', comments='#')
xp_ss   = ATES[:, 0]
geom    = ATES[:, 1]
flux_ss = ATES[:, 5]

shock3  = np.loadtxt('lightcurve_triplet_model3_2.dat', comments='#')
xp_s3   = shock3[:, 0]
flux_s3 = shock3[:, 5]

# --- Máximas absorciones ---
abs_geom = (1 - geom.min())    * 100
abs_ss   = (1 - flux_ss.min()) * 100
abs_s3   = (1 - flux_s3.min()) * 100

# --- Figura principal ---
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(xp_ss, geom,    color='grey', lw=1.5,
        label=f'Geométrico  max={abs_geom:.5f} %')
ax.plot(xp_ss, flux_ss, color='C0',   lw=1.5,
        label=f'Sin Choque  max={abs_ss:.5f} %')
ax.plot(xp_s3, flux_s3, color='C1',   lw=1.5,
        label=f'Choque 3    max={abs_s3:.5f} %')

ax.set_xlabel(r'$x_p\ (R_p)$',         fontsize=12)
ax.set_ylabel('Flujo normalizado',      fontsize=12)
ax.set_title('Curva de luz — He I 10830 Å', fontsize=13)
ax.set_ylim(0.91, 0.92)
ax.legend(fontsize=9, loc='lower left',
          handlelength=2, framealpha=0.9,
          prop={'family': 'monospace'})
ax.grid(True, alpha=0.3)

# --- Inset: zoom del fondo del tránsito (Sin Choque vs Choque 3) ---
# Región de zoom: alrededor del mínimo de flux_s3
idx_min = np.argmin(flux_s3)
xc      = xp_s3[idx_min]
dx_zoom = 200   # semiancho en xp — ajustá si necesitás más/menos rango

mask_ss = (xp_ss >= xc - dx_zoom) & (xp_ss <= xc + dx_zoom)
mask_s3 = (xp_s3 >= xc - dx_zoom) & (xp_s3 <= xc + dx_zoom)

y_vals = np.concatenate([flux_ss[mask_ss], flux_s3[mask_s3]])
ylo = y_vals.min() - 0.00003
yhi = y_vals.max() + 0.00003

# Crear inset axes dentro del plot principal (esquina superior derecha)
ax_ins = inset_axes(ax, width='42%', height='45%', loc='upper right',
                    bbox_to_anchor=ax.bbox, bbox_transform=ax.transData,
                    borderpad=0)

ax_ins.plot(xp_ss[mask_ss], flux_ss[mask_ss], color='C0', lw=1.2)
ax_ins.plot(xp_s3[mask_s3], flux_s3[mask_s3], color='C1', lw=1.2)

ax_ins.set_xlim(xc - dx_zoom, xc + dx_zoom)
ax_ins.set_ylim(ylo, yhi)
ax_ins.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%.5f'))
ax_ins.tick_params(labelsize=7)
ax_ins.grid(True, alpha=0.3)
ax_ins.set_title('Zoom fondo tránsito', fontsize=7, pad=3)

# Líneas de conexión entre inset y región del plot principal
#mark_inset(ax, ax_ins, loc1=3, loc2=4, fc='none', ec='0.5', lw=0.8, ls='--')

#plt.tight_layout()
plt.savefig('comparacion2_CLs_He10830.png', dpi=150)
plt.show()