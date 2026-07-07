import numpy as np
import matplotlib.pyplot as plt
import os
import re

# ------------------------------------------------------------
# Modelo sin choque
# ------------------------------------------------------------
# Cambia este nombre si tu archivo se llama distinto
# --- Datos ---
ATES    = np.loadtxt('lightcurve_triplet_sinchoque2.dat', comments='#')
xp_ss   = ATES[:, 0]
geom    = ATES[:, 1]
flux_ss = ATES[:, 5]

# --- Maximas absorciones ---
abs_geom = (1 - geom.min())    * 100
abs_ss   = (1 - flux_ss.min()) * 100

# ------------------------------------------------------------
# Crear figura
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 8))

# Modelo geometrico
ax.plot(
    xp_ss,
    geom,
    color='grey',
    ls='--',
    lw=3,
    label=f'Geom max={abs_geom:.5f} %'
)
# Modelo sin choque
ax.plot(
    xp_ss,
    flux_ss,
    color='k',
    lw=3,
    label=f'Sin choque max={abs_ss:.5f} %'
)

# ------------------------------------------------------------
# Modelos con choque
# ------------------------------------------------------------
base_path = "/home/mbaracchi/guacho/quimica_choques/HP11_shock3D_model3"

folders = [
    "phi0_alpha45",
    "phi0_alpha60",
    "phi0_alpha75",
    "phi0_alpha90",
    "phi45_alpha45",
    "phi45_alpha60",
    "phi45_alpha75",
    "phi45_alpha90",
    "phi77_alpha60",
    "phi77_alpha75",
    "phi77_alpha90",
    "phi90_alpha45",
    "phi90_alpha60",
    "phi90_alpha75",
    "phi90_alpha90",
    "shock_fiducial",
]
#folders = sorted([d for d in os.listdir(".") if os.path.isdir(d) and d.startswith("phi")])

# agregar fiducial al final
#folders.append("shock_fiducial")

for folder in folders:

    filename = os.path.join(base_path, folder, "lightcurve_triplet_model3.dat")

    if not os.path.exists(filename):
        print(f"No existe {filename}")
        continue

    data = np.loadtxt(filename, comments="#")

    phase = data[:,0]
    flux  = data[:,5]
    abs_   = (1 - flux.min()) * 100

    if folder == "shock_fiducial":

        label = r"$\phi=77^\circ,\ \alpha=45^\circ$ (fiducial)"

    else:

        m = re.match(r"phi(\d+)_alpha(\d+)", folder)

        if m is None:
            label = folder
        else:
            phi = m.group(1)
            alpha = m.group(2)
            label = rf"$\phi={phi}^\circ,\ \alpha={alpha}^\circ$, max={abs_:.5f} %"

    ax.plot(phase, flux, lw=1.8, label=label)

# ------------------------------------------------------------
# Estetica
# ------------------------------------------------------------
ax.set_xlabel("Orbital distance", fontsize=14)
ax.set_ylabel("Absorcion", fontsize=14)
print(phase.min(), phase.max())
print(xp_ss.min(), xp_ss.max())
#ax.set_xlim = (-30.0,-15.0)

ax.grid(alpha=0.3)

ax.legend(
    loc='upper center',
    bbox_to_anchor=(0.5,-0.18),
    ncol=4,
    fontsize=9
)


plt.tight_layout()

plt.savefig("comparison_lightcurves_model3_pretans.png", dpi=300, bbox_inches='tight')
plt.show()