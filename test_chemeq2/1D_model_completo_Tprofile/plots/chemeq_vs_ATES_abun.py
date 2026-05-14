import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# =========================
# 1. Leer archivos cell_*.dat
# =========================

path = "../output/species/"
files = sorted(glob.glob(os.path.join(path, "cell_*.dat")),
               key=lambda x: int(x.split("_")[-1].split(".")[0]))

nr = len(files)
ns = 6  # número de especies

# arrays
r = np.zeros(nr)

X_initial = np.zeros((ns, nr))
X_final   = np.zeros((ns, nr))

# =========================
# 2. Leer cada archivo
# =========================

for i, f in enumerate(files):
    data = np.loadtxt(f, comments='#')

    # radio (columna 2)
    r[i] = data[0, 1]

    # fila inicial
    X_initial[:, i] = data[0, 2:]

    # fila final
    X_final[:, i] = data[-1, 2:]

# =========================
# 3. Ordenar por radio creciente
# =========================

sort_idx = np.argsort(r)

r_sorted = r[sort_idx]
X_i_sorted = X_initial[:, sort_idx]
X_f_sorted = X_final[:, sort_idx]

# =============================
# Leer archivo de abundancias de ATES
# =============================

file = 'Ion_species_adv.txt'

r, nhi, nhii, nhei, nheii, nheiii, nheiTR = np.loadtxt(file, unpack=True)

nh  = nhi[:] + nhii[:]      # Total hydrogen density
fhi    = nhi[:]/nh[:]       # HI fraction
fhii   = nhii[:]/nh[:]      # HII fraction

# ===================================
# Perfil de Termepratura de ATES
# ===================================
data_T = np.loadtxt('../inputs/temperature_profile_uniform.dat', skiprows=1)

r_T = data_T[:, 0]
T = data_T[:, 1]
# =========================
# 4. Graficar
# =========================

plt.figure(figsize=(8,6))

Species_names = ["HI" , "HII", "HeIS", "HeIM", "HeII", "e-"]

colors = {
    "HI":   "blue",
    "HII":  "red",
    "HeIS": "green",
    "HeIM": "orange",
    "HeII": "cyan",
    "e-":   "black"
}

ax1 = plt.gca()

# ---- perfiles finales ----
for s in range(ns):
    sp = Species_names[s]
    ax1.plot(
        r_sorted,
        X_f_sorted[s, :],
        '-',
        color=colors[sp],
        lw=2,
        label=f"{sp} final"
    )

# ATES
ax1.plot(r, nhi, '--', color=colors["HI"], label=r'$n_{HI}$')
ax1.plot(r, nhii, '--', color=colors["HII"], label=r'$n_{HII}$')
ax1.plot(r, nhei, '--', color=colors["HeIS"], label=r'$n_{HeI}$')
ax1.plot(r, nheii, '--', color=colors["HeII"], label=r'$n_{HeII}$')
ax1.plot(r, nhii + nheii, '--', color=colors["e-"], label=r'$n_e$')
ax1.plot(r, nheiTR, '--', color=colors["HeIM"], label=r'$n_{HeIM}$')

# ---- temperatura ----
ax2 = ax1.twinx()

ax2.plot(r_T, T, color='k', lw=2, ls=':', label='T')

ax2.set_ylabel('Temperatura (K)', color='k')
ax2.tick_params(axis='y', labelcolor='k')

# formato
ax1.set_xlabel("Radio (Rp)")
ax1.set_ylabel("Densidad numerica")
ax1.set_yscale("log")
ax1.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
ax1.set_xlim(0.95, r_sorted[-1])
ax1.set_ylim(1e-1, 1e10)

# leyenda conjunta
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, ncol=2)

plt.tight_layout()
plt.savefig("./perfiles/chemeq_vs_ATES_abun_adv.png", dpi=300)
plt.show()

