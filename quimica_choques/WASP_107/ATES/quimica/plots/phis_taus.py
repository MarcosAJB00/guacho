import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# =========================
# 1. Leer archivos cell_*.dat (phis)
# =========================

path = "../output/phis/"
files = sorted(glob.glob(os.path.join(path, "cell_*.dat")),
               key=lambda x: int(x.split("_")[-1].split(".")[0]))

nr = len(files)

# arrays
r = np.zeros(nr)

# phis
phi_initial = np.zeros((3, nr))
phi_final   = np.zeros((3, nr))

# taus
tau_initial = np.zeros((3, nr))
tau_final   = np.zeros((3, nr))

# =========================
# 2. Leer cada archivo
# =========================

for i, f in enumerate(files):
    data = np.loadtxt(f, comments='#')

    # radio
    r[i] = data[0, 1]

    # inicial
    phi_initial[:, i] = data[0, 2:5]
    tau_initial[:, i] = data[0, 5:8]

    # final
    phi_final[:, i] = data[-1, 2:5]
    tau_final[:, i] = data[-1, 5:8]

# =========================
# 3. Ordenar por radio
# =========================

sort_idx = np.argsort(r)

r_sorted = r[sort_idx]

phi_i = phi_initial[:, sort_idx]
phi_f = phi_final[:, sort_idx]

tau_i = tau_initial[:, sort_idx]
tau_f = tau_final[:, sort_idx]

# =========================
# 4. GRAFICO PHIS
# =========================

plt.figure(figsize=(8,6))

phi_names = ["phiHI", "phiHeIS", "phiHeIM"]

colors = {
    "phiHI":   "blue",
    "phiHeIS": "green",
    "phiHeIM": "red"
}

# inicial (punteado)
for s in range(3):
    sp = phi_names[s]
    plt.plot(
        r_sorted,
        phi_i[s, :],
        '--',
        color=colors[sp],
        lw=1.5,
        label=f"{sp} initial"
    )

# final (sólido)
for s in range(3):
    sp = phi_names[s]
    plt.plot(
        r_sorted,
        phi_f[s, :],
        '-',
        color=colors[sp],
        lw=2,
        label=f"{sp} final"
    )

plt.xlabel("Radio (Rp)")
plt.ylabel("Photoionization rate (phi)")
plt.yscale("log")
plt.xlim(0.95, r_sorted[-1])
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.legend()
plt.tight_layout()

plt.savefig("./perfiles/phis.png", dpi=300)
plt.show()

# =========================
# 5. GRAFICO TAUS
# =========================

plt.figure(figsize=(8,6))

tau_names = ["tauHI", "tauHeIS", "tauHeIM"]

colors = {
    "tauHI":   "blue",
    "tauHeIS": "green",
    "tauHeIM": "red"
}

# inicial (punteado)
for s in range(3):
    sp = tau_names[s]
    plt.plot(
        r_sorted,
        tau_i[s, :],
        '--',
        color=colors[sp],
        lw=1.5,
        label=f"{sp} initial"
    )

# final (sólido)
for s in range(3):
    sp = tau_names[s]
    plt.plot(
        r_sorted,
        tau_f[s, :],
        '-',
        color=colors[sp],
        lw=2,
        label=f"{sp} final"
    )

plt.xlabel("Radio (Rp)")
plt.ylabel("Profundidad óptica")
plt.xlim(0.95, r_sorted[-1])
plt.yscale("log")
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig("./perfiles/taus.png", dpi=300)
plt.show()


output_file = "../output/phi_tau_finales.dat"
output_data = np.column_stack([
    r_sorted,
    phi_f[0, :],  # phiHI
    phi_f[1, :],  # phiHeIS
    phi_f[2, :],  # phiHeIM
    tau_f[0, :],  # tauHI
    tau_f[1, :],  # tauHeIS
    tau_f[2, :]   # tauHeIM
    ])

header = ("r(Rp)   phiHI   phiHeIS   phiHeIM   tauHI   tauHeIS   tauHeIM")
np.savetxt(output_file, output_data, fmt="%.8e", header=header)
print(f"Archivo guardado: {output_file}")
