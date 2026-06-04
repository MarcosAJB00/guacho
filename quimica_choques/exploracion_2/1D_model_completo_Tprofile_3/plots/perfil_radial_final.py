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

output_file = "../output/perfiles_finales.dat"
output_data = np.column_stack([
    r_sorted,
    X_f_sorted[0, :],  # HI
    X_f_sorted[1, :],  # HII
    X_f_sorted[2, :],  # HeIS
    X_f_sorted[3, :],  # HeIM
    X_f_sorted[4, :],  # HeII
    X_f_sorted[5, :]   # e-
    ])

header = ("r(Rp)   HI   HII   HeIS   HeIM   HeII   e-")
np.savetxt(output_file, output_data, fmt="%.8e", header=header)
print(f"Archivo guardado: {output_file}")

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

# ---- perfiles iniciales (punteado) ----
for s in range(ns):
    sp = Species_names[s]
    plt.plot(
        r_sorted,
        X_i_sorted[s, :],
        '--',
        color=colors[sp],
        lw=1.5,
        label=f"{sp} initial"
    )

# ---- perfiles finales (sólido) ----
for s in range(ns):
    sp = Species_names[s]
    plt.plot(
        r_sorted,
        X_f_sorted[s, :],
        '-',
        color=colors[sp],
        lw=2,
        label=f"{sp} final"
    )

plt.xlabel("Radio (Rp)")
plt.ylabel("Densidad numerica")
plt.yscale("log")
plt.xlim(0.95, r_sorted[-1])
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.legend(ncol=2)
plt.tight_layout()

plt.savefig("./perfiles/inicial_final.png", dpi=300)
plt.show()

plt.figure(figsize=(8,6))
# ---- perfiles finales (sólido) ----
for s in range(ns):
    sp = Species_names[s]
    plt.plot(
        r_sorted,
        X_f_sorted[s, :],
        '-',
        color=colors[sp],
        lw=2,
        label=f"{sp} final"
    )

plt.xlabel("Radio (Rp)")
plt.ylabel("Densidad numerica")
plt.yscale("log")
#plt.xlim(0.95, 10.0)

plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.legend(ncol=2)
plt.tight_layout()
plt.xlim(0.95, r_sorted[-1])
plt.savefig("./perfiles/final.png", dpi=300)
plt.show()

plt.figure(figsize=(8,6))
h_total = X_f_sorted[0, :] + X_f_sorted[1, :]
plt.plot(r_sorted, X_f_sorted[0, :]/h_total, "-", color='royalblue', lw=2, label="HI fraction")
plt.plot(r_sorted, X_f_sorted[1, :]/h_total, "--", color='royalblue', lw=2, label="HII fraction")
plt.xlim(0.95, r_sorted[-1])
plt.ylim(-0.05, 1.05)
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("Radio (Rp)")
plt.ylabel("fracción de hidrógeno")
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("./perfiles/fraccion_h_final.png", dpi=300)
plt.show()

plt.figure(figsize=(8,6))
he_total = X_f_sorted[2, :] + X_f_sorted[3, :] + X_f_sorted[4, :]
plt.plot(r_sorted, X_f_sorted[0, :]/h_total, "-", color='royalblue', lw=2, label="HI fraction")
plt.plot(r_sorted, X_f_sorted[1, :]/h_total, "--", color='royalblue', lw=2, label="HII fraction")
plt.xlim(0.95, r_sorted[-1])
plt.ylim(-0.05, 1.05)
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("Radio (Rp)")
plt.ylabel("fracción de hidrógeno")
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("./perfiles/fraccion_h_final.png", dpi=300)
plt.show()

plt.figure(figsize=(8,6))
#dest_total_final= np.sum(X_f_sorted[:, :], axis=0)
#dest_total_ini = np.sum(X_i_sorted[:, :], axis=0)

r, rho = np.loadtxt("../inputs/initial_density_profile.dat", comments="#", unpack=True)
mp = 1.67e-24  # masa de un protón en gramos
me = 9.11e-28  # masa de un electrón en gramos
mhe = 4*mp     # masa de un átomo de helio en gramos
dest_total_final = (X_f_sorted[0, :]*mp + X_f_sorted[1, :]*mp + 
                   X_f_sorted[2, :]*mhe + X_f_sorted[3, :]*mhe + 
                   X_f_sorted[4, :]*mhe )#+ X_f_sorted[5, :]*me )
dest_total_ini = (X_i_sorted[0, :]*mp + X_i_sorted[1, :]*mp + 
                 X_i_sorted[2, :]*mhe + X_i_sorted[3, :]*mhe + 
                 X_i_sorted[4, :]*mhe )#+ X_i_sorted[5, :]*me )
#diff = np.abs(dest_total_final - dest_total_ini)
#plt.plot(r_sorted, diff, "-", color='royalblue', lw=2, label="dest total |final - initial|")
plt.plot(r_sorted, dest_total_final, "-", color='blue', lw=2, label="dest total final")
plt.plot(r_sorted, dest_total_ini, "--", color='red', lw=2, label="dest total initial")
plt.plot(r, rho, ":", color='gray', lw=2, label="densidad total input")
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("Radio (Rp)")
plt.ylabel("densidad total (g/cm$^3$)")
plt.yscale("log")
plt.xlim(0.95, r_sorted[-1])
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("./perfiles/dest_total.png", dpi=300)
plt.show()

plt.figure(figsize=(8,6))
dest_H_final= np.sum(X_f_sorted[:2, :], axis=0)
dest_H_ini = np.sum(X_i_sorted[:2, :], axis=0)
diff_H = np.abs(dest_H_final - dest_H_ini)
#plt.plot(r_sorted, diff_H, "-", color='royalblue', lw=2, label="dest H |final - initial|")
plt.plot(r_sorted, dest_H_final, "-", color='blue', lw=2, label="dest H final")
plt.plot(r_sorted, dest_H_ini, "--", color='red', lw=2, label="dest H initial")
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("Radio (Rp)")
plt.ylabel("densidad numerica total")
plt.xlim(0.95, r_sorted[-1])
plt.legend(ncol=2)
plt.yscale("log")
plt.tight_layout()
plt.savefig("./perfiles/dest_H.png", dpi=300)
plt.show()

plt.figure(figsize=(8,6))
dest_He_final= np.sum(X_f_sorted[2:5, :], axis=0) 
dest_He_ini = np.sum(X_i_sorted[2:5, :], axis=0) 
diff_He = np.abs(dest_He_final - dest_He_ini)
#plt.plot(r_sorted, diff_He, "-", color='royalblue', lw=2, label="dest He |final - initial|")
plt.plot(r_sorted, dest_He_final, "-", color='blue', lw=2, label="dest He final")
plt.plot(r_sorted, dest_He_ini, "--", color='red', lw=2, label="dest He initial")
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("Radio (Rp)")
plt.ylabel("densidad numerica total")
plt.xlim(0.95, r_sorted[-1])
plt.legend(ncol=2)
plt.yscale("log")
plt.tight_layout()
plt.savefig("./perfiles/dest_He.png", dpi=300)
plt.show()

plt.figure(figsize=(8,6))
nheim_f = X_f_sorted[3, :]
nheis_f = X_f_sorted[2, :]
nheim_i = X_i_sorted[3, :]
nheis_i = X_i_sorted[2, :]
plt.plot(r_sorted, nheim_f/nheis_f, "-", color='blue', lw=2, label="HeIM/HeIS final")
plt.plot(r_sorted, nheim_i/nheis_i, "--", color='red', lw=2, label="HeIM/HeIS initial")
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("Radio (Rp)")
plt.ylabel("HeIM/HeIS")
plt.xlim(0.95, r_sorted[-1])
plt.legend(ncol=2)
plt.yscale("log")
plt.tight_layout()
plt.savefig("./perfiles/alfa.png", dpi=300)
plt.show()