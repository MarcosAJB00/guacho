import numpy as np
import matplotlib.pyplot as plt
import os

# =========================
# Paths
# =========================

path = "../chemeq2_1D/output/"
models_0 = np.loadtxt(
    "../shock_1D/output/model_list.dat", comments="#"
)

path_1 = "../chemeq2_1D_1/output/"

models_1 = np.loadtxt(
    "../shock_1D/output_1/model_list.dat", comments="#"
)

model = np.concatenate([models_0[:, 0], models_1[:, 0]])
T0 = np.concatenate([models_0[:, 1], models_1[:, 1]])
n0 = np.concatenate([models_0[:, 2], models_1[:, 2]])
u0 = np.concatenate([models_0[:, 3], models_1[:, 3]])
y0 = np.concatenate([models_0[:, 4], models_1[:, 4]])
models_num = len(model)

# =========================
# Config
# =========================

ns = 6

Species_names = [
    "HI", "HII",
    "HeIS", "HeIM",
    "HeII", "e-"
]

phi_tau_names = [
    "HI", "HeIS", "HeIM"
]

# Crear carpetas si no existen
os.makedirs("./perfiles/species", exist_ok=True)
os.makedirs("./perfiles/phis", exist_ok=True)
os.makedirs("./perfiles/taus", exist_ok=True)

# =========================
# Loop sobre modelos
# =========================

for i in range(models_num):

    if i < len(models_0):
        filename = f"{path}final_model_{int(model[i])}.dat"
    else:
        filename = f"{path_1}final_model_{int(model[i])}.dat"

    print(f"Procesando archivo: {filename}")

    data = np.loadtxt(filename, comments="#")

    # Columnas
    r = data[:, 0]
    X = data[:, 1:1+ns]
    phis = data[:, 1+ns:1+ns+3]
    taus = data[:, 1+ns+3:1+ns+6]

    u0[i] = u0[i] / 1.0e5  # Convertir de cm/s a km/s

    # =====================
    # 1) Species
    # =====================

    fig, ax = plt.subplots(figsize=(8,5))

    for j in range(ns):
        ax.plot(r, X[:, j], label=Species_names[j])

    ax.set_xlabel("Distancia [cm]")
    ax.set_ylabel("Densidad numerica")
    ax.set_yscale("log")
    ax.set_title(f"Modelo {int(model[i])}: T0={T0[i]:.1e} K, n0={n0[i]:.1e} cm^-3, u0={u0[i]:.1e} km/s, y0={y0[i]:.1e}")
    ax.grid(True)

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        f"./perfiles/species/model_{int(model[i])}_species.png",
        dpi=300
    )

    plt.close(fig)

    # =====================
    # 2) Phis
    # =====================

    fig, ax = plt.subplots(figsize=(8,5))

    for j in range(3):
        ax.plot(r, phis[:, j], label=phi_tau_names[j])

    ax.set_xlabel("Distancia [cm]")
    ax.set_ylabel("Phi")
    ax.set_yscale("log")
    ax.set_title(f"Modelo {int(model[i])}: T0={T0[i]:.1e} K, n0={n0[i]:.1e} cm^-3, u0={u0[i]:.1e} km/s, y0={y0[i]:.1e}")

    ax.grid(True)

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        f"./perfiles/phis/model_{int(model[i])}_phis.png",
        dpi=300
    )

    plt.close(fig)

    # =====================
    # 3) Taus
    # =====================

    fig, ax = plt.subplots(figsize=(8,5))

    for j in range(3):
        ax.plot(r, taus[:, j], label=phi_tau_names[j])

    ax.set_xlabel("Distancia [cm]")
    ax.set_ylabel("Tau")
    ax.set_yscale("log")
    ax.set_title(f"Modelo {int(model[i])}: T0={T0[i]:.1e} K, n0={n0[i]:.1e} cm^-3, u0={u0[i]:.1e} km/s, y0={y0[i]:.1e}")
    ax.grid(True)

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        f"./perfiles/taus/model_{int(model[i])}_taus.png",
        dpi=300
    )

    plt.close(fig)