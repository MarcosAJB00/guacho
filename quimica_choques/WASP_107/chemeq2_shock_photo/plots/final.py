import numpy as np
import matplotlib.pyplot as plt
import os

# =========================
# Paths
# =========================

path = "../output/"
models = np.loadtxt(
    "../../shock_1D/output/model_list.dat", comments="#"
)

if len(models.shape) == 1:
    model = np.array([models[0]])
    T0 = np.array([models[1]])
    n0 = np.array([models[2]])
    u0 = np.array([models[3]])
    y0 = np.array([models[4]])
else:
    model = models[:, 0]
    T0 = models[:, 1]
    n0 = models[:, 2]
    u0 = models[:, 3]  
    y0 = models[:, 4]

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

    filename = f"{path}final_model_{int(model[i])}.dat"

    print(f"Procesando archivo: {filename}")

    data = np.loadtxt(filename, comments="#")

    # Columnas
    Rjup = 7.1492e9  # cm
    r = data[:, 0]/Rjup
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

    ax.set_xlabel("Distancia [Rjup]")
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

    ax.set_xlabel("Distancia [Rjup]")
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

    ax.set_xlabel("Distancia [Rjup]")
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

    # =====================
    # 4) fHe
    # =====================

    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8,5))
   
    fHeIM_IS = X[:, 3] / X[:, 2]
    fHeIM_tot = X[:, 3] / (X[:, 2] + X[:, 3] + X[:, 4])
    fHeIS_tot = X[:, 2] / (X[:, 2] + X[:, 3] + X[:, 4])
    fHeII_tot =  X[:, 4] / (X[:, 2] + X[:, 3] + X[:, 4])
    fHeI_tot = (X[:, 2] + X[:, 3]) / (X[:, 2] + X[:, 3] + X[:, 4])
         
    ax.plot(r, fHeIM_IS, label="fHeIM/IS")
    ax.plot(r, fHeIM_tot, label="fHeIM/He_tot")
    ax.plot(r, fHeIS_tot, label="fHeIS/He_tot")
    ax.plot(r, fHeII_tot, label="fHeII/He_tot")
    ax.plot(r, fHeI_tot, label="fHeI/He_tot", linestyle='--')

    ax.set_xlabel("Distancia [Rjup]")
    ax.set_ylabel("Fracción de Helio")
    ax.set_yscale("log")
    ax.set_title(f"Modelo {int(model[i])}: T0={T0[i]:.1e} K, n0={n0[i]:.1e} cm^-3, u0={u0[i]:.1e} km/s, y0={y0[i]:.1e}")
    ax.grid(True)

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        f"./perfiles/frac_He/model_{int(model[i])}_fHe.png",
        dpi=300
    )

    plt.close(fig)