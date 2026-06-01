import numpy as np
import matplotlib.pyplot as plt
import os
import glob

# =========================
# Paths
# =========================

path = "../output/"

# Buscar todos los modelos existentes
files = sorted(glob.glob(path + "final_model_*.dat"))

print(f"Se encontraron {len(files)} modelos")

# =========================
# Config
# =========================

ns = 6

Species_names = [
    "HI",
    "HII",
    "HeIS",
    "HeIM",
    "HeII",
    "e-"
]

phi_tau_names = [
    "HI",
    "HeIS",
    "HeIM"
]

# Crear carpetas
os.makedirs("./perfiles/species", exist_ok=True)
os.makedirs("./perfiles/phis", exist_ok=True)
os.makedirs("./perfiles/taus", exist_ok=True)

# =========================
# Loop sobre archivos
# =========================

for filename in files:

    model_id = os.path.basename(filename)
    model_id = model_id.replace("final_model_", "")
    model_id = model_id.replace(".dat", "")

    print(f"Procesando archivo: {filename}")

    data = np.loadtxt(filename, comments="#")

    # Verificación mínima
    if data.ndim == 1:
        print(f"Archivo {filename} vacío o con una sola fila. Saltando.")
        continue

    # =====================
    # Columnas
    # =====================

    r = data[:, 0]

    X = data[:, 1:1+ns]

    phis = data[:, 1+ns:1+ns+3]

    taus = data[:, 1+ns+3:1+ns+6]

    # =====================
    # 1) Species
    # =====================

    fig, ax = plt.subplots(figsize=(8,5))

    for j in range(ns):
        ax.plot(r, X[:, j], label=Species_names[j])

    ax.set_xlabel("Distancia [cm]")
    ax.set_ylabel("Densidad numérica")
    ax.set_yscale("log")

    ax.set_title(f"Modelo {model_id}")

    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        f"./perfiles/species/model_{model_id}_species.png",
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

    ax.set_title(f"Modelo {model_id}")

    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        f"./perfiles/phis/model_{model_id}_phis.png",
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

    ax.set_title(f"Modelo {model_id}")

    ax.grid(True)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        f"./perfiles/taus/model_{model_id}_taus.png",
        dpi=300
    )

    plt.close(fig)

print("Listo.")