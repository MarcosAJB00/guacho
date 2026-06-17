import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.ticker import MultipleLocator

# =========================
# Paths
# =========================

path = "../output/"

files = sorted(glob.glob(path + "final_model_*.dat"))

print(f"Se encontraron {len(files)} modelos")

# =========================
# Config
# =========================

HeIM_index = 4   # r HI HII HeIS HeIM HeII e-

model_ids = []
HeIM_dens_tot = []

# =========================
# Leer modelos existentes
# =========================

for filename in files:

    model_id = os.path.basename(filename)
    model_id = model_id.replace("final_model_", "")
    model_id = model_id.replace(".dat", "")

    print(f"Procesando archivo: {filename}")

    try:
        data = np.loadtxt(filename, comments="#")

        if data.ndim == 1:
            print(f"Archivo {filename} inválido. Saltando.")
            continue

        HeIM = data[:, HeIM_index]

        model_ids.append(int(model_id))
        #HeIM_dens_tot.append(np.sum(HeIM))
        HeIM_dens_tot.append(np.trapezoid(HeIM))

    except Exception as e:
        print(f"Error leyendo {filename}: {e}")

model_ids = np.array(model_ids)
HeIM_dens_tot = np.array(HeIM_dens_tot)

# =========================
# Graficar
# =========================

dens_tot_norm = HeIM_dens_tot / np.sum(HeIM_dens_tot)

plt.figure(figsize=(8,6))

# Índice del mejor modelo
index_max = np.argmax(HeIM_dens_tot)

# Todos los modelos
plt.scatter(
    model_ids,
    dens_tot_norm,
    s=35
)

# Líneas verticales
plt.vlines(
    model_ids,
    0,
    dens_tot_norm,
    linestyles='dashed',
    alpha=0.5,
    linewidth=0.8
)

# Resaltar mejor modelo
plt.scatter(
    model_ids[index_max],
    dens_tot_norm[index_max],
    s=80,
    marker='*',
    color='red',
    zorder=10
)

plt.xlabel("Modelo")
plt.ylabel("HeIM total normalizado")

plt.ylim(0, np.max(dens_tot_norm)*1.1)

ax = plt.gca()

# Ticks mayores cada 10 modelos
ax.xaxis.set_major_locator(MultipleLocator(10))

# Ticks menores cada modelo
ax.xaxis.set_minor_locator(MultipleLocator(1))

# Grid
ax.grid(True, which='major', alpha=0.5)
ax.grid(True, which='minor', alpha=0.2)

# Longitud de ticks
ax.tick_params(axis='x', which='major', length=8)
ax.tick_params(axis='x', which='minor', length=4)

plt.title(
    f"Best model = {model_ids[index_max]}"
)

plt.savefig(
    "dens_HeIM.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print(
    f"Mejor modelo: {model_ids[index_max]}"
)