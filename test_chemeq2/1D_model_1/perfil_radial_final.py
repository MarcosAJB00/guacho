import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# =========================
# 1. Leer radios desde .log
# =========================

log_file = "./output/status.log"  

# salteamos header
log_data = np.loadtxt(log_file, skiprows=8)

# ahora los radios están en la PRIMERA COLUMNA
r = log_data[:, 0]

# =========================
# 2. Leer archivos salida_*.txt
# =========================

path = "./output/"
files = sorted(glob.glob(os.path.join(path, "salida_*.txt")),
               key=lambda x: int(x.split("_")[-1].split(".")[0]))

nr = len(files)
ns = 6  # número de especies

# array: especies x radio
X_final = np.zeros((ns, nr))

# =========================
# 3. Extraer última línea de cada archivo
# =========================

for i, f in enumerate(files):
    data = np.loadtxt(f)
    
    last_row = data[-1]  # última fila
    
    abundances = last_row[1:]  # ignoramos tiempo
    
    X_final[:, i] = abundances

# =========================
# 4. Chequeo de consistencia
# =========================

if len(r) != nr:
    raise ValueError(f"Cantidad de radios ({len(r)}) != cantidad de archivos ({nr})")

# =========================
# 5. Ordenar por radio creciente 
# =========================

sort_idx = np.argsort(r)

r_sorted = r[sort_idx]
X_sorted = X_final[:, sort_idx]

# =========================
# 6. Graficar
# =========================

plt.figure(figsize=(8,6))
Species_names = ["HI", "HII", "HeIS", "HeIM", "HeII", "e-"]
for s in range(ns):
    plt.plot(r_sorted, X_sorted[s, :], label=f"{Species_names[s]}")

plt.xlabel("Radio")
plt.ylabel("Abundancia final")
plt.yscale("log")  # clave
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig("perfil_radial_final.png", dpi=300)
plt.show()