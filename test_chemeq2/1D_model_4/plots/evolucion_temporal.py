import numpy as np
import matplotlib.pyplot as plt
import sys

# ==================================================
# Leer archivo: tiempo + 6 especies químicas
# ==================================================

# nombre del archivo
n = int(sys.argv[1])
salida = n
archivo_sp = f"../output/species/cell_{salida}.dat"  # modificar según tu caso
archivo_ph = f"../output/phis/cell_{salida}.dat"  # modificar según tu caso

# cargar datos
data_sp = np.loadtxt(archivo_sp, comments="#")  # omitir líneas de comentario
data_ph = np.loadtxt(archivo_ph, comments="#")  # omitir líneas de comentario

# columnas
t = data_sp[:, 0]          # tiempo
r = data_sp[:, 1]          # posición (si es relevante)
sp1 = data_sp[:, 2]
sp2 = data_sp[:, 3]
sp3 = data_sp[:, 4]
sp4 = data_sp[:, 5]
sp5 = data_sp[:, 6]
sp6 = data_sp[:, 7]

t_phi = data_ph[:, 0]          # tiempo
phi_HI = data_ph[:, 2]
phi_HeIS = data_ph[:, 3]
phi_HeIM = data_ph[:, 4]
tau_HI = data_ph[:, 5]
tau_HeIS = data_ph[:, 6]
tau_HeIM = data_ph[:, 7]

# ==================================================
# Gráfico de evolución temporal de especies químicas
# ==================================================
labels = [
    "HI",
    "HII",
    "HeI",
    "He(2³S)",
    "HeII",
    "e⁻"
]

plt.figure(figsize=(10,6))

plt.plot(t, sp1, label=labels[0])
plt.plot(t, sp2, label=labels[1])
plt.plot(t, sp3, label=labels[2])
plt.plot(t, sp4, label=labels[3])
plt.plot(t, sp5, label=labels[4])
plt.plot(t, sp6, label=labels[5])

plt.xlabel("Tiempo [s]")
plt.ylabel("Densidad numérica [cm⁻³]")
plt.title(f"Evolución temporal de especies químicas en R={r[0]:.2e} Rp")
plt.yscale("log")     # opcional
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"./evolucion/species_{salida}.png", dpi=300)  # guardar figura
plt.show()


# ======================================
# Gráfico de evolución temporal de phis 
# ======================================
plt.figure(figsize=(10,6))
plt.plot(t_phi, phi_HI, label='phi_HI')
plt.plot(t_phi, phi_HeIS, label='phi_HeIS')
plt.plot(t_phi, phi_HeIM, label='phi_HeIM')
plt.xlabel("Tiempo [s]")
plt.ylabel("Seccion eficaz de fotoionizacion [cm²]")
plt.title(f"Evolución temporal de phis en R={r[0]:.2e} Rp")
plt.yscale("log")     # opcional
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"./evolucion/phis_{salida}.png", dpi=300)  # guardar figura
plt.show()

plt.figure(figsize=(10,6))
plt.plot(t_phi, tau_HI, label='tau_HI')
plt.plot(t_phi, tau_HeIS, label='tau_HeIS')
plt.plot(t_phi, tau_HeIM, label='tau_HeIM')
plt.xlabel("Tiempo [s]")
plt.ylabel("Profundidad óptica")
plt.title(f"Evolución temporal de taus en R={r[0]:.2e} Rp")
plt.yscale("log")     # opcional
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"./evolucion/taus_{salida}.png", dpi=300)  # guardar figura
plt.show()