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

model = models[:, 0]
T0 = models[:, 1]
n0 = models[:, 2]
u0 = models[:, 3]  
y0 = models[:, 4]
models_num = len(model)

HeIM_index = 4 # r HI HII HeIS HeIM HeII e- 
HeIM_dens_tot = []

for i in range(models_num):

    filename = f"{path}final_model_{int(model[i])}.dat"
    print(f"Procesando archivo: {filename}")
    data = np.loadtxt(filename, comments="#")
    HeIM_dens_tot.append(np.sum(data[:, HeIM_index]))

plt.figure(figsize=(8, 6))
plt.plot(model, HeIM_dens_tot, marker='o')
plt.xlabel('Modelo')
plt.ylabel('Densidad total de HeIM')
plt.title('Densidad total de HeIM para cada modelo')
plt.grid()
index_max = np.argmax(HeIM_dens_tot)
T0_max = T0[index_max]
n0_max = n0[index_max]
u0_max = u0[index_max]/1e5 # Convertir a km/s
y0_max = y0[index_max]
plt.title(f'Best: Model {int(model[index_max])}, T0={T0_max:.1e}, n0={n0_max:.1e}, u0={u0_max:.1e}, y0={y0_max:.1e}')
plt.savefig("dens_HeIM.png", dpi=300, bbox_inches='tight')
plt.show()

