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
plt.savefig("dens_HeIM.png", dpi=300, bbox_inches='tight')
plt.show()