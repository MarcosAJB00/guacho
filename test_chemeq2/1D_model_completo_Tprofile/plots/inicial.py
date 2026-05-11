import numpy as np
import matplotlib.pyplot as plt
import glob
import os

init_data = np.loadtxt("../inputs/neutral_density_profile.dat", comments="#")

r_init = init_data[:, 0]
rho_init = init_data[:, 1]

init_data_test = np.loadtxt("../output/perfiles_iniciales.txt", comments="#")
HI_test = init_data_test[:, 1]
HII_test = init_data_test[:, 2]
HeIS_test = init_data_test[:, 3]
HeIM_test = init_data_test[:, 4]
HeII_test = init_data_test[:, 5]
e_test = init_data_test[:, 6]

#total = HI_test + HeIS_test + HeIM_test + HII_test + HeII_test
mp = 1.67e-24  # masa de un protón en gramos

plt.figure(figsize=(9,6))
plt.plot(r_init, rho_init/mp, label="Densidad total")
plt.plot(r_init, HI_test, label="HI")
plt.plot(r_init, HII_test, label="HII")
plt.plot(r_init, HeIS_test, label="HeIS")
plt.plot(r_init, HeIM_test, label="HeIM")
plt.plot(r_init, HeII_test, label="HeII")
plt.plot(r_init, e_test, label="e-")

plt.xlabel("Radio (Rp)")
plt.ylabel("Densidad numérica [cm$^{-3}$]")
plt.yscale("log")     # opcional
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("./perfiles/perfil_inicial.png", dpi=300)
plt.show()