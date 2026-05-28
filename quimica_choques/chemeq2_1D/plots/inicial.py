import numpy as np
import matplotlib.pyplot as plt
import glob
import os

init_data = np.loadtxt("../inputs/density_profile_Roman.dat", comments="#")

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

plt.figure(figsize=(8,6))
mp = 1.67e-24  # masa de un protón en gramos
mhe = 4*mp     # masa de un átomo de helio en gramos
me = 9.11e-28  # masa de un electrón en gramos

rho_total_ini = rho_init
density_total_ini = (mp*(HI_test + HII_test) + 
                    mhe*(HeIS_test + HeIM_test + HeII_test) 
                    + me*e_test)

plt.plot(r_init, density_total_ini, "-", color='blue', lw=2, label="densidad total inicial code")
plt.plot(r_init, rho_total_ini, ":", color='red', lw=2, label="densidad total input")
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("Radio (Rp)")
plt.ylabel("densidad total (g/cm$^3$)")
plt.yscale("log")
plt.xlim(0.95, r_init[-1])
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("./perfiles/rho_input_vs_code.png", dpi=300)