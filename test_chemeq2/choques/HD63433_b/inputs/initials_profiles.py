import numpy as np
import matplotlib.pyplot as plt

data_rho = np.loadtxt('density_profile_ATES.dat', skiprows=1)

r_rho = data_rho[:, 0]
rho = data_rho[:, 1]

data_T = np.loadtxt('temperature_profile_uniform.dat', skiprows=1)

r_T = data_T[:, 0]
T = data_T[:, 1]

plt.figure(figsize=(8, 6))
plt.plot(r_rho, rho, label='Density Profile')
plt.xlabel('Radius (Rp)')
plt.ylabel('Density (g/cm^3)')
plt.title('Initial Profiles density')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_density_profile.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(r_T, T, label='Temperature Profile', color='red')
plt.xlabel('Radius (Rp)')
plt.ylabel('Temperature (K)')
plt.title('Initial Profiles temperature')
plt.legend()
#plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_temperature_profile.png', dpi=300, bbox_inches='tight')
plt.show()