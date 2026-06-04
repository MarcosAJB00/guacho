import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('model_3.dat', skiprows=1)

r = data[:, 0]
rho = data[:, 1]
T = data[:, 2]


plt.figure(figsize=(8, 6))
plt.plot(r, rho, label='Density Profile')
plt.xlabel('Radius (cm)')
plt.ylabel('Density (g/cm^3)')
plt.title('Initial Profiles density')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_density_profile.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(r, T, label='Temperature Profile', color='red')
plt.xlabel('Radius (cm)')
plt.ylabel('Temperature (K)')
plt.title('Initial Profiles temperature')
plt.legend()
#plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_temperature_profile.png', dpi=300, bbox_inches='tight')
plt.show()

np.savetxt('initial_density_profile.dat', 
            np.column_stack((r, rho)),
            header='Radius (cm) Density (g/cm^3)'
            , fmt='%.6e'
            )
np.savetxt('initial_temperature_profile.dat',
            np.column_stack((r, T)),
            header='Radius (cm) Temperature (K)',
            fmt='%.6e'
            )