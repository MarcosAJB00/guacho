import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('profiles_Z_physical_HD63433b.csv', delimiter=',', skiprows=1)

r = data[:, 0] #en Rp 
density = data[:, 1] # en g/cm^3
pressure = data[:, 2] # en dyn/cm^2
temperature = data[:, 3] # en K
tau = data[:, 4]
velocity = data[:, 5] # en Km/s

mask = (r >= -5.0) & (r <= -0.95)
r = r[mask]
density = density[mask]
temperature = temperature[mask]
velocity = velocity[mask]

plt.figure(figsize=(8, 6))
plt.plot(r, density, label='Density Profile')
plt.xlabel('Radius (Rp)')
plt.ylabel('Density (g/cm^3)')
plt.title('Initial Profiles density')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_density_profile.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(r, temperature, label='Temperature Profile', color='red')
plt.xlabel('Radius (Rp)')
plt.ylabel('Temperature (K)')
plt.title('Initial Profiles temperature')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_temperature_profile.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(r, velocity, label='Velocity Profile', color='blue')
plt.xlabel('Radius (Rp)')
plt.ylabel('Velocity (Km/s)')
plt.title('Initial Profiles velocity')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_velocity_profile.png', dpi=300, bbox_inches='tight')
plt.show()