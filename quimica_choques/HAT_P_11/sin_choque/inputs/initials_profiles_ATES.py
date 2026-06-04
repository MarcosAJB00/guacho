import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator


r,rho,v,p,T,heat,cool = \
      np.loadtxt('Hydro_ioniz.txt',unpack = True)

mp = 1.6726219e-24 # g
rho = rho*mp
v = v/1e5 # cm/s a km/s
mask = (r >= 1.001) & (r <= 9.0)
r = r[mask]
density = rho[mask]
temperature = T[mask]
velocity = v[mask]

# Ordenar por radio creciente
sort_idx = np.argsort(r)

r = r[sort_idx]
density = density[sort_idx]
temperature = temperature[sort_idx]
velocity = velocity[sort_idx]

# Grid uniforme nuevo
N_new = 1000
r_uniform = np.linspace(r.min(), r.max(), N_new)

# Interpoladores
dens_interp = PchipInterpolator(r, np.log10(density))
temp_interp = PchipInterpolator(r, temperature)
vel_interp = PchipInterpolator(r, velocity)

# Evaluar
density_uniform = 10**dens_interp(r_uniform)
temperature_uniform = temp_interp(r_uniform)
velocity_uniform = vel_interp(r_uniform)

np.savetxt('density_profile_ATES.dat',
    np.column_stack((r_uniform, density_uniform)),
    header='r[R_p] density[g/cm^3]',
    fmt='%.10e'
)

np.savetxt('temperature_profile_ATES.dat',
    np.column_stack((r_uniform, temperature_uniform)),
    header='r[R_p] temperature[K]',
    fmt='%.10e'
)

np.savetxt('velocity_profile_ATES.dat',
    np.column_stack((r_uniform, velocity_uniform)),
    header='r[R_p] velocity[Km/s]',
    fmt='%.10e'
)

plt.figure(figsize=(8, 6))
plt.plot(r, density, label='Density Profile')
plt.plot(r_uniform, density_uniform,color='green', label='Interpolated Density Profile', linestyle='--')
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
plt.plot(r_uniform, temperature_uniform,color='orange', label='Interpolated Temperature Profile', linestyle='--')
plt.xlabel('Radius (Rp)')
plt.ylabel('Temperature (K)')
plt.title('Initial Profiles temperature')
plt.legend()
#plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_temperature_profile.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(r, velocity, label='Velocity Profile', color='blue')
plt.plot(r_uniform, velocity_uniform,color='purple', label='Interpolated Velocity Profile', linestyle='--')
plt.xlabel('Radius (Rp)')
plt.ylabel('Velocity (Km/s)')
plt.title('Initial Profiles velocity')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('initial_velocity_profile.png', dpi=300, bbox_inches='tight')
plt.show()