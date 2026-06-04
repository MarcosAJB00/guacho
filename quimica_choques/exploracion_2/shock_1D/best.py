import numpy as np 
import matplotlib.pyplot as plt

data = np.loadtxt("./output/best_model.dat", comments="#")

Rjup = 7.1492e9  # cm
x = data[:, 0]/Rjup
rho = data[:, 1]
T = data[:, 2]
n = data[:, 3]
u = data[:, 4]
y = data[:, 5]
P = data[:, 6]
h = data[:, 7]
L = data[:, 8]

model_name = "dx1e8"

plt.figure(figsize=(8,6))
plt.plot(x, rho, label='Density [g/cm^3]')
plt.xlabel('x [Rjup]')
plt.ylabel('Density [g/cm^3]')
plt.yscale('log')
plt.title('Density profile for the best model')
plt.grid()
plt.savefig(f"{model_name}_density_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, T, label='Temperature [K]')
plt.xlabel('x [Rjup]')
plt.ylabel('Temperature [K]')
plt.yscale('log')
plt.title('Temperature profile for the best model')
plt.grid()
plt.savefig(f"{model_name}_temperature_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, n, label='Number density [cm^-3]')
plt.xlabel('x [Rjup]')
plt.ylabel('Number density [cm^-3]')
plt.yscale('log')
plt.title('Number density profile for the best model')
plt.grid()
plt.savefig(f"{model_name}_number_density_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, u, label='Velocity [cm/s]')
plt.xlabel('x [Rjup]')
plt.ylabel('Velocity [cm/s]')
plt.yscale('log')
plt.title('Velocity profile for the best model')
plt.grid()
plt.savefig(f"{model_name}_velocity_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, y, label='Ionization fraction')
plt.xlabel('x [Rjup]')
plt.ylabel('Ionization fraction')
plt.yscale('log')
plt.title('Ionization fraction profile for the best model')
plt.grid()
plt.savefig(f"{model_name}_ionization_fraction_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, P, label='Pressure [dyn/cm^2]') 
plt.xlabel('x [Rjup]')
plt.ylabel('Pressure [dyn/cm^2]')
plt.yscale('log')
plt.title('Pressure profile for the best model')
plt.grid()
plt.savefig(f"{model_name}_pressure_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, h, label='Enthalpy [erg/g]')
plt.xlabel('x [Rjup]')
plt.ylabel('Enthalpy [erg/g]')
plt.yscale('log')
plt.title('Enthalpy profile for the best model')
plt.grid()
plt.savefig(f"{model_name}_enthalpy_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
L = np.maximum(L, 1e-30)  # Evitar valores cero o negativos para el logaritmo
plt.plot(x[::100], L[::100], label='Radiative loss [erg/g]')
plt.xlabel('x [Rjup]')
plt.ylabel('Radiative loss [erg/g]')
plt.yscale('log')
plt.xscale('log')
plt.title('Radiative loss profile for the best model')
plt.grid()

print("N puntos =", len(x))
print("min =", np.nanmin(L))
print("max =", np.nanmax(L))
print("NaN =", np.any(np.isnan(L)))
print("Inf =", np.any(np.isinf(L)))

plt.savefig(f"{model_name}_radiative_loss_profile.png", dpi=300, bbox_inches='tight')
