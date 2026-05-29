import numpy as np 
import matplotlib.pyplot as plt

data = np.loadtxt("./output/best_model.dat", comments="#")

x = data[1:, 0]
rho = data[1:, 1]
T = data[1:, 2]
n = data[1:, 3]
u = data[1:, 4]
y = data[1:, 5]
P = data[1:, 6]
h = data[1:, 7]
L = data[1:, 8]

plt.figure(figsize=(8,6))
plt.plot(x, rho, label='Density [g/cm^3]')
plt.xlabel('x [cm]')
plt.ylabel('Density [g/cm^3]')
plt.yscale('log')
plt.title('Density profile for the best model')
plt.grid()
plt.savefig("best_density_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, T, label='Temperature [K]')
plt.xlabel('x [cm]')
plt.ylabel('Temperature [K]')
plt.yscale('log')
plt.title('Temperature profile for the best model')
plt.grid()
plt.savefig("best_temperature_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, n, label='Number density [cm^-3]')
plt.xlabel('x [cm]')
plt.ylabel('Number density [cm^-3]')
plt.yscale('log')
plt.title('Number density profile for the best model')
plt.grid()
plt.savefig("best_number_density_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, u, label='Velocity [cm/s]')
plt.xlabel('x [cm]')
plt.ylabel('Velocity [cm/s]')
plt.yscale('log')
plt.title('Velocity profile for the best model')
plt.grid()
plt.savefig("best_velocity_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, y, label='Ionization fraction')
plt.xlabel('x [cm]')
plt.ylabel('Ionization fraction')
plt.yscale('log')
plt.title('Ionization fraction profile for the best model')
plt.grid()
plt.savefig("best_ionization_fraction_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, P, label='Pressure [dyn/cm^2]') 
plt.xlabel('x [cm]')
plt.ylabel('Pressure [dyn/cm^2]')
plt.yscale('log')
plt.title('Pressure profile for the best model')
plt.grid()
plt.savefig("best_pressure_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, h, label='Enthalpy [erg/g]')
plt.xlabel('x [cm]')
plt.ylabel('Enthalpy [erg/g]')
plt.yscale('log')
plt.title('Enthalpy profile for the best model')
plt.grid()
plt.savefig("best_enthalpy_profile.png", dpi=300, bbox_inches='tight')

plt.figure(figsize=(8,6))
plt.plot(x, L, label='Radiative loss [erg/g]')
plt.xlabel('x [cm]')
plt.ylabel('Radiative loss [erg/g]')
plt.yscale('log')
plt.xscale('log')
plt.title('Radiative loss profile for the best model')
plt.grid()
plt.savefig("best_radiative_loss_profile.png", dpi=300, bbox_inches='tight')
