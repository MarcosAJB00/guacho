import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('output.dat')

ci_importante = 'u50kms'

x = data[1:, 0]
h = data[1:, 1]
y = data[1:, 2]
rho = data[1:, 3]
u = data[1:, 4]
T = data[1:, 5]
P = data[1:, 6]
n = data[1:, 7]

plt.figure(figsize=(10, 6))
plt.plot(x, P, label='Presión (P)', color='blue')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Distancia (x)')
plt.ylabel('Presión (P)')
plt.title('Perfil de Presión a lo largo de la Distancia')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig(f'./plots/presion_perfil.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, T, label='Temperatura (T)', color='red')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Distancia (x)')
plt.ylabel('Temperatura (T)')
plt.title('Perfil de Temperatura a lo largo de la Distancia')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('./plots/temperatura_perfil.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, n, label='Densidad (1/cm^3)', color='green')
#plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distancia (x)')
plt.ylabel('Densidad (n)')
plt.title('Perfil de Densidad a lo largo de la Distancia')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('./plots/dens_num_perfil.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, rho, label='Densidad (rho)', color='green')
#plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distancia (x)')
plt.ylabel('Densidad (g/cm^3)')
plt.title('Perfil de Densidad a lo largo de la Distancia')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('./plots/densidad_perfil.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, u/1e5, label='Velocidad (u)', color='blue')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Distancia (x)')
plt.ylabel('Velocidad (km/s)')
plt.title('Perfil de Velocidad a lo largo de la Distancia')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('./plots/velocidad_perfil.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, h, label='Entalpia (h)', color='blue')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Distancia (x)')
plt.ylabel('Entalpia (erg/g)')
plt.title('Perfil de Entalpia a lo largo de la Distancia')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('./plots/entalpia_perfil.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, y, label=' fraccion de ionizacion (y)', color='blue')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Distancia (x)')
plt.ylabel(' fraccion de ionizacion H')
plt.title('Perfil de fraccion de ionizacion a lo largo de la Distancia')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('./plots/frac_ioni_perfil.png',dpi=300,bbox_inches='tight')
plt.show()