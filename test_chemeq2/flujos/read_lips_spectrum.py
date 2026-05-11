import numpy as np
import matplotlib.pyplot as plt

# leer ignorando header
data = np.loadtxt('sorce_ssi_l3.csv', delimiter=',', skiprows=1)

# extraer columnas
wvl_nm = data[:,1]
F_earth = data[:,2]

# distancia
a = 0.04747  # AU

# factor geometrico
geom = (1.0 / a)**2

# conversion unidades
conv = 100.0

# flujo en el planeta
F_planet = F_earth * conv * geom  # erg/s/cm2/A

wvl = wvl_nm * 10.0  # nm → Angstrom
# guardar
np.savetxt('spectrum_lips.dat', np.column_stack([wvl, F_planet]))

plt.plot(wvl, F_planet)
plt.scatter(wvl, F_planet/10.0, s=10, color='red', label='LIPS data')
plt.xlim(10, 2700)
plt.yscale('log')
plt.xlabel('Wavelength [Angstrom]')
plt.ylabel('Flux [erg s$^{-1}$ cm$^{-2}$ Angstrom$^{-1}$]')
plt.title('Flux at HD 209458b')
plt.savefig('spectrum_lips.png', dpi=300)
plt.show()
