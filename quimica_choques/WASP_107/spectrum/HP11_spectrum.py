import numpy as np 
import matplotlib.pyplot as plt

spectrum = np.loadtxt('HAT-P-11-stellar-spectrum-1AU-REG.dat', skiprows=1)

wvl = spectrum[:,0] # en Angstroms
flux = spectrum[:,1] # en erg/s/cm^2/Angstrom

mask = (flux > 0)
flux = flux[mask] # Filtrar valores no positivos
wvl = wvl[mask] # Filtrar longitudes de onda correspondientes

au = 1.496e13 # cm
a_spectrum = 1.0*au # cm
a_p = 0.0532*au # cm

flux_p = flux * (a_spectrum/a_p)**2

np.savetxt('HP11_spectrum_beni_scaled.dat',
    np.column_stack((wvl, flux_p)),
    header='wavelength (A)   Flux (erg/s/cm^2/Angstrom) ',
    fmt='%.10e'
)

plt.figure(figsize=(10,6))
plt.plot(wvl, flux, label=f'HAT-P-11 at {a_spectrum/au:.3f} AU')
plt.plot(wvl, flux_p, label=f'HAT-P-11 at {a_p/au:.3f} AU')
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (erg/s/cm^2/Angstrom)')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1, 10000)
plt.legend()    
plt.grid()
plt.title('Stellar Spectrum at HAT-P-11b')
plt.savefig('HAT-P-11b_spectrum.png', dpi=300, bbox_inches='tight')