import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

# leer ignorando header
data = np.loadtxt('sorce_ssi_l3.csv', delimiter=',', skiprows=1)

# extraer columnas
wvl_nm = data[:,1]
F_earth = data[:,2]

# distancia
a = 0.04747  # AU

# conversion unidades, de W/m2/nm a erg/s/cm2/A
# 1 W = 10^7 erg/s
# 1 m2 = 10^4 cm2
# 1 nm = 10 A
conv = 100.0 # 10^7 / 10^4/ 10

# flujo en el planeta
F_planet = F_earth * conv * (1.0 / a)**2  # erg/s/cm2/A

wvl = wvl_nm * 10.0  # nm → Angstrom

# =========================
# INTEGRACIONES
# =========================

# Ly-alpha [121.5–121.7 nm]
lya_min = 115.0 * 10.0  # Angstrom
lya_max = 130.0 * 10.0  # Angstrom
mask_lya = (wvl > lya_min) & (wvl < lya_max)
F_Lya = integrate.simpson(F_planet[mask_lya], wvl[mask_lya])
print(f'Flux Ly-alpha: {F_Lya:.2e} erg/s/cm2')

# ==========================
# ESTIMACION EUV A PARTIR DE Ly-alpha
# Referencias: Linsky+2014
# ==========================
# flujo Ly-alpha (erg/s/cm^2)
logF = np.log10(F_Lya)

# 40–50 nm
F_40_50 = 10**(-2.294 + 0.258*logF)

# 50–60 nm
F_50_60 = 10**(-2.098 + 0.572*logF)

# 60–70 nm
F_60_70 = 10**(-1.920 + 0.240*logF)

# 70–80 nm
F_70_80 = 10**(-1.894 + 0.518*logF)

# 80–91.2 nm
F_80_91 = 10**(-1.811 + 0.764*logF)

# 91.2–117 nm
F_91_117 = 10**(-1.004 + 0.065*logF)

print("EUV fluxes:")
print(f"40–50 nm:   {F_40_50:.3e}")
print(f"50–60 nm:   {F_50_60:.3e}")
print(f"60–70 nm:   {F_60_70:.3e}")
print(f"70–80 nm:   {F_70_80:.3e}")
print(f"80–91.2 nm: {F_80_91:.3e}")
print(f"91.2–117 nm:{F_91_117:.3e}")

# centros (en Angstrom)
wvl_euv = np.array([
    450,   # 40–50 nm
    550,   # 50–60 nm
    650,   # 60–70 nm
    750,   # 70–80 nm
    856,   # 80–91.2 nm
    1041   # 91.2–117 nm
])

F_bands = np.array([
    F_40_50,
    F_50_60,
    F_60_70,
    F_70_80,
    F_80_91,
    F_91_117
])

# guardar
#np.savetxt('spectrum_lips.dat', np.column_stack([wvl, F_planet]))


plt.plot(wvl, F_planet, label='Original spectrum')
plt.scatter(wvl, F_planet, s=10, color='red', label='LIPS data')
# puntos EUV
plt.scatter(wvl_euv, F_bands, color='green', s=40, label='EUV (Linsky+2014)')
# opcional: unirlos
plt.plot(wvl_euv, F_bands, linestyle='--', color='green', alpha=0.7)
#Ly-alpha
plt.axvline(1215.67, color='blue', linestyle=':', label='Ly-alpha')
plt.axvline(lya_min, color='red', linestyle='--', label='Ly-alpha min')
plt.axvline(lya_max, color='red', linestyle='--', label='Ly-alpha max')

plt.legend()
plt.xlim(10, 2700)
plt.yscale('log')
plt.xlabel('Wavelength [Angstrom]')
plt.ylabel('Flux [erg s$^{-1}$ cm$^{-2}$ Angstrom$^{-1}$]')
plt.title('Flux at HD 209458b')
plt.savefig('spectrum_lips.png', dpi=300)
plt.show()


wvl_full = np.concatenate([wvl_euv, wvl])
F_full = np.concatenate([F_bands, F_planet])

# ordenar por longitud de onda
idx = np.argsort(wvl_full)
wvl_full = wvl_full[idx]
F_full = F_full[idx]

np.savetxt('spectrum_lips_euv.dat',
           np.column_stack([wvl_full, F_full]),
           header='wavelength[A] flux[erg/s/cm2/A]')