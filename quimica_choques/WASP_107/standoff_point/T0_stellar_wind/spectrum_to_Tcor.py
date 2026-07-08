import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

# =========================
# CONSTANTES
# =========================
Rsun = 6.955e10     # cm
au   = 1.496e13     # cm
h    = 6.626e-27    # erg s
c    = 2.998e10     # cm/s
evtoerg = 1.602e-12
Lsun = 3.826e33     # erg/s

# =========================
# LEER ESPECTRO
# =========================
# wvl en Angstrom, F_lambda en erg/s/cm2/A
wvl, F_lambda = np.loadtxt('HD85512_spectrum_scaled_to_WASP107_Antonia.dat', unpack=True, skiprows=1)

# Flujo escalado a la superficie de la estrella
a_espec = 0.055*au  # semieje mayor de WASP-107b
Rstar = 0.66*Rsun  # radio de WASP-107
F_lambda_surface = F_lambda * (a_espec/Rstar)**2

mask_x = (wvl > 10.0) & (wvl < 100)
Fx_trap = integrate.trapezoid(F_lambda_surface[mask_x], wvl[mask_x])
Fx_simpson = integrate.simpson(F_lambda_surface[mask_x], wvl[mask_x])


print(f"Flujo X-ray (10.0-100 A) en la superficie de la estrella: {Fx_trap:.3e} erg/s/cm2 (trapecio)")
#print(f"Flujo X-ray (10.0-100 A) en la superficie de la estrella: {Fx_simpson:.3e} erg/s/cm2 (Simpson)")

# Temperatura coronal empírica
def T_cor(F):
    return 0.11 * (F**0.26)  # en MK

T_cor_trap = T_cor(Fx_trap)
T_cor_simpson = T_cor(Fx_simpson)

print(f"Temperatura coronal en la superficie de la estrella: {T_cor_trap:.3e} MK (trapecio)")
#print(f"Temperatura coronal (10.0-100 A) en la superficie de la estrella: {T_cor_simpson:.3e} MK (Simpson)")

T0 = T_cor_trap/1.36 #Relacion entre T_cor y T0 (Aline 2018)

print(f"Temperatura coronal T0 en la superficie de la estrella: {T0:.3e} MK (trapecio)")