import numpy as np
import matplotlib.pyplot as plt

# ======================================
# Temperatura
# ======================================

T = np.logspace(2, 6, 1000)  # 10^2 - 10^6 K

kB_eV = 8.617333262145e-5
logT = np.log10(T)

# ======================================
# Gammas Bray+2000
# ======================================

gamma13 = np.zeros_like(T)
gamma31a = np.zeros_like(T)
gamma31b = np.zeros_like(T)

mask = (logT >= 3.75) & (logT <= 5.75)

gamma13[mask] = (
    -1.5699
    + 1.303 * logT[mask]
    - 0.388 * logT[mask]**2
    + 0.0516 * logT[mask]**3
    - 0.0026 * logT[mask]**4
)

gamma31a[mask] = (
    -209.42
    + 170.03 * logT[mask]
    - 50.073 * logT[mask]**2
    + 6.42 * logT[mask]**3
    - 0.3045 * logT[mask]**4
)

gamma31b[mask] = (
    -36.477
    + 22.988 * logT[mask]
    - 4.5935 * logT[mask]**2
    + 0.2968 * logT[mask]**3
)

# ======================================
# Rates asociados a HeIM
# ======================================

# ---------- PRODUCCIÓN ----------

# HeIS -> HeIM
iche13a = (
    2.10e-8
    * np.sqrt(13.6 / (kB_eV * T))
    * np.exp(-19.81 / (kB_eV * T))
    * gamma13
)

# HeII -> HeIM
iaheiim_b = 2.10e-13 * (1e4 / T)**0.778


# ---------- DESTRUCCIÓN ----------

# Decaimiento radiativo
irheii = np.full_like(T, 1.272e-4)

# Colisión con HI
iche31 = np.full_like(T, 5.e-10)

# Colisiones electrónicas
iche31a = (
    2.10e-8
    * np.sqrt(13.6 / (kB_eV * T))
    * np.exp(-0.80 / (kB_eV * T))
    * gamma31a / 3
)

iche31b = (
    2.10e-8
    * np.sqrt(13.6 / (kB_eV * T))
    * np.exp(-1.40 / (kB_eV * T))
    * gamma31b / 3
)

# Ionización colisional
icheIM = (
    8.335e-10
    * np.sqrt(T)
    * np.exp(-55338.0 / T)
)

# ======================================
# Plot
# ======================================

plt.figure(figsize=(10,7))

# Producción
plt.loglog(
    T, iche13a, ls='--',
    label='Production: iche13a (HeIS → HeIM)'
)

plt.loglog(
    T, iaheiim_b, ls='--',
    label='Production: iaheiim_b (HeII → HeIM)'
)

# Destrucción
plt.loglog(
    T, irheii,
    label='Destruction: irheii (radiative decay)'
)

plt.loglog(
    T, iche31,
    label='Destruction: iche31 (HI collisions)'
)

plt.loglog(
    T, iche31a,
    label='Destruction: iche31a (e⁻ collisions)'
)

plt.loglog(
    T, iche31b,
    label='Destruction: iche31b (e⁻ collisions)'
)

plt.loglog(
    T, icheIM,
    label='Destruction: icheIM (collisional ionization)'
)

plt.xlabel('Temperature [K]')
plt.ylabel(r'Rate coefficient [cm$^3$ s$^{-1}$]')

plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-30, 1e-6)

plt.grid(True, which='both', alpha=0.3)

plt.legend(fontsize=9)

plt.tight_layout()

plt.savefig(
    'HeIM_rates_vs_T.png',
    dpi=300,
    bbox_inches='tight'
)

plt.show()