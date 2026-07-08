# %%
import numpy as np

# Parámetros de la estrella
logLx_HD = 28.88
Lx_HD = 10**logLx_HD  # erg/s
R_HD = 0.912 * 6.957e10  # cm


# Flujo de rayos X en superficie
Fx_HD = Lx_HD / (4 * np.pi * (R_HD**2))  # erg/cm^2/s


# Temperatura coronal empírica
def T_cor(F):
    return 0.11 * (F**0.26)  # en MK

T_HD = T_cor(Fx_HD)

# Constantes físicas
M_HD = 0.99 * 1.989e30  # kg
k_B = 1.3807e-23        # J/K
m_H = 1.673e-27         # kg
mu = 0.602              # peso molecular medio
G = 6.672e-11           # N m^2/kg^2

# Velocidad del sonido (m/s)
a_HD = np.sqrt(k_B * T_HD * 1e6 / (mu * m_H))

# Radio crítico en metros
r_c = G * M_HD / (2 * a_HD**2)

# Conversión a radios solares y UA
r_c_Rsol = r_c / 6.957e8
r_c_UA = r_c / (1.496e11)

print(f"Fx = {Fx_HD:.2e} erg/cm^2/s")
print(f"T_HD = {T_HD:.2f} MK")
print(f"a_HD = {a_HD:.2e} m/s")
print(f"r_c = {r_c_Rsol:.2f} R_sun")
print(f"r_c = {r_c_UA:.4f} AU")
