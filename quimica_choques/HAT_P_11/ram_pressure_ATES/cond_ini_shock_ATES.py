import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

data = np.loadtxt('Hydro_ioniz.txt', comments='#')

mu = 1.3
m_p = 1.6726219e-24  # g
mp = mu * m_p  # masa media por partícula (g)
rp    = data[30:, 0]            # [Rp] origen en el planeta
rho = data[30:, 1]         # [1/cm³]
v   = data[30:, 2]            # [cm/s]
T   = data[30:, 4]            # [K]

r_standoff = 4.558 # en Rp

# Interpoladores
dens_interp = PchipInterpolator(rp, np.log10(rho))
temp_interp = PchipInterpolator(rp, T)
vel_interp = PchipInterpolator(rp, v)

radios_grid = [r_standoff - 2.0, r_standoff - 1.0, r_standoff, r_standoff + 1.0, r_standoff + 2.0]
rho_grid = []
T_grid = []
v_grid = []
for r in radios_grid:
    rho_r = 10**dens_interp(r)
    T_r = temp_interp(r)
    v_r = vel_interp(r)

    rho_grid.append(rho_r)
    T_grid.append(T_r)
    v_grid.append(v_r)

    print(f'En r = {r:.2f} Rp: ρ = {rho_r:.3e} g/cm³, T = {T_r:.3e} K, v = {v_r/1e5:.3e} km/s')

np.savetxt('cond_ini_shock_ATES.dat',
    np.column_stack((radios_grid, rho_grid, T_grid, v_grid)),
    header='r[Rp] density[g/cm^3] temperature[K] velocity[cm/s]',
    fmt='%.10e'
)
