import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# ============================================================
# CONSTANTES
# ============================================================
G   = const.G.value
k_B = const.k_B.value
m_p = const.m_p.value

# ============================================================
# VIENTO DE PARKER
# ============================================================

def parker_eq(r, v, M, a):
    return v * (2*a**2/r - G*M/r**2) / (v**2 - a**2)

def parker(M, R, T, rmax_factor=30, mu=0.62):
    T_K  = T.to(u.K).value
    R_m  = R.to(u.m).value
    M_kg = M.to(u.kg).value

    cs = np.sqrt(k_B * T_K / (mu * m_p))
    rc = G * M_kg / (2 * cs**2)
    delta = 1e-6

    sol_in = solve_ivp(
        parker_eq, [rc*(1-delta), R_m], [cs*(1-delta)],
        args=(M_kg, cs), rtol=1e-8, atol=1e-10
    )
    sol_out = solve_ivp(
        parker_eq, [rc*(1+delta), rmax_factor*R_m], [cs*(1+delta)],
        args=(M_kg, cs), rtol=1e-8, atol=1e-10
    )

    r = np.concatenate([sol_in.t[::-1], sol_out.t])
    v = np.concatenate([sol_in.y[0][::-1], sol_out.y[0]])
    return r, v, rc, cs 

def density_profile(r_cm, v_cms, mdot):
    return mdot / (4 * np.pi * r_cm**2 * v_cms)


# ============================================================
# PARÁMETROS
# ============================================================
T_star    = 1.0 * u.MK
R_star    = 0.586 * u.R_sun
M_star    = 0.577 * u.M_sun
mdot_star = 0.45e12              # g/s 
T_star_K  = T_star.to(u.K).value

a_orb     = (0.0344 * u.au).to(u.cm).value   # cm

# ============================================================
# VIENTO ESTELAR
# rmax_factor solo necesita cubrir a_orb:
#   a_orb ~ 7e11 cm, R_sun ~ 6.96e10 cm → factor ~10, usamos 30
# ============================================================
r_star_m, v_star_ms, r_sonic_m, cs = parker(M_star, R_star, T_star, rmax_factor=25)

r_star_cm  = r_star_m  * 100.0
v_star_cms = v_star_ms * 100.0
r_sonic_cm = r_sonic_m * 100.0

rho_star  = density_profile(r_star_cm, v_star_cms, mdot_star)

au_cm = (1.0 * u.au).to(u.cm).value

plt.figure(figsize=(8, 5))
plt.plot(r_star_cm/au_cm, rho_star, label='Densidad')
plt.yscale('log')
plt.xlabel('Distancia desde la estrella [au]')
plt.ylabel('Densidad [g/cm$^3$]')
plt.title('Perfil de densidad del viento estelar')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('stellar_density_wind.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(8, 5))
plt.plot(r_star_cm/au_cm, v_star_cms/1e5, label='Velocidad')
plt.axvline(r_sonic_cm/au_cm,color='r',ls='--',label='Punto sónico')
plt.yscale('log')
plt.xlabel('Distancia desde la estrella [au]')
plt.ylabel('Velocidad [Km/s]')
plt.title('Perfil de velocidad del viento estelar')
plt.legend()
plt.grid(alpha=0.3, axis = 'both')
plt.tight_layout()
plt.savefig('stellar_velocity_wind.png', dpi=300, bbox_inches='tight')


np.savetxt('star_wind.dat',
    np.column_stack((r_star_cm, rho_star, v_star_cms)),
    header=' r_s (cm)  density (g/cm3)  velocity (cm/s)',
    fmt='%.10e'
)