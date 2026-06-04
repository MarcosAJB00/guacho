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
    return r, v

def density_profile(r_cm, v_cms, mdot):
    return mdot / (4 * np.pi * r_cm**2 * v_cms)

# ============================================================
# PARÁMETROS
# ============================================================
T_star    = 1.0 * u.MK
R_star    = 1.0 * u.R_sun
M_star    = 1.0 * u.M_sun
mdot_star = 1.0e10              # g/s
T_star_K  = T_star.to(u.K).value

Rp        = 1.38 * u.R_jup
Rp_cm     = Rp.to(u.cm).value

a_orb     = (0.047 * u.au).to(u.cm).value   # cm

# ============================================================
# VIENTO ESTELAR
# rmax_factor solo necesita cubrir a_orb:
#   a_orb ~ 7e11 cm, R_sun ~ 6.96e10 cm → factor ~10, usamos 30
# ============================================================
r_star_m, v_star_ms = parker(M_star, R_star, T_star, rmax_factor=15)

r_star_cm  = r_star_m  * 100.0
v_star_cms = v_star_ms * 100.0

rho_star  = density_profile(r_star_cm, v_star_cms, mdot_star)
n_star    = rho_star / m_p
Pdyn_star = rho_star * v_star_cms**2
Pth_star  = n_star * k_B * T_star_K
Ptot_star = Pdyn_star + Pth_star

au_cm = (1.0 * u.au).to(u.cm).value

plt.figure(figsize=(8, 5))
plt.plot(r_star_cm/au_cm, Ptot_star, label='Presión total')
plt.plot(r_star_cm/au_cm, Pdyn_star, label='Presión dinámica')
plt.plot(r_star_cm/au_cm, Pth_star, label='Presión térmica')
plt.axvline(a_orb/au_cm, ls='--', color='k', label=f'Orbita planetaria = {a_orb/au_cm:.3f} au')
plt.yscale('log')
plt.xlabel('Distancia desde la estrella [au]')
plt.ylabel('Presión [dyn cm$^{-2}$]')
plt.title('Perfil de presión del viento estelar')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('stellar_pressure_wind.png', dpi=300, bbox_inches='tight')

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
plt.yscale('log')
plt.xlabel('Distancia desde la estrella [au]')
plt.ylabel('Velocidad [Km/s]')
plt.title('Perfil de velocidad del viento estelar')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('stellar_velocity_wind.png', dpi=300, bbox_inches='tight')


f_Ptot_star = interp1d(r_star_cm, Ptot_star,
                        bounds_error=False, fill_value="extrapolate")
f_Pdyn_star = interp1d(r_star_cm, Pdyn_star,
                        bounds_error=False, fill_value="extrapolate")
f_Pth_star  = interp1d(r_star_cm, Pth_star,
                        bounds_error=False, fill_value="extrapolate")

# ============================================================
# PERFIL PLANETARIO
# ============================================================
data_p = np.loadtxt('Hydro_ioniz.txt', comments='#')

mp =1.3 * m_p  # masa media por partícula (g)
rp    = data_p[30:, 0] * Rp_cm   # [Rp] → [cm], origen en el planeta
rho_p = data_p[30:, 1]*mp            # [g/cm³]
v_p   = data_p[30:, 2]            # [cm/s]
T_p   = data_p[30:, 4]            # [K]

print(f'rp.min(): {rp.min():.3e}, rp.max(): {rp.max():.3e}')
print(f'T_p.min(): {T_p.min():.3e}, T_p.max(): {T_p.max():.3e}')
print(f'v_p.min(): {v_p.min():.3e}, v_p.max(): {v_p.max():.3e}')
print(f'rho_p.min(): {rho_p.min():.3e}, rho_p.max(): {rho_p.max():.3e}')

mu_p   = 1.3
n_p    = rho_p / (mu_p * m_p)
Pdyn_p = rho_p * v_p**2
Pth_p  = n_p * k_B * T_p
Ptot_p = Pdyn_p + Pth_p

plt.figure(figsize=(8, 5))
plt.plot(rp/Rp_cm, rho_p, label='Densidad')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [R_p]')
plt.ylabel('Densidad [g/cm$^3$]')
plt.title('Perfil de densidad del viento planetario')
plt.savefig('planetary_density_profile.png', dpi=300, bbox_inches='tight')
plt.legend()

plt.figure(figsize=(8, 5))
plt.plot(rp/Rp_cm, v_p/1e5, label='Velocidad')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [R_p]')
plt.ylabel('Velocidad [Km/s]')
plt.title('Perfil de velocidad del viento planetario')
plt.savefig('planetary_velocity_profile.png', dpi=300, bbox_inches='tight')
plt.legend()

plt.figure(figsize=(8, 5))
plt.plot(rp/Rp_cm, T_p/1e5, label='Temperatura')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [R_p]')
plt.ylabel('Temperatura [K]')
plt.title('Perfil de temperatura del viento planetario')
plt.savefig('planetary_temperature_profile.png', dpi=300, bbox_inches='tight')
plt.legend()

plt.figure(figsize=(8, 5))
plt.plot(rp/Rp_cm, Ptot_p, label='Presión total')
plt.plot(rp/Rp_cm, Pdyn_p, label='Presión dinámica')
plt.plot(rp/Rp_cm, Pth_p, label='Presión térmica')
#plt.xscale('log')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [R_p]')
plt.ylabel('Presión [dyn cm$^{-2}$]')
plt.title('Perfil de presión del viento planetario')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('planetary_pressure_profile.png', dpi=300, bbox_inches='tight')   

f_Ptot_p = interp1d(rp, Ptot_p, bounds_error=False, fill_value="extrapolate")
f_Pdyn_p = interp1d(rp, Pdyn_p, bounds_error=False, fill_value="extrapolate")
f_Pth_p  = interp1d(rp, Pth_p,  bounds_error=False, fill_value="extrapolate")

# ============================================================
# PUNTO DE ESTANCAMIENTO
#
# Sobre la línea estrella-planeta, un punto a distancia r del planeta
# está a (a_orb - r) de la estrella.
# Buscamos r tal que:  P_planeta(r) = P_estrella(a_orb - r)
# ============================================================

# Grilla fina en r (desde el planeta), limitada a donde a_orb-r > 0
r_grid = np.linspace(rp.min(), min(rp.max(), a_orb * 0.999), 10000)

P_p_grid = f_Ptot_p(r_grid)
P_s_grid = f_Ptot_star(a_orb - r_grid)   # distancia desde la estrella

diff = P_p_grid - P_s_grid
sign_changes = np.where(np.diff(np.sign(diff)))[0]

print('Sign_canges indices:', sign_changes)

print()
print("=== Diagnóstico ===")
print(f"rp range      : {rp.min():.3e} – {rp.max():.3e} cm")
print(f"P_planeta     : {P_p_grid.min():.3e} – {P_p_grid.max():.3e} dyn/cm²")
print(f"P_estrella    : {P_s_grid.min():.3e} – {P_s_grid.max():.3e} dyn/cm²")
print(f"a_orb         : {a_orb:.3e} cm")
print()

if len(sign_changes) > 0:
    i0 = sign_changes[-1] # último cambio de signo → punto de stand-off más cercano al planeta
    # interpolación lineal para afinar
    r_standoff = r_grid[i0] - diff[i0] * (r_grid[i0+1] - r_grid[i0]) / (diff[i0+1] - diff[i0])
    print("=== Stand-off encontrado ===")
    print(f"  r_standoff = {r_standoff:.3e} cm")
    print(f"  r_standoff = {r_standoff/Rp_cm:.3f} Rp")
    print(f"  Distancia desde la estrella = {(a_orb - r_standoff):.3e} cm")
else:
    r_standoff = None
    print("=== Sin intersección en el rango del perfil planetario ===")
    if P_p_grid[0] > P_s_grid[0]:
        print("  P_planeta > P_estrella en todo el dominio")
        print("  → el stand-off está más allá del borde exterior del perfil")
    else:
        print("  P_estrella > P_planeta en todo el dominio")
        print("  → el stand-off está dentro de Rp (viento estelar aplasta al planetario)")
print()

# ============================================================
# GRÁFICO
# ============================================================
plt.figure(figsize=(8, 5))
rp_Rp   = r_grid / Rp_cm

ax = plt.gca() 
ax.plot(rp_Rp, P_p_grid, lw=2, label='Viento planetario')
ax.plot(rp_Rp, P_s_grid, lw=2, label='Viento estelar')
if r_standoff is not None:
    ax.axvline(r_standoff/Rp_cm, ls='--', color='k',
               label=f'Stand-off = {r_standoff/Rp_cm:.2f} $R_p$')
ax.set_yscale('log')
ax.set_xlabel(r'$r$ [$R_p$] (desde el planeta)')
ax.set_ylabel(r'$P_{tot}$ [dyn cm$^{-2}$]')
ax.set_title('Presión total')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('pressure_balance.png', dpi=300, bbox_inches='tight')
plt.show()