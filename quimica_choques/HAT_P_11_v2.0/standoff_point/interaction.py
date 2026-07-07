import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from scipy.interpolate import interp1d

# Constantes (cgs)
m_p = 1.67e-24  # g
mp_s = 1.0 * m_p  # masa media por partícula estrella (g), H puro
mp_p = 1.3 * m_p  # masa media por partícula planeta (g), H y He
k_B = 1.38e-16  # erg/K
au_cm = (1.0 * u.au).to(u.cm).value

# ===================================================
# ESTRELLA
# ===================================================
# Viento de la estrella (Parker)
star = np.loadtxt('HP11_star_wind.dat', comments='#')

r_s = star[:,0] #en cm
rho_s = star[:,1] #en g/cm3
vel_s = star[:,2] #en cm/s

# Presion viento estelar
T_star = 2.5 * u.MK
n_star    = rho_s / mp_s

Pdyn_s = rho_s * vel_s**2 # en dyn/cm2
Pth_s  = n_star * k_B * T_star.value # en dyn/cm2
Ptot_s = Pdyn_s + Pth_s  # en dyn/cm2

# ===================================================
# PLANETA
# ===================================================

# Hydro profiles
#r,rho,v,p,T,heat,cool = \
#      np.loadtxt('../ATES/Hydro_ioniz.txt',unpack = True)
#rho = rho*mu
planet = np.loadtxt('../ATES/Hydro_ioniz.txt')

r_p = planet[:,0] # en Rp
rho_p = planet[:,1]*mp_p # en g/cm3
vel_p = planet[:,2] # en cm/s
T_p = planet[:,4] #en K

Rp = (0.4466*u.Rjup).to(u.cm).value

n_p    = rho_p / mp_p

Pdyn_p = rho_p * vel_p**2 # en dyn/cm2
Pth_p  = n_p * k_B * T_p # en dyn/cm2
Ptot_p = Pdyn_p + Pth_p # en dyn/cm2

cs_p = np.sqrt(k_B * T_p / mp_p) # velocidad del sonido (cm/s)
r_sonic = r_p[np.argmin(np.abs(vel_p - cs_p))]
print(f"Radio sónica: {r_sonic:.3f} Rp")

a_orb = (0.0532*u.au).to(u.cm).value #semi eje mayor del planeta

cs_p 


#Grafico presion estrella
plt.figure(figsize=(8, 5))
plt.plot(r_s/au_cm, Ptot_s, label='Presión total')
plt.plot(r_s/au_cm, Pdyn_s, ls='--', label='Presión dinámica')
#plt.plot(r_s/au_cm, Pth_s, label='Presión térmica')
plt.axvline(a_orb/au_cm, ls='--', color='k', label=f'Orbita planetaria = {a_orb/au_cm:.3f} au')
plt.yscale('log')
#plt.xlim((a_orb/au_cm) -0.005, (a_orb/au_cm) +0.005)
plt.xlabel('Distancia desde la estrella [au]')
plt.ylabel('Presión [dyn cm$^{-2}$]')
plt.title('Perfil de presión del viento estelar')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('stellar_pressure_wind.png', dpi=300, bbox_inches='tight')

#Grafico presion planetaria
plt.figure(figsize=(8, 5))
plt.plot(r_p, Ptot_p, label='Presión total')
plt.plot(r_p, Pdyn_p, label='Presión dinámica')
plt.plot(r_p, Pth_p, label='Presión térmica')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [Rp]')
plt.ylabel('Presión [dyn cm$^{-2}$]')
plt.title('Perfil de presión del viento planetario')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('planet_pressure_wind.png', dpi=300, bbox_inches='tight')

#Grafico densidad planetaria
plt.figure(figsize=(8, 5))
plt.plot(r_p, rho_p, label='Densidad total')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [Rp]')
plt.ylabel('Densidad [g/cm3]')
plt.title('Perfil de densidad del viento planetario')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('planet_density_wind.png', dpi=300, bbox_inches='tight')

#Grafico velocidad planetaria
plt.figure(figsize=(8, 5))
plt.plot(r_p, vel_p/1e5, label='Velocidad total')
plt.scatter(r_sonic, vel_p[np.argmin(np.abs(vel_p - cs_p))]/1e5, color='red', zorder=5, label=f'Punto sónico = {r_sonic:.3f} Rp')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [Rp]')
plt.ylabel('Velocidad [Km/s]')
plt.title('Perfil de velocidad del viento planetario')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('planet_velocity_wind.png', dpi=300, bbox_inches='tight')

# ===================================================================
# INTERACCION ESTRELLA-PLANETA --> STAND-OFF POINT
# ===================================================================


r_p_cm = r_p*Rp #de Rp a cm
#interpolo la presion estelar y planetario
f_Ptot_star = interp1d(r_s, Ptot_s,bounds_error=False, fill_value="extrapolate")
f_Ptot_p = interp1d(r_p_cm, Ptot_p, bounds_error=False, fill_value="extrapolate")

# Sobre la línea estrella-planeta, un punto a distancia r del planeta
# está a (a_orb - r) de la estrella.
# Buscamos r tal que:  P_planeta(r) = P_estrella(a_orb - r)

# Grilla fina en r (desde el planeta)
r_grid = np.linspace(r_p_cm.min(), r_p_cm.max(), 10000)

P_p_grid = f_Ptot_p(r_grid)
P_s_grid = f_Ptot_star(a_orb - r_grid)   # distancia desde la estrella

diff = P_p_grid - P_s_grid # diferencia de presiones
sign_changes = np.where(np.diff(np.sign(diff)))[0] #donde cambia de signo las diferencias

f_rho_star = interp1d(r_s, rho_s,
                      bounds_error=False,
                      fill_value="extrapolate")

f_vel_star = interp1d(r_s, vel_s,
                      bounds_error=False,
                      fill_value="extrapolate")
print(sign_changes)

if len(sign_changes) > 0:
    i = sign_changes[-1]
    r_standoff = r_grid[i]

    # Distancia desde la estrella al punto de stand-off
    r_star_off = a_orb - r_standoff

    # Condiciones del viento estelar en el stand-off
    rho_star_off = f_rho_star(r_star_off)
    vel_star_off = f_vel_star(r_star_off)

    print("Stand-off point:")
    print(f"{r_standoff/Rp:.3f} Rp")
    print(f"Distancia a la estrella = {r_star_off/au_cm:.6f} au")
    print(f"Densidad estelar = {rho_star_off:.5e} g/cm^3")
    print(f"Velocidad estelar = {vel_star_off/1e5:.3f} km/s")


#Extraer las condiciones inicales para los choques

if r_standoff is not None:

    r_off = r_standoff / Rp

    r_list = [
        r_off - 2.0,
        r_off - 1.0,
        r_off,
        r_off + 1.0,
        r_off + 2.0
    ]

    f_n = interp1d(r_p, n_p)
    f_T = interp1d(r_p, T_p)
    f_u = interp1d(r_p, vel_p)

    n0_list = []
    T0_list = []
    u0_list = []

    for r in r_list:
        n0_list.append(f_n(r))
        T0_list.append(f_T(r))
        u0_list.append(f_u(r))

    np.savetxt(
        'shocks_ci.dat',
        np.column_stack((r_list, n0_list, T0_list, u0_list)),
        header='r(Rp) n0(cm^-3) T0(K) u0(cm/s)',
        fmt='%.5e'
    )

      
#Grafico presion planetaria y estelar
plt.figure(figsize=(8, 5))
plt.plot(r_grid/Rp, P_p_grid, label='Presión total planetaria')
plt.plot(r_grid/Rp, P_s_grid, label='Presión total estelar')
f_Pp = interp1d(r_grid, P_p_grid)
plt.scatter(r_sonic, f_Pp(r_sonic*Rp), color='red', zorder=5, label=f'Punto sónico = {r_sonic:.3f} Rp')
if r_standoff is not None:
    plt.axvline(r_standoff/Rp, ls='-', color='k',
               label=f'Stand-off = {r_standoff/Rp:.2f} $R_p$')
    for r in r_list:
        plt.axvline(r, ls='--', color='red', alpha=0.6)

plt.yscale('log')
plt.xlabel('Distancia desde el planeta [Rp]')
plt.ylabel('Presión [dyn cm$^{-2}$]')
plt.title('Perfil de presión del viento planetario')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('interaction_pressure_wind.png', dpi=300, bbox_inches='tight')



