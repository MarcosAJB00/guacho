import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import	astropy.units as u
import	astropy.constants as const
from scipy.interpolate import interp1d

# Constantes, por default en MKS, les pongo .value para no renegar en las funciones
G = const.G.value
k_B = const.k_B.value
#max_rang = (10.0*1.38*u.R_jup).to(u.m).value
max_rang = 4.5 #in Rp
m_p	= const.m_p.value

# =========================
def parker_eq(r, v, M, a):
    ecu = v*(2 * a**2 / r - G * M / r**2) / (v**2 - a**2 )
    return ecu

def parker(M, R, T, mu=0.602):

    R_sc = R
    T = (T.to(u.K)).value
    R = (R.to(u.m)).value
    M = (M.to(u.kg)).value


    a = np.sqrt(k_B * T / (mu * m_p))  # Velocidad del sonido

    r_c = G * M / (2 * a**2)
    v_c = a

    # Radio en el que se resuelve la ecuación
    r_min = R
    r_max = max_rang * R

    # Condiciones iniciales cerca del punto crítico (cerca para que no explote)
    delta = 1e-6
    v_izq = v_c * (1 - delta)
    v_der = v_c * (1 + delta)

    # Integrar hacia adentro
    sol_dentro = solve_ivp(parker_eq, [r_c * (1 - delta), r_min], [v_izq], args=(M, a), rtol=1e-8)

    # Integrar hacia afuera
    sol_afuera = solve_ivp(parker_eq, [r_c * (1 + delta), r_max], [v_der], args=(M, a), rtol=1e-8)

    # Juntar soluciones
    r_combinado = np.concatenate((sol_dentro.t[::-1], sol_afuera.t))
    v_combinado = np.concatenate((sol_dentro.y[0][::-1], sol_afuera.y[0]))

    # Unidades
    r_combinado = r_combinado * u.m
    v_combinado = v_combinado * u.m / u.s
    r_c	= r_c * u.m
    a = a * u.m / u.s

    return (r_combinado).to(u.R_jup)/R_sc, v_combinado.to(u.km/u.s), (r_c).to(u.R_jup)/R_sc, a.to(u.km/u.s)

def density_profile(r, v, mdot):
    """
    Calcula el perfil de densidad para un viento de Parker.
    """
    return mdot / (4.0 * np.pi * r**2 * v)


def resample_profile(r, rho, n_points=200):
    """
    Reinterpola el perfil en radios equiespaciados.
    """
    # asegurar arrays numpy
    r = np.array(r, dtype=float)
    rho = np.array(rho, dtype=float)

    # radios equiespaciados
    r_new = np.linspace(r.min(), r.max(), n_points)

    # interpolación lineal
    f = interp1d(r, rho, kind='linear')
    rho_new = f(r_new)

    return r_new, rho_new
# =========================

#===================================
# Viento de Parker - VELOCIDAD
#===================================
# Parámetros
T_star = 9.0e3*u.K # Temperatura de la corona
r_star = 1.38*u.R_jup # Radio estelar
M_star = 0.69*u.M_jup # Masa estelar
r, v, r_c, a = parker(M_star, r_star, T_star)

# Graficar
plt.plot(r, v, label=f'T = {T_star.value} K')
plt.xlabel(r'$r$')
plt.ylabel('v')
plt.title('Solución de Parker para el viento lento')
plt.legend()
plt.grid(True)
plt.yscale('log')
#plt.savefig('parker_viento.png', dpi=300, bbox_inches='tight')
#plt.show()
#===================================

#===================================
# Viento de Parker - DENSIDAD
#===================================
mdot = 4.0e10 # en g/s
Rp = 1.38 * u.R_jup
# perfil original
rho = density_profile((r * Rp).to(u.cm).value,
                      v.to(u.cm/u.s).value,
                      mdot)

# nuevos puntos equiespaciados
n_points = 100
r_new, rho_new = resample_profile(r, rho, n_points=n_points)

# guardar archivo
np.savetxt(
    "neutral_density_profile.dat",
    np.column_stack([r_new, rho_new]),
    header="r(Rp)   rho(g/cm^3)"
)

# gráfico
plt.plot(r, rho, 'o', label="Original")
plt.plot(r_new, rho_new, '-', label="Interpolado")
plt.yscale('log')
plt.xlabel("r / Rp")
plt.ylabel("rho [g/cm^3]")
plt.legend()
plt.grid(True)
plt.title("Perfil de densidad del viento de Parker")
#plt.savefig('parker_density.png', dpi=300, bbox_inches='tight')
#plt.show()