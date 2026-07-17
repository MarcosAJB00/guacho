import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d

# ============================================================
# Parametros generales
# ============================================================
Rjup = 7.1492e9  # cm
Rp = 0.4466 * Rjup
r_standoff = 4.77          # en Rp

n_planeta = 1e8             # densidad dentro del planeta

# ============================================================
# Perfil 1D del viento planetario
# ============================================================
planeta = np.loadtxt('perfiles_finales.dat', comments='#')
r_p = planeta[:, 0]          # en Rp
HeIM_p = planeta[:, 4]       # en 1/cm3

interp_he_p = interp1d(r_p, HeIM_p, bounds_error=False, fill_value=0.0)

# ============================================================
# Perfil 1D del choque (shock)
# ============================================================
shock = np.loadtxt('final_model_1.dat', comments='#')

r_s = shock[:, 0]                      # en cm
shock_size = r_s[-1] / Rp              # tamano del choque en Rp
r_s = r_s / Rp + r_standoff - shock_size   # en Rp, corrido para arrancar en r_standoff - shock_size
HeIM_s = shock[:, 4]                   # en 1/cm3

interp_he_s = interp1d(r_s, HeIM_s, bounds_error=False, fill_value=0.0)

# ------------------------------------------------------------
# Parche del choque
# ------------------------------------------------------------
# theta_c = 90 grados en el codigo original (parche centrado en el ecuador,
# es decir, contenido en este mismo plano z=0). Por eso ac치 la distancia
# angular al centro del parche se reduce simplemente a una diferencia de phi.
phi_c = np.deg2rad(77.0)
alpha = np.deg2rad(88.0)   # tamano angular (radio) del parche

# ============================================================
# Grilla 2D (plano XY, z = 0 -> el mismo plano que graficaba el codigo original)
# ============================================================
N = 512
L = r_p[-1]

x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
X, Y = np.meshgrid(x, y, indexing='ij')   # X[i,j] = x[i], Y[i,j] = y[j]
R = np.sqrt(X**2 + Y**2)

# ============================================================
# Coordenadas polares
# ============================================================
phi = np.arctan2(Y, X)
# ============================================================
# Angulo respecto del eje del choque
# ============================================================

theta = np.abs(np.angle(np.exp(1j*(phi - phi_c))))
# ============================================================
# Radio del choque de Wilkin
# ============================================================

eps = 1e-8

theta_safe = np.maximum(theta, eps)

cot = np.cos(theta_safe)/np.sin(theta_safe)

Rshock = (
    r_standoff
    / np.sin(theta_safe)
    * np.sqrt(
        3.0*(1.0-theta_safe*cot)
    )
)

Rshock[theta < 1e-4] = r_standoff

rho2D = np.zeros_like(R)

# interior del planeta
mask_planet = (R < 1.0)
rho2D[mask_planet] = n_planeta

mask_wind = (
    (R >= 1.0)
    &
    (R < (Rshock - shock_size))
)

rho2D[mask_wind] = interp_he_p(R[mask_wind])
# ============================================================
# Máscara angular del parche
# ============================================================
dphi = np.angle(np.exp(1j*(phi - phi_c)))
mask_angular = np.abs(dphi) <= alpha

# ============================================================
# Máscara radial del choque
# ============================================================
mask_radial = (
    (R >= (Rshock - shock_size))
    &
    (R <= Rshock)
)

mask_shock = mask_angular & mask_radial

# ============================================================
# Radio equivalente dentro del perfil 1D del choque
# ============================================================

r_local = (
    R[mask_shock]
    -
    (Rshock[mask_shock] - shock_size)
    +
    (r_standoff - shock_size)
)

rho2D[mask_shock] = interp_he_s(r_local)

# resto (fuera de planeta, viento y choque)
mask_resto = ~(mask_planet | mask_wind | mask_shock)
rho2D[mask_resto] = 0.0


# ============================================================
# Curva analítica de Wilkin (solo para graficar)
# ============================================================

theta_plot = np.linspace(1e-4, np.pi-1e-4, 2000)

cot = np.cos(theta_plot)/np.sin(theta_plot)

Rplot = (
    r_standoff
    / np.sin(theta_plot)
    * np.sqrt(
        3.0*(1.0-theta_plot*cot)
    )
)

# Rama superior
phi_plot = phi_c + theta_plot

x1 = Rplot*np.cos(-phi_plot)
y1 = Rplot*np.sin(-phi_plot)

# Rama inferior
phi_plot = phi_c - theta_plot

x2 = Rplot*np.cos(-phi_plot)
y2 = Rplot*np.sin(-phi_plot)

# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(8, 7))

plt.imshow(
    rho2D.T,
    origin='upper',              # importante: coherente con como armamos X,Y
    extent=[-L, L, -L, L],
    norm=LogNorm(vmin=1.0e0, vmax=1.0e6)
)
plt.colorbar(label=r'$n_{\rm HeI^*}$ [cm$^{-3}$]')

plt.xlabel(r'$x$ [R$_p$]')
plt.ylabel(r'$y$ [R$_p$]')
plt.title('Perfil 1D proyectado en el plano XY (z=0)')

# planeta
circ1 = plt.Circle((0, 0), 1.0, color='white', fill=False, lw=2, label='Planeta')
plt.gca().add_artist(circ1)

# standoff
circ2 = plt.Circle((0, 0), r_standoff, color='k', fill=False, ls='--', lw=1.2, label=r'$r_{standoff}$')
plt.gca().add_artist(circ2)

plt.plot(x1, y1, 'r--', lw=1.5, label='Choque analítico')
plt.plot(x2, y2, 'r--', lw=1.5)
#plt.xlim(-10*L, 10*L)
#plt.ylim(-5*L, 10*L)

plt.xlim(-1.2*L, 1.2*L)
plt.ylim(-0.8*L, 1.2*L)

plt.legend(loc='upper right', fontsize=8)
plt.gca().set_aspect('equal')

plt.tight_layout()
plt.savefig('2D_xy_shock_analitico.png', dpi=300, bbox_inches='tight')
plt.show()