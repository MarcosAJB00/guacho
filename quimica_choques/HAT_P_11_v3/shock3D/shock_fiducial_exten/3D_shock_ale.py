import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d, RegularGridInterpolator
import os

# ============================================================
# Parametros generales
# ============================================================
Rjup = 7.1492e9  # cm
Rp = 0.4466 * Rjup
r_standoff = 4.77          # en Rp

n_planeta = 1e8
T_planeta = 847.0
v_planeta = 0.0

# ============================================================
# Perfil 1D del viento planetario
# ============================================================
planeta = np.loadtxt('perfiles_finales.dat', comments='#')
planet_T_profile = np.loadtxt('temperature_profile_ATES.dat', comments='#')
planet_v_profile = np.loadtxt('velocity_profile_ATES.dat', comments='#')

r_p = planeta[:, 0]          # en Rp
HeIM_p = planeta[:, 4]       # en 1/cm3
T_p = planet_T_profile[:, 1]        # en K
v_p = planet_v_profile[:, 1] * 1e5  # de Km/s a cm/s

interp_he_p = interp1d(r_p, HeIM_p, bounds_error=False, fill_value=0.0)
interp_T_p  = interp1d(r_p, T_p, bounds_error=False, fill_value=0.0)
interp_v_p  = interp1d(r_p, v_p, bounds_error=False, fill_value=0.0)

# ============================================================
# Perfil 1D del choque (shock)
# ============================================================
shock = np.loadtxt('final_model_1.dat', comments='#')
shock_T_v = np.loadtxt('model_1.dat', comments='#')

r_s = shock[:, 0]                      # en cm
shock_size = r_s[-1] / Rp              # tamano del choque en Rp
r_s = r_s / Rp + r_standoff - shock_size   # en Rp, corrido para arrancar en r_standoff - shock_size
HeIM_s = shock[:, 4]                   # en 1/cm3
T_s = shock_T_v[:, 2]                  # en K
v_s = shock_T_v[:, 3]                  # en cm/s

interp_he_s = interp1d(r_s, HeIM_s, bounds_error=False, fill_value=0.0)
interp_T_s  = interp1d(r_s, T_s, bounds_error=False, fill_value=0.0)
interp_v_s  = interp1d(r_s, v_s, bounds_error=False, fill_value=0.0)

# ============================================================
# Geometria del choque: eje de simetria (superficie de revolucion)
# ============================================================
# theta_c = 90 grados -> el eje del choque queda contenido en el plano
# orbital XY (igual que en plot.py). phi_c es la direccion azimutal del
# eje dentro de ese plano (misma convencion que shock_ale_plot.py).
theta_c = np.deg2rad(90.0)
phi_c   = np.deg2rad(77.0)
alpha   = np.deg2rad(88.0)   # extension angular del choque (parche), medida desde el eje

u_c = np.array([
    np.sin(theta_c) * np.cos(phi_c),
    np.sin(theta_c) * np.sin(phi_c),
    np.cos(theta_c)
])
print('u_c = ', u_c)

# ============================================================
# Grilla 3D
# ============================================================
N = 256
L = r_p[-1]

x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
z = np.linspace(-L, L, N)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
R = np.sqrt(X**2 + Y**2 + Z**2)

# versor radial de cada celda
ux = np.zeros_like(R)
uy = np.zeros_like(R)
uz = np.zeros_like(R)
mask_nonzero = (R > 0)
ux[mask_nonzero] = X[mask_nonzero] / R[mask_nonzero]
uy[mask_nonzero] = Y[mask_nonzero] / R[mask_nonzero]
uz[mask_nonzero] = Z[mask_nonzero] / R[mask_nonzero]

# angulo respecto del eje del choque: este "gamma" es el analogo 3D
# del "theta" usado en shock_ale_plot.py (alli theta era, por construccion,
# la distancia angular a phi_c dentro del plano XY). Al ser el choque una
# superficie de revolucion en torno a u_c, la formula de Wilkin se evalua
# con este mismo angulo, sin importar en que plano caiga la celda.
cos_gamma = ux * u_c[0] + uy * u_c[1] + uz * u_c[2]
cos_gamma = np.clip(cos_gamma, -1.0, 1.0)
gamma = np.arccos(cos_gamma)

# ============================================================
# Radio del choque de Wilkin (superficie de revolucion en torno a u_c)
# ============================================================
eps = 1e-8
gamma_safe = np.maximum(gamma, eps)
cot = np.cos(gamma_safe) / np.sin(gamma_safe)

arg = 3.0 * (1.0 - gamma_safe * cot)
arg = np.clip(arg, 0.0, None)   # evita NaN mas alla del angulo maximo valido de Wilkin

Rshock = r_standoff / np.sin(gamma_safe) * np.sqrt(arg)
Rshock[gamma < 1e-4] = r_standoff

# ============================================================
# Clasificacion de regiones
# ============================================================
rho3D = np.zeros_like(R)
T3D = np.zeros_like(R)
v3D = np.zeros_like(R)

# interior del planeta
mask_planet = (R < 1.0)
rho3D[mask_planet] = n_planeta
T3D[mask_planet] = T_planeta
v3D[mask_planet] = v_planeta

# viento planetario: llena la cavidad hasta el frente del choque de Wilkin,
# en cualquier direccion (igual criterio que shock_ale_plot.py, donde
# mask_wind no depende del parche angular)
mask_wind = (R >= 1.0) & (R < (Rshock - shock_size))
rho3D[mask_wind] = interp_he_p(R[mask_wind])
T3D[mask_wind]   = interp_T_p(R[mask_wind])
v3D[mask_wind]   = interp_v_p(R[mask_wind])

# parche angular del choque: cono de semi-angulo alpha en torno a u_c
mask_angular = (gamma <= alpha)

# cascaron radial del choque (entre el frente y el fondo del choque)
mask_radial = (R >= (Rshock - shock_size)) & (R <= Rshock)

mask_shock = mask_angular & mask_radial

# radio equivalente dentro del perfil 1D del choque
r_local = (
    R[mask_shock]
    - (Rshock[mask_shock] - shock_size)
    + (r_standoff - shock_size)
)

rho3D[mask_shock] = interp_he_s(r_local)
T3D[mask_shock]   = interp_T_s(r_local)
v3D[mask_shock]   = interp_v_s(r_local)

# resto (fuera de planeta, viento y choque): vacio
mask_resto = ~(mask_planet | mask_wind | mask_shock)
rho3D[mask_resto] = 0.0
T3D[mask_resto]   = 10.0   # floor del Fortran
v3D[mask_resto]   = 0.0

# ============================================================
# Rotacion de ejes: de la convencion Python (plano orbital XY,
# observador en -Y) a la convencion de Guacho
# (plano orbital XZ, observador en +Z = LOS)
# ============================================================

def python_to_guacho_axes(arr):
    """
    Transformacion de ejes:
        x_guacho = x_python
        y_guacho = z_python
        z_guacho = -y_python
    Es una rotacion propia (ortogonal, det=+1), asi que el remapeo
    es exacto, sin interpolar.
    """
    arr_perm = np.transpose(arr, axes=(0, 2, 1))   # (x,y,z) -> (x,z,y)
    arr_guacho = np.flip(arr_perm, axis=2)          # invierte el eje que ahora es -y_python
    return np.ascontiguousarray(arr_guacho)


def write_fortran_binary(array, filename, dtype=np.float64):
    """
    Escribe un array 3D en binario stream, en orden Fortran (column-major),
    compatible con real*8 dens(nx,ny,nz), access='stream'.
    """
    arr = np.asarray(array, dtype=dtype)
    with open(filename, 'wb') as f:
        f.write(arr.tobytes(order='F'))


def read_fortran_binary(filename, nx, ny, nz, dtype=np.float64):
    """Lectura de verificacion, simulando como lo leeria Fortran."""
    data = np.fromfile(filename, dtype=dtype)
    return data.reshape((nx, ny, nz), order='F')


# --- Aplicar la rotacion a los 3 cubos ---
rho3D_g = python_to_guacho_axes(rho3D)
T3D_g   = python_to_guacho_axes(T3D)
v3D_g   = python_to_guacho_axes(v3D)

# --- Escribir los binarios ya en convencion Guacho ---
write_fortran_binary(rho3D_g, 'density.bin')
write_fortran_binary(T3D_g,   'temperature.bin')
write_fortran_binary(v3D_g,   'velocity.bin')

print("Archivos binarios escritos (convencion Guacho: plano orbital XZ, LOS=+Z):")
for fname in ['density.bin', 'temperature.bin', 'velocity.bin']:
    size_mb = os.path.getsize(fname) / 1e6
    print(f"  {fname}: {size_mb:.1f} MB")

# --- Verificacion roundtrip ---
rho_check = read_fortran_binary('density.bin', N, N, N)
print("Verificacion roundtrip (density):",
      "OK" if np.allclose(rho_check, rho3D_g) else "FALLO")
v_check = read_fortran_binary('velocity.bin', N, N, N)
print("Verificacion roundtrip (velocity):",
      "OK" if np.allclose(v_check, v3D_g) else "FALLO")
T_check = read_fortran_binary('temperature.bin', N, N, N)
print("Verificacion roundtrip (temperature):",
      "OK" if np.allclose(T_check, T3D_g) else "FALLO")

# ============================================================
# Graficos de verificacion 2D
# ============================================================
# Interpolador sobre la grilla 3D (en Python-space, antes de rotar a Guacho)
# para poder extraer cortes en planos arbitrarios.
rho_interp3D = RegularGridInterpolator((x, y, z), rho3D, bounds_error=False, fill_value=0.0)

Ngrid = 256
Lplot = L

# ------------------------------------------------------------
# Corte 1: plano XY (z=0). Debe reproducir el mapa de shock_ale_plot.py,
# ya que el eje del choque esta contenido en este plano (theta_c = 90).
# ------------------------------------------------------------
iz0 = N // 2
rho_xy = rho3D[:, :, iz0]

plt.figure(figsize=(8, 7))
plt.imshow(
    rho_xy.T,
    #origin='lower',
    extent=[-L, L, -L, L],
    norm=LogNorm(vmin=1.0e0, vmax=1.0e6)
)
plt.colorbar(label=r'$n_{\rm HeI^*}$ [cm$^{-3}$]')
plt.xlabel(r'$x$ [R$_p$]')
plt.ylabel(r'$y$ [R$_p$]')
plt.title('Corte z=0 (plano XY): debe coincidir con shock_ale_plot.py')

circ1 = plt.Circle((0, 0), 1.0, color='white', fill=False, lw=2, label='Planeta')
plt.gca().add_artist(circ1)
circ2 = plt.Circle((0, 0), r_standoff, color='k', fill=False, ls='--', lw=1.2, label=r'$r_{standoff}$')
plt.gca().add_artist(circ2)

plt.legend(loc='upper right', fontsize=8)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig('3D_check_xy_z0.png', dpi=300, bbox_inches='tight')
plt.show()

# ------------------------------------------------------------
# Corte 2: plano meridiano, que contiene al eje u_c y al eje Z.
# Al ser el choque una superficie de revolucion, este corte debe
# tener EXACTAMENTE la misma forma que el corte anterior (rotado
# a un sistema (s, z) donde s es la coordenada a lo largo de phi_c).
# ------------------------------------------------------------
s_arr = np.linspace(-Lplot, Lplot, Ngrid)
z_arr = np.linspace(-Lplot, Lplot, Ngrid)
S, Zm = np.meshgrid(s_arr, z_arr, indexing='ij')

Xm = S * np.cos(phi_c)
Ym = S * np.sin(phi_c)

pts_meridian = np.stack([Xm.ravel(), Ym.ravel(), Zm.ravel()], axis=-1)
rho_meridian = rho_interp3D(pts_meridian).reshape(Xm.shape)

plt.figure(figsize=(8, 7))
plt.imshow(
    rho_meridian.T,
    #origin='lower',
    extent=[-Lplot, Lplot, -Lplot, Lplot],
    norm=LogNorm(vmin=1.0e0, vmax=1.0e6)
)
plt.colorbar(label=r'$n_{\rm HeI^*}$ [cm$^{-3}$]')
plt.xlabel(r'$s$ [R$_p$] (a lo largo de $\phi_c$)')
plt.ylabel(r'$z$ [R$_p$]')
plt.title('Corte meridiano (plano que contiene al eje $u_c$ y al eje Z)')

circ1 = plt.Circle((0, 0), 1.0, color='white', fill=False, lw=2, label='Planeta')
plt.gca().add_artist(circ1)
circ2 = plt.Circle((0, 0), r_standoff, color='k', fill=False, ls='--', lw=1.2, label=r'$r_{standoff}$')
plt.gca().add_artist(circ2)

plt.legend(loc='upper right', fontsize=8)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig('3D_check_meridian.png', dpi=300, bbox_inches='tight')
plt.show()

# ------------------------------------------------------------
# Corte 3: plano perpendicular al eje u_c, a una distancia s0 a lo
# largo del eje. Al ser una superficie de revolucion, este corte debe
# verse como un anillo (circulo) centrado en el eje, confirmando la
# simetria azimutal alrededor de u_c.
# ------------------------------------------------------------
s0 = r_standoff  # distancia a lo largo del eje donde se corta

# base ortonormal del plano perpendicular a u_c
e1 = np.array([0.0, 0.0, 1.0])                                   # eje Z (perpendicular a u_c porque theta_c=90)
e2 = np.array([-np.sin(phi_c), np.cos(phi_c), 0.0])              # perpendicular a u_c, dentro del plano XY

a_arr = np.linspace(-Lplot, Lplot, Ngrid)
b_arr = np.linspace(-Lplot, Lplot, Ngrid)
A, B = np.meshgrid(a_arr, b_arr, indexing='ij')

Xp = s0 * u_c[0] + A * e1[0] + B * e2[0]
Yp = s0 * u_c[1] + A * e1[1] + B * e2[1]
Zp = s0 * u_c[2] + A * e1[2] + B * e2[2]

pts_perp = np.stack([Xp.ravel(), Yp.ravel(), Zp.ravel()], axis=-1)
rho_perp = rho_interp3D(pts_perp).reshape(A.shape)

plt.figure(figsize=(7, 7))
plt.imshow(
    rho_perp.T,
    #origin='lower',
    extent=[-Lplot, Lplot, -Lplot, Lplot],
    norm=LogNorm(vmin=1.0e0, vmax=1.0e6)
)
plt.colorbar(label=r'$n_{\rm HeI^*}$ [cm$^{-3}$]')
plt.xlabel(r'$a$ [R$_p$] (a lo largo de $\hat z$)')
plt.ylabel(r'$b$ [R$_p$] (perpendicular a $u_c$ en el plano XY)')
plt.title(f'Corte perpendicular a $u_c$ en s={s0:.2f} $R_p$\n(debe verse como un anillo)')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig('3D_check_perp_ring.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(rho3D_g[N // 2, N // 2, :])
plt.xlabel('x')
plt.ylabel('rho')
plt.yscale('log')
plt.savefig('rho_1D.png')
plt.show()