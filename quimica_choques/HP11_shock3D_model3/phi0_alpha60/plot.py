import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d, RegularGridInterpolator
import os

Rjup = 7.1492e9  # cm
Rp = 0.4466*Rjup
r_standoff = 4.77       # en Rp

n_planeta = 1e8       # densidad dentro del planeta
T_planeta = 847.0  #
v_planeta = 0.0

# PERFIL 1D PLANETA
planeta = np.loadtxt('perfiles_finales.dat',comments='#')
planet_T_profile = np.loadtxt('temperature_profile_ATES.dat',comments='#')
planet_v_profile = np.loadtxt('velocity_profile_ATES.dat',comments='#')

r_p = planeta[:,0] # en Rp
HeIM_p = planeta[:,4] # en 1/cm3
T_p = planet_T_profile[:,1] # en K
v_p = planet_v_profile[:,1] # en Km/s
v_p = v_p*1e5 #de Km/s a cm/s

# interpolador
interp_he_p = interp1d(r_p,HeIM_p, bounds_error=False,fill_value=0.0)
interp_T_p = interp1d(r_p,T_p, bounds_error=False,fill_value=0.0)
interp_v_p = interp1d(r_p,v_p, bounds_error=False,fill_value=0.0)

# PERFIL 1D SHOCK
shock = np.loadtxt('final_model_3.dat',comments='#')
shock_T_v = np.loadtxt('model_3.dat',comments='#')

r_s = shock[:,0] # en cm
shock_size = r_s[-1]/Rp #tamaño del choque en Rp
r_s = r_s/Rp + r_standoff - shock_size  # en Rp, le sumo el r_standoff para que no arranque en 0
HeIM_s = shock[:,4] #en 1/cm3
T_s = shock_T_v[:,2] #en K
v_s = shock_T_v[:,3] #en cm/s

# interpolador
interp_he_s = interp1d(r_s,HeIM_s,bounds_error=False,fill_value=0.0)
interp_T_s = interp1d(r_s,T_s,bounds_error=False,fill_value=0.0)
interp_v_s = interp1d(r_s,v_s,bounds_error=False,fill_value=0.0)

# Grilla 3D
N = 256
L = r_s[-1]

x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
z = np.linspace(-L, L, N)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
R = np.sqrt(X**2 + Y**2 + Z**2)

rho3D = np.zeros_like(R)
T3D = np.zeros_like(R)
v3D = np.zeros_like(R)

# interior del planeta
mask_planet = (R < 1.0)

rho3D[mask_planet] = n_planeta
T3D[mask_planet] = T_planeta
v3D[mask_planet] = v_planeta

# viento planetario
mask_wind = (R >= 1.0) & (R <= (r_standoff - shock_size))
rho3D[mask_wind] = interp_he_p(R[mask_wind])
T3D[mask_wind] = interp_T_p(R[mask_wind])
v3D[mask_wind] = interp_v_p(R[mask_wind])


# PARCHE DE SHOCK
#coordenadas del centro del parche
theta_c = np.deg2rad(90.0) #coordenada polar (90° es el ecuador)
phi_c   = np.deg2rad(0.0) # coordenda azimutal (0° apunta hacia x positivos)

alpha   = np.deg2rad(60.0) #tamaño angular radio del parche del parche

#versor radial en la direccion del parche
u_c = np.array([
    np.sin(theta_c)*np.cos(phi_c),
    np.sin(theta_c)*np.sin(phi_c),
    np.cos(theta_c)
])

print('u_c = ',u_c)
# versor de todas las celdas
ux = np.zeros_like(R)
uy = np.zeros_like(R)
uz = np.zeros_like(R)
mask_nonzero = (R > 0)
ux[mask_nonzero] = X[mask_nonzero]/R[mask_nonzero]
uy[mask_nonzero] = Y[mask_nonzero]/R[mask_nonzero]
uz[mask_nonzero] = Z[mask_nonzero]/R[mask_nonzero]

#distancia angular al centro del parche
cos_gamma = (ux*u_c[0]+ uy*u_c[1]+ uz*u_c[2])
cos_gamma = np.clip(cos_gamma,-1.0,1.0)
gamma = np.arccos(cos_gamma)

#mascara angular del shock
mask_angular = (gamma <= alpha) #me da un parche circular

#mascara radial del shock
mask_radial = ((R >= (r_standoff - shock_size)) & (R <= L))

#mascara final del shock
mask_shock = (mask_angular & mask_radial)

#rellenar perfil del shock 
rho3D[mask_shock] = interp_he_s(R[mask_shock])
T3D[mask_shock] = interp_T_s(R[mask_shock])
v3D[mask_shock] = interp_v_s(R[mask_shock])

# todo lo que no es planeta, ni viento, ni shock
mask_resto = ~(mask_planet | mask_wind | mask_shock)

rho3D[mask_resto] = 0.0   # densidad insignificante
T3D[mask_resto]   = 10.0     # temperatura mínima (floor del Fortran es 10 K)
v3D[mask_resto]   = 0.0     # sin velocidad
# ============================================================
# Rotación de ejes: de tu convención (plano orbital XY, 
# observador en -Y) a la convención de Guacho 
# (plano orbital XZ, observador en +Z = LOS)
# ============================================================

def python_to_guacho_axes(arr):
    """
    Transformación de ejes:
        x_guacho = x_python
        y_guacho = z_python
        z_guacho = -y_python
    Es una rotación propia (ortogonal, det=+1), así que el remapeo
    es exacto, sin interpolar.
    """
    arr_perm = np.transpose(arr, axes=(0, 2, 1))   # (x,y,z) -> (x,z,y)
    arr_guacho = np.flip(arr_perm, axis=2)          # invierte el eje que ahora es -y_python
    return np.ascontiguousarray(arr_guacho)

def write_fortran_binary(array, filename, dtype=np.float64):
    """
    Escribe un array 3D en binario stream, en orden Fortran (column-major),
    compatible con real*8 dens(nx,ny,nz), access='stream'.
    OJO: tofile() de numpy SIEMPRE escribe en C-order, por eso usamos
    tobytes(order='F') explícitamente.
    """
    arr = np.asarray(array, dtype=dtype)
    with open(filename, 'wb') as f:
        f.write(arr.tobytes(order='F'))

def read_fortran_binary(filename, nx, ny, nz, dtype=np.float64):
    """Lectura de verificación, simulando cómo lo leería Fortran."""
    data = np.fromfile(filename, dtype=dtype)
    return data.reshape((nx, ny, nz), order='F')


# --- Aplicar la rotación a los 3 cubos ---
rho3D_g = python_to_guacho_axes(rho3D)
T3D_g   = python_to_guacho_axes(T3D)
v3D_g   = python_to_guacho_axes(v3D)

# --- Escribir los binarios ya en convención Guacho ---
write_fortran_binary(rho3D_g, 'density.bin')
write_fortran_binary(T3D_g,   'temperature.bin')
write_fortran_binary(v3D_g,   'velocity.bin')

print("Archivos binarios escritos (convención Guacho: plano orbital XZ, LOS=+Z):")
for fname in ['density.bin', 'temperature.bin', 'velocity.bin']:
    size_mb = os.path.getsize(fname) / 1e6
    print(f"  {fname}: {size_mb:.1f} MB")

# --- Verificación roundtrip ---
rho_check = read_fortran_binary('density.bin', N, N, N)
print("Verificación roundtrip (density):",
      "OK" if np.allclose(rho_check, rho3D_g) else "FALLO")
# --- Verificación roundtrip ---
rho_check = read_fortran_binary('velocity.bin', N, N, N)
print("Verificación roundtrip (velocity):",
      "OK" if np.allclose(rho_check, v3D_g) else "FALLO")
# --- Verificación roundtrip ---
rho_check = read_fortran_binary('temperature.bin', N, N, N)
print("Verificación roundtrip (temperature):",
      "OK" if np.allclose(rho_check, T3D_g) else "FALLO")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# índice del plano central
iz = N// 2
iy = N // 2
ix = N // 2

plt.figure(figsize=(8,7))

plt.imshow(
    rho3D[:, :, iy].T,
    #origin='lower',
    extent=[-L, L, -L, L],
    norm=LogNorm(
        vmin=1e-25,#rho3D.min(),
        vmax= rho3D.max()
    )
)

plt.colorbar(label=r'$n_{\rm HeI^*}$ [cm$^{-3}$]')

plt.xlabel(r'$x$ [R$_p$]')
plt.ylabel(r'$z$ [R$_p$]')
plt.title('Perfil 1D proyectado a esfera 3D')

# planeta
circ1 = plt.Circle((0,0),1.0,color='white',
                   fill=False,lw=2, label='Planet')
plt.gca().add_artist(circ1)

# standoff
circ2 = plt.Circle((0,0),r_standoff,color='cyan',
                   fill=False,ls='--',lw=2, label='r_standoff')
plt.gca().add_artist(circ2)

# shock
circ3 = plt.Circle((0,0),r_s[-1],color='k',
                   fill=False,ls='--',lw=2,label='Shock')
plt.gca().add_artist(circ3)

plt.tight_layout()
#plt.legend()
plt.savefig('2D_xy.png',dpi=300,bbox_inches='tight')
plt.show()


plt.figure()

plt.plot(rho3D_g[N//2,N//2,:])
plt.xlabel('x')
plt.ylabel('rho')
plt.yscale('log')
plt.savefig('rho_1D.png')
