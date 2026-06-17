import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ====================================================================
# 1. PARÁMETROS DEL SISTEMA (Podés modificarlos)
# ====================================================================
R_planet = 1.0        # Radio del planeta (unidades arbitrarias, ej. radios de Júpiter)
R_hill = 5.0          # Radio de Hill
R_shock_end = 7.0     # Radio donde termina el parche de choque
delta_phi = 30.0      # Semi-ángulo del parche de choque (en grados)

# Orientación del centro del parche (vector unitario)
# Supongamos que el choque apunta hacia la estrella (dirección X)
shock_dir = np.array([1.0, 0.0, 0.0]) 

# Parámetros de la grilla 3D
N_grid = 150          # Resolución de la grilla (N_grid x N_grid x N_grid)
L_box = R_shock_end * 1.1 # Tamaño de la caja (un poco más grande que el choque)

# ====================================================================
# 2. CARGA DE TUS PERFILES 1D
# ====================================================================
# Aquí es donde cargarías tus archivos. Ejemplo:
# r_wind, rho_wind = np.loadtxt('perfil_viento.txt', unpack=True)
# r_shock, rho_shock = np.loadtxt('perfil_choque.txt', unpack=True)

# ---> PARA ESTE EJEMPLO, GENERAMOS PERFILES DE PRUEBA <---
r_wind = np.linspace(R_planet, R_hill * 1.5, 100)
# Perfil de viento que decae con r^-2
rho_wind = 1.0 * (R_planet / r_wind)**2  

r_shock = np.linspace(R_hill, R_shock_end * 1.2, 100)
# Perfil del choque (ejemplo: densidad alta en el radio de Hill que decae)
rho_shock = 3.0 * (R_hill / r_shock)**3  

# Creamos funciones de interpolación
# bounds_error=False y fill_value=0 asume densidad 0 fuera del rango definido
interp_wind = interp1d(r_wind, rho_wind, bounds_error=False, fill_value=0.0)
interp_shock = interp1d(r_shock, rho_shock, bounds_error=False, fill_value=0.0)

# ====================================================================
# 3. CREACIÓN DE LA GRILLA 3D
# ====================================================================
print("Generando la grilla 3D...")
x = np.linspace(-L_box, L_box, N_grid)
y = np.linspace(-L_box, L_box, N_grid)
z = np.linspace(-L_box, L_box, N_grid)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Calculamos la distancia radial (r) para cada punto de la grilla
R_grid = np.sqrt(X**2 + Y**2 + Z**2)

# Calculamos el ángulo entre cada punto y la dirección del choque
# Producto punto: A dot B = |A||B| cos(theta)
dot_product = X*shock_dir[0] + Y*shock_dir[1] + Z*shock_dir[2]
cos_alpha = dot_product / (R_grid + 1e-10) # Evitamos división por cero
alpha_grid = np.degrees(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))

# ====================================================================
# 4. ASIGNACIÓN DE DENSIDADES (VIENTO Y CHOQUE)
# ====================================================================
print("Calculando densidades...")
Density_3D = np.zeros_like(R_grid)

# A. Llenamos con el viento planetario (entre R_planet y R_hill)
mask_wind = (R_grid >= R_planet) & (R_grid <= R_hill)
Density_3D[mask_wind] = interp_wind(R_grid[mask_wind])

# B. Agregamos el parche de choque
# Condición: Dentro del ángulo delta_phi Y entre R_hill y R_shock_end
mask_shock = (alpha_grid <= delta_phi) & (R_grid >= R_hill) & (R_grid <= R_shock_end)

# Pisamos (o sumamos, según la física de tu modelo) la densidad en el parche
Density_3D[mask_shock] = interp_shock(R_grid[mask_shock])

# Opcional: Anulamos la densidad dentro del planeta sólido
Density_3D[R_grid < R_planet] = 0.0

# ====================================================================
# 5. PROYECCIÓN EN EL PLANO DEL CIELO (TRÁNSITO)
# ====================================================================
print("Proyectando a lo largo de la línea de visión...")
# Asumimos que observamos desde el eje Y (mirando hacia el plano X-Z)
# Integramos la densidad a lo largo del eje Y (índice 1 de la matriz)
dy = y[1] - y[0]
Column_Density = np.sum(Density_3D, axis=1) * dy

# Creamos una máscara circular para simular el cuerpo sólido y opaco del planeta
X_2D, Z_2D = np.meshgrid(x, z, indexing='ij')
R_2D = np.sqrt(X_2D**2 + Z_2D**2)
planet_mask = R_2D <= R_planet

# ====================================================================
# 6. VISUALIZACIÓN
# ====================================================================
plt.figure(figsize=(8, 6))

# Dibujamos la densidad de columna (escala logarítmica suele ser útil, acá usamos lineal)
# Rotamos la matriz para que coincidan los ejes X y Z con la imagen
img = plt.imshow(Column_Density.T, extent=[-L_box, L_box, -L_box, L_box], 
                 origin='lower', cmap='viridis')

# Superponemos el planeta sólido en color negro
plt.contourf(X_2D, Z_2D, planet_mask, levels=[0.5, 1.5], colors=['black'])

plt.colorbar(img, label='Densidad de Columna (Unidades arbitrarias)')
plt.xlabel('Distancia en X (Ej. Dirección a la estrella)')
plt.ylabel('Distancia en Z')
plt.title('Proyección de Viento Planetario + Parche de Choque')

# Marcamos el contorno del Radio de Hill como referencia
circle_hill = plt.Circle((0, 0), R_hill, color='white', fill=False, linestyle='--', alpha=0.5, label='Radio de Hill')
plt.gca().add_patch(circle_hill)
plt.legend()

plt.show()