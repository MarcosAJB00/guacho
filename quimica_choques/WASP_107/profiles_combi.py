import numpy as np
import matplotlib.pyplot as plt

m_p = 1.67e-24  # g
mp =1.3 * m_p  # masa media por partícula (g)
# ============================================================
# STAND-OFF Y CONDICIONES INICIALES DEL CHOQUE
# shape: (n_puntos, ncols) con unpack=False (default)
# ============================================================
data_stand_off = np.loadtxt('./standoff_point/shocks_ci.dat')
r_stand_off  = data_stand_off[:, 0]   # [Rp]
rho_stand_off = data_stand_off[:, 1]*mp  # [g/cm³]
temp_stand_off = data_stand_off[:, 2] # [K]
v_stand_off  = data_stand_off[:, 3]   # [cm/s]

# ============================================================
# LISTA DE MODELOS DE CHOQUE
# ============================================================
shock_models = np.loadtxt('./shock_1D/output/model_list.dat')
model_number = shock_models[:, 0]
shock_temp   = shock_models[:, 1]   # [K]
n0_shock     = shock_models[:, 2]   # [cm⁻³]
v_shock      = shock_models[:, 3]   # [cm/s]
y0_shock     = shock_models[:, 4]

n_models = len(model_number)

# ============================================================
# QUÍMICA DE LOS CHOQUES
# Cada archivo puede tener distinto número de filas,
# así que guardamos listas en vez de arrays 2D fijos.
# ============================================================
path = './chemeq2_shock/output/'

r_HeIM  = []   # lista de arrays, uno por modelo
X_HeIM  = []
Rjup = 7.1492e9  # cm
Rp = 0.4466*Rjup
for i in range(n_models):
    fname = path + 'final_model_' + str(int(model_number[i])) + '.dat'
    data_chem = np.loadtxt(fname)          # shape (nrows, ncols), sin unpack
    r_HeIM.append(data_chem[:, 0]/Rp)         # [Rp] — radio desde el frente de choque
    X_HeIM.append(data_chem[:, 4])         # [cm⁻³]

# ============================================================
# PERFIL SIN CHOQUE
# ============================================================
sin_choque = np.loadtxt('./ATES/quimica/output/perfiles_finales.dat')
r_wos    = sin_choque[:, 0]   # [Rp]
HeIM_wos = sin_choque[:, 4]   # [cm⁻³]

data_p = np.loadtxt('./ATES/Hydro_ioniz.txt', comments='#')


k_B = 1.38e-16  # erg/K

rp    = data_p[:, 0] * Rp   # [Rp] → [cm], origen en el planeta
rho_p = data_p[:, 1]*mp            # [g/cm³]
v_p   = data_p[:, 2]            # [cm/s]
T_p   = data_p[:, 4]            # [K]

cs_p = np.sqrt(k_B * T_p / mp) # velocidad del sonido (cm/s)
r_sonic = rp[np.argmin(np.abs(v_p - cs_p))]
print(f"Radio sónica: {r_sonic/Rp:.3f} Rp")
plt.figure(figsize=(8, 5))
plt.plot(rp/Rp, v_p/1e5, label='Velocidad')
plt.plot(rp/Rp, cs_p/1e5, label='Velocidad del sonido')
plt.scatter(r_sonic/Rp, v_p[np.argmin(np.abs(v_p - cs_p))]/1e5, color='red', zorder=5, label=f'Punto sónica = {r_sonic/Rp:.3f} Rp')
plt.yscale('log')
plt.xlabel('Distancia desde el planeta [R_p]')
plt.ylabel('Velocidad [Km/s]')
plt.title('Perfil de velocidad del viento planetario')
plt.legend()
plt.savefig('planetary_velocity_profile.png', dpi=300, bbox_inches='tight')

# Integral del perfil sin choque (referencia, modelo 0)
HeIM_trap = []
HeIM_trap.append(np.trapezoid(HeIM_wos, r_wos))

# ============================================================
# COMBINAR PERFILES Y CALCULAR INTEGRALES
#
# Para cada modelo i:
#   - tomamos el perfil sin choque desde r=1.05 Rp hasta r_stand_off[i]
#   - a partir de r_stand_off[i] pegamos el perfil del choque,
#     desplazando su eje r para que arranque en r_stand_off[i]
# ============================================================
fig1, ax1 = plt.subplots(figsize=(9, 6))
ax1.plot(r_wos, HeIM_wos, color='k', lw=2, label='Sin choque')
ax1.scatter(r_sonic/Rp, HeIM_wos[np.argmin(np.abs(r_wos - r_sonic/Rp))], color='red', zorder=5, label=f'Punto sónico (sin choque) = {r_sonic/Rp:.3f} Rp')

r_comb_list    = []
HeIM_comb_list = []

for i in range(n_models):

    rs = r_stand_off[i]   # radio de stand-off del modelo i [Rp]

    # Parte del perfil sin choque: desde 1.05 Rp hasta el stand-off
    mask = (r_wos >= 1.0001) & (r_wos <= rs)
    r_part1    = r_wos[mask]
    HeIM_part1 = HeIM_wos[mask]

    # Parte del choque: el radio del archivo empieza en 0 (desde el frente),
    # lo desplazamos para que en el gráfico arranque en rs
    r_part2    = rs + r_HeIM[i]
    HeIM_part2 = X_HeIM[i]

    r_comb    = np.concatenate([r_part1, r_part2])
    HeIM_comb = np.concatenate([HeIM_part1, HeIM_part2])

    r_comb_list.append(r_comb)
    HeIM_comb_list.append(HeIM_comb)

    HeIM_trap.append(np.trapezoid(HeIM_comb, r_comb))

    ax1.plot(r_comb, HeIM_comb,
             label=f'Modelo {int(model_number[i])} (rs={rs:.2f} Rp)',
             alpha=0.7)

ax1.set_xlabel(r'Radio [$R_p$]')
ax1.set_ylabel(r'Densidad HeIM [cm$^{-3}$]')
ax1.set_yscale('log')
ax1.set_ylim(10.0,5e4)
ax1.set_title('Perfiles de HeIM con y sin choque')
ax1.grid(alpha=0.3, ls='--', lw=0.5, which='both')
ax1.legend(fontsize=9,loc='lower left')
plt.tight_layout()
plt.savefig('comparacion_HeIM_profiles.png', dpi=300, bbox_inches='tight')
plt.show()

# ============================================================
# GRÁFICO DE INTEGRALES
# índice 0 = sin choque, 1..n = modelos con choque
# ============================================================
HeIM_trap = np.array(HeIM_trap)
labels    = ['Sin choque'] + [f'M{int(m)}' for m in model_number]
x         = np.arange(len(HeIM_trap))

fig2, ax2 = plt.subplots(figsize=(9, 5))
ax2.scatter(x, HeIM_trap, color='red', zorder=3)
ax2.set_xticks(x)
ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
ax2.set_ylabel('HeIM')
ax2.set_title('Trapezoidal HeIM por modelo')
ax2.grid(alpha=0.3, ls='--', lw=0.5, which='both')
plt.tight_layout()
plt.savefig('comparacion_HeIM_trap.png', dpi=300, bbox_inches='tight')
plt.show()