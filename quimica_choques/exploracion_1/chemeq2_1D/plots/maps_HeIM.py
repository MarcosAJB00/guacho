import numpy as np
import matplotlib.pyplot as plt

# =========================
# Leer lista de modelos
# =========================

path = "../output/"

models = np.loadtxt(
    "../../shock_1D/output/model_list.dat",
    comments="#"
)

model = models[:, 0]
T0 = models[:, 1]
n0 = models[:, 2]
u0 = models[:, 3] / 1e5   # km/s
y0 = models[:, 4]

models_num = len(model)

# Columna de HeIM
# r HI HII HeIS HeIM HeII e-
HeIM_index = 4

HeIM_sum = []
HeIM_trapz = []

# =========================
# Leer outputs
# =========================

for i in range(models_num):

    filename = f"{path}final_model_{int(model[i])}.dat"

    print(f"Procesando archivo: {filename}")

    data = np.loadtxt(filename, comments="#")

    r = data[:, 0]
    HeIM = data[:, HeIM_index]

    # 1) suma simple
    HeIM_sum.append(np.sum(HeIM))

    # 2) integral trapezoidal
    HeIM_trapz.append(np.trapezoid(HeIM, r))

HeIM_sum = np.array(HeIM_sum)
HeIM_trapz = np.array(HeIM_trapz)

# =====================================================
# 1) Gráfico original (modelo vs HeIM total)
# =====================================================

plt.figure(figsize=(8,6))

norm_sum = np.sum(HeIM_sum)
norm_trapz = np.sum(HeIM_trapz)
plt.plot(model, HeIM_sum/norm_sum, marker='o', color='blue', label='HeIM total (sum)')
plt.plot(model, HeIM_trapz/norm_trapz, marker='s', color='red', label='HeIM total (trapz)')

plt.xlabel('Modelo')
plt.ylabel('HeIM total')
plt.grid()
plt.legend()

index_max_sum = np.argmax(HeIM_sum)
index_max_trapz = np.argmax(HeIM_trapz)

T0_max_sum = T0[index_max_sum]
n0_max_sum = n0[index_max_sum]
u0_max_sum = u0[index_max_sum]
y0_max_sum = y0[index_max_sum]

T0_max_trapz = T0[index_max_trapz]
n0_max_trapz = n0[index_max_trapz]
u0_max_trapz = u0[index_max_trapz]
y0_max_trapz = y0[index_max_trapz]

plt.title(
    f'Best (sum): Model {int(model[index_max_sum])}, '
    f'T0={T0_max_sum:.1e}, '
    f'n0={n0_max_sum:.1e}, '
    f'u0={u0_max_sum:.1f} km/s, '
    f'y0={y0_max_sum:.1e}\n'
    f'Best (trapz): Model {int(model[index_max_trapz])}, '
    f'T0={T0_max_trapz:.1e}, '
    f'n0={n0_max_trapz:.1e}, '
    f'u0={u0_max_trapz:.1f} km/s, '
    f'y0={y0_max_trapz:.1e}' 
)

plt.savefig(
    "dens_HeIM.png",
    dpi=300,
    bbox_inches='tight'
)

plt.show()

# =====================================================
# Construcción de grilla
# =====================================================

u0_unique = np.sort(np.unique(u0))
n0_unique = np.sort(np.unique(n0))

Z_sum = np.zeros((len(n0_unique), len(u0_unique)))
Z_trapz = np.zeros((len(n0_unique), len(u0_unique)))

for i in range(models_num):

    ix = np.where(u0_unique == u0[i])[0][0]
    iy = np.where(n0_unique == n0[i])[0][0]

    Z_sum[iy, ix] = HeIM_sum[i]
    Z_trapz[iy, ix] = HeIM_trapz[i]

# =====================================================
# Heatmap SUM
# =====================================================

plt.figure(figsize=(8,6))

im = plt.imshow(
    Z_sum/norm_sum,
    origin='lower',
    aspect='auto'
)

plt.colorbar(im, label='HeIM total (sum)')

plt.xticks(
    np.arange(len(u0_unique)),
    u0_unique.astype(int)
)

plt.yticks(
    np.arange(len(n0_unique)),
    [f'{x:.0e}' for x in n0_unique]
)

plt.xlabel(r'$u_0$ [km/s]')
plt.ylabel(r'$n_0$ [cm$^{-3}$]')

plt.title(
    f'Best model: {int(model[index_max_sum])}, '
    f'n0={n0[index_max_sum]:.1e}, '
    f'u0={u0[index_max_sum]:.1f}'
)

plt.savefig(
    "grid_HeIM_sum.png",
    dpi=300,
    bbox_inches='tight'
)

plt.show()
# =====================================================
# 3) Heatmap usando TRAPEZOID
# =====================================================

plt.figure(figsize=(8,6))

im = plt.imshow(
    Z_trapz/norm_trapz,
    origin='lower',
    aspect='auto'
)

plt.colorbar(im, label='HeIM total (trapz)')

plt.xticks(
    np.arange(len(u0_unique)),
    u0_unique.astype(int)
)

plt.yticks(
    np.arange(len(n0_unique)),
    [f'{x:.0e}' for x in n0_unique]
)

plt.xlabel(r'$u_0$ [km/s]')
plt.ylabel(r'$n_0$ [cm$^{-3}$]')

plt.title(
    f'Best model: {int(model[index_max_trapz])}, '
    f'n0={n0[index_max_trapz]:.1e}, '
    f'u0={u0[index_max_trapz]:.1f}'
)

plt.savefig(
    "grid_HeIM_trapz.png",
    dpi=300,
    bbox_inches='tight'
)

plt.show()