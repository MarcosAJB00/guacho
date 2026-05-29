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
    HeIM_trapz.append(np.trapz(HeIM, r))

HeIM_sum = np.array(HeIM_sum)
HeIM_trapz = np.array(HeIM_trapz)

# =====================================================
# 1) Gráfico original (modelo vs HeIM total)
# =====================================================

plt.figure(figsize=(8,6))

plt.plot(model, HeIM_sum, marker='o', color='blue', label='HeIM total (sum)')
plt.plot(model, HeIM_trapz, marker='s', color='red', label='HeIM total (trapz)')

plt.xlabel('Modelo')
plt.ylabel('HeIM total')
plt.grid()
plt.legend()

index_max = np.argmax(HeIM_sum)

T0_max = T0[index_max]
n0_max = n0[index_max]
u0_max = u0[index_max]
y0_max = y0[index_max]

plt.title(
    f'Best: Model {int(model[index_max])}, '
    f'T0={T0_max:.1e}, '
    f'n0={n0_max:.1e}, '
    f'u0={u0_max:.1f} km/s, '
    f'y0={y0_max:.1e}'
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
# 2) Heatmap usando SUM
# =====================================================

plt.figure(figsize=(8,6))

im = plt.imshow(
    Z_sum,
    origin='lower',
    aspect='auto',
    extent=[
        u0_unique.min(),
        u0_unique.max(),
        n0_unique.min(),
        n0_unique.max()
    ]
)

plt.colorbar(im, label='HeIM total (sum)')

plt.xlabel(r'$u_0$ [km/s]')
plt.ylabel(r'$n_0$ [cm$^{-3}$]')

plt.title('Grid de modelos - SUM')

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
    Z_trapz,
    origin='lower',
    aspect='auto',
    extent=[
        u0_unique.min(),
        u0_unique.max(),
        n0_unique.min(),
        n0_unique.max()
    ]
)

plt.colorbar(im, label='HeIM total (trapz)')

plt.xlabel(r'$u_0$ [km/s]')
plt.ylabel(r'$n_0$ [cm$^{-3}$]')

plt.title('Grid de modelos - TRAPEZOID')

plt.savefig(
    "grid_HeIM_trapz.png",
    dpi=300,
    bbox_inches='tight'
)

plt.show()