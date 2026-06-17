import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.ticker import LogFormatterSciNotation

# =========================
# Leer lista de modelos
# =========================

path = "../output/"

models_data = np.loadtxt(
    "../../shock_1D/output_1/model_list.dat",
    comments="#"
)

model = models_data[:, 0]
T0    = models_data[:, 1]
n0    = models_data[:, 2]
u0    = models_data[:, 3] / 1e5   # km/s
y0    = models_data[:, 4]

models_num = len(model)
HeIM_index = 4   # columna: r HI HII HeIS HeIM HeII e-

HeIM_trapz = []

# =========================
# Leer outputs
# =========================

for i in range(models_num):
    filename = f"{path}final_model_{int(model[i])}.dat"
    print(f"Procesando: {filename}")

    data = np.loadtxt(filename, comments="#")
    r    = data[:, 0]
    HeIM = data[:, HeIM_index]

    # trapezoid pesa cada punto por el espaciado real dr -> físicamente correcto
    HeIM_trapz.append(np.trapezoid(HeIM, r))

HeIM_trapz = np.array(HeIM_trapz)

# =========================
# Parámetros y etiquetas
# =========================

PARAMS = {
    'T0': {'values': T0,  'label': r'$T_0$ [K]',          'log': True},
    'n0': {'values': n0,  'label': r'$n_0$ [1/cm3]',  'log': True},
    'u0': {'values': u0,  'label': r'$u_0$ [km/s]',        'log': False},
    'y0': {'values': y0,  'label': r'$y_0$',               'log': True},
}
param_keys = list(PARAMS.keys())
N = len(param_keys)

# =========================
# Corner plot
# =========================

def make_corner_plot(metric_values, metric_name, cmap='viridis', output_file='corner_HeIM.png'):

    norm_factor = np.max(metric_values)
    z = metric_values / norm_factor

    best_idx = np.argmax(metric_values)
    best_params = {k: PARAMS[k]['values'][best_idx] for k in param_keys}

    fig = plt.figure(figsize=(12, 12))
    gs  = gridspec.GridSpec(N, N, figure=fig, hspace=0.05, wspace=0.05)

    axes = {}

    for row in range(N):
        for col in range(N):
            if col > row:
                continue   # solo triángulo inferior + diagonal

            ax = fig.add_subplot(gs[row, col])
            axes[(row, col)] = ax
            prow = param_keys[row]
            pcol = param_keys[col]

            prow_vals = PARAMS[prow]['values']
            pcol_vals = PARAMS[pcol]['values']

            # ── Diagonal: distribución marginal ──────────────────────────────
            if row == col:
                unique_vals = np.sort(np.unique(pcol_vals))
                margin = np.array([
                    np.sum(z[pcol_vals == v]) for v in unique_vals
                ])
                margin /= margin.max()

                ax.bar(range(len(unique_vals)), margin, color=plt.colormaps[cmap](margin))
                ax.set_xticks(range(len(unique_vals)))
                ax.set_xticklabels(
                    [f'{v:.1e}' if PARAMS[pcol]['log'] else f'{v:.0f}' for v in unique_vals],
                    fontsize=13, rotation=45
                )
                ax.set_ylim(0, 1.1)
                ax.yaxis.set_visible(False)
                ax.set_facecolor('#f0f0f0')

                # Marcar el mejor valor
                best_v = best_params[pcol]
                best_xi = np.where(unique_vals == best_v)[0]
                if len(best_xi):
                    ax.axvline(best_xi[0], color='crimson', lw=1.5, ls='--', alpha=0.8)

            # ── Off-diagonal: heatmap 2D ──────────────────────────────────────
            else:
                x_unique = np.sort(np.unique(pcol_vals))
                y_unique = np.sort(np.unique(prow_vals))
                nx, ny = len(x_unique), len(y_unique)

                Z = np.zeros((ny, nx))
                x_idx = {v: i for i, v in enumerate(x_unique)}
                y_idx = {v: i for i, v in enumerate(y_unique)}

                for i in range(models_num):
                    xi = x_idx[pcol_vals[i]]
                    yi = y_idx[prow_vals[i]]
                    Z[yi, xi] += z[i]

                Z /= (Z.max() + 1e-30)

                im = ax.imshow(
                    Z, origin='lower', aspect='auto',
                    cmap=cmap, vmin=0, vmax=1,
                    extent=[-0.5, nx-0.5, -0.5, ny-0.5]
                )

                # Contornos (requiere al menos (2,2))
                #if Z.max() > 0 and nx >= 2 and ny >= 2:
                #    ax.contour(Z, levels=[0.5, 0.9], colors='white',
                #               linewidths=0.7, alpha=0.6,
                #               extent=[-0.5, nx-0.5, -0.5, ny-0.5])

                # Marcar el mejor modelo
                bxi = x_idx.get(best_params[pcol])
                byi = y_idx.get(best_params[prow])
                if bxi is not None and byi is not None:
                    ax.plot(bxi, byi, '*', color='crimson', ms=8, zorder=5,
                            markeredgecolor='white', markeredgewidth=0.5)

                ax.set_xticks(range(nx))
                ax.set_yticks(range(ny))

                x_labels = [f'{v:.0e}' if PARAMS[pcol]['log'] else f'{v:.0f}' for v in x_unique]
                y_labels = [f'{v:.0e}' if PARAMS[prow]['log'] else f'{v:.0f}' for v in y_unique]

                if row == N - 1:
                    ax.set_xticklabels(x_labels, fontsize=13, rotation=45)
                else:
                    ax.set_xticklabels([])

                if col == 0:
                    ax.set_yticklabels(y_labels, fontsize=13)
                else:
                    ax.set_yticklabels([])

            # Etiquetas de ejes exteriores
            if col == 0 and row != col:
                ax.set_ylabel(PARAMS[prow]['label'], fontsize=14, labelpad=4)
            if row == N - 1:
                ax.set_xlabel(PARAMS[pcol]['label'], fontsize=14, labelpad=4)

    # Colorbar global
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.5])
    sm  = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(f'HeIM {metric_name} (norm.)', fontsize=14)

    # Título con mejor modelo
    title_parts = [f"{PARAMS[k]['label'].replace('$','').strip()} = {best_params[k]:.2e}" for k in param_keys]
    fig.suptitle(
        f'Corner Plot — HeIM ({metric_name})\n'
        f'Mejor modelo: {", ".join(title_parts)}',
        fontsize=15, y=0.98
    )

    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"Guardado: {output_file}")
    plt.show()

# =========================
# Generar ambos gráficos
# =========================

make_corner_plot(HeIM_trapz, metric_name='trapz', cmap='viridis', output_file='corner_HeIM_trapz.png')