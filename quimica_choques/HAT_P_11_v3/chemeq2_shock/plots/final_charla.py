import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (
    AutoMinorLocator,
    LogLocator,
    MaxNLocator,
    NullFormatter
)

# =====================================================
# Constants
# =====================================================

Rjup = 7.1492e9  # cm

# =====================================================
# Load model 3
# =====================================================

model = 3
ns = 6

Species_names = [
    r"$\mathrm{H\,I}$",
    r"$\mathrm{H\,II}$",
    r"$\mathrm{He\,I}\;(1^{1}S)$",
    r"$\mathrm{He\,I}\;(2^{3}S)$",
    r"$\mathrm{He\,II}$",
    r"$e^{-}$"
]

data = np.loadtxt(
    f"../output/final_model_{model}.dat",
    comments="#"
)

r = data[:,0]/Rjup
X = data[:,1:1+ns]

# =====================================================
# Plot style
# =====================================================

plt.rcParams.update({

    "font.size": 15,

    "axes.linewidth": 1.4,

    "xtick.direction": "in",
    "ytick.direction": "in",

    "xtick.top": True,
    "ytick.right": True,

    "xtick.major.size": 8,
    "ytick.major.size": 8,

    "xtick.major.width": 1.3,
    "ytick.major.width": 1.3,

    "xtick.minor.size": 4,
    "ytick.minor.size": 4,

    "xtick.minor.width": 1.0,
    "ytick.minor.width": 1.0,
})

# =====================================================
# Colors (tab10 es muy sobrio)
# =====================================================

colors = plt.get_cmap("tab10").colors

fig, ax = plt.subplots(figsize=(7.2,5.2))

for j in range(ns):

    ax.plot(
        r,
        X[:,j],
        lw=2.2,
        color=colors[j],
        label=Species_names[j]
    )

# -----------------------------------------------------

ax.set_yscale("log")

ax.set_xlabel(r"Distance [$R_{\rm J}$]", fontsize=20)
ax.set_ylabel(r"Number density [cm$^{-3}$]",fontsize=20)

# -----------------------------------------------------
# Ticks
# -----------------------------------------------------

ax.xaxis.set_major_locator(MaxNLocator(6))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax.yaxis.set_major_locator(LogLocator(base=10))
ax.yaxis.set_minor_locator(
    LogLocator(base=10, subs=np.arange(2,10)*0.1)
)
ax.yaxis.set_minor_formatter(NullFormatter())

ax.tick_params(axis='both', which='major', pad=8)

# -----------------------------------------------------
# Grid
# -----------------------------------------------------

ax.grid(True, which='major', alpha=0.30, linewidth=0.7)
ax.grid(True, which='minor', alpha=0.10, linewidth=0.5)

# -----------------------------------------------------
# Legend
# -----------------------------------------------------

leg = ax.legend(
    loc="best",
    frameon=False,
    fontsize=12,
    ncol=2
)

for line in leg.get_lines():
    line.set_linewidth(3)

fig.tight_layout()

fig.savefig(
    "species_model3.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()