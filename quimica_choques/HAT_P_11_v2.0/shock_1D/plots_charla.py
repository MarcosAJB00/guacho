import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (
    AutoMinorLocator,
    LogLocator,
    MaxNLocator,
    NullFormatter
)

# ----------------------------
# Constants
# ----------------------------
Rjup = 7.1492e9      # cm
mp = 1.67e-24        # g

# ----------------------------
# Load model 3
# ----------------------------
model = 3

data = np.loadtxt(
    f'./output/output_{model}.dat',
    skiprows=1
)

x = data[1:,0] / Rjup
rho = data[1:,1] / mp
T = data[1:,2]

# ----------------------------
# Plot style
# ----------------------------
plt.rcParams.update({

    "font.size": 15,

    "axes.linewidth": 1.4,

    "xtick.direction": "in",
    "ytick.direction": "in",

    "xtick.top": True,
    "ytick.right": True,

    # Major ticks
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.major.width": 1.3,
    "ytick.major.width": 1.3,

    # Minor ticks
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.width": 1.0,
    "ytick.minor.width": 1.0,
})

# ======================================================
# Density
# ======================================================

fig, ax = plt.subplots(figsize=(6.5,5.0))

ax.plot(
    x,
    rho,
    color='black',
    lw=2.2
)

ax.set_yscale('log')

ax.set_xlabel(r'Distance [$R_{\rm J}$]',fontsize=25)
ax.set_ylabel(r'n [cm$^{-3}$]',fontsize=25)

# Major ticks
ax.xaxis.set_major_locator(MaxNLocator(6))
ax.yaxis.set_major_locator(LogLocator(base=10))

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(
    LogLocator(base=10, subs=np.arange(2,10)*0.1)
)
ax.yaxis.set_minor_formatter(NullFormatter())

# Tick appearance
ax.tick_params(axis='both', which='major', pad=8)
ax.tick_params(axis='both', which='minor')

# Grid
ax.grid(True, which='major', alpha=0.30, linewidth=0.7)
ax.grid(True, which='minor', alpha=0.10, linewidth=0.5)

fig.tight_layout()

fig.savefig(
    "density_profile_model3.png",
    dpi=300,
    bbox_inches="tight"
)

# ======================================================
# Temperature
# ======================================================

fig, ax = plt.subplots(figsize=(6.5,5.0))

ax.plot(
    x,
    T,
    color='black',
    lw=2.2
)

ax.set_xlabel(r'Distance [$R_{\rm J}$]',fontsize=25)
ax.set_ylabel(r'T [K]',fontsize=25)

# Major ticks
ax.xaxis.set_major_locator(MaxNLocator(6))
ax.yaxis.set_major_locator(MaxNLocator(6))

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

# Tick appearance
ax.tick_params(axis='both', which='major', pad=8)
ax.tick_params(axis='both', which='minor')

# Grid
ax.grid(True, which='major', alpha=0.30, linewidth=0.7)
ax.grid(True, which='minor', alpha=0.10, linewidth=0.5)

fig.tight_layout()

fig.savefig(
    "temperature_profile_model3.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()