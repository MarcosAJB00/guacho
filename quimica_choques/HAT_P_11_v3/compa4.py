import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Load light curves
# ------------------------------------------------------------
path1 = '/home/mbaracchi/guacho/quimica_choques/HAT_P_11_v3/shock3D/sinchoque_exten/'
path2 = '/home/mbaracchi/guacho/quimica_choques/HAT_P_11_v3/shock3D/shock_fiducial_exten/'

curve1 = np.loadtxt(path1 + 'lightcurve_triplet_sinchoque2.dat')
curve2 = np.loadtxt(path2 + 'lightcurve_triplet_model3.dat')

phase = curve1[:,0]

geom       = curve1[:,1]
flux        = curve1[:,5]

geom_photo = curve2[:,1]
flux_photo = curve2[:,5]

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 15,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

fig, ax = plt.subplots(figsize=(7,5))

ax.plot(
    phase,
    geom,
    lw=2.5,
    color="black",
    label="Geometric"
)

ax.plot(
    phase,
    flux,
    lw=1.5,
    color="tab:red",
    label=r"Sin Choque"
)


ax.plot(
    phase,
    flux_photo,
    lw=1.5,
    ls="--",
    color="tab:red",
    label=r"Con Choque Ale"
)

ax.set_xlabel("Orbital phase")
ax.set_ylabel("Normalized stellar flux")

ax.set_xlim(phase.min(), phase.max())
ax.grid(alpha=0.25, linestyle="--")
ax.legend(frameon=False)

plt.tight_layout()

plt.savefig("lightcurve_comparison.jpg", dpi=300, bbox_inches="tight")

plt.show()