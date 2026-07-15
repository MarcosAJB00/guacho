import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator, MaxNLocator, NullFormatter

# =====================================================
# Constants
# =====================================================

Rjup = 7.1492e9
kB = 1.380649e-16
eV_to_erg = 1.602176634e-12

# =====================================================
# Load hydrodynamic model (T, density, velocity)
# =====================================================

hydro = np.loadtxt("model_3.dat", comments="#")

r = hydro[:,0] / Rjup
rho = hydro[:,1]
T = hydro[:,2]

# =====================================================
# Load chemistry model
# =====================================================

chem = np.loadtxt("final_model_3.dat", comments="#")

# columns:
# HI HII HeIS HeIM HeII e phiHI phiHeIS phiHeIM tauHI tauHeIS tauHeIM

nHI   = chem[:,0]
nHII  = chem[:,1]
nHeIS = chem[:,2]
nHeIM = chem[:,3]
nHeII = chem[:,4]
ne    = chem[:,5]

phiHeIM = chem[:,8]

# =====================================================
# Reaction rates (same as Fortran)
# =====================================================

def rates(T):

    logT = np.log10(T)

    gamma13 = np.zeros_like(T)
    gamma31a = np.zeros_like(T)
    gamma31b = np.zeros_like(T)

    mask = (logT >= 3.75) & (logT <= 5.75)

    gamma13[mask] = (
        -1.5699
        + 1.303*logT[mask]
        - 0.388*logT[mask]**2
        + 0.0516*logT[mask]**3
        - 0.0026*logT[mask]**4
    )

    gamma31a[mask] = (
        -209.42
        + 170.03*logT[mask]
        - 50.073*logT[mask]**2
        + 6.42*logT[mask]**3
        - 0.3045*logT[mask]**4
    )

    gamma31b[mask] = (
        -36.477
        + 22.988*logT[mask]
        - 4.5935*logT[mask]**2
        + 0.2968*logT[mask]**3
    )

    rate = {}

    # HI ionization
    rate["chi_HI"] = 5.83e-11*np.sqrt(T)*np.exp(-157828.0/T)

    # He I ionization
    rate["chi_HeIS"] = 2.379e-11*np.sqrt(T)*np.exp(-285335.4/T)

    # He metastable ionization
    rate["chi_HeIM"] = 8.335e-10*np.sqrt(T)*np.exp(-55338.0/T)

    # HeII ionization channels
    rate["chi_HeII_a"] = (
        2.10e-8*np.sqrt(13.6/(8.617e-5*T))
        * np.exp(-0.80/(8.617e-5*T))
        * gamma31a/3
    )

    rate["chi_HeII_b"] = (
        2.10e-8*np.sqrt(13.6/(8.617e-5*T))
        * np.exp(-1.40/(8.617e-5*T))
        * gamma31b/3
    )

    # recombination
    rate["rec_HeII"] = 2.10e-13*(1e4/T)**0.778

    # radiative decay
    rate["A31"] = 1.272e-4

    # charge exchange
    rate["ce1"] = 1.25e-15*(300./T)**(-0.25)
    rate["ce2"] = 1.75e-11*(300./T)**0.75*np.exp(-128000.0/T)

    return rate

# evaluate rates
rt = rates(T)

# =====================================================
# Production terms
# =====================================================

P1 = rt["chi_HeIS"] * ne * nHeIS
P2 = rt["rec_HeII"] * ne * nHeII

# =====================================================
# Destruction terms
# =====================================================

D1 = rt["A31"] * nHeIM
D2 = rt["ce1"] * nHI * nHeIM
D3 = rt["chi_HeIM"] * ne * nHeIM
D4 = rt["ce2"] * ne * nHeIM
D5 = phiHeIM * nHeIM

# =====================================================
# Total
# =====================================================

Ptot = P1 + P2
Dtot = D1 + D2 + D3 + D4 + D5

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
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
})

fig, ax = plt.subplots(figsize=(7.5,5.5))

# =====================================================
# Plot all rates
# =====================================================

#ax.plot(r, P1, lw=2,ls="--", label="Collisional ionization (He I)", color="tab:blue")
ax.plot(r, P2, lw=2,ls="--", label=r"Recombination to He I(2$^3$S)", color="tab:green")

#ax.plot(r, D1, lw=2, label="Radiative decay", color="black")
#ax.plot(r, D2, lw=2, label="Charge exchange (H I)", color="tab:red")
#ax.plot(r, D3, lw=2, label="Electron impact ionization", color="tab:orange")
#ax.plot(r, D4, lw=2, label="Electron processes", color="tab:purple")
ax.plot(r, D5, lw=2, label="Photoionization", color="tab:brown")

# totals
#ax.plot(r, Ptot, lw=3.5, ls="--", color="navy", label="Total production")
#ax.plot(r, Dtot, lw=3.5, ls="--", color="darkred", label="Total destruction")

# =====================================================
# Axes
# =====================================================

ax.set_yscale("log")

ax.set_xlabel(r"Distance [$R_J$]")
ax.set_ylabel(r"Rate [cm$^{-3}$ s$^{-1}$]")
#ax.set_ylim(1e-50,1e10)

ax.xaxis.set_major_locator(MaxNLocator(6))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax.yaxis.set_major_locator(LogLocator(base=10))
ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2,10)*0.1))
ax.yaxis.set_minor_formatter(NullFormatter())

ax.grid(True, which="major", alpha=0.3)
ax.grid(True, which="minor", alpha=0.1)

ax.legend(frameon=False, fontsize=11, ncol=2)

fig.tight_layout()
fig.savefig("HeIM_reaction_budget_model3.png", dpi=300, bbox_inches="tight")

plt.show()