import numpy as np
import matplotlib.pyplot as plt

# ---------- leer datos ----------
#r_m, dens, T, vel = np.loadtxt("model_3.dat", unpack=True)
r_m, T = np.loadtxt("temperature_profile_ATES.dat", unpack=True)
#data = np.loadtxt("final_model_3.dat")
data = np.loadtxt("perfiles_finales.dat")
#r_f, nHI, nHII, nHeIS, nHeIM, nHeII, ne, phiHI, phiHeIS, phiHeIM, tauHI, tauHeIS, tauHeIM = data.T
r_f, nHI, nHII, nHeIS, nHeIM, nHeII, ne = data.T

data_phi = np.loadtxt("phi_tau_finales.dat")
r_2, phiHI, phiHeIS, phiHeIM, tauHI, tauHeIS, tauHeIM = data_phi.T

# mismo grid en r -> uso directamente (si no, habria que interpolar T(r) a r_f)
assert len(r_m) == len(r_f), "r de model_3.dat y final_model_3.dat no coinciden en longitud"
r = r_f

# ---------- gamma_ij (Bray 2000), validos solo 3.75<logT<5.75 ----------
logT = np.log10(T)
gamma13  = np.zeros_like(T)
gamma31a = np.zeros_like(T)
gamma31b = np.zeros_like(T)
mask = (logT >= 3.75) & (logT <= 5.75)
lt = logT[mask]
gamma13[mask]  = -1.5699 + 1.303*lt - 0.388*lt**2 + 0.0516*lt**3 - 0.0026*lt**4
gamma31a[mask] = -209.42 + 170.03*lt - 50.073*lt**2 + 6.42*lt**3 - 0.3045*lt**4
gamma31b[mask] = -36.477 + 22.988*lt - 4.5935*lt**2 + 0.2968*lt**3

kB_eV = 8.617333262145e-5

# rate coefficients que dependen de T (cm^3/s)
che13a = 2.10e-8*np.sqrt(13.6/(kB_eV*T))*np.exp(-19.81/(kB_eV*T))*gamma13 
che31a = 2.10e-8*np.sqrt(13.6/(kB_eV*T))*np.exp(-0.80/(kB_eV*T))*gamma31a/3
che31b = 2.10e-8*np.sqrt(13.6/(kB_eV*T))*np.exp(-1.40/(kB_eV*T))*gamma31b/3
che31  = 5.0e-10                                   # HeI(triplet) - HI collisions
aheiim_b = 2.10e-13*(1.0e4/T)**0.778                # recombinacion a 2^3S (case B)

# ---------- terminos de la ecuacion del HeIM (s^-1, "por particula") ----------
rate_photo   = phiHeIM                       # fotoionizacion del metaestable
rate_ecoll   = ne*(che31a + che31b + che13a)          # colisiones con electrones (perdida 3->1)
rate_hcoll   = nHI*che31                     # colisiones con atomos de H (perdida 3->1)
rate_recomb  = ne*aheiim_b                   # recombinacion HeII -> HeIM

# altitud normalizada (no se conoce Rp explicito: se usa r/r[0] como proxy)
alt = r/r[0]

# ---------- plot ----------
fig, ax = plt.subplots(figsize=(6,4.5))
ax.plot(alt, rate_photo,  color="orange",     lw=2, label="photoionization")
ax.plot(alt, rate_ecoll,  color="firebrick",  lw=2, label="electron collisions")
ax.plot(alt, rate_hcoll,  color="teal",       lw=2, label="H-atom collisions")
ax.plot(alt, rate_recomb, color="steelblue",  lw=2, ls="--", label="recombination")

ax.set_yscale("log")
ax.set_xlabel("altitude [r/r$_0$]")  # cambiar a r/Rp si tenes Rp en cm
ax.set_ylabel("rate [s$^{-1}$]")
#ax.set_ylim(1e-12, 1e-1)
ax.legend(frameon=False, fontsize=9)
ax.set_title("HeIM (triplet) rate terms")
plt.tight_layout()
plt.savefig("heim_rates_HP11.png", dpi=150, bbox_inches='tight')
plt.show()
print("listo")