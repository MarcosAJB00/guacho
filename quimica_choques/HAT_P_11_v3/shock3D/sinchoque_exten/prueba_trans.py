import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

# Parámetros (igual que tu código original)
nxmap      = 256
nymap      = 256
nvmap      = 250
velinit    = -150e5
velfinal   =  150e5
rorb       = 0.0532 * 1.49e13
Rjup       = 7.1492E9
Rsun       = 6.955e10
Rstar      = 0.76 * Rsun
Rp         = 0.4466 * Rjup
dx         = 2.0*5.09421794574017 * Rp / float(nxmap)
rstar      = Rstar / dx
rplan      = Rp / dx
miss_pix   = 144044
path       = './helines/'
kind       = 'd'

Emission_missing_pixel    = np.ones(nvmap)
x_vel = np.linspace(velinit, velfinal, nvmap)

lightcurve_triplet = []
phases = np.arange(0, 81)
theta_deg = -5.0 + (5.0 - (-5.0)) * phases / 80.0   # -5° a +5°

fout = open(path + 'dat/HE_tau_000_triplete.dat', "w")

for nout in range(0, 81):

    # --- leer los 3 cubos de tau ---
    tau_total = np.zeros((nvmap, nymap, nxmap))
    for typeout in ['3pj0', '3pj1', '3pj2']:
        filein = path + f'HE_tau_{typeout}_{nout:03d}.bin'
        f      = FortranFile(filein, 'r')
        data   = f.read_reals(kind).reshape(nxmap, nymap, nvmap, order='F').T
        tau_total += data   # suma de tau ANTES de exp()

    # --- transmision combinada del triplete ---
    emtaunu = np.exp(-tau_total)

    # --- integrar sobre el disco estelar (igual que tu código original) ---
    EmissionVel   = np.zeros(nvmap)
    TotalEmission = 0.
    Ref           = 0.

    theta_plan = (nout - 40) * np.pi * 0.125 / 180
    xplan      = rorb / dx * np.tan(theta_plan)
    yplan      = 0.

    for i in range(nxmap):
        for j in range(nymap):
            x  = (float(i - nxmap/2) + 0.5)
            y  = (float(j - nymap/2) + 0.5)
            r1 = np.sqrt(x**2 + y**2)

            if r1 <= rstar:
                r2 = np.sqrt((x - xplan)**2 + (y - yplan)**2)
                if r2 < rplan:
                    emtaunu[:, j, i] = 0.
                    #EmissionVel[:] = EmissionVel[:] + emtaunu[:,j,i]
                    #TotalEmission  = TotalEmission + np.sum(emtaunu[:,j,i])
                else:
                    EmissionVel   += emtaunu[:, j, i]
                    TotalEmission += np.sum(emtaunu[:, j, i])
                Ref += 1.

    # --- pixeles faltantes ---
    Ref            = Ref + miss_pix
    EmissionVel    = (EmissionVel + miss_pix * Emission_missing_pixel) / Ref
    TotalEmission  = (TotalEmission + miss_pix * 250) / Ref / float(nvmap)

    lightcurve_triplet.append(np.mean(EmissionVel))
    fout.write(f"{nout} {np.mean(EmissionVel)}\n")
    print(f"fase {nout:02d} | absorción triplete: {(1-np.mean(EmissionVel))*100:.3f}%")

fout.close()

# --- plot curva de luz ---
lightcurve_triplet = np.array(lightcurve_triplet)

geom_name = path+ 'dat/geom_trans_phase.dat'
geom_data = np.loadtxt(geom_name)
geom_fase = geom_data[:,0]
geom_abs = geom_data[:,1]

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(theta_deg, geom_abs, label='Geometric',color='grey')
ax.plot(theta_deg, lightcurve_triplet, 'o-', markersize=4, label='He I triplete')
#ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='continuo')
ax.set_xlabel('Fase orbital (grados)')
ax.set_ylabel('Flujo normalizado')
ax.set_title('Curva de luz - He I 10830 Å (triplete completo)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('lightcurve_He10830_triplete.png', dpi=150)
plt.show()