import numpy as np
import matplotlib.pyplot as plt

path = './helines/'
typeouts = ['3pj0', '3pj1', '3pj2']

geom_name = path+ 'dat/geom_trans_phase.dat'
geom_data = np.loadtxt(geom_name)
geom_fase = geom_data[:,0]
geom_abs = geom_data[:,1]

# Ángulo real en grados según tu código (paso de 0.125°)
phases = np.arange(0, 81)
theta_degrees = (phases - 40) * 0.125

fig, ax = plt.subplots(figsize=(10, 5))

abs_total = np.zeros(len(phases))

P = 4.88 # periodo en días
time_hours = theta_degrees * P * 24.0 / 360.0

for typeout in typeouts:
    fname = path + f'dat/HE_tau_000_{typeout}.dat'
    data  = np.loadtxt(fname)
    
    nout       = data[:, 0]
    mean_flux  = data[:, 1]
    absorption = (1. - mean_flux) * 100   # en porcentaje
    abs_total += (1.0 - mean_flux)

    #ax.plot(theta_degrees, mean_flux, label=typeout, marker='o', markersize=3)
    ax.plot(time_hours, mean_flux, label=typeout, marker='o', markersize=3)

flux_total = 1.0 - abs_total

np.savetxt('CL_He_Hp11_choque3.dat',
    np.column_stack((theta_degrees, geom_abs, flux_total)),
    header=' theta (degrees)  Flux norm',
    fmt='%.10e'
)

#ax.plot(theta_degrees, flux_total, label='Total Abs',color='k',marker='o', markersize=3)
#ax.plot(theta_degrees, geom_abs, label='Geometric',color='grey')
ax.plot(time_hours, geom_abs, label='Geometric',color='grey')
ax.plot(time_hours, flux_total, label='Total Abs',color='k',marker='o', markersize=3)
#ax.set_xlabel('Fase orbital (grados)')
ax.set_xlabel('Tiempo de transito (horas)')
ax.set_ylabel('Absorción (%)')
ax.set_title('Curva de luz - He I 10830 Å')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('lightcurve_He10830.png', dpi=150)
plt.show()