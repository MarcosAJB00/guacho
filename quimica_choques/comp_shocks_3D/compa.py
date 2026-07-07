import numpy as np
import matplotlib.pyplot as plt

#sin choques
ATES = np.loadtxt('lightcurve_triplet_sinchoque2.dat',comments='#')

xp_ss = ATES[:,0]
geom = ATES[:,1]
flux_ss = ATES[:,5]


shock3 = np.loadtxt('lightcurve_triplet_model3_2.dat',comments='#')

xp_s3 = shock3[:,0]
flux_s3 = shock3[:,5]

#shock5 = np.loadtxt('lightcurve_triplet_model5.dat',comments='#')

#xp_s5 = shock5[:,0]
#flux_s5 = shock5[:,5]

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(xp_ss, geom, label='Geometrico',color='grey')
ax.plot(xp_ss, flux_ss, '-', label='Sin Choque',lw=1)
ax.plot(xp_s3, flux_s3, '-', label='Choque 3',lw=1)
#ax.plot(xp_s5, flux_s5, '-', label='Choque 5',lw=3)

ax.set_xlabel('xp (Rp)')
ax.set_ylabel('Absorcion (%)')
ax.set_title('Curva de luz - He I 10830 A')
ax.legend()
#ax.set_yscale('log')
ax.set_ylim(0.91,.92)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('comparacion_CLs_He10830.png', dpi=150)
plt.show()