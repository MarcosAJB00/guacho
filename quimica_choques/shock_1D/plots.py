import numpy as np
import matplotlib.pyplot as plt

models = np.loadtxt('./output/model_list.dat', skiprows=1)
model_numbers = models[:,0]

fig_rho, ax_rho = plt.subplots()
fig_T, ax_T = plt.subplots()

for i in range(len(model_numbers)):

    data = np.loadtxt(
        f'./output/output_{int(model_numbers[i])}.dat',
        skiprows=1
    )

    x = data[:,0]
    rho = data[:,1]
    T = data[:,2]
    
    mp = 1.67e-24 # g
    ax_rho.plot(x, rho/mp, label=f'Model {int(model_numbers[i])}')
    ax_T.plot(x, T, label=f'Model {int(model_numbers[i])}')

ax_rho.set_xlabel('x [cm]')
ax_rho.set_ylabel('rho [1/cm^3]')
ax_rho.set_yscale('log')
ax_rho.set_xscale('log')
ax_rho.grid(True, which='both', ls='--')
ax_rho.set_title('Density profiles')
#ax_rho.legend()

ax_T.set_xlabel('x [cm]')
ax_T.set_ylabel('T [K]')
ax_T.set_yscale('log')
ax_T.set_xscale('log')
ax_T.grid(True, which='both', ls='--')
ax_T.set_title('Temperature profiles')
#ax_T.legend()

fig_rho.savefig('density_profiles.png', dpi=300, bbox_inches='tight')
fig_T.savefig('temperature_profiles.png', dpi=300, bbox_inches='tight')


plt.show()