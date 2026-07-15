import numpy as np
import matplotlib.pyplot as plt

models = np.loadtxt('./output/model_list.dat', skiprows=1)
if len(models.shape) == 1:
    model_numbers = np.array([models[0]])
    T0_values = np.array([models[1]])
    n0_values = np.array([models[2]])
    u0_values = np.array([models[3]])
    y0_values = np.array([models[4]])
else:
    model_numbers = models[:,0]
    T0_values = models[:,1]
    n0_values = models[:,2]
    u0_values = models[:,3]
    y0_values = models[:,4]

fig_rho, ax_rho = plt.subplots()
fig_T, ax_T = plt.subplots()

Rjup = 7.1492e9 # cm

for i in range(len(model_numbers)):

    data = np.loadtxt(
        f'./output/output_{int(model_numbers[i])}.dat',
        skiprows=1
    )
    
    x = data[:,0]/Rjup
    rho = data[:,1]
    T = data[:,2]
    
    mp = 1.67e-24 # g

    #if model_numbers[i] != 3 and model_numbers[i] != 15:
    #    ax_rho.plot(x, rho/mp, linestyle='--', alpha=0.6)
    #    ax_T.plot(x, T, linestyle='--', alpha=0.6)
    #elif model_numbers[i] == 3 or model_numbers[i] == 15:
    ax_rho.plot(x, rho/mp, lw =1.5, label=f'Model {int(model_numbers[i])}')
    ax_T.plot(x, T, lw =1.5, label=f'Model {int(model_numbers[i])}')

ax_rho.set_xlabel('x [Rjup]')
ax_rho.set_ylabel('rho [1/cm^3]')
ax_rho.set_yscale('log')
ax_rho.set_xscale('log')
#ax_rho.set_xlim(0.0, 20.0)
ax_rho.grid(True, which='both', ls='--')
ax_rho.set_title('Density profiles')
ax_rho.legend()

ax_T.set_xlabel('x [Rjup]')
ax_T.set_ylabel('T [K]')
ax_T.set_yscale('log')
ax_T.set_xscale('log')
ax_T.grid(True, which='both', ls='--')
ax_T.set_title('Temperature profiles')
ax_T.legend()

fig_rho.savefig('density_profiles.png', dpi=300, bbox_inches='tight')
fig_T.savefig('temperature_profiles.png', dpi=300, bbox_inches='tight')


plt.show()