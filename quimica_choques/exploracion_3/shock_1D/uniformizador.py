import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

models = np.loadtxt('./output/model_list.dat', skiprows=1)
model_numbers = models[:,0]

fig_rho, ax_rho = plt.subplots()
fig_T, ax_T = plt.subplots()

for i in range(len(model_numbers)):

    print(f'Processing model {int(model_numbers[i])} of {len(model_numbers)}')
    data = np.loadtxt(
        f'./output/output_{int(model_numbers[i])}.dat',
        skiprows=1
    )
    
    x = data[1:,0] # en cm
    rho = data[1:,1] # en g/cm^3
    T = data[1:,2]  # en K

    # Grid uniforme nuevo
    N_new = 1000
    x_uniform = np.linspace(x.min(), x.max(), N_new)

    # Interpoladores
    dens_interp = PchipInterpolator(x, np.log10(rho))
    temp_interp = PchipInterpolator(x, T)

    # Evaluar
    density_uniform = 10**dens_interp(x_uniform)
    temperature_uniform = temp_interp(x_uniform)

    #print(f'Model {int(model_numbers[i])} x.min(): {x.min()}')
    #print(f'Model {int(model_numbers[i])} x_uniform.min(): {x_uniform.min()}')


    np.savetxt(f'./uniform_output/model_{int(model_numbers[i])}.dat',
        np.column_stack((x_uniform, density_uniform, temperature_uniform)),
        header='r[R_p] density[g/cm^3] temperature[K]',
        fmt='%.10e'
    )


    ax_rho.plot(x, rho, label='Density Profile')
    ax_rho.plot(x_uniform, density_uniform, linestyle='--')

    ax_T.plot(x, T, label='Temperature Profile')
    ax_T.plot(x_uniform, temperature_uniform, linestyle='--')
   
ax_rho.set_xlabel('Radius (Rp)')
ax_rho.set_ylabel('Density (g/cm^3)')
ax_rho.set_title('Initial Profiles density')
#ax_rho.legend()
ax_rho.set_xscale('log')
ax_rho.set_yscale('log')
ax_rho.grid(True, which='both', linestyle='--', linewidth=0.5)
fig_rho.savefig('uniform_density_profile.png', dpi=300, bbox_inches='tight')


ax_T.set_xlabel('Radius (Rp)')
ax_T.set_ylabel('Temperature (K)')
ax_T.set_title('Initial Profiles temperature')
#ax_T.legend()
ax_T.set_xscale('log')
ax_T.set_yscale('log')
ax_T.grid(True, which='both', linestyle='--', linewidth=0.5)
fig_T.savefig('uniform_temperature_profile.png', dpi=300, bbox_inches='tight')

plt.show()

