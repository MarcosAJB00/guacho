import numpy as np
import matplotlib.pyplot as plt

comp_quimica = True

models = np.loadtxt('./shock_1D/output/model_list.dat', skiprows=1)
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

models_photo = np.loadtxt('./shock_1D_photo/output/model_list.dat', skiprows=1)
if len(models_photo.shape) == 1:
    model_numbers = np.array([models_photo[0]])
    T0_values = np.array([models_photo[1]])
    n0_values = np.array([models_photo[2]])
    u0_values = np.array([models_photo[3]])
    y0_values = np.array([models_photo[4]])
else:
    model_numbers = models_photo[:,0]
    T0_values = models_photo[:,1]
    n0_values = models_photo[:,2]
    u0_values = models_photo[:,3]
    y0_values = models_photo[:,4]


fig_rho, ax_rho = plt.subplots()
fig_T, ax_T = plt.subplots()
fig_HeIM, ax_HeIM = plt.subplots()

Rjup = 7.1492e9 # cm

for i in range(len(model_numbers)):

    data = np.loadtxt(
        f'./shock_1D/uniform_output/model_{int(model_numbers[i])}.dat',
        skiprows=1
    )

    if comp_quimica:
        data_chem = np.loadtxt(
            f'./chemeq2_shock/output/final_model_{int(model_numbers[i])}.dat',
            skiprows=1
        )   

    data_photo = np.loadtxt(
        f'./shock_1D_photo/uniform_output/model_{int(model_numbers[i])}.dat',
        skiprows=1
    )

    if comp_quimica:
        data_chem_photo = np.loadtxt(
        f'./chemeq2_shock_photo/output/final_model_{int(model_numbers[i])}.dat',
        skiprows=1
        )   
    
    x = data[:,0]/Rjup
    rho = data[:,1]
    T = data[:,2]
    HeIM = data_chem[:,4]

    x_photo = data_photo[:,0]/Rjup
    rho_photo = data_photo[:,1]
    T_photo = data_photo[:,2]
    HeIM_photo = data_chem_photo[:,4]

    mp = 1.67e-24 # g
    if model_numbers[i] == 1:
        ax_rho.plot(x, rho/mp, lw =1.5, label=f'Model')
        ax_rho.plot(x_photo, rho_photo/mp, ls = '--', lw =1.5, label=f'Model photo ')
        ax_T.plot(x, T, lw =1.5, label=f'Model ')
        ax_T.plot(x_photo, T_photo, ls = '--', lw =1.5, label=f'Model photo ')
        ax_HeIM.plot(x, HeIM, lw =1.5, label=f'Model ')
        ax_HeIM.plot(x_photo, HeIM_photo, ls = '--', lw =1.5, label=f'Model photo ')
    else:
        ax_rho.plot(x, rho/mp, lw =1.5)
        ax_rho.plot(x_photo, rho_photo/mp, ls = '--', lw =1.5)
        ax_T.plot(x, T, lw =1.5)
        ax_T.plot(x_photo, T_photo, ls = '--', lw =1.5)
        ax_HeIM.plot(x, HeIM, lw =1.5)
        ax_HeIM.plot(x_photo, HeIM_photo, ls = '--', lw =1.5)

ax_rho.set_xlabel('x [Rjup]')
ax_rho.set_ylabel('rho [1/cm^3]')
ax_rho.set_yscale('log')
#ax_rho.set_xscale('log')
#ax_rho.set_xlim(0.0, 20.0)
ax_rho.grid(True, which='both', ls='--')
ax_rho.set_title('Density profiles')
ax_rho.legend()

ax_T.set_xlabel('x [Rjup]')
ax_T.set_ylabel('T [K]')
#ax_T.set_yscale('log')
#ax_T.set_xscale('log')
ax_T.grid(True, which='both', ls='--')
ax_T.set_title('Temperature profiles')
ax_T.legend()

ax_HeIM.set_xlabel('x [Rjup]')
ax_HeIM.set_ylabel('HeIM')
ax_HeIM.set_yscale('log')
#ax_HeIM.set_xscale('log')
ax_HeIM.grid(True, which='both', ls='--')
ax_HeIM.set_title('HeIM profiles')
ax_HeIM.legend()

fig_rho.savefig('density_profiles.png', dpi=300, bbox_inches='tight')
fig_T.savefig('temperature_profiles.png', dpi=300, bbox_inches='tight')
fig_HeIM.savefig('HeIM_profiles.png', dpi=300, bbox_inches='tight')

plt.show()