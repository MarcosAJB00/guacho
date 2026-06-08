import numpy as np 
import matplotlib.pyplot as plt


# Stand-off distance
data_stand_off = np.loadtxt('../ram_pressure_ATES/condi_ini_shock_ATES.dat', unpack=True)

r_stand_off = data_stand_off[:,0] #en Rp
rho_stand_off = data_stand_off[:,1] #en g/cm^3
temp_stand_off = data_stand_off[:,2] #en K
v_stand_off = data_stand_off[:,3] #en cm/s

# Condiciones iniciales modelos shocks

shock_models = np.loadtxt('./shock_1D/output/model_list.dat', unpack=True)

model_number = shock_models[:,0]
shock_temp = shock_models[:,1] #en K
n0_shock = shock_models[:,2] #en cm^-3
v_shock = shock_models[:,3] #en cm/s
y0_shock = shock_models[:,4] 

# Quimca de los choques
path = './chemeq2_1D/output/'
n_steps = 1000
X_HeIM = np.zeros(len(model_number), n_steps)
r_HeIM = np.zeros(len(model_number), n_steps)
for i in range(len(model_number)):
    data_chem = np.loadtxt(path + 'final_model_' + str(int(model_number[i])) + '.dat', unpack=True)
    r_HeIM[i,:] = data_chem[:,0] #en Rp
    X_HeIM[i,:] = data_chem[:,4] #en cm^-3

# Quimica sin choques

sin_choque = np.loadtxt('./sin_choque/output/perfiles_finales.dat', unpack=True)

r_wos = sin_choque[:,0] #en Rp
HeIM_wos = sin_choque[:,4] #en cm^-3
HeIM_trap =[]
HeIM_trap.append(np.trapz(HeIM_wos, r_wos))
# Combinar perfiles de HeIM

plt.figure(figsize=(8, 6))
plt.plot(r_wos, HeIM_wos, label='HeIM Profile without Shock')

for i in range(len(model_number)):
    mask = (r >= 1.05) & (r <= r_stand_off[i])
    r_comb = np.concatenate((r[mask], r_stand_off[i] + r_HeIM[i,:]))
    HeIM_combi = np.concatenate((HeIM_wos[mask], X_HeIM[i,:]))
    HeIM_trap.append(np.trapz(HeIM_combi, r_comb))
    
    plt.plot(r_comb, HeIM_combi, label='HeIM Profile with Shock Model ' + str(int(model_number[i])))

plt.xlabel('Radius [Rp]')
plt.ylabel('HeIM Density [cm^-3]')
plt.title('Comparison of HeIM Profiles with and without Shock')
plt.grid(alpha=0.3, linestyle='--', linewidth=0.5, which='both')
plt.legend()
plt.savefig('comparacion_HeIM_profiles.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(8, 6))
plt.scatter(HeIM_trap, color='red')
plt.xlabel('Model Number')
plt.ylabel('Integrated HeIM Density')
plt.title('Comparison of HeIM Profiles with and without Shock')
plt.grid(alpha=0.3, linestyle='--', linewidth=0.5, which='both')
plt.savefig('comparacion_HeIM_trap.png', dpi=300, bbox_inches='tight')
plt.show()



