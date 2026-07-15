import numpy as np
import matplotlib.pyplot as plt

FluxHeIM = 4.101e+15  
FluxHeIM_2 = 4.101e+13 
FluxHeIM_3 = 4.101e+11

taus = np.loadtxt('./chemeq2_shock/output/final_model_3.dat', comments='#')

taus = taus[:, -1]
  
#a0_H    = 3.87e-18 #seccion eficaz HI
#a0_HeIS = 4.03e-18
a0_HeIM = 3.59e-18

#phiHI   = a0_H*FluxHI*np.exp(-taus(1))
#phiHeIS = a0_HeIS*FluxHeIS*np.exp(-taus(2))
#phiHeIM = a0_HeIM*FluxHeIM*np.exp(-taus)
phiHeIM = np.exp(-taus)
#phiHeIM_2 = a0_HeIM*FluxHeIM_2*np.exp(-taus)
#phiHeIM_3 = a0_HeIM*FluxHeIM_3*np.exp(-taus)

plt.plot(taus, label = f'taus, flux = {FluxHeIM:.2e}')
plt.xlabel('Index')
plt.ylabel('taus HeIM')
plt.title('taus HeIM vs Index')
plt.legend()
plt.grid()
plt.savefig('tauHeIM_plot.png')
plt.close()

plt.plot(phiHeIM, label = f'phiHEIM, flux = {FluxHeIM:.2e}')
#plt.plot(phiHeIM_2, label = f'phiHEIM, flux = {FluxHeIM_2:.2e}')
#plt.plot(phiHeIM_3, label = f'phiHEIM, flux = {FluxHeIM_3:.2e}')
plt.xlabel('Index')
plt.ylabel('phiHeIM')
plt.title('phiHeIM vs Index')
#plt.yscale('log')
plt.legend()
plt.grid()
plt.savefig('phiHeIM_plot.png')