import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('salida_datos.txt')  

time = data[:, 0]
O2_minus = data[:, 1]
Cs_plus = data[:, 2]
Cs = data[:, 3]
CsO2 = data[:, 4]
O2 = data[:, 5]
N2 = data[:, 6]
e_minus = data[:, 7]

plt.figure(figsize=(10, 6))
plt.plot(time, O2_minus, label='O2-', color='blue')
plt.plot(time, Cs_plus, label='Cs+', color='orange')
plt.plot(time, Cs, label='Cs', color='green')
plt.plot(time, CsO2, label='CsO2', color='red')
plt.plot(time, O2, label='O2', color='purple')
plt.plot(time, N2, label='N2', color='brown')
plt.plot(time, e_minus, label='e-', color='black')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Time')
plt.ylabel('Density')
plt.title('Density vs Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('density_vs_time.png', dpi=300)
