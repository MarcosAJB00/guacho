import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#---------------------------------------
# Leer archivo
#---------------------------------------
# Ajustá el separador si no es coma
df = pd.read_csv('sorce_ssi_l3.csv')

# Asegurate de que las columnas tengan estos nombres
# (modificá si hace falta)
# fecha | lambda | irradiancia

df.columns = ['fecha', 'lambda', 'irr']

#---------------------------------------
# Agrupar por fecha
#---------------------------------------
fechas = df['fecha'].unique()

plt.figure(figsize=(8,6))

for f in fechas:  # Limitar a las primeras 10 fechas para no saturar el gráfico
    sub = df[df['fecha'] == 2455008.0]

    # ordenar por longitud de onda (importante)
    sub = sub.sort_values('lambda')

    plt.plot(sub['lambda'][:10], sub['irr'][:10], label=f'{f}')
    #np.savetxt(f'espectro_{f}.txt', sub[['lambda', 'irr']].values, header='lambda\tirr', fmt='%.6e', comments='')

#---------------------------------------
# Plot
#---------------------------------------
plt.yscale('log') 
plt.xscale('log') 
#plt.xlim(1, 1000)  # Ajustá según el rango de longitudes de onda
plt.xlabel('Longitud de onda')
plt.ylabel('Irradiancia')
plt.title('Comparación de espectros')

# si hay muchas fechas esto puede saturar
# podés comentar esta línea si molesta
#plt.legend(fontsize=6)

plt.tight_layout()
plt.savefig('espectros_multiepoca.png', dpi=300)
plt.show()