#!/usr/bin/env python3
from scipy.io import FortranFile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ===========================================================================
# PARÁMETROS
# ===========================================================================
path    = './helines/'
nout    = 30        # fase a visualizar (0-80)
typeout = '3pj0'    # componente
vel_kms = -10.0       # canal de velocidad a graficar [km/s]

nxmap   = 256
nymap   = 256
nvmap   = 250
velinit  = -150e5
velfinal =  150e5
kind    = 'd'
# ===========================================================================

# índice del canal de velocidad
iv = int((vel_kms*1e5 - velinit) / (velfinal - velinit) * nvmap)
iv = np.clip(iv, 0, nvmap-1)
print(iv)

filein = path + f'HE_tau_{typeout}_{nout:03d}.bin'
f    = FortranFile(filein, 'r')
data = f.read_reals(kind).reshape(nxmap, nymap, nvmap, order='F').T
f.close()
# data shape: (nvmap, nymap, nxmap)

plt.figure(figsize=(6, 5))
plt.imshow(data[iv,:, :], origin='lower', cmap='inferno',
           extent=[-nxmap/2, nxmap/2, -nymap/2, nymap/2], norm=LogNorm())
plt.colorbar(label='τ')
plt.title(f'Mapa de τ — {typeout} — nout={nout} — v={vel_kms:+.0f} km/s')
plt.xlabel('x [pixel]')
plt.ylabel('y [pixel]')

plt.tight_layout()
plt.savefig(f'tau_map_{typeout}_{nout:03d}_v{vel_kms:+.0f}.png', dpi=150)
plt.show()


print(f'shape : {data.shape}')
print(f'min   : {data.min():.6f}')
print(f'max   : {data.max():.6f}')
print(f'mean  : {data.mean():.6f}')
print(f'celdas con tau > 0 : {(data > 0).sum()}')
print(f'celdas con tau > 1 : {(data > 1).sum()}')
print(f'NaN : {np.isnan(data).sum()}')
print(f'Inf : {np.isinf(data).sum()}')
print(f'negativos : {(data < 0).sum()}')