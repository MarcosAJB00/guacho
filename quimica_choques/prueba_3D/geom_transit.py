#!/usr/bin/python
#import fortranfile
#from scipy.io import FortranFile
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable

velinit=-100e5
velfinal=100e5
nvmap=250
dvel=(velfinal-velinit)/float(nvmap)
nxmap=50 #1032
nymap=1500 # ficticius size to include the star diameter 
nymap_real= 50
rorb = 0.03397*1.49e13
torb = 2.15*86400.

Rjup=7.1492E9
Rsun=6.955e10
Rstar=0.76*Rsun
Rp = 0.4466*Rjup

dx=5.09421794574017*Rp/float(nxmap) #  dx in cm ( 1AU = 1.49e13 cm)
rstar= Rstar/dx     #  star radius in in pixel x
rplan= Rp/dx
print(f'rstar: {rstar}, rplan: {rplan}')
print('DF_geom=', (Rp/Rstar)**2)

path ='./helines/'

'''
this corresponds to the averaged (over the stellar disk)  absortion 
as function of tinme (phase)
'''
fout1 = open(path+'dat/geom_trans_phase.dat',"w")

#  nout corresponds to the phase
for nout in range(0,81):

    miss_pix = 0
    plan_pix = 0
    star_pix = 0

    EmissionVel=np.zeros(shape=(nvmap))
    TotalEmission=0.
    Ref = 0.
    
    theta_plan= (nout-40)*np.pi/180
    yplan= 0.#nymap/2
    xplan = rorb/dx*np.tan(theta_plan) 

    for i in range(nxmap):
        for j in range(nymap):
            x = (float(i-nxmap/2)+0.5)
            y = (float(j-nymap/2)+0.5)
            r1=np.sqrt( x**2 + y**2 )
            if (r1 <= rstar ):
                r2 = np.sqrt((x - xplan)**2 + (y-yplan)**2)
                star_pix += 1
                if (r2 < rplan):
                        plan_pix += 1
                        # inside the planet we dont add emission
                else:
                        if (abs(y) > float(nymap_real/2) ):
                          miss_pix += 1
                          #print(j)
                        EmissionVel += 1.0
                        TotalEmission += 1.0
                Ref = Ref + 1.

    EmissionVel=EmissionVel[:] / Ref
    TotalEmission = TotalEmission /Ref
    x= np.linspace(velinit,velfinal,nvmap)

    fout2= open(path+f'dat/geom_trans_vel-{nout:03d}.dat',"w")
    for l in range(nvmap):
        fout2.write(str(x[l])+' '+str(EmissionVel[l])+ "\n" )
    fout2.close()
    print (f'total absorption for phase {nout:0d}, is :, {(1.-TotalEmission):.3%}')
    fout1.write(str(nout)+" "+str(np.mean(EmissionVel))+ "\n")

    print(f'number of star/planet pixels: {star_pix:0d}/{plan_pix:0d}')
    print(f'number of outside pixels: {miss_pix:0d}\n')
fout1.close()

