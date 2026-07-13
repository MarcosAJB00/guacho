#!/usr/bin/python
#import fortranfile
from scipy.io import FortranFile
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable

endian='<'      # little endian
kind='d'        # double precision

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    return fig.colorbar(mappable, cax=cax,label='', extend='max')

velinit  = -150e5
velfinal =  150e5
nvmap    =  250
dvel=(velfinal-velinit)/float(nvmap)

#nxmap = 1032
#nymap = 258

nxmap = 256
nymap = 256

rorb  = 0.0532*1.49e13
torb  = 4.887802443*86400.
Rjup  = 7.1492E9
Rsun  = 6.955e10
Rstar = 0.76 * Rsun
Rp    = 0.4466  * Rjup
dx    = 2.0*4.77*Rp/float(nxmap)# * Rp/float(nxmap) #  dx in cm ( 1AU = 1.49e13 cm)
rstar = Rstar/dx     #  star radius in in pixel x
rplan = Rp/dx

print('DF_geom=', (Rp/Rstar)**2)
runs     = [0 ] #Numero de salida de GUACHO
typeouts = [ '3pj0', '3pj1', '3pj2' ] #Lineas espectrales

#  number of missing pixels of stellar surface |z_max|<Rst
miss_pix = 144044#1316

for run in runs:
    for typeout in typeouts:
        print (run, typeout)

        path ='./helines/'
        fout1 = open(path+'dat/HE_tau_'+str(run).zfill(3)+'_'+typeout+'.dat',"w")
        Emission_missing_pixel = np.zeros(nvmap)
        Emission_missing_pixel[:] = 1.

        #  nout = phase (los)
        for nout in range(0,81):

            #filein=path+'HE_tau_'+typeout+'_'+str(run).zfill(3)+'_'+str(nout).zfill(3)+'.bin'
            filein=path+'HE_tau_'+typeout+'_'+str(nout).zfill(3)+'.bin'
            f=FortranFile(filein,'r')
            data=f.read_reals(kind).reshape(nxmap,nymap,nvmap,order='F').T
            emtaunu=np.exp(-data)

            EmissionVel=np.zeros(shape=(nvmap))
            TotalEmission=0.
            Ref = 0.
    
            theta_plan= (nout-40)*np.pi*0.125/180
            yplan= 0.#nymap/2
            xplan = rorb/dx*np.tan(theta_plan) 
            #print(xplan,yplan)
            #print(rstar,rplan)

            for i in range(nxmap):
                for j in range(nymap):
                    x = (float(i-nxmap/2)+0.5)
                    y = (float(j-nymap/2)+0.5)
                    r1=np.sqrt( x**2 + y**2 ) #del centro del mapa

                    if (r1 <= rstar ):
                        r2 = np.sqrt((x - xplan)**2 + (y-yplan)**2)

                        if (r2 < rplan):
                                emtaunu[:,j,i]=0.
                                #EmissionVel[:] = EmissionVel[:] + emtaunu[:,j,i]
                                #TotalEmission  = TotalEmission + np.sum(emtaunu[:,j,i])
                        else:
                                EmissionVel=EmissionVel[:] + emtaunu[:,j,i]
                                TotalEmission=TotalEmission + np.sum(emtaunu[:,j,i])
                        Ref = Ref + 1.
            #adding missing pixels from not accounting for the upper and lower edges of the star
            Ref = Ref + miss_pix
            EmissionVel=(EmissionVel[:] + miss_pix*Emission_missing_pixel[:])/Ref
            TotalEmission = (TotalEmission + miss_pix*250)/Ref/float(nvmap)
            x= np.linspace(velinit,velfinal,nvmap)

            fout2= open(path+'dat/HE_tau-'+str(run).zfill(3)+'-abs-'+str(nout).zfill(3)+'_'+typeout+'.dat',"w")
            for l in range(nvmap):
                fout2.write(str(x[l])+' '+str(EmissionVel[l])+ "\n" )
            fout2.close()
            print ("total absorption for output ",nout,", is :", (1.-TotalEmission)*100)

            fout1.write(str(nout)+" "+str(np.mean(EmissionVel))+ "\n")

            fig = plt.figure(2)
            plt.clf()
            ax = fig.add_subplot(1, 1, 1)
            #circ2=plt.Circle((xp-100,yp),radius=5*rplan,color='k',fill=False)
            cax=ax.imshow(np.sum(emtaunu[:,:,:],0)/250,origin='lower',norm=LogNorm(vmin=9.0e-1,vmax=1.0))
            #ax.add_patch(circ2)
            cbar = fig.colorbar(cax)
            plt.savefig(path+'dat/'+typeout+'_tau_map-'+str(run).zfill(3)+'-abs-'+str(nout).zfill(3)+'.png')


        fout1.close()
