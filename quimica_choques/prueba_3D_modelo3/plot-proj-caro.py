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

velinit  = -100e5
velfinal = 100e5
nvmap    = 250
dvel=(velfinal-velinit)/float(nvmap)

nxmap = 50
nymap = 206

rorb  = 0.0532*1.49e13
torb  = 4.887802443*86400.
Rjup  = 7.1492E9
Rsun  = 6.955e10
Rstar = 0.76 * Rsun
Rp    = 0.4466  * Rjup
dx    = 5.09421794574017  * Rp/float(nxmap)
rstar = Rstar/dx
rplan = Rp/dx

print('DF_geom=', (Rp/Rstar)**2)

runs  = [ 0 ]
typeouts = [ '3pj0', '3pj1', '3pj2']
#line='LA'

#  number of missing pixels of stellar surface |z_max|<Rst
miss_pix = 13672

for run in runs:
    for typeout in typeouts:
        print (run, typeout)

        path = '/home/cvillarreal/guacho-he/hatp32/hatp32/hlines/'
        fout1 = open(path+'dat/'+typeout+'_tau_'+str(run).zfill(3)+'.dat',"w")
        if(typeout=='LA'):
            fout1.write('k'+" "+'blue_wing'+" "+'red_wing'+"\n")

        Emission_missing_pixel = np.zeros(nvmap)
        Emission_missing_pixel[:] = 1.

        for nout in range(0,81):

            filein= path+typeout+'_tau-'+str(run).zfill(3)+'_'+str(nout).zfill(3)+'.bin'
            f=FortranFile(filein,'r')
            data=f.read_reals(kind).reshape(nxmap,nymap,nvmap,order='F').T
            emtaunu=np.exp(-data)

            EmissionVel=np.zeros(shape=(nvmap))
            TotalEmission=0.
            Ref = 0.
            
            theta_plan= (nout-40)*np.pi/180
            yplan=0.
            xplan= rorb/dx*np.tan(theta_plan)

          
            for i in range(nxmap):
                for j in range(nymap):
                    x = (float(i-nxmap/2)+0.5)
                    y = (float(j-nymap/2)+0.5)
                    r1=np.sqrt(x**2+y**2)
                    
                    if (r1<=rstar):
                        r2= np.sqrt((x - xplan)**2 + (y-yplan)**2)

                        if (r2 < rplan):
                            emtaunu[:,i,j]=0.
                        else:
                            EmissionVel=EmissionVel[:] + emtaunu[:,j,i]
                            TotalEmission=TotalEmission + np.sum(emtaunu[:,j,i])
                        Ref = Ref + 1.
            #adding missing pixels from not accounting for the upper and lower edges of the star
            Ref = Ref + miss_pix
            EmissionVel=(EmissionVel[:] + miss_pix*Emission_missing_pixel[:])/Ref
            TotalEmission = (TotalEmission + miss_pix*250)/Ref/float(nvmap)
            x= np.linspace(velinit,velfinal,nvmap)

            fout2= open(path+'dat/'+typeout+'_tau-'+str(run).zfill(3)+'-abs-'+str(nout).zfill(3)+'.dat',"w")
            for l in range(nvmap):
                fout2.write(str(x[l])+' '+str(EmissionVel[l])+ "\n")
            fout2.close()
            print ('total absorption for output', nout,'is :', (1.-TotalEmission)*100)
            
           if(typeout=='LA'):
               fout1.write(str(nout)+" "+str(np.mean(EmissionVel[np.where(np.logical_and(x<= -40.e5, x>=-120.e5))]))+" "+
               str(np.mean(EmissionVel[np.where(np.logical_and(x<=120.e5, x>=30.e5))])) +"\n")
           else:
               fout1.write(str(nout)+" "+str(np.mean(EmissionVel))+ "\n")
   
         fout1.close()
