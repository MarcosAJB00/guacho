'''
 ! ================ HeIM =================
    ! Producción
    q(iHeIM) = + rate(iche13a)  * y(ie )* y(iHeIS)   &
               + rate(iaheiim_b)* y(ie )* y(iHeII)


    ! Pérdida
    d(iHeIM) = + rate(irheii)   * y(iHeIM)           &
               + rate(iche31)   * y(iHI)* y(iHeIM)   &
               + rate(iche31a)  * y(ie )* y(iHeIM)   &
               + rate(iche31b)  * y(ie )* y(iHeIM)   &
               + rate(iphiHeIM) * y(iHeIM )          &
               + rate(icheIM)   * y(ie) *y(iHeIM)

! fitted values for gammaij from bray2000 only valid 3.75<logT<5.75 [K]
     logT = log10(T)
     if (logT >= 3.75 .and. logT<=5.75) then
       gamma13  = - 1.5699 + 1.303*log10(T) - 0.388*log10(T)**2 &
                    + 0.0516*log10(T)**3 - 0.0026*log10(T)**4
       gamma31a = - 209.42 + 170.03*log10(T)- 50.073*log10(T)**2  &
                    + 6.42*log10(T)**3 - 0.3045*log10(T)**4
       gamma31b = - 36.477 + 22.988*log10(T) - 4.5935*log10(T)**2 &
                    + 0.2968*log10(T)**3
     else
       gamma13  = 0.
       gamma31a = 0.
       gamma31b = 0.
     endif
     ! collisional ionisation rates for H (raga) and He (Lampon 2020)
     rate(ichi   ) = 5.830e-11*sqrt(T)*exp(-157828.0/T)    !colissional ioniz. of HI [cm3/s]
     rate(iche13a) = 2.10e-8*sqrt(13.6/(8.617333262145e-5*T))*exp(-19.81/(8.617333262145e-5*T))*gamma13    !colissional ioniz. of HeI [cm3/s]
     rate(iche31a) = 2.10e-8*sqrt(13.6/(8.617333262145e-5*T))*exp( -0.80/(8.617333262145e-5*T))*gamma31a/3 !colissional ioniz. of HeII [cm3/s]
     rate(iche31b) = 2.10e-8*sqrt(13.6/(8.617333262145e-5*T))*exp( -1.40/(8.617333262145e-5*T))*gamma31b/3 ![cm3/s]
     rate(iche31)  = 5.e-10                         !HeI colisions with HI
     rate(icheIS)  = 2.379e-11*sqrt(T)*exp(-285335.4/T)     !colisional ionisation of helium 11S (cen 1992)
     rate(icheIM)  = 8.335e-10*sqrt(T)*exp(-55338.0/T)      !colisional ionisation of helium 23S (Black 1981)

     ! recombination rates for He case B from Benjamin+1999
     rate(iahii_b)   = 2.59e-13*(1.0e4/T)**0.700     !recombination of HII [cm3/s]
     rate(iaheii_a)  = 4.27e-13*(1.0e4/T)**0.678     !recombination of HeII case A tot [cm3/s]
     rate(iaheii_b)  = 2.72e-13*(1.0e4/T)**0.789     !recombination of HeII case B tot [cm3/s]
     rate(iaheii_1)  = 1.54e-13*(1.0e4/T)**0.486     !recombination of HeII to ground only [cm3/s]
     rate(iaheiim_b) = 2.10e-13*(1.0e4/T)**0.778     !recombination of HeII to 2s3 case B [cm3/s]

     rate(irheii ) = 1.272e-4           !radiative transition He(23s) to He(11s)

     !charge exchange (okclopick, lampon)
     rate(iceheIS) = 1.25e-15*(300./T)**(-0.25)
     rate(iceheII) = 1.75e-11*(300./T)**0.75*exp(-128000.0/T)

     rate(iphiHI  ) = phHI             !photoionization of HI
     rate(iphiHeIS) = phHeIS           !photoionization of HeI from ground
     rate(iphiHeIM) = phHeIM           !photoionization of HeI from triplet

'''
import numpy as np
import matplotlib.pyplot as plt


T_profile = np.loadtxt('./inputs/temperature_profile_ATES.dat', unpack=True, skiprows=1)
T = T_profile[1,:]

data = np.loadtxt('./output/perfiles_finales.dat', unpack=True, skiprows=1)
r, HI, HII, HeIS, HeIM, HeII, e_ = data

phi_profile = np.loadtxt('./output/phi_tau_finales.dat', unpack=True, skiprows=1)
phiHeIM = phi_profile[3,:]

# Rates del Metaestable
gamma13  = np.zeros_like(T)
gamma31a = np.zeros_like(T)
gamma31b = np.zeros_like(T)

mask = (np.log10(T) >= 3.75) & (np.log10(T)<=5.75)
gamma13[mask]  = - 1.5699 + 1.303*np.log10(T[mask]) - 0.388*np.log10(T[mask])**2 + 0.0516*np.log10(T[mask])**3 - 0.0026*np.log10(T[mask])**4
gamma31a[mask] = - 209.42 + 170.03*np.log10(T[mask])- 50.073*np.log10(T[mask])**2 + 6.42*np.log10(T[mask])**3 - 0.3045*np.log10(T[mask])**4
gamma31b[mask] = - 36.477 + 22.988*np.log10(T[mask]) - 4.5935*np.log10(T[mask])**2

# collisional ionisation rates for H (raga) and He (Lampon 2020)
K = 8.617333262145e-5 # eV/K
rate_iche13a = 2.10e-8*np.sqrt(13.6/(K*T))*np.exp(-19.81/(K*T))*gamma13    #colissional ioniz. of HeI [cm3/s]
rate_iche31a = 2.10e-8*np.sqrt(13.6/(K*T))*np.exp( -0.80/(K*T))*gamma31a/3 #colissional ioniz. of HeII [cm3/s]
rate_iche31b = 2.10e-8*np.sqrt(13.6/(K*T))*np.exp( -1.40/(K*T))*gamma31b/3 #[cm3/s]
rate_iche31 = 5.e-10                         #HeI colisions with HI
rate_icheIM = 8.335e-10*np.sqrt(T)*np.exp(-55338.0/T)      #colisional ionisation of helium 23S (Black 1981)

rate_iaheiim_b = 2.10e-13*(1.0e4/T)**0.778     #recombination of HeII to 2s3 case B [cm3/s]

rate_irheii = 1.272e-4           #radiative transition He(23s) to He(11s)
rate_iphiHeIM = phiHeIM           #photoionization of HeI from triplet

#tasas de recombinacion (1/s)
rate_recom = rate_iaheiim_b*e_

#tasa de decaimiento radiativo (1/s)
#rate_decay = rate_irheii

#tasas de fotoionizacion (1/s)
rate_photo = rate_iphiHeIM

#tasas de colisiones con electrones (1/s)
rate_coll = rate_iche13a*e_ + rate_iche31a*e_ + rate_iche31b*e_ + rate_iche31*HI + rate_icheIM*e_

#tasas de colisiones con HI (1/s)
rate_coll_HI = rate_iche31*HI

plt.figure(figsize=(8, 6))
plt.plot(r, rate_recom, label='Recombination Rate')
#plt.plot(r, rate_decay, label='Radiative Decay Rate')
plt.plot(r, rate_photo, label='Photoionization Rate')
plt.plot(r, rate_coll, label='Collisional Rate (e-)')
plt.plot(r, rate_coll_HI, label='Collisional Rate (HI)')
plt.xlabel('Radius (Rp)')
plt.ylabel('Rates (1/s)')
plt.title('Rates Profile for HeIM')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('rates_profile_HeIM.png', dpi=300,bbox_inches='tight')