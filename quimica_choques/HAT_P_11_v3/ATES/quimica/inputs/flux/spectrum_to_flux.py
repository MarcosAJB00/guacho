import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

# =========================
# CONSTANTES
# =========================
Rsun = 6.955e10     # cm
au   = 1.496e13     # cm
h    = 6.626e-27    # erg s
c    = 2.998e10     # cm/s
evtoerg = 1.602e-12
Lsun = 3.826e33     # erg/s

# =========================
# PLANETA: HD 209458b
# =========================
a = 0.04747 * au   # distancia orbital

# =========================
# LEER ESPECTRO
# =========================
# archivo: wavelength[nm], L_lambda[Lsun]
#wvl_nm, L_lambda_Lsun = np.loadtxt('lxuvdata_min_fullspectra.dat', unpack=True, skiprows=1)
wvl, F_lambda = np.loadtxt('../../../../spectrum/HP11_spectrum_beni_scaled.dat', unpack=True, skiprows=1)
# conversiones
#wvl = wvl_nm * 10.0              # nm → Angstrom
#L_lambda = L_lambda_Lsun * Lsun  # → erg/s

# flujo en el planeta
#F_lambda = L_lambda / (4 * np.pi * a**2)   # erg/s/cm2/A

# =========================
# PLOT ESPECTRO
# =========================
plt.figure()
plt.yscale('log')
plt.plot(wvl, F_lambda)
plt.xlabel('Wavelength [A]')
plt.ylabel('Flux [erg s$^{-1}$ cm$^{-2}$ A$^{-1}$]')
plt.title('Flux at planet')
plt.savefig('spectrum.png', dpi=300)
plt.show()

# =========================
# INTEGRACIONES
# =========================

# X-ray [0.1–100 Å]
mask_x = (wvl > 0.1) & (wvl < 100)
Fx = integrate.simpson(F_lambda[mask_x], wvl[mask_x])
Lx = 4*np.pi*a**2 * Fx

# EUV [100–912 Å]
mask_euv = (wvl > 100) & (wvl < 912)
Feuv = integrate.simpson(F_lambda[mask_euv], wvl[mask_euv])
Leuv = 4*np.pi*a**2 * Feuv

# XEUV [0.1–912 Å]
mask_xeuv = (wvl > 0.1) & (wvl < 912)
Fxeuv = integrate.simpson(F_lambda[mask_xeuv], wvl[mask_xeuv])
Lxeuv = 4*np.pi*a**2 * Fxeuv

print('Fx   =', Fx, 'erg/s/cm2')
print('Lx   =', Lx, 'erg/s')
print('Feuv =', Feuv, 'erg/s/cm2')
print('Leuv =', Leuv, 'erg/s')
print('Fxeuv=', Fxeuv, 'erg/s/cm2')
print('Lxeuv=', Lxeuv, 'erg/s')

# =========================
#  HIDROGENO (504–912 Å)
# =========================
nu_h = 13.6 * evtoerg / h #nu treshold for H
lambda_h = c / nu_h * 1e8  # # treshold for H [Å]

def sigma_h(l):
    return 6.3e-18 * (l / lambda_h)**3 #[cm^2]

h_range= np.where((wvl < 912) & (wvl>504))
soft_euv = integrate.trapezoid(F_lambda[h_range], wvl[h_range])

q_h_nu = F_lambda[h_range]*(lambda_h-wvl[h_range])/lambda_h*sigma_h(wvl[h_range])
j_h_nu = F_lambda[h_range]*wvl[h_range]*sigma_h(wvl[h_range])/h/c*1e-8

q_h    = integrate.trapezoid(q_h_nu, wvl[h_range])
j_h    = integrate.trapezoid(j_h_nu, wvl[h_range])

#wavelenght associated to photon energy in the wavelenght bin
lbda_soft_euv= h*c*1e8/(q_h/j_h + 13.6*evtoerg)
#cross section of H11S at that lambda
a0_h = sigma_h(lbda_soft_euv)
a0_h_threshold = sigma_h(lambda_h)

print(' ======== HIDROGENO (504–912 Å)========')
print('lmbda_h=',lambda_h)
print( 'soft euv [erg/s/cm2]= ', soft_euv)
print('Excess E per phot H (ev)=',q_h/j_h/evtoerg)
print('Excess E per phot H (erg)=',q_h/j_h)
print('wvl for photon total energy [ang]=', lbda_soft_euv)
print('a0 h [cm2]= ', a0_h)
print('a0 h threshold [cm2]= ', a0_h_threshold)
#print('average heating (e_nu - 13.6) (erg)=',(13.6*evtoerg - q_h/j_h))
print(f'# sof_euv phot H (1/cm2/s)={(soft_euv/(q_h/j_h + 13.6*evtoerg)):.3e}')
print('')



# =========================
# CARGAR SECCIONES EFICACES He
# =========================
lda_he21, sigma_he21, lda_he23, sigma_he23, lda_he11, sigma_he11 = \
    np.loadtxt('he_photoionzation_crossection.txt', unpack=True, skiprows=2)

f_sigma_he11 = interp1d(lda_he11, sigma_he11, bounds_error=False, fill_value=0)
f_sigma_he23 = interp1d(lda_he23, sigma_he23, bounds_error=False, fill_value=0)

# =========================
# HE singlete (200–504 Å)
# =========================
he11_range= np.where((wvl < 504) & (wvl >200))

nu_he11 = 24.6*evtoerg/h
wvl_he11= c/nu_he11*1e8

he11_range= np.where((wvl < 504) & (wvl >200))
hard_euv = integrate.trapezoid(F_lambda[he11_range], wvl[he11_range])

#the following calculation assumes that the heating is important at tau=1
q_he11_nu = F_lambda[he11_range]*(wvl_he11-wvl[he11_range])/wvl_he11*f_sigma_he11(wvl[he11_range])
j_he11_nu = F_lambda[he11_range]*wvl[he11_range]*f_sigma_he11(wvl[he11_range])/h/c*1e-8

q_he11    = integrate.trapezoid(q_he11_nu, wvl[he11_range])
j_he11    = integrate.trapezoid(j_he11_nu, wvl[he11_range])
#wavelenght associated to photon energy in the wavelenght bin
lbda_hard_euv= h*c*1e8/(q_he11/j_he11 + 24.6*evtoerg)
#cross section of H11S at that lambda
a0_he11 = f_sigma_he11(lbda_hard_euv)
a0_he11_threshold = f_sigma_he11(502.66)

print(' ======== HELIO S(200–504 Å)========')
print('wvl_he11=',wvl_he11)
print( 'hard euv [erg/s/cm2]= ', hard_euv)
print('Excess E per phot He11 (ev)=',q_he11/j_he11/evtoerg)
print('Excess E per phot He11 (erg)=',q_he11/j_he11)
print('wvl for photon total energy [ang]=', lbda_hard_euv)
print('a0 he11 [cm2]= ', a0_he11)
print('a0 he11 threshold [cm2]= ', a0_he11_threshold)
#print('average heating (e_nu - 24.6) (erg)=',(24.6*evtoerg - q_he11/j_he11 ))
print(f'hard_euv # phot He11 (1/cm2/s)={(hard_euv/(q_he11/j_he11 + 24.6*evtoerg)):.3e}')
print('')


# =========================
# HE triplete (912–2600 Å)
# =========================
mask_He3 = (wvl > 912) & (wvl < 2600)

nu_he23 =  4.8*evtoerg/h
wvl_he23= c/nu_he23*1e8

he23_range = np.where((wvl < 2600) & (wvl>912))
lw_euv = integrate.trapezoid(F_lambda[he23_range], wvl[he23_range])

q_he23_nu = F_lambda[he23_range]*(wvl_he23-wvl[he23_range])/wvl_he23*f_sigma_he23(wvl[he23_range])
j_he23_nu = F_lambda[he23_range]*wvl[he23_range]*f_sigma_he23(wvl[he23_range])/h/c*1e-8

q_he23    = integrate.trapezoid(q_he23_nu, wvl[he23_range])
j_he23    = integrate.trapezoid(j_he23_nu, wvl[he23_range])

#wavelenght associated to photon energy in the wavelenght bin
lbda_mid_euv= h*c*1e8/(q_h/j_h + 4.8*evtoerg)
#cross section of H11S at that lambda
a0_h23 = f_sigma_he23(lbda_mid_euv)
a0_h23_threshold = f_sigma_he23(wvl_he23)

print(' ======== HELIO T(912–2600 Å)========')
print('wvl_he23=',wvl_he23)
print( 'lw euv [erg/s/cm2]= ', lw_euv)
print('Excess E per phot He23 (ev)=',q_he23/j_he23/evtoerg)
print('Excess E per phot He23 (erg)=',q_he23/j_he23)
print('wvl for photon total energy [ang]=', lbda_mid_euv)
print('a0 he23 [cm2]= ', a0_h23)
print('a0 he23 threshold [cm2]= ', a0_h23_threshold)
#print('average heating (e_nu - 4.6) (erg)=',(4.6*evtoerg - q_he23/j_he23 ))
print(f'lw # phot He23 (1/cm2/s)={(lw_euv/(q_he23/j_he23 + 4.6*evtoerg)):.3e}')


#==========================
# PLOT FOR PAPER
#==========================
wvl_h= np.arange(0,913)

fig,ax1 = plt.subplots(figsize=(10,8))
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Flux [erg/s/cm$^2$/\AA]')
ax1.set_xlabel('Wavelength [\AA]')
ax1.set_xlim(10e1,3e3)
ax1.plot(wvl,F_lambda, c='red')
ax1.axvspan(200,504,facecolor='g', alpha=0.3)
ax1.axvspan(504,912,facecolor='b', alpha=0.3)
ax1.axvspan(912,2600,facecolor='orange', alpha=0.3)
ax1.text(220,1e-3,'Hard-EUV')
ax1.text(520,1e-3,'Soft-EUV')
ax1.text(1400,1e-3,'Mid-UV')
ax2=ax1.twinx()
ax2.set_yscale('log')
ax2.set_ylim(1e-20,1e-17)
ax2.set_ylabel('Cros-section [cm$^2$]')
ax2.plot(wvl_h, sigma_h(wvl_h), label='$a_{HI}$ ')
#plt.loglog(lda_he21, sigma_he21,label='sigma He21 ')
ax2.plot(lda_he23, sigma_he23,label='$a_{He2^3S} $')
ax2.plot(lda_he11, sigma_he11,label='$a_{He1^1S} $')

plt.legend(loc='upper left')
fig.tight_layout()
plt.savefig('spectrum_wvlbands.png',dpi=300,transparent=False,bbox_inches='tight')
plt.show()