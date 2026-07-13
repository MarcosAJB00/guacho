import numpy as np

#----------------------------------------------------
# Constantes (CGS)
#----------------------------------------------------
G    = 6.67430e-8          # cm3 g-1 s-2
MSUN = 1.98847e33          # g
AU   = 1.495978707e13      # cm

#----------------------------------------------------
# Parametros de entrada
#----------------------------------------------------
ur = 2.874e7          # velocidad radial del viento [cm/s]
Mstar = 0.69        # masa estelar [M_sun]
a = 0.055            # semieje mayor [AU]


# Si quieres incluir corrotacion (sub-Alfvenico)
Prot = 	17.1                     # dias
omega = 2*np.pi/(Prot*86400.0)  # rad/s
uphi = omega*a*AU               # cm/s

#----------------------------------------------------
# Calculo
#----------------------------------------------------
ukep = np.sqrt(G*(Mstar*MSUN)/(a*AU))

theta0 = np.degrees(np.arctan2(ur, np.abs(ukep-uphi)))

print(f"u_Kep  = {ukep/1e5:.2f} km/s")
print(f"u_phi  = {uphi/1e5:.2f} km/s")
print(f"Theta0 = {theta0:.2f} grados")
