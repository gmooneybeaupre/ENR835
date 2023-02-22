#coding: utf-8
#
# exemple solaire 4.1
#
# Exemple du calcul de la radiation absorbée par un capteur vitre (1 vitre)
#
#
from solar_mod import *
import numpy as np
from properties_mod import *
#
# données du problème
#
It = 2.8    # MJ/m2 M2
Itb = 2.0    # MJ/m2 M2
Itd = 0.6    # MJ/m2 M2
beta = 50
theb = 45
Tin = 30.0
Tamb = 0.0
# donnes capteur
etao = 0.717
a1 = 4.0141   # W/m2 K
a2 = 0.01872
Y = 0.726
m = -5.513
Ac = 2.874      # surface totale
Aa = 2.69       # surface d'ouverture
pourcentage = 0.2
Tref = Tin  + 273.15
debit_vol_test = 19.5/1000*Ac/1000
Cp = propylene_glycol_prop('Cp',Tref,pourcentage)
rho = propylene_glycol_prop('rho',Tref,pourcentage)
mp =debit_vol_test*rho
#
# Calcul du (tau*alpha) pour la
# radiation directe
#
Itr = It - Itb - Itd
# debut
mod = 2  #
if mod ==2:
    bo = 0.11
    b1 = 0.051
elif mod ==1:
    bo = 0.16
    b1 = 0

def  kta(the):
    if the > 60:
        kta = 1 - bo - b1
        kta = kta*(1 - (the - 60.0)/30.0)
    else:
        kta = 1 - bo*(1/cosd(the)-1) - b1*(1/cosd(the)-1)**2
    return kta
ktab = kta(theb)
thed = angle_diffus(beta)
ktad = kta(thed)
ther = angle_reflechi(beta)
ktar = kta(ther)
ktam = (ktab*Itb + ktad*Itd + ktar*Itr)/It
print ('ktam = ',ktam)
Itw = It*1000.0/3.6         # donnes en Watts/m2
etaoc = etao*ktam
eta = etaoc - a1*(Tin - Tamb)/Itw - a2*(Tin - Tamb)**2/Itw
print ('rendement  = ',eta)
qu = eta*Itw*Ac
Tout = Tin + qu/(mp*Cp)
Tm = (Tin + Tout)/2
print ('Tout   = ',Tout,' C')
print ('Tmoyenne   = ',Tm,' C')
print ('qutile (a) = ',qu,' Watts')
print ('qutile (a) = ',qu*3.6,' kJ/hr')
