#coding: utf-8
#
# exemple 4.2
#
#
from solar_mod import *
import numpy as np
#from matplotlib.pyplot import *
#
# données du problème
#
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
DT = Tin - Tamb
# donnes capteur
etao = 0.717
a1 = 4.0141   # W/m2 K
a2 = 0.01872
Y = 0.726
m = -5.513
Ac = 2.874      # surface totale
Aa = 2.69       # surface d'ouverture
pourcentage = 0.2
Tref = Tin + 273.15
debit_vol_test = 19.5/1000*Ac/1000
Cptest = propylene_glycol_prop('Cp',Tref,pourcentage)
rho = propylene_glycol_prop('rho',Tref,pourcentage)
mp_test =debit_vol_test*rho
#
mp_reel = 80/3600  # kg/s (eau)
Cpreel = eau_prop('Cp',Tref)
mptestAc = mp_test/Ac
mpreelAc = mp_reel/Ac
#
# debut
Itr = It - Itb - Itd
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

CCtestAc = mptestAc*Cptest
CCreelAc = mpreelAc*Cpreel
cas = 2
if cas ==2:
    FrUL = a1 + a2*DT
else:
    FrUL = -m
FpUL = -CCtestAc*np.log(1-FrUL/(CCtestAc))
xtest = CCtestAc/FpUL
xreel = CCreelAc/FpUL
r = xreel*(1-np.exp(-1/xreel))/(xtest*(1-np.exp(-1/xtest)))
etaoc = etao*ktam
etaoc = r*etao*ktam
a1c = a1*r
a2c = a2*r
print ('r = ',r)
print ('etao = ',etaoc)
print ('a1 = ',a1c,' W/m2 K')
print ('a2 = ',a2c,' W/m2 K2')
Itw = It*1000.0/3.6         # donnes en Watts/m2
etac = etaoc - a1c*(Tin - Tamb)/Itw - a2c*(Tin - Tamb)**2/Itw
qu = etac*Itw*Ac
print ('eta  = ',etac)
print ('qutile  = ',qu,' Watts')
print ('qutile  = ',qu*3.6,' KJ/hr')
Tout = Tin + qu/(mpreelAc*Ac*Cpreel)
print ('Tout   = ',Tout,' C')

