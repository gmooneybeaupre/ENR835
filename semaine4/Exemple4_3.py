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
mp =debit_vol_test*rho
#
mp_reel = mp  # kg/s
Cpreel = Cptest
mptestAc = mp/Ac
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

cas = 2
if cas ==2:
    FrUL = a1 + a2*DT
else:
    FrUL = -m


CCtestAc = mptestAc*Cptest
CCreelAc = mpreelAc*Cpreel
z = 1-FrUL/(2.0*CCtestAc)
etaos = etao/z
a1s = a1/z
a2s = a2/z
print ('a) etaos = ',etaos)
print ('a) a1s = ',a1s,' W/m2 K')
print ('a) a2s = ',a2s,' W/m2 K2')
print ('a) a1s = ',a1s*3.6,' kJ/hr-m2 K')
print ('a) a2s = ',a2s*3.6,' kJ/hr-m2 K2')
etaosc = etaos*ktam
#
# a
#
Tm = 32.33
Itw = It*1000.0/3.6         # donnes en Watts/m2
etac = etaosc - a1s*(Tm - Tamb)/Itw - a2s*(Tm - Tamb)**2/Itw
qu = etac*Itw*Ac
Tout = Tin + qu/(mpreelAc*Ac*Cpreel)
print ('eta  = ',etac)
print ('qutile  = ',qu,' Watts')
print ('qutile  = ',qu*3.6,' KJ/hr')
print ('Tout   = ',Tout,' C')
#
ok = False
compt = 0
compt_max = 10
Tout = Tin + 4
delta = 0.01
while not ok:
    Tm = (Tin + Tout)/2
    etac = etaosc - a1s*(Tm - Tamb)/Itw - a2s*(Tm - Tamb)**2/Itw
    qu = etac*Itw*Ac
    Toutn = Tin + qu/(mpreelAc*Ac*Cpreel)
    diff = abs(Tout - Toutn)
    if diff < delta:
        ok = True
    else:
        compt = compt + 1
        Tout = Toutn
        if compt > compt_max:
            print('Erreur non convergence')
            err = 1
            ok = True
print ('eta  = ',etac)
print ('qutile  = ',qu,' Watts')
print ('qutile  = ',qu*3.6,' KJ/hr')
Tout = Tin + qu/(mpreelAc*Ac*Cpreel)
print ('Tout   = ',Tout,' C')
#
# B effet de la surface
#
etaob = etao*Ac/Aa
a1b = a1*Ac/Aa
a2b = a2*Ac/Aa
print ('b) etao = ',etaob)
print ('b) a1 = ',a1b,' W/m2 K')
print ('b) a2 = ',a2b,' W/m2 K2')
print ('b) a1 = ',a1b*3.6,' kJ/hr-m2 K')
print ('b) a2 = ',a2b*3.6,' kJ/hr-m2 K2')
#
# c effet des 2
#
etaosc = etaos*Ac/Aa
a1sc = a1s*Ac/Aa
a2sc = a2s*Ac/Aa
print ('c) etaos = ',etaosc)
print ('c) a1s = ',a1sc,' W/m2 K')
print ('c) a2s = ',a2sc,' W/m2 K2')
print ('c) a1s = ',a1sc*3.6,' kJ/hr-m2 K')
print ('c) a2s = ',a2sc*3.6,' kJ/hr-m2 K2')