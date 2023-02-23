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


def Calcul_the(thez,gams,theb,beta):
    theT = abs(arctand(sind(thez)*sind(abs(gams-gam))/cosd(theb)))
    theL = abs(arctand(tand(thez)*cosd(abs(gams-gam)))-beta)
    return theT,theL

theTv = np.array([0,10,20,30,40,50,60,90])
theLv = theTv
KtheTv  = np.array([1.0,1.01,1.02,1.04,1.04,0.99,0.90,0])
KtheLv  = np.array([1.0,1.0,0.99,0.97,0.95,0.91,0.83,0])
#
# données du problème
#
beta = 50
phi = 45.5
n = jour_mois_jour_annee(9,'juil')
delt  = decl_solaire(n)      # calcul de la declinaison solaire
omen = 45+7.5   # 3: 30
gam = 0
theb = normale_solaire(delt,phi,omen,beta,gam)   # derection de la radiation directe
thez = zenith_solaire(delt,phi,omen)     # zenith solaire
gams = azimuth_solaire(thez,delt,phi,omen)
print(theb,thez,gams)

theT,theL = Calcul_the(thez,gams,theb,beta)
print(theT,theL)

It = 2.8    # MJ/m2 M2
Itb = 2.8    # MJ/m2 M2
Itd = 0.0    # MJ/m2 M2
Tin = 30.0
Tamb = 0.0
# donnes capteur
etao = 0.761
a1 = 1.36   # W/m2 K
a2 = 0.0074
Aa = 2.16       # surface d'ouverture
Tref = Tin  + 273.15
mp = 0.02*Aa
Cp = 4180
rho = 1000
#
KthT = np.interp(theT,theTv,KtheTv)
KthL = np.interp(theL,theLv,KtheLv)
Kt = KthT*KthL
print(KthL,KthT,Kt)

# radiation directe
#
Itw = It*1000/3.6
etaoc = etao*Kt
ok = False
compt = 0
compt_max = 10
Tout = Tin + 4
delta = 0.01
Tm = (Tin + Tout)/2
while not ok:
    eta = etaoc - a1*(Tm - Tamb)/Itw - a2*(Tm - Tamb)**2/Itw
    qu = eta*Itw*Aa
    Toutn = Tin + qu/(mp*Cp)
    Tm = (Tin + Toutn)/2
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

print ('rendement  = ',eta)
qu = eta*Itw*Aa
Tout = Tin + qu/(mp*Cp)
Tm = (Tin + Tout)/2
print ('Tout   = ',Tout,' C')
print ('Tmoyenne   = ',Tm,' C')
print ('qutile (a) = ',qu,' Watts')
print ('qutile (a) = ',qu*3.6,' kJ/hr')
