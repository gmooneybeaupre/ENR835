#coding: utf-8
#
#
# exemple solaire 6.3
#
# Exemple de calcul des pertes thermiques d'un capteur plan vers le haut
#
#

from solar_mod import *
import numpy as np
from matplotlib.pyplot import *
#
# donnees du probleme
#
modele = 'iso'
#
# données du problème
#
gam = 0.0        # azimuth du capteur
beta = 45.0
H = 1.0
Y = 2.0
Ac = H*Y            # surface du capteur
P = 2*(H+Y)         # périmètre du capteur
uinf = 4.0          # vitesse du vent
Lair  = 0.025       # epaisseur d'air
e1 = 0.95           # emissivité de l'absorbeur
e2 = 0.88         # emissivité de la vitre
e3 = 0.88
Tinfc = 10.0        # température de l'air ambiant (Celsius)
Tinfk = Tinfc + 273.0       # température de l'air ambiant (Kelvin)
Tskyc = 0.0          # température du ciel   (Celsius)
Tskyk = Tskyc+273.0  # température du ciel   (Kelvin)
Liso = 0.1
kiso = 0.04
Acote = 0.6
W = 0.166
D = 0.01
Rpjoint = 0
deltaa = 0.0005
ka = 385
#ka = 15
N = int(round(H/W))
mp1 = 0.008
mpt = mp1*N
Ubas = kiso/Liso
Ucote = Ubas*Acote/Ac
It = 2.82  # MJ/m2 K
S = 2.25  # MJ/m2 K
tau_al_moy = S/It
print ('tau_al_moy = ',tau_al_moy)
Sw = S*1000.0/3.6  # ew W/m2
Itw = It*1000.0/3.6 # ew W/m2

#
# Calcul du rendement pour différentes températures d'entrée
#
Tfiv = np.arange(Tinfc+1,85,1)
nTfi = len(Tfiv)
rendv = np.zeros(nTfi)
ULv = np.zeros(nTfi)
for i in range(0,nTfi):
    Tfi = Tfiv[i]
    rendv[i],qupp2,Tpcn2,Tfon2 = Calcul_rendement(Tfi,Itw,Sw,N,ka,deltaa,D,Rpjoint,mp1,Ucote,Ubas,beta,H,Y,uinf,Tinfc,Tskyc,Lair,e1,e2,e3 = e3,flag_turb = 0)
    ULv[i] = (Sw - qupp2)/(Tpcn2 - Tinfc)
x = (Tfiv - Tinfc)
#
# évaluation des paramètres de capteurs par régression
#
p = np.polyfit(x,rendv,2)
a2 = -p[0]*Itw
a1 = - p[1]*Itw
etao = p[2]
print (etao,a1,a2)
Tfiv2 = np.arange(Tinfc+1,85,4)
x2 = Tfiv2 - Tinfc
# donnes capteur   Enerwors ( cours 5)
etaob = 0.717
a1b = 4.0141   # W/m2 K
a2b = 0.01872
renv2 = etaob - a1b*x2/Itw - a2b*x2**2/Itw
renv3 = etao - a1*x2/Itw - a2*x2**2/Itw
#plot(Tfiv,rendv,Tfiv2,renv2,Tfiv2,renv3,'*')
Itcr = ULv*(Tfiv-Tinfc)/tau_al_moy
plot(Tfiv,rendv,Tfiv2,renv2)
legend(('capteur theorique','Capteur Enerworks'))
show()
