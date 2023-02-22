#coding: utf-8
#
#
# exemple solaire 6.2
#
# Exemple de calcul des pertes thermiques d'un capteur plan vers le haut
#
#
from sys import *
from solar_mod import *
import numpy as np

#
# donnees du probleme
#
It = 2.82           # MJ/m2/hr
S = 2.22            # MJ/m2/hr
Sw = S*1000.0/3.6                # S en Watts/m2
Itw = It*1000.0/3.6                # S en Watts/m2
g = 9.8
sig = 5.67e-8
beta = 45.0
H = 1.0
Y = 2.0
Ac = H*Y            # surface du capteur
P = 2*(H+Y)         # périmètre du capteur
uinf = 4.0          # vitesse du vent
Lair  = 0.025       # epaisseur d'air
e1 = 0.95           # emissivité de l'absorbeur
e2 = 0.88           # emissivité de la vitre
Tinfc = 10.0;         # température de l'air ambiant (Celsius)
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
N = int(round(H/W))
mp1 = 0.008
mpt = mp1*N
Tfi = 80.0
Tf = 90.0
UL = 7.1   # vient de l'exemple 5.4
# Calcul des coefficeint F
m = np.sqrt(UL/(ka*deltaa))
x = m*(W-D)/2
# Calcul du rendement d'ailette F
F = np.tanh(x)/x
# Propriétées de l'eau
Tfk = Tf + 273.0
mu = eau_prop('muf',Tfk)
Pr = eau_prop('Prf',Tfk)
kf = eau_prop('kf',Tfk)
Cp = eau_prop('Cpf',Tfk)
Re = 4.0*mp1/(pi*D*mu)
if Re<2300:      # laminaire
    Nud = 3.66
else:
    f = (0.790*np.log(Re)-1.64)**(-2)
    Nud = (f/8.0)*(Re-1000.0)*Pr/(1.0+12.7*np.sqrt(f/8.0)*(Pr**(2.0/3.0)-1.0))
hf = Nud*kf/D
Rpconv = 1.0/(hf*pi*D)
den = W*(1.0/(UL*(D+(W-D)*F))+Rpjoint+Rpconv)
# Calcul du rendement d'absorbeur F'
Fp = 1.0/(UL*den)
print ('F = ',F)
print ('Fp = ',Fp)
Fp = 1.0/(UL*den)
X = UL*Ac*Fp/(mpt*Cp)
# Calcul du facteur de récupération FR
Fpp = 1/X*(1-np.exp(-X))
Fr = Fp*Fpp
print ('Fpp = ',Fpp)
print ('Fr = ',Fr)
qupp = Fr*(Sw - UL*(Tfi-Tinfc))
qu = qupp*Ac
print(qupp)
eta = qupp/Itw
Tfon = Tfi + qu/(mpt*Cp)
Tfn = Tfi+qupp/(Fr*UL)*(1-Fpp)
Tpn = Tfi+qupp/(Fr*UL)*(1-Fr)
Tfl = (Tfi + Tfon)/2
print('rend = ',eta)
print('T fluide out = ',Tfon)
print('T plaque  moy = ',Tpn)
print('T fluide moy = ',Tfn)
print('T fluide moy = ',Tfl)



#
# 1ère hypothèse
#
Tpc = 100
Tfo =  100

# Calcul du Uhaut
T2,Uhaut =  Calcul_pertes(Tpc,beta,H,Y,uinf,Tinfc,Tskyc,Lair,e1,e2)
Ubas = kiso/Liso
Ucote = Ubas*Acote/Ac
# Calcul du Utotal
UL = Uhaut+Ubas+Ucote
# Calcul des coefficeint F
m = np.sqrt(UL/(ka*deltaa))
x = m*(W-D)/2
# Calcul du rendement d'ailette F
F = np.tanh(x)/x
# Propriétées de l'eau
Tfk = (Tfi+Tfo)/2.0 + 273.0
mu = eau_prop('muf',Tfk)
Pr = eau_prop('Prf',Tfk)
kf = eau_prop('kf',Tfk)
Cp = eau_prop('Cpf',Tfk)

Re = 4.0*mp1/(pi*D*mu)
if Re<2300:      # laminaire
    Nud = 3.66
else:
    f = (0.790*np.log(Re)-1.64)**(-2)
    Nud = (f/8.0)*(Re-1000.0)*Pr/(1.0+12.7*np.sqrt(f/8.0)*(Pr**(2.0/3.0)-1.0))
hf = Nud*kf/D
Rpconv = 1.0/(hf*pi*D)
den = W*(1.0/(UL*(D+(W-D)*F))+Rpjoint+Rpconv)
# Calcul du rendement d'absorbeur F'
Fp = 1.0/(UL*den)
zz = UL*Ac*Fp/(mpt*Cp)
# Calcul du facteur de récupération FR
Fpp = 1/zz*(1-np.exp(-zz))
Fr = Fp*Fpp
UDT = UL*(Tfi-Tinfc)

quppa= Fr*(Sw-UL*(Tfi-Tinfc))     # pertes en W/m2
print ('puissance utile (Eq 6.9) = ',quppa)
quppb =(Sw-UL*(Tpc-Tinfc))
print ('puissance utile (Eq 5.30) = ',quppb)
qu = quppa*Ac
Itw = It*1000.0/3.6
rend = quppa/Itw
# Calcul des nouvelles températures
Tfon = Tfi + qu/(mpt*Cp)
Tfn = Tfi+quppa/(Fr*UL)*(1-Fpp)
Tpcn = Tfi+quppa/(Fr*UL)*(1-Fr)
print('UL = ',UL)
print('rend = ',rend)
print('T fluide out = ',Tfon)
print('T plaque  moy = ',Tpcn)
print('T fluide moy = ',Tfn)
quppb =(Sw-UL*(Tpcn-Tinfc))
print ('puissance utile (Eq 5.30) = ',quppb)

#
#
# solution avec la fonction Calcul_rendement

rend2,qupp2,Tpcn2,Tfon2 = Calcul_rendement(Tfi,Itw,Sw,N,ka,deltaa,D,Rpjoint,mp1,Ucote,Ubas,beta,H,Y,uinf,Tinfc,Tskyc,Lair,e1,e2)
UL2 = (Sw - qupp2)/(Tpcn2 - Tinfc)
print('UL = ',UL2)
print('rend = ',rend2)
print('q utile = ',qupp2,' W/m2')
print('q utile = ',qupp2*Ac*3.6,' kJ/hr')
print('T plaque  moy = ',Tpcn2)
print('T fluide out = ',Tfon2)
