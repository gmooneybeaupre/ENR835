#coding: utf-8
#
# exemple solaire 6.1
#
# Exemple de calcul des pertes thermiques d'un capteur plan
#
#
from sys import *
from solar_mod import *
import numpy as np
from matplotlib.pyplot import *

#
# donnees du probleme
#
Y = 2
H = 1
beta = 45.0
W = 0.166          # distance entre les tubes
D = 0.01
Rpjoint = 0
deltaa = 0.0005
ka = 385
N = int(round(H/W))
mp1 = 0.008
mp = mp1*N
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

