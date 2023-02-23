#coding: utf-8
#
# exemple solaire 5.2
#
# Exemple de calcul des pertes thermiques d'un capteur plan vers le haut
#
#
from solar_mod import *
from properties_mod import *
import numpy as np

#
# donnees du probleme
#
g = 9.8
sig = 5.67e-8
beta = 45.0
H = 1.0
Y = 2.0
Ac = H*Y            # surface du capteur
Lair  = 0.025       # epaisseur d'air
e1 = 0.95           # emissivité de l'absorbeur
e2 = 0.88           # emissivité de la vitre
T1c = 100.0           # Temperature de la plaque (Celsius)
T1k = T1c+273        # Temperature de la plaque (Kelvin)
#
T2c =  60
Dt = T1c-T2c
T2k = T2c + 273.0
Tairk = (T1k+T2k)/2
# Calcul du coef de convection interne
betas = 1/Tairk
nui = air_prop('nu',Tairk)
ali = air_prop('al',Tairk)
ki = air_prop('k',Tairk)
Ra = g*betas*abs(Dt)*Lair**3/(nui*ali)
Ra_c = 1708.0/cosd(beta)
ct = cosd(beta)
f1 = max(0,1.0-1708.0/(Ra*ct))
f2 = 1.0-1708*(sind(1.8*beta))**1.6/(Ra*ct)
f3 = max(0,(Ra*ct/5830)**(1.0/3.0)-1.0)
Nui = 1.0 + 1.44*f1*f2+f3
hconvi = Nui*ki/Lair
# Calcul du coefficient de radiation interne
hradi = sig*(T1k+T2k)*(T1k**2+T2k**2)/(1.0/e1+1.0/e2-1.0)
Ripp = 1.0/(hconvi+hradi)    # résistance équivalente interne
# Calcul du coef de convection externe
qppi = (T1c - T2c)/Ripp
qconvi = hconvi*(T1c-T2c)
qradi = hradi*(T1c-T2c)
qppi2 = qconvi+qradi
print ('q" int radiation ',qradi,' W/m2')
print ('q" int convection ',qconvi,' W/m2')
print ('q" int totales ',qppi,' W/m2',qppi2,' W/m2')
