#coding: utf-8
#
# exemple solaire 7.3
#
#
#
from solar_mod import *
import numpy as np

#
# données du problème
#
Ti = 10.0
Ac = 10.0
Vc = 0.5*Ac     # litres
Vtcond = 20.0
Vtot = Vtcond + Vc
H = 15.0
Pmin = 50.0
Pmax = 500.0    # kPa relatif
Patm = 100.0
DTnom = 130.0
n = 0.085
g = 9.8
gam = n/DTnom
Ve = n*Vtot
Tgel  = propylene_glycol_prop('temp',300,0.40)
rho  = propylene_glycol_prop('rho',300,0.40)
DT2 = Ti - Tgel
n2 = gam*DT2
Vre = n2*Vtot
Vre = max(Vre,3)
Vmax = Vc + Ve + Vre
Pgaz = Pmin + rho*g*H/1000.0
Vnom = Vmax*(Pmax+Patm)/(Pmax - Pgaz)
print ('Vnominal = ',Vnom)
Prea = (Pmax+Patm)*(Vnom - Vmax)/(Vnom  - Vre)
Prer = Prea - Patm
print ('P reservoir = ',Prer)
Prea2 = (Pgaz + Patm)*(Vnom)/(Vnom  - Vre)
Prer2 = Prea2 - Patm
print ('P reservoir = ',Prer2)