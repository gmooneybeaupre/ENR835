#coding: utf-8
#
# exemple solaire 7.2
#
#
#
from solar_mod import *
import  numpy as np

#
# données du problème
#
Trefech = 300.0
Trefco = 310.0
mpcoll  = 0.08
mpech = 0.1
Cpech  = eau_prop('Cp',Trefech)
Cpcoll  = propylene_glycol_prop('Cp',Trefco,0.30)
Ccoll = mpcoll*Cpcoll
Cech = mpech*Cpech
Cmin = min(Ccoll,Cech)
Cmax = max(Ccoll,Cech)
Cr = Cmin/Cmax
UA = 500.0
Fr = 0.9
Ac = 8.0
UL = 7.0
NTU = UA/Cmin
#
# efficacité échangeur contre-courant
#
eff = (1-np.exp(-NTU*(1-Cr)))/(1-Cr*np.exp(-NTU*(1-Cr)))
z = 1 + Ac*Fr*UL/Ccoll*(Ccoll/(eff*Cmin) - 1.0)
Frp = Fr/z
print ('FRp = ',Frp)


