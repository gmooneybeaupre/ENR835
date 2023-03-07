#coding: latin-1
#
# exemple 9.3
#
# Exemple de la méthode f-chart
#
#
from solar_mod import *
import numpy as np

j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])      # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
#
# paramètres du capteur lors de l'essai
#

phi = 43.0
Fr_tan = 0.79
Fr_UL = 3.8
Htm = 10.7*1e6  # Radiation solaire mensuelle plan incliné
UA = 145.0      # coefficient de pertes
Tinf = -8.0
Lecd = 2.06*1e9     # charge mensuelle ECD en Joules
DJ = 828
Lchauf = UA*DJ*24.0*3600 # charge mensuelle chauffage en Joules
L = Lecd + Lchauf           # charge mensuelle totale en Joules
Vref = 75.0
Ac = 22.0
V1 = 1650/Ac
ta_tan = 0.96
N = 31
Dt = N*24*3600
Tref = 100
epCmin_UA = 2
#
# a)
#
X = Fr_UL*(Tref-Tinf)*Dt*Ac/L
Xa = X*(V1/75.0)**(-0.25)
Ya = Fr_tan*ta_tan*Htm*N*Ac/L
fa = fchart(Xa,Ya)
print ('f a) = ',fa)
# b)
#
Xb = Xa*(0.5)**(-0.25)
fb = fchart(Xb,Ya)
print ('f b) = ',fb)
#
# c
#
epCmin_UAc = 0.5
fcY = 0.39+0.65*np.exp(-0.139/epCmin_UAc)
Yc=Ya*fcY
fc1 = fchart(Xa,Yc)
print ('f c)= ',fc1)
fc2 = fchart(Xb,Yc)
print ('f c)= ',fc2)



