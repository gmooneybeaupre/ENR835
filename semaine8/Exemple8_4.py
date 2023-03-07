#coding: utf-8
#
# exemple 8.4
#
# Exemple de la méthode f-chart
#
#
from solar_mod import *
import numpy as np
from scipy.optimize import newton,brentq
from matplotlib.pyplot import *

j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])      # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
#
# paramètres du capteur lors de l'essai
#

phi = 43.0
Fr_tan = 0.60
Fr_UL = 4.0
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
fvoulue = 0.5

#
# a)
#
def fct1(A):
    X = Fr_UL*(Tref-Tinf)*Dt*A/L
    Y = Fr_tan*ta_tan*Htm*N*A/L
    f = fchart_air(X,Y)
    y = f - fvoulue
    return y
Ac = newton(fct1,45)
print ('Ac = ','%.2f' % Ac, 'm2')
#
# vérification
#
X = Fr_UL*(Tref-Tinf)*Dt*Ac/L
Y = Fr_tan*ta_tan*Htm*N*Ac/L
f = fchart_air(X,Y)
print ('f janvier = ','%.2f' % f)