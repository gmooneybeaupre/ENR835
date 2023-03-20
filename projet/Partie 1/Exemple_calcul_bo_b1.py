#coding: utf-8
#
# exemple calcul bo b1 à partir d'un tableua de données
#
#
from  solar_mod import *
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.pyplot import *

the = np.arange(10.0,90,10.0)
the2 = np.arange(1.0,90,2.0)
the  = np.array([0,10,20,30,40,50,60,70,80,89])
iam  = np.array([1.,1.0,0.99,0.98,0.96,0.92,0.86,0.71,0.3,0.0])
#
# on ne garde que pour theta < 60
x = 1.0/cosd(the[0:7]) - 1.0
y = iam[0:7]
#
# calcul en prenant une seule valeur
#
#
#
def  kta(the,bo,b1):
    if the > 60:
        kta = 1 - bo - b1  # kta à 60 degres
        kta = kta*(1 - (the - 60.0)/30.0)
    else:
        kta = 1 - bo*(1/cosd(the)-1) - b1*(1/cosd(the)-1)**2
    return kta

def calcul_bo_b1(x,bo,b1):
    yn  = 1.0 - bo*x  -b1*x**2
    return  yn
#
#
# régression polynomiale
#
param,resn = curve_fit(calcul_bo_b1,x,y)
bo = param[0]
b1 = param[1]
n = len(the2)
iam2 = np.zeros(n)
for i in range(0,n):
    iam2[i] = kta(the2[i],bo,b1)
print ('bo ( regression 2nd)  = ', bo)
print ('b1 ( regression 2nd)  = ', b1)
plot(the2,iam2,the,iam,'*')
legend(('régression','expérimental'),loc =3)
show()
