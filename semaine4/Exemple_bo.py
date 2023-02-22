#coding: utf-8
#
# exemple calcul bo b1
#
#
from solar_mod import *
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.pyplot import *

the = np.arange(0.0,90,10.0)
the2 = np.arange(1.0,82,2.0)
y  = np.array([1.0,1.0,0.99,0.98,0.96,0.92,0.86,0.71,0.3])
x = 1.0/cosd(the) - 1.0
x2 = 1.0/cosd(the2) - 1.0
#
# calcul en prenant une seule valeur

#
# exemple  de régression polynomiale simple d'ordre 1
def fct(x,bo,b1):
    return  1.0 - bo*x - b1*x**2
#
# valeur initiale
#

#
po = np.array([0.1,0])
params,resn = curve_fit(fct,x,y,po)
bo = params[0]
b1 = params[1]
y2  = 1.0 - bo*x2 - b1*x2**2
print ('bo ( regression 2nd)  = ', bo)
print ('b1 ( regression 2nd)  = ', b1)
plot(the2,y2,the,y,'*')
legend(('bo 0','bo 1','ordre 2','exp'),loc =3)
show()
