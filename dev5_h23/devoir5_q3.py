#Devoir 5
#Question 3

from solar_mod import *
import numpy as np
from matplotlib.pyplot import *

#  constantes
k = 1.381e-23
q = 1.602e-19
Tref = 25.0 + 273.15

# valeurs du capteur

N = 60
A = 1.64
Isc = 8.73
Voc = 37.5
Im = 8.04
Vm = 30.5
muI = 0.00463 
muV = -0.130
Pmax = Im*Vm

#Question a
print ('Question a ')
# Calcul des valeurs initales pour la recherche des parametres
Rsh1 = 100
a1 = 1.5*k*Tref*N/q
IL1 = Isc
Io1 = (IL1-Voc/Rsh1)/(np.exp(Voc/a1)-1)
Rs1 = (a1*np.log((IL1-Im)/Io1+1) - Vm)/Im

if Rs1 < 0:
    Rs1 =0.1
# On met ces valeurs dans un vecteur
xi = np.zeros(5)
xi[0] = IL1
xi[1] = Io1
xi[2] = a1
xi[3] = Rsh1
xi[4] = Rs1
# On met les valeurs des parametres dans le vecteur  param
param = np.array([Isc,Voc,Im,Vm,muV,muI])
# Appel de la fonction pv_module qui trouve les5 valeurs du modele à 5 parametres
# param : liste de paramètres du panneau
# xi : vecteur des valeurs initiales
xf = pv_module(xi,param)

# Le vecteur xf retourne les 5 parametres nécessaire à la modélistaion du panneau PV
# resultats
IL = xf[0]
Io = xf[1]
a = xf[2]
Rsh = xf[3]
Rs = xf[4]
print ('Le courant IL est = ', IL)
print ('Le courant Io est = ', Io)
print ('Le parametre a est = ', a )
print ('La résistance Rs est = ', Rs)
print ('La résistance Rsh est = ', Rsh)


#Question b
print ('Question b ')

Rload = 3
Gb = 900
Tb = 30 + 273.15

Ib = I_pvR(xf,Rload,G = Gb,T = Tb)
Vb = Ib*Rload
Pb = Vb*Ib

print ('La puissance délivrée par le panneau est ',Pb,' Watts')

#Question c
print ('Question c ')
Vc = 16
Pc = 90
Tc = 25 + 273.15
Ic = Pc/Vc
Gc = G_pvI(xf,Ic,Vc,T = Tc)

print ('Lirradiation incidente est ',Gc,' W/m2')


