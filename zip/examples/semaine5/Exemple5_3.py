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
P = 2*(H+Y)
Lair  = 0.025       # epaisseur d'air
Lv = 0.002          # épaisseur de la vitre
kv = 1.4            # conductivité de la vitre
e1 = 0.95           # emissivité de l'absorbeur
e2 = 0.88           # emissivité de la vitre
T1c = 100.0           # Temperature de la plaque (Celsius)
T1k = T1c+273        # Temperature de la plaque (Kelvin)
uinf = 4
Tinfc = 10
Tinfk = Tinfc + 273
Tskyc = 0
Tskyk = Tskyc + 273
#
ok = False
compt = 1
compt_max = 30
delta = 0.01
T2c = 60
while not ok:
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
    #
    # pertes extérieures
    #
    # Calcul du coef de convection externe
    Rpp_vitre = Lv/kv
    T3c = T2c - qppi*Rpp_vitre
    T3k = T3c + 273
    Recr = 500000
    Textk = (T3k+Tinfk)/2.0
    nue = air_prop('nu',Textk)
    ale = air_prop('al',Textk)
    ke = air_prop('k',Textk)
    Pr = air_prop('Pr',Textk)
    Lc = 4*Ac/P
    #
    Re = uinf*Lc/nue
    if Re< Recr:
        Nu = 0.86*np.sqrt(Re)*Pr**(1.0/3.0)     # laminaire
    else:
        Nu = (0.037*Re**(4.0/5.0)-871.0)*Pr**(1.0/3.0)     # turbulent
    hconve = Nu*ke/Lc
    #
    # Calcul du coefficient de radiation externe
    #
    hrade = e2*sig*(T3k**2+Tskyk**2)*(T3k+Tskyk)
    Rcpp = 1/hconve
    Rrpp = 1/hrade
    #
    # fin du calcul des résistances thermiques
    # vérification des bilans thermiques avec hypothese
    #
    qconve = hconve*(T3c-Tinfc)
    qrade = hrade*(T3c-Tskyc)
    qppe = qconve+qrade
    print ('q" ext radiation ',qrade,' W/m2')
    print ('q" ext convection ',qconve,' W/m2')
    print ('q" ext totales ',qppe,' W/m2')

    #
    # Bilan radiatif à la surface externe de la vitre
    #
    Rintpp = Ripp + Rpp_vitre
    C1 = 1/Rintpp +  1/Rcpp + 1/Rrpp
    T3n = (T1c/Rintpp + Tinfc/Rcpp + Tskyc/Rrpp)/C1
    # vérification des bilans thermiques avec nouvelle valeur de T2
    qconve = hconve*(T3n-Tinfc)
    qrade = hrade*(T3n-Tskyc)
    qppen = qconve+qrade
    qppin = (T1c - T2c)/Ripp
    qppv= (T2c - T3n)/Rpp_vitre

    diff = abs(T3n - T3c)
    if diff < delta:
        ok = True
    else:
        T2c = T3n + qppin*Rpp_vitre
        compt = compt + 1
        if compt > compt_max:
            ok = True
            erreur = 1
print(qppin,qppe,qppv)
Uhaut = qppin/(T1c-Tinfc)
print ('T vitre ',T2c,T3c,T3n)
print ('U haut ',Uhaut,' W/m2-k')
print ('q" int correct ',qppin,' W/m2')
print ('q" ext correct ',qppen,' W/m2')
T2m,T3m,Uhautm =  Calcul_pertes(T1c,beta,H,Y,uinf,Tinfc,Tskyc,Lair,e1,e2,Lv,kv)
print ('T vitre ',T2m,T3m)
print ('U haut ',Uhautm,' W/m2-k')
