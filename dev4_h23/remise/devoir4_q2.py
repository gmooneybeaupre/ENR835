#coding: utf-8
#
from solar_mod import *
import numpy as np

j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])      # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
#
# paramètres du capteur lors de l'essai
Tinf = np.array([-10.3,-8.8,-2.4,5.7,12.9,18.0,20.8,19.4,14.5,8.3,1.6,-6.9])
HmKWh = np.array([1.5817,2.5884,3.8525,4.5808,5.6846,6.2114,6.0228,4.8676,4.0267,2.5372,1.2235,1.1742])  # Radiation plan horizontal kWh /jr/m2
Hm = HmKWh * 3.6e6     # Radiation plan horizontal Joules /jr/m2

phi = 45.5
gam = 0
beta = 60
debitactuel = 12 
capaciteactuelle = 0.15
Fr_tan = 0.76
Fr_UL = 4.85
UA = 300.0      # coefficient de pertes
Lecd = 0     # charge mensuelle ECD en Joules
DJ = np.array([874,747])
#Vref = 75.0
Ac = 40.0
#V1 = 1650/Ac
ta_tan = 0.88
Tref = 100
rhog = 0.7

# Fr'/Fr = 1 car il n'y a pas d'échangeur dans les capteurs à air
Frp_Fr = 1

Lchauf = np.array([0.0,0.0])
L = np.array([0.0,0.0])
Hmt = np.array([0.0,0.0])
f = np.array([0.0,0.0])
for m in range(0,2) :
    Dt = nm[m]*24*3600
    N = nm[m]
    Lchauf[m] = UA*DJ[m]*24*3600   # charge mensuelle chauffage en Joules
    L[m] = Lecd + Lchauf[m]        # charge mensuelle totale en Joules
    
    njt = j_type[m]                # njt est le jour type du mois
    #
    # Radiation plan incliné Joules /jr/m2
    Hom = irradiation_extraterrestre_jour(njt,phi)
    KTm = Hm[m]/Hom
    delt = decl_solaire(njt)
    omesm = angle_sunset(phi,delt)
    Hdtm = Erbs_mois(KTm,omesm)
    HmdE = Hm[m]*Hdtm          # Contribution diffus donnée par Erbs
    HmbE = Hm[m] - HmdE        # Contribution directe donnée par Erbs
    omeb = 37.5                # angle horaire correspondant à 2:30 du midi
    theb = normale_solaire(delt,phi,omeb,beta,gam)
    Rbb = Calcul_Rb_mois(phi,njt,beta,gam)
    HmtE,HmtbE,HmtdE,HmtrE = modele_isotropique(Hm[m],HmbE,HmdE,beta,Rbb,rhog)
    Hmt[m] = HmtE
    #
    # Calcul des paramètres adimensionnels X et Y
    X = Fr_UL*Frp_Fr*(Tref-Tinf[m])*Dt*Ac/L[m]
    Y = Fr_tan*Frp_Fr*ta_tan*Hmt[m]*N*Ac/L[m]
    #
    #correction du débit 
    Xc = X*(debitactuel/10)**0.28  #12/10 = 1.2 compris entre 0.5 et 2.0
    #
    #correction du stockage 
    Xcc = Xc*(capaciteactuelle/0.25)**(-0.30)  #0.15/0.25 = 0.6 compris entre 0.5 et 4.0
    #
    # Calcul de la fraction solaire
    f[m] = fchart_air(Xcc,Y)

print ('En janvier, la fraction solaire est de', f[0])
print ('En février, la fraction solaire est de', f[1])