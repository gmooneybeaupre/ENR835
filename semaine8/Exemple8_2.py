#coding: utf-8
#
# exemple 8.1
#
# Exemple de la méthode f-chart
#
#
from sys import *
from solar_mod import *
import numpy as np
from matplotlib.pyplot import *

j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])      # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
#
# données du problème
#
L = np.array([36,30.4,26.7,15.7,9.2,4.1,2.9,3.4,6.3,13.2,22.8,32.5])      #  charge en GJoules
L = L*1.0e9                                                              # charge en Joules
Tinf = np.array([-10.3,-8.8,-2.4,5.7,12.9,18.0,20.8,19.4,14.5,8.3,1.6,-6.9])
HmKWh = np.array([1.5817,2.5884,3.8525,4.5808,5.6846,6.2114,6.0228,4.8676,4.0267,2.5372,1.2235,1.1742])  # Radiation plan horizontal kWh /jr/m2
Hm = HmKWh * 3.6e6     # Radiation plan horizontal Joukles /jr/m2

phi = 45.5
beta = 45.0
gam = 0.0
rhog = 0.2
a1 = 4.18           # parametre du capteur
a2 = 0.02           # parametre du capteur
etao = 0.674        # parametre du capteur
bo = 0.21           # parametre du capteur
b1 = 0.0
Fr_tan = etao
DeltaT = 0.0
Fr_UL = a1 + a2*DeltaT
ep = 0.6
Cp = 4190
mpp = 0.034         # debit nominal par m2 de capteur


thed = angle_diffus(beta)
ther = angle_reflechi(beta)
def  kta(the):
    if the > 60:
        kta = 1 - bo - b1  # kta à 60 degres
        kta = kta*(1 - (the - 60.0)/30.0)
    else:
        kta = 1 - bo*(1/cosd(the)-1) - b1*(1/cosd(the)-1)**2
    return kta
Kd = kta(thed)
Kr = kta(ther)

ta_tan1 =  np.zeros(12)
ta_tan2 =  np.zeros(12)
ta_tan3 =  np.zeros(12)
f1 = np.zeros(12)      # fraction solaire mensuelle approche 2
f2 = np.zeros(12)      # fraction solaire mensuelle approche 3
f3 = np.zeros(12)      # fraction solaire mensuelle approche 4

Hmt = np.zeros(12)
    #

for im in range(0,12):     # im index des moix
    njt = j_type[im]                 # njt est le jour type du mois
    # valeur de référence ( Eq. 8.13)
    #
    # calcul à partir des valeurs mensuelles uniquement
    #
    Hom = irradiation_extraterrestre_jour(njt,phi)
    KTm = Hm[im]/Hom
    delt = decl_solaire(njt)
    omesm = angle_sunset(phi,delt)
    Hdtm = Erbs_mois(KTm,omesm)
    HmdE = Hm[im]*Hdtm              # Contribution diffus donnée par Erbs
    HmbE = Hm[im] - HmdE        # Contribution directe donnée par Erbs
    omeb = 37.5                          # angle horaire correspondant à 2:30 du midi
    theb = normale_solaire(delt,phi,omeb,beta,gam)
    Kb = kta(theb)
    Rbb = Calcul_Rb_mois(phi,njt,beta,gam)
    HmtE,HmtbE,HmtdE,HmtrE = modele_isotropique(Hm[im],HmbE,HmdE,beta,Rbb,rhog)
    Hmt[im] = HmtE
    ta_tan1[im]  = (Kb*HmtbE +  Kd*HmtdE + Kr**HmtrE)/HmtE    # Eq 8.15 - 8.16
    ta_tan2[im] = Calcul_Kta_moyen(phi,beta,gam,rhog,bo,b1,njt,Hm[im],KTm)
    if ((im < 9) & (im > 3)):
        ta_tan3[im] = 0.92  # été
    else:
        ta_tan3[im] = 0.96   # hiver


# Calcul de la fraction solaire par la méthode f-chart
Accv = np.arange(10,55,5)
nn = len(Accv)
Frp_Fr = (1.0+a1/(Cp*mpp)*(1.0/ep-1.0))**(-1)
Tref  = 100.0
dhw = 0.0
Fra1 = np.zeros(nn)
Fra2 = np.zeros(nn)
Fra3 = np.zeros(nn)

for ip in range(0,nn):
    Acc = Accv[ip]
    for i in range(0,12):
        Dt = nm[i]*24*3600
        N = nm[i]
        X = Fr_UL*Frp_Fr*(Tref-Tinf[i])*Dt*Acc/L[i]
        if dhw == 1:
            fc = (11.6 + 1.18*Tw+3.86*Tm-2.32*Tinf[i])/(100-Tinf[i])
            X = X*fc
        Y1 = Fr_tan*Frp_Fr*ta_tan1[i]*Hmt[i]*N*Acc/L[i]
        Y2 = Fr_tan*Frp_Fr*ta_tan2[i]*Hmt[i]*N*Acc/L[i]
        Y3 = Fr_tan*Frp_Fr*ta_tan3[i]*Hmt[i]*N*Acc/L[i]
        f1[i] = fchart(X,Y1)
        f2[i] = fchart(X,Y2)
        f3[i] = fchart(X,Y3)
    Charge = sum(L)
    Fra1[ip] = sum(f1*L)/Charge
    Fra2[ip] = sum(f2*L)/Charge
    Fra3[ip] = sum(f3*L)/Charge
plot(Accv,Fra1,Accv,Fra2,Accv,Fra3)
legend(('Eq 8.15b','Eq 8.14 ','regle pouce'))
show()
