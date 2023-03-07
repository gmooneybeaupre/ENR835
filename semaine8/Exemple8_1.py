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
Hmj = 5.694       # MJ/m2-jr
Hmbj = 2.977       # MJ/m2-jr
Hmdj = Hmj - Hmbj
phi = 45.5
beta = 45.0
gam = 0.0
rhog = 0.2
a1 = 4.18           # parametre du capteur
a2 = 0.02           # parametre du capteur
etao = 0.674        # parametre du capteur
bo = 0.21           # parametre du capteur
b1 = 0.0
thed = angle_diffus(beta)
ther = angle_reflechi(beta)
def  kta(the):
    if the > 60:
        kta = 1 - bo - b1  # kta à 60 degres
        kta = kta*(1 - (the - 60.0)/30.0)
    else:
        kta = 1 - bo*(1/cosd(the)-1) - b1*(1/cosd(the)-1)**2
    return kta


def Calcul_Kta(phi,beta,gam,rhog,njt,Hm,Hmd):


    # Eq 8.14 connaisssant Hm , Hmd

    delt = decl_solaire(njt)
    omesm = angle_sunset(phi,delt)
    om1 = -omesm - np.mod(-omesm,-15)
    om2 = omesm - np.mod(omesm,15)
    nn = int((om2-om1)/15.0)   # nombre d'heures d'ensoleillement
    s1 = 0
    s2 = 0
    thed = angle_diffus(beta)
    Kd = kta(thed)
    ther = angle_reflechi(beta)
    Kr = kta(ther)
    for ih in range(0,nn):
        ome = om1 + 15*ih + 7.5
        Rb = Calcul_Rb(phi,njt,ome,beta,gam)
        rt = Collares_total(ome,omesm)
        rd = Collares_diffus(ome,omesm)
        I = rt*Hm
        Id =  rd*Hmd
        Ib = I - Id
        theb = normale_solaire(delt,phi,ome,beta,gam)
        Kb = kta(theb)
        IT = Ib*Rb + Id*(1+cosd(beta))/2 + I*rhog*(1-cosd(beta))/2
        S = Kb*Ib*Rb + Kd*Id*(1+cosd(beta))/2 + Kr*I*rhog*(1-cosd(beta))/2
        s1 = s1 + S
        s2 = s2 + IT
    ta =  s1/s2
    return ta


Kd = kta(thed)
Kr = kta(ther)

njt = j_type[0]   # njt est le jour type du mois de janvier
delt = decl_solaire(njt)
omesm = angle_sunset(phi,delt)
omeb = 37.5                          # angle horaire correspondant à 2:30 du midi
theb = normale_solaire(delt,phi,omeb,beta,gam)
Kb = kta(theb)
Rbb = Calcul_Rb_mois(phi,njt,beta,gam)
HmtA,HmtbA,HmtdA,HmtrA = modele_isotropique(Hmj,Hmbj,Hmdj,beta,Rbb,rhog)
ta_tan_janv1  = (Kb*HmtbA +  Kd*HmtdA + Kr*HmtrA)/HmtA
Hom = irradiation_extraterrestre_jour(njt,phi)      # J/jr-m2
HomMJ = Hom/1e6
KTm = Hmj/HomMJ
#
#  Fonction du module solaire Calcul_Kta_moyen où on suppose qu'uniquement  la radiation totale
#       sur plan horizontal est connu, la fonction utilise les correlation de Erbs pour estimer les autres
ta_tan_janv2 = Calcul_Kta_moyen(phi,beta,gam,rhog,bo,b1,njt,Hmj,KTm)
#
#  Fonction Calcul_Kta (voir plus haut) où on suppose que  la radiation totale, diffus et direct
#       sur plan horizontal sopnt connues
ta_tan_janv3=  Calcul_Kta(phi,beta,gam,rhog,njt,Hmj,Hmdj)
print('rapport ta / tan janvier est de (Eq 8.15b) :',ta_tan_janv1)
print('rapport ta / tan janvier est de (Eq 8.14) :',ta_tan_janv2,ta_tan_janv3)

#
#  Analyse à partir des valeurs  horaires
# importer le fichier meteo
#
donnees =  np.loadtxt("Meteo_horaire_Montreal.txt",skiprows = 1)   # donnes climatiques horaire Montreal  % valeurs en kJ/hr/m2
hr = donnees[:,0]   # heures
Tamb = donnees[:,1]   # valeurs Tinf
hr = donnees[:,0]   # heures
I =  1000*donnees[:,3]   # valeurs de I total donné par TRNSYS en Joules/m2/ hr
Ib =  1000*donnees[:,4]   # valeurs de I direct donné par TRNSYS en Joules/m2/ hr
Id  = 1000*donnees[:,5]   # valeurs de I diffus donné par TRNSYS en Joules/m2/ hr
It = np.zeros(8760)         # radiation horaire plan incliné
Ktav = np.zeros(8760)
Hm = np.zeros(12)      # radiation horaire
Hmb = np.zeros(12)
Hmd = np.zeros(12)
Hmt = np.zeros(12)
Hmtb = np.zeros(12)
Hmtd = np.zeros(12)
Hmtr = np.zeros(12)
ta_tan1 =  np.zeros(12)
ta_tan2 =  np.zeros(12)
ta_tan3 =  np.zeros(12)
ta_tan4 =  np.zeros(12)
ta_tan5 =  np.zeros(12)
i_jour = 0      # index des jours (0 à 365)
i_hr = 0        # index des heures (0 à 8760)


for im in range(0,12):     # im index des moix
    som1 = 0
    som2 = 0
    som3 = 0
    som4 = 0
    som5 = 0
    som6 = 0
    som7 = 0
    som8  =0
    for ij in range(0,nm[im]):        # ij : jour du mois
        for j in range(0,24):
            ome1 = -180.0 + 15.0*j
            ome2 = ome1 + 15.0
            omem = (ome1+ome2)/2.0
            Rb = Calcul_Rb(phi,i_jour+1,omem,beta,gam)
            delt  = decl_solaire(i_jour+1)
            theb = normale_solaire(delt,phi,omem,beta,gam)
            thez = zenith_solaire(delt,phi,omem)
            Itt,Itb,Itd,Itr = modele_isotropique(I[i_hr],Ib[i_hr],Id[i_hr],beta,Rb,rhog)
            It[i_hr] = Itt
            if (Itt > 0):
                Kb = kta(theb)
                Ka = (Kb*Itb+Kd*Itd+Kr*Itr)/Itt
            else:
                Kb = 0
                Ka = 0
            Ktav[i_hr] = Ka
            som1 = som1 + Kb*Itb
            som2 = som2 + It[i_hr]
            som3 = som3 + I[i_hr]
            som4 = som4 + Ib[i_hr]
            som5 = som5 + Id[i_hr]
            som6 = som6 + Itb
            som7 = som7 + Itd
            som8 = som8 + Itr
            i_hr = i_hr+1
        i_jour = i_jour+1
    Hm[im] = som3/nm[im]    # Joules/m2/jr
    Hmt[im] = som2/nm[im]    # Joules/m2/jr
    Hmb[im] = som4/nm[im]    # Joules/m2/jr
    Hmd[im] = som5/nm[im]    # Joules/m2/jr
    Hmtb[im] = som6/nm[im]    # Joules/m2/jr
    Hmtd[im] = som7/nm[im]    # Joules/m2/jr
    Hmtr[im] = som8/nm[im]    # Joules/m2/jr
    ta_tan1[im]  = (som1/nm[im] +  Kd*Hmtd[im] + Kr*Hmtr[im])/Hmt[im]    # Eq 8.13
#
    # calcul  à partir des données menulelles sur un plan horizontal
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
    ta_tan2[im] = Calcul_Kta_moyen(phi,beta,gam,rhog,bo,b1,njt,Hm[im],KTm)
    HmtA,HmtbA,HmtdA,HmtrA = modele_isotropique(Hm[im],Hmb[im],Hmd[im],beta,Rbb,rhog)
    ta_tan3[im]  = (Kb*HmtbA +  Kd*HmtdA + Kr**HmtrA)/HmtA    # Eq 8.15 - 8.16
    HmtE,HmtbE,HmtdE,HmtrE = modele_isotropique(Hm[im],HmbE,HmdE,beta,Rbb,rhog)
    ta_tan4[im]  = (Kb*HmtbE +  Kd*HmtdE + Kr**HmtrE)/HmtE    # Eq 8.15 - 8.16
    if ((im < 9) & (im > 3)):
        ta_tan5[im] = 0.92  # été
    else:
        ta_tan5[im] = 0.96   # hiver
imm = np.arange(0,12)
plot(imm,ta_tan1,imm,ta_tan2,'-*',imm,ta_tan3,'x',imm,ta_tan4,imm,ta_tan5)
legend(('reference Eq 8.13','Eq 8.14','Eq 8.15b','Eq 8.15b (Erbs)','Eq 8.16'))
show()