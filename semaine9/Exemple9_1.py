#coding: utf-8
#
#
#
from solar_mod import *
import numpy as np
from tabulate import tabulate


j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])  # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures Ã©coulÃ©es aprÃ¨s chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
Gsc = 1367.0
rhogv = np.array([0.7, 0.7,0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4 ])
Tinfv = np.array([-10.3,-8.8,-2.4,5.7,12.9,18.0,20.8,19.4,14.5,8.3,1.6,-6.9])
ta_talv = np.array([0.96,0.96,0.96,0.92,0.92,0.92,0.92, 0.92, 0.92, 0.92, 0.96, 0.96 ])
HmKWhv = np.array([1.5817,2.5884,3.8525,4.5808,5.6846,6.2114,6.0228,4.8676,4.0267,2.5372,1.2235,1.1742])  # Radiation plan horizontal kWh /jr/m2
Hmv = HmKWhv*3.6     # Radiation plan horizontal MJoukles /jr/m2

#
#
def Calcul_HT(imois,Hm,rhog,bet,gamma):
    njt = j_type[imois]                 # njt est le jour type du mois
    Hom = irradiation_extraterrestre_jour(njt,phi)
    HomMJ = Hom/1e6
    KTm = Hm/HomMJ
    delt = decl_solaire(njt)
    omesm = angle_sunset(phi,delt)
    Hdtm = Erbs_mois(KTm,omesm)
    Hdm = Hm*Hdtm              # Contribution diffus donnÃ©e par Erbs
    Hbm = Hm - Hdm        # Contribution directe donnÃ©e par Erbs
    Rbb = Calcul_Rb_mois(phi,njt,bet,gamma)
    HT = Rbb*Hbm + Hdm*(1+cosd(bet))/2 + Hm*rhog*(1 - cosd(bet))/2
    return HT


phi = 45.5
Ac = 50.0
gam = 0
beta = 60.0
Vol = 250.0    # m3
rhol = 1000.0
Cp = 4190.0
CC = rhol*Vol*Cp/1e9        # GJ/K
Frp_tan = 0.78
Frp_UL = 4.55
ta_tan = 0.96
DeltaT = 2

#
# calcul iteratif pour le mois de Mars
#
#
Tini = 25
im = 2
#
rhog = rhogv[im]
Tinf = Tinfv[im] + DeltaT       # +2 tel que suggéré
Hm = Hmv[im]  # Hm est en MJ/m2 jr
HmkWh = HmKWhv[im]  # Hm est en kWh/m2 jr
Frp_ta = Frp_tan*ta_tan
#
Tmoyh = 25
ITcW = Frp_UL*(Tmoyh - Tinf)/Frp_ta # W/m2
ITcMJ = ITcW*3.6/1e3
njt = j_type[im]
Hom = irradiation_extraterrestre_jour(njt,phi)
Hom = Hom/1e6  # en MJ/m2-jr
KTm = Hm/Hom
# Utilisation de la fonction Calcul_phi_bar
phib,Xc = Calcul_phi_bar(njt,phi,beta,gam,Hm,ITcMJ,rhog,KTm)
print('Coefficient Xc, mois de mars = ',Xc)
print('Utilisatbilité mois mars = ',phib)

#
# calcul  au long
#
delt = decl_solaire(njt)  # declinaison solaire
omesm = angle_sunset(phi,delt) # angle du coucher du soleil le 16 mars
HT = Calcul_HT(im,Hm,rhog,beta,gam)   # calcul de HT à partiur des corrélations de Erbs
Rbar = HT/Hm
A = 2.943 - 9.271*KTm + 4.031*KTm**2
B = -4.345 + 8.853*KTm - 3.602*KTm**2
C = -0.170 - 0.306*KTm + 2.936*KTm**2
ome = 0 # midi
rtn = Collares_total(ome,omesm)
rdn = Collares_diffus(ome,omesm)
Hdt = Erbs_jour(KTm,omesm)
the = normale_solaire(delt,phi,ome,beta,gam)
thez = zenith_solaire(phi,delt,ome)
Rbn = cosd(the)/cosd(thez)
Rn = (1 - rdn*Hdt/rtn)*Rbn + rdn*Hdt/rtn*(1+cosd(beta))/2 + rhog*(1 - cosd(beta))/2
Xc1 = ITcW/(rtn*Rn*HmkWh*1000)
Xc2 = ITcMJ/(rtn*Rn*Hm)
phib_mars = np.exp((A + B*Rn/Rbar)*(Xc1 + C*Xc1**2))
print('Coefficient Xc, mois de mars = ',Xc1,Xc2)
print('Utilisatbilité mois mars = ',phib_mars)
Eutile_mars = nm[im]*phib*HT*Frp_ta*Ac/1000  # en GJ
print('Énergie récupérée mois mars = ',Eutile_mars,' GJ')
#
# Calcul pour l'année
#

Euv = np.zeros(12)
Xcv = np.zeros(12)
phiv = np.zeros(12)
for im in range(0,12):
    ta_tan = ta_talv[im]
    rhog = rhogv[im]
    Tinf = Tinfv[im] + DeltaT       # +2 tel que suggéré
    Hm = Hmv[im]  # Hm est en MJ/m2 jr
    HmkWh = HmKWhv[im]  # Hm est en kWh/m2 jr
    Frp_ta = Frp_tan*ta_tan
    ITcW = Frp_UL*(Tmoyh - Tinf)/Frp_ta # W/m2
    ITcMJ = ITcW*3.6/1e3
    njt = j_type[im]
    Hom = irradiation_extraterrestre_jour(njt,phi)
    Hom = Hom/1e6  # en MJ/m2-jr
    KTm = Hm/Hom
    # Utilisation de la fonction Calcul_phi_bar
    phib,Xc = Calcul_phi_bar(njt,phi,beta,gam,Hm,ITcMJ,rhog,KTm)
    delt = decl_solaire(njt)  # declinaison solaire
    omesm = angle_sunset(phi,delt) # angle du coucher du soleil le 16 mars
    HT = Calcul_HT(im,Hm,rhog,beta,gam)   # calcul de HT à partiur des corrélations de Erbs
    Xcv[im] = Xc
    Euv[im] = nm[im]*phib*HT*Frp_ta*Ac/1000  # en GJ
    phiv[im] = phib
head = ['mois','Energie \nrécupérée (GJ)','Xc','phi']
my_list = []
for im in range(0,12):
        pt = [nom[im],'%.2f' % Euv[im],'%.2f' % Xcv[im],'%.2f' % phiv[im]]
        my_list.append(pt)
print(tabulate(my_list,headers= head,numalign="left"))
print('Énergie totale récupérée = ','%.2f' % sum(Euv),' GJ')