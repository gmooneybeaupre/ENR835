#coding: utf-8
#
#
# Exemple de la méthode f-chart air chaud
#
#
from time import time
from solar_mod import *
from numpy import *
from scipy.optimize import newton
#
# donnees due l'Exemple 8.2
j_type = array([17,47,75,105,135,162,198,228,258,288,318,344])      # jour type de chaque mois
nm = array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
Tinfv = np.array([-10.3,-8.8,-2.4,5.7,12.9,18.0,20.8,19.4,14.5,8.3,1.6,-6.9])
HmKWh = np.array([1.5817,2.5884,3.8525,4.5808,5.6846,6.2114,6.0228,4.8676,4.0267,2.5372,1.2235,1.1742])  # Radiation plan horizontal kWh /jr/m2
nom = ['janvier','février','mars','avril','mai','juin','juillet','aout','septembre','octobre','novembre','decembre']
Hm = HmKWh*3.6*1e6     # Radiation plan horizontal MJ /jr/m2


#
# paramètres du capteur lors de l'essai
#

phi = 45.5
Fr_tan = 0.76
Fr_UL = 4.85
UA = 300.0      # coefficient de pertes
DJv = np.array([874,747])
Tref = 100
Volume = 0.15  #m3/m2
Debit = 12.0  # l/s / m2
beta = 60
Ac = 40
gam = 0
rhog = 0.7
ta_tan = 0.88

#
nn = 0
txx = ['a)','b)']
fv = np.zeros(2)
Lv = np.zeros(2)
for im  in range(0,2):
#
    print(txx[im],'\n')
    Tinf = Tinfv[im]
    njt = j_type[im]                 # njt est le jour type du mois
    Hom = irradiation_extraterrestre_jour(njt,phi)
    KTm = Hm[im]/Hom
    delt = decl_solaire(njt)
    omesm = angle_sunset(phi,delt)
    Hdtm = Erbs_mois(KTm,omesm)
    print('\tDéclinaison solaire pour le ',njt -nn,nom[im],' est :','%0.2f' % delt)
    print('\tangle du lever-coucher solaire  pour le ',njt -nn,nom[im],' est :','%0.2f' % omesm)
    print('\tindice de clarté mensuel pour le mois de ',nom[im],' est :','%0.2f' % KTm)
    print('\tLe rapport diffus-total donnée par les relations de Erbs pour le mois de ',nom[im],' est :','%0.2f' % Hdtm)
    Hdm = Hm[im]*Hdtm              # Contribution diffus donnée par Erbs
    Hbm = Hm[im] - Hdm        # Contribution directe donnée par Erbs
    Rbb = Calcul_Rb_mois(phi,njt,beta,gam)
    print('\tLa radiation diffuse estimée  pour le mois de ',nom[im],' est :','%0.2f' % (Hdm/1e6),' MJ/jr-m2')
    print('\tLa radiation directe estimée  pour le mois de ',nom[im],' est :','%0.2f' % (Hbm/1e6),' MJ/jr-m2')
    print('\tRb moyen mensuel pour le mois de ',nom[im],' est :','%0.2f' % Rbb)
    #
    # Calcul avec les valeurs du fichier meteo
    Ht,Htb,Htd,Htr = modele_isotropique(Hm[im],Hbm,Hdm,beta,Rbb,rhog)
    print('\tLa radiation totale sur plan incliné  pour le mois de ',nom[im],' est :','%0.2f' % (Ht/1e6),' MJ/jr-m2')
    Dt = nm[im]*24*3600
    N = nm[im]
    DJ = DJv[im]
    Lm = UA*DJ*24.0*3600 # charge mensuelle chauffage en Joules
    Lv[im] = Lm
    print('\tLa charge de chauffage totale pour le mois de ',nom[im],' est :','%0.2f' % (Lm/1e9),' GJ')
    Fr_ta = Fr_tan*ta_tan
    Xn = Fr_UL*(Tref-Tinf)*Dt*Ac/Lm
    print('\tLe coefficient X non corrigé pour le mois de ',nom[im],' est :','%0.2f' % Xn)
    Xc1 = (Debit/10)**(0.28)
    Xc2 = (Volume/0.25)**(-0.3)
    X = Xn*Xc1*Xc2
    print('\tLe coefficient X corrigé pour le mois de ',nom[im],' est :','%0.2f' % X)
    Y = Fr_tan*ta_tan*Ht*N*Ac/Lm
    print('\tLe coefficient Y  pour le mois de ',nom[im],' est :','%0.2f' % Y)
    f = fchart_air(X,Y)
    print('\tLa fraction solaire prédite par al méthode f-chart pour le mois de ',nom[im],' est :','%0.2f' % f)
    nn = nn + nm[im]
    fv[im] = f
F = sum(Lv*fv)/sum(Lv)
print('\n')
print('La fraction solaire moyenne pour les 2 mois est de :','%0.2f' % F)

