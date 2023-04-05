#coding: utf-8
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
Tinfv = np.array([-10.3,-8.8,-2.4,5.7,12.9,18.0,20.8,19.4,14.5,8.3,1.6,-6.9])
Tground = np.array([6,4.5,4,5,7.5,10,13.5,15.5,16,14.5,12,10])
Lv = np.array([15,13,9,6,3,2,2,2,3,4,8,12])
rhogv = np.array([0.7, 0.7,0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4 ])
HmKWhv = np.array([1.5817,2.5884,3.8525,4.5808,5.6846,6.2114,6.0228,4.8676,4.0267,2.5372,1.2235,1.1742])  # Radiation plan horizontal kWh /jr/m2
Hmv = HmKWhv*3.6     # Radiation plan horizontal MJoukles /jr/m2
ta_talv = np.array([0.96,0.96,0.96,0.92,0.92,0.92,0.92, 0.92, 0.92, 0.92, 0.96, 0.96 ])
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
DeltaT = 2
UAsol = 30

#

#
nannees = 1
nn = 12*nannees
KTv = np.zeros(nn)
phiv = np.zeros(nn)
Xcv = np.zeros(nn)
KTv = np.zeros(nn)
Euv = np.zeros(nn)
Pertesv = np.zeros(nn)
Tfinalv = np.zeros(nn)
Tiniv = np.zeros(nn)
mois  = []
Energie_annuelle = np.zeros(nannees)
#
#
compt_max = 20
delta = 0.01
Tini = 24
na = -1
#
for ix in range(0,nn):
    im = np.mod(2 + ix,12)   # on débute la simulation au mois de mars selon l'énoncé
    ij = int(ix/12)
    imm = ij*12 + im
    imn = np.mod(ix,12)
    if imn == 0:
        som1 = 0
        na = na + 1
    rhog = rhogv[im]
    Tinf = Tinfv[im] + DeltaT
    Hm = Hmv[im]    #
    HT = Calcul_HT(im,Hm,rhog,beta,gam)
    Rm = HT/Hm
    ta_tan = ta_talv[im]
    Frp_ta = Frp_tan*ta_tan
    njt = j_type[im]
    Hom = irradiation_extraterrestre_jour(njt,phi)  # J/n2-jr
    Hom = Hom/1e6    # on met en MJ/m2-jr
    KTm = Hm/Hom
    KTv[imm] = KTm
    Lm = Lv[im]
    Tg = Tground[im]
    mois.append(nom[imn])


    #
    # debut du processus iteratif
    #
    ok = False
    compt = 0
    Tfin = 26.0 # hypothese sur la valeur finale du mois
    while not ok:
        Tmoy = (Tfin + Tini)/2
        Pertes = UAsol*(Tmoy - Tg)*3600*24*nm[im]/1e9  # GJ
        ITc = Frp_UL*(Tmoy - Tinf)/Frp_ta # W/m2
        ITcMJ = ITc*3600/1e6    # MJ/m2
        phib,Xc = Calcul_phi_bar(njt,phi,beta,gam,Hm,ITcMJ,rhog,KTm)
        Eutile = nm[im]*phib*HT*Frp_ta*Ac/1000  # en GJ
        Tfin_n = Tini + (Eutile - Lm -  Pertes)/CC
        diff = abs(Tfin_n - Tfin)
        if diff > delta:
            Tfin = Tfin_n
            compt = compt + 1
            if compt > compt_max:
                ok = True
                print('Erreur non convergence')
        else:
            ok = True
            Tfinalv[imm] = Tfin
            Tiniv[imm] = Tini
            Euv[imm] = Eutile
            Pertesv[imm] = Pertes
            Tini = Tfin  # ceci sera la température initiale du prochain mois
            phiv[imm] = phib
            Xcv[imm] = Xc
    som1 = som1 + Eutile
    if imn == 11:
        Energie_annuelle[na] = som1


head = ['mois','Energie \nrécupérée (GJ)','Pertes\nsol(GJ)','Température \ndébut mois','Température \nfin mois','Xc','phi']
my_list = []
for ix in range(0,nn):
        imm = np.mod(2 + ix,12)   # on débute la simulation au mois de mars selon l'énoncé
        imn = np.mod(ix,12)
        ij = int(ix/12)
        im = ij*12 + imm
        pt = [mois[im],'%.2f' % Euv[im],'%.2f' % Pertesv[im],'%.2f' % Tiniv[im],'%.2f' % Tfinalv[im],'%.3f' % Xcv[im],'%.3f' % phiv[im]]
        my_list.append(pt)
print(tabulate(my_list,headers= head,numalign="left"))
for ia in range(0,nannees):
    print('Énergie totale récupérée pour l''année', '%d' % (ia +1) , ' = ','%.2f' % Energie_annuelle[ia],' GJ')
print('Énergie totale récupérée  = ','%.2f' % sum(Euv),' GJ')