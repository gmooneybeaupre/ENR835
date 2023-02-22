#coding: utf-8
#
# exemple solaire 3.2
#
# Comparaison entre les valeurs de radiation directe et diffuse venant d'un fichier meteo
# et celles évaluées par les corrélations de Erbs et de Corrales
#
from solar_mod import *
import numpy as np
from matplotlib.pyplot import *

j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])  # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
Gsc = 1367.0
#
# donnees du probleme
#
phi = 45.5  # latitude  Montreal
lon = -73.0 -35.0/60.0  # longitude  Montreal (73 deg 35 min ouest)
Lst = -75        # longitude  méridien HNE -5
beta = 45
gam = 0.0
rhog = 0.3


#
# importer le fichier meteo
#
donnees = np.loadtxt("Meteo_horaire_Montreal.txt",skiprows = 1)   # donnes climatiques horaire Montreal  % valeurs en kJ/hr/m2
t = donnees[:,0]   # heures
I =  donnees[:,3]/1000   # valeurs de I total en MJoules/m2/ hr
Ib =  donnees[:,4]/1000   # valeurs de I direct  en MJoules/m2/ hr
Id = donnees[:,5]/1000   # valeurs de I diffus en MJoules/m2/ hr
It_trnsys = donnees[:,6]/1000   # valeurs de I diffus en MJoules/m2/ hr
#
# debut du programme
#
n =jour_mois_jour_annee(5,'mai')
print ('n= ',n )
#
# Question a
#
n_hr = (n-1)*24+12  # nombre correspondant aux nombres d'heures depuis le début de l'année
Ia = I[n_hr]
Iba = Ib[n_hr]
Ida = Id[n_hr]
Ita = It_trnsys[n_hr]
print ('Radiation le 5 mai de midi à 13:00 hrs (fichier meteo)')
print ('I total  = ', '%.2f'% Ia,' MJ/m2 hr')
print ('I direct  = ', '%.2f'% Iba,' MJ/m2 hr')
print ('I diffus  = ', '%.2f'% Ida,' MJ/m2 hr')

#
# Calcul entre midi et 13:00 hrs
#
ome1 = 0.0
ome2 = 15.0
omen = (ome1+ome2)/2.0
# calcul de thetaz
# calcul de theta
Rb = Calcul_Rb(phi,n,omen,beta,gam)
Io = irradiation_extraterrestre_horaire(n,phi,ome1,ome2)
IoM = Io/1e6 # transforme en MJ
Itb1 = Iba*Rb
Itd1 = Ida*(1+cosd(beta))/2.0
Itr1 = Ia*rhog*(1-cosd(beta))/2.0
It1 =  Itb1+Itd1+Itr1
It1b,Itb1,Itd1,Irt1 = modele_isotropique(Ia,Iba,Ida,beta,Rb,rhog)
It2b,Itb2,Itd2,Itr2 = modele_hay_davis(Ia,Iba,Ida,beta,Rb,rhog,IoM)
It3b,Itb3,Itd3,Itr3 = modele_hdkr(Ia,Iba,Ida,beta,Rb,rhog,IoM)
print ('modele isotropique',  '%.2f'% It1,' MJ/hr-m2')
print ('modele Hay davis',  '%.2f'% It2b,' MJ/hr-m2')
print ('modele HDKR',  '%.2f'% It3b,' MJ/hr-m2')
# model Perez plus complexe
delt  = decl_solaire(n)
thez = zenith_solaire(phi,delt,omen)
Gon = irradiation_extraterrestre_normale(n)
Ion = Gon*3600.0
IonM = Ion/1e6
the = normale_solaire(delt,phi,omen,beta,gam)
It4,Itb4,Itd4,Itr4 = modele_perez(Ia,Iba,Ida,beta,Rb,rhog,IoM,IonM,thez,the)
print ('modele Perez',  '%.2f'% It4,' MJ/hr-m2')
print ('TRNSYS',  '%.2f'% Ita,' MJ/hr-m2' )

#
#
# Calcul des valeurs journalieres et mensuelles à partir des valeurs horaires
#
# Initialisations des vecteurs
# valeurs journalieres

#
# Initialisations des vecteurs
# valeurs journalieres

It_iso = np.zeros(8760)
It_HD = np.zeros(8760)
It_HDKR = np.zeros(8760)
It_Perez = np.zeros(8760)
Ht_iso = np.zeros(365)      # totales
Ht_HD = np.zeros(365)      # totales
Ht_HDKR = np.zeros(365)      # totales
Ht_Perez = np.zeros(365)      # totales
Ht_trnsys = np.zeros(365)      # totales
# valeurs mensuelles
Hm = np.zeros(12)
Hmb = np.zeros(12)
Hmd = np.zeros(12)
Htm_iso = np.zeros(12)      # totales
Htm_HD = np.zeros(12)      # totales
Htm_HDKR = np.zeros(12)      # totales
Htm_Perez = np.zeros(12)      # totales
Htm_trnsys = np.zeros(12)      # totales
#
kj = 0  # index des jours
kh = 0  # index des heures
for im in range(0,12):     # boucle allant de 0 à 11
    # boucle sur les 12 mois im = 0 pour janvier, etc..
    somm1 = 0
    somm2 = 0
    somm3 = 0
    somm4 = 0
    somm5 = 0
    somm6 = 0
    somm7 = 0
    somm8 = 0
    for ij  in range(0,nm[im]):
        # boucle sur tous les jours du mois im, ij va de 0 à 27,29 ou 30
        som1 = 0
        som2 = 0
        som3 = 0
        som4 = 0
        som5 = 0
        som6 = 0
        som7 = 0
        som8 = 0
        delt  = decl_solaire(kj+1)
        for ih in range(0,24):
            som1 = som1 + I[kh]   # Joules/m2 /hr
            som2 = som2 + Ib[kh]  # Joules/m2 /hr
            som3 = som3 + Id[kh]  # Joules/m2 /hr
            # boucle sur les 24 heures du jour
            ome1 = -180.0+ih*15.0
            ome2 = ome1+15.0
            omen = (ome1+ome2)/2.0
            nj = kj+1
            Rb = Calcul_Rb(phi,nj,omen,beta,gam)
            Io = irradiation_extraterrestre_horaire(nj,phi,ome1,ome2)
            Iom = Io/1e6
            It1,Itb1,Itd1,Itr1 = modele_isotropique(I[kh],Ib[kh],Id[kh],beta,Rb,rhog)
            It2,Itb2,Itd2,Itr2 = modele_hay_davis(I[kh],Ib[kh],Id[kh],beta,Rb,rhog,Iom)
            It3,Itb3,Itd3,Itr3 = modele_hdkr(I[kh],Ib[kh],Id[kh],beta,Rb,rhog,Iom)
            # model Perez plus complexe
            thez = zenith_solaire(phi,delt,omen)
            Gon = irradiation_extraterrestre_normale(n)
            Ion = Gon*3600.0
            Ionm = Ion/1e6
            the = normale_solaire(delt,phi,omen,beta,gam)
            It4,Itb4,Itd4,Itr4 = modele_perez(I[kh],Ib[kh],Id[kh],beta,Rb,rhog,Iom,Ionm,thez,the)
            It_iso[kh] = It1
            It_HD[kh] = It2
            It_HDKR[kh] = It3
            It_Perez[kh] = It4
            som4 = som4 + It1
            som5 = som5 + It2
            som6 = som6 + It3
            som7 = som7 + It4
            som8 = som8 + It_trnsys[kh]
            kh = kh + 1
        Ht_iso[kj] = som4      # Joules/m2/jr
        Ht_HD[kj] = som5      # Joules/m2/jr
        Ht_HDKR[kj] = som6      # Joules/m2/jr
        Ht_Perez[kj] = som7     # Joules/m2/jr
        Ht_trnsys[kj] = som8     # Joules/m2/jr
        somm1 = somm1 + som1
        somm2 = somm2 + som2
        somm3 = somm3 + som3
        somm4 = somm4 + som4
        somm5 = somm5 + som5
        somm6 = somm6 + som6
        somm7 = somm7 + som7
        somm8 = somm8 + som8
        kj = kj + 1
    Hm[im] = somm1/nm[im]      # Joules/m2/jr
    Hmb[im] = somm2/nm[im]       # Joules/m2/jr
    Hmd[im] = somm3/nm[im]       # Joules/m2/jr
    Htm_iso[im] = somm4/nm[im]       # Joules/m2/jr
    Htm_HD[im] = somm5/nm[im]       # Joules/m2/jr
    Htm_HDKR[im] = somm6/nm[im]       # Joules/m2/jr
    Htm_Perez[im] = somm7/nm[im]      # Joules/m2/jr
    Htm_trnsys[im] = somm8/nm[im]      # Joules/m2/jr
#
# exemple de l'irradiation sur plan incliné
#
ij = n - 1
print ('Radiation plan incliné le 5 mai modele isotrope  ')
print ('HT total  = ', '%.2f'% Ht_iso[ij],' MJ/jr-m2' )
im = 4
print ('Radiation plan incliné pour le mois de mai')
print ('HTm total  = ', '%.2f'% Htm_iso[im],' MJ/jr-m2')
Rbb = Calcul_Rb_mois(phi,j_type[im],beta,gam)
Hom = irradiation_extraterrestre_jour_moyen(im+1,phi)
Hom = Hom/1e6
print ('Radiation totale plan horizontal pour le mois de mai')
print ('HTm total  = ', '%.2f'% Hm[im],' MJ/jr-m2')
print ('Radiation directe plan horizontal pour le mois de mai')
print ('HTm total  = ', '%.2f'% Hmb[im],' MJ/jr-m2')
print ('Radiation diffuse plan horizontal pour le mois de mai')
print ('HTm total  = ', '%.2f'% Hmd[im],' MJ/jr-m2')

KTm = Hm[im]/Hom
HTm2 = Hmb[im]*Rbb+Hmd[im]*(1+cosd(beta))/2.0+Hm[im]*rhog*(1.0-cosd(beta))/2.0
print ('En mai,   la radiation totale incliné (approchée) est = ',  '%.2f'% HTm2,' MJ/jr-m2')
#
# exemple de l'irradiation sur plan incliné
#
err_iso = abs(Ht_iso - Ht_trnsys)
err_HD = abs(Ht_HD - Ht_trnsys)
err_HDKR = abs(Ht_HDKR - Ht_trnsys)
err_Perez = abs(Ht_Perez - Ht_trnsys)
ww = sum(Ht_trnsys)
print('Difference moldele ISO = ','%.1f'% (sum(err_iso)/ww*100), ' %')
print('Difference moldele HD = ','%.1f'% (sum(err_HD)/ww*100), ' %')
print('Difference moldele HDKR = ','%.1f'% (sum(err_HDKR)/ww*100), ' %')
print('Difference moldele Perez = ','%.1f'%  (sum(err_Perez)/ww*100), ' %')
figure(1)
plot(t,It_HDKR,t,It_trnsys,'x')
legend(('HDKR','TRNSYS'))
title('Valeurs horaire')
figure(2)
imj = np.arange(0,365)
plot(imj,Ht_iso,imj,Ht_HD,imj,Ht_HDKR,imj,Ht_Perez,imj,Ht_trnsys,'x')
legend(('iso','HD','HDKR','Perez','TRNSYS'))
title('Valeurs journalières')
figure(3)
imm = np.arange(0,12)
plot(imm,Htm_iso,imm,Htm_HD,imm,Htm_HDKR,imm,Htm_Perez,imm,Htm_trnsys,'x')
legend(('iso','HD','HDKR','Perez','TRNSYS'))
title('Valeurs mensuelles')
show()

