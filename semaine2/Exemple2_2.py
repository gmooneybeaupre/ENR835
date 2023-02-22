#coding: utf-8
#
# exemple solaire 2.2
#
# Comparaison entre les valeurs de radiation directe et diffuse venant d'un fichier meteo
# et celles évaluées par les corrélations de Erbs et de Corrales
#
import numpy as np
from solar_mod import *
pi = np.pi

j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])  # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
Gsc = 1367.0
phi = 45.5  # latitude  Montreal
lon = -73.0 -35.0/60.0  # longitude  Montreal (73 deg 35 min ouest)
Lst = -75        # longitude  méridien HNE -5
#
# importer le fichier meteo
#
donnees = np.loadtxt("Meteo_horaire_Montreal.txt",skiprows = 1)   # donnes climatiques horaire Montreal  % valeurs en kJ/hr/m2
temps = donnees[:,0]   # heures
I =  donnees[:,3]/1000   # valeurs de I total en MJoules/m2/ hr
Ib =  donnees[:,4]/1000   # valeurs de I direct  en MJoules/m2/ hr
Id = donnees[:,5]/1000   # valeurs de I diffus en MJoules/m2/ hr

#
# debut du programme
#
n =jour_mois_jour_annee(5,'mai')
print ('n= ',n)
#
# Question a
#
n_hr = (n-1)*24+12  # nombre correspondant aux nombres d'heures depuis le début de l'année
Ia = I[n_hr]
Iba = Ib[n_hr]
Ida = Id[n_hr]
print ('a) Radiation le 5 mai de midi à 13:00 hrs (fichier meteo)')
print ('\t I total (fichier meteo)  = ','%.2f' %Ia,' MJ/hr-m2')
print ('\t I direct (fichier meteo) = ','%.2f' %Iba,' MJ/hr-m2')
print ('\t I diffus (fichier meteo) = ','%.2f' % Ida,' MJ/hr-m2')
#
# b correlation de Erbs
#
Io = irradiation_extraterrestre_horaire(n,phi,0,15)
Io = Io/1e6      # en MJ/hr-m2
kT = Ia/Io           # indice de clarté horaire

Idt = Erbs_horaire(kT)
IdE = Idt*Ia
IbE = Ia - IdE
print ('b) Radiation le 5 mai de midi à 13:00 hrs (Erbs)')
print ('\t I direct (Erbs) = ','%.2f' %IbE,' MJ/hr-m2')
print ('\t I diffus (Erbs) = ','%.2f' %IdE,' MJ/hr-m2')
#
#
# Calcul des valeurs journalieres et mensuelles à partir des valeurs horaires
#
# Initialisations des vecteurs
# valeurs journalieres
H = np.zeros(365)      # totales
Hb = np.zeros(365)     # directes
Hd = np.zeros(365)     # diffus
# valeurs mensuelles
Hm = np.zeros(12)
Hmb = np.zeros(12)
Hmd = np.zeros(12)
#
kj = 0  # index des jours
kh = 0  # index des heures
for im in range(0,12):     # boucle sur les mois allant de 0 à 11
    somm1 = 0           # sommation sur les mois pour la radiation totale
    somm2 = 0           # sommation sur les mois pour la radiation directe
    somm3 = 0
    for ij  in range(0,nm[im]):
        # boucle sur tous les jours du mois im, ij va de 0 à 27,29 ou 30
        som1 = 0            # sommation sur les jours pour la radiation totale
        som2 = 0
        som3 = 0
        for ih in range(0,24):  # boucle sur les 24 heures du jour
            som1 = som1 + I[kh]   # Joules/m2 /hr
            som2 = som2 + Ib[kh]  # Joules/m2 /hr
            som3 = som3 + Id[kh]  # Joules/m2 /hr
            kh = kh + 1
        H[kj] = som1      # Joules/m2/jr
        Hb[kj] = som2     # Joules/m2/jr
        Hd[kj] = som3     # Joules/m2/jr
        somm1 = somm1 + som1
        somm2 = somm2 + som2
        somm3 = somm3 + som3
        kj = kj + 1
    Hm[im] = somm1/nm[im]       # MJoules/m2/jr
    Hmb[im] = somm2/nm[im]      # MJoules/m2/jr
    Hmd[im] = somm3/nm[im]      # MJoules/m2/jr
#
# exemple de l'irradiation extra-terrestre sur le plan horizontal pour le 5 mai
#
ij = n - 1
print ('b ii) Radiation le 5 mai (fichier meteo)')
print ('\t H total  = ','%.2f' % H[ij],' MJ/m2 jr')
print ('\t H direct (fichier meteo)  = ','%.2f' % Hb[ij],' J/m2 jr')
print ('\t H diffus (fichier meteo)  = ','%.2f' % Hd[ij],' J/m2 jr')
#
# Comparaison avec les correlations de Erbs
#
Ho = irradiation_extraterrestre_jour(n,phi)
Ho = Ho/1e6     # transforme en MJ
print ('\t Ho = ','%.2f' % Ho,' MJ/jr-m2')
KT = H[ij]/Ho
delt = decl_solaire(n)
omes = angle_sunset(phi,delt)
Hdt = Erbs_jour(KT,omes)
HdE = H[ij]*Hdt
HbE = H[ij] - HdE
print ('\t H direct (Erbs)  = ','%.2f' % HbE,' MJ/jr-m2')
print ('\t H diffus (Erbs) = ','%.2f' % HdE,' MJ/jr-m2')
#
# exemple de l'irradiation extra-terrestre sur le plan horizontal pour le mois de mai
#
mois = 5 # mois d emai
im = mois - 1
print ('b iii) Radiation le mois de mai (fichier meteo)')
print ('\t Hm total (fichier meteo)  = ','%.2f' % Hm[im],' MJ/jr-m2')
print ('\t Hm direct (fichier meteo)  = ','%.2f' % Hmb[im],' MJ/jr-m2')
print ('\t Hm diffus (fichier meteo)   = ','%.2f' % Hmd[im],' MJ/jr-m2')
#
# Comparaison avec les correlations de Erbs
#
Hom = irradiation_extraterrestre_jour_moyen(im+1,phi)
Hom = Hom/1e6     # transforme en MJ
KTm = Hm[im]/Hom
nm = j_type[im]  # jour type du mois de mai
delt = decl_solaire(nm)
omesm = angle_sunset(phi,delt)
Hdtm = Erbs_mois(KTm,omesm)
HdmE = Hm[im]*Hdtm
HbmE = Hm[im] - HdmE
print ('\t Hm direct (Erbs)  = ','%.2f' % HbmE,' MJ/jr-m2')
print ('\t Hm diffus (Erbs)  = ','%.2f' % HdmE,' MJ/jr-m2')
print ('c ) Exemple relation Collares')

#
# exemple d'utilisation des formules de Collares
#
# evaluation des valeurs horaires à partir des valeurs journalieres
# pour la tranche horaire 12:00 13:00 heure solaire
#
ome = 7.5   # 12:30 milieu de la tranche horaire
a = 0.409 + 0.5016*sind(omes - 60)
b = 0.6609 - 0.4767*sind(omes - 60)

rt1 = pi/24*(a + b*cosd(ome))*(cosd(ome)-cosd(omes))/(sind(omes)-pi*omes/180*cosd(omes))
rt2 = Collares_total(ome,omes)
I_Corales = rt2*H[ij]
print ('\t irradiation horaire totale fichier meteo', '%.2f' % Ia,' MJ/hr-m2')
print ('\t irradiation horaire totale  estimée par Corrales','%.2f' % I_Corales,' MJ/hr-m2')
rd1 = pi/24*(cosd(ome)-cosd(omes))/(sind(omes)-pi*omes/180*cosd(omes))
rd2 = Collares_diffus(ome,omes)
Id_Corales = rd2*Hd[ij]
Ib_Corales = I_Corales - Id_Corales
print ('\t irradiation horaire diffuse fichier meteo', Ida)
print ('\t irradiation horaire diffuse Erbs', IdE)
print ('\t irradiation horaire diffuse  estimée par Corrales', Id_Corales)
print ('\t irradiation horaire directe fichier meteo', Iba)
print ('\t irradiation horaire directe Erbs', IbE)
print ('\t irradiation horaire directe  estimée par Corrales', Ib_Corales)

