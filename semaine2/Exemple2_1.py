#coding: utf-8
#
# importation des modules utilisés
import numpy as np
from solar_mod import *
pi = np.pi
#
# Vecteurs utilisés dans certains calculs
j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])
jm = np.array([17,47,75,105,135,162,198,228,258,288,318,344])      # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])               # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
Gsc = 1367.0
#
def calcul_n(j,m):
    if m<3:
        n = j + 31*(m-1)
    else:
        n = j + 31*(m-1) - np.floor(0.4*m + 2.3)
    return int(n)



#
phi = 45.5  # latitude  Montreal
lon = -73.0 -35.0/60.0  # longitude  Montreal (73 deg 35 min ouest)
Lst = -75        # longitude  méridien HNE -5
#
n = jour_mois_jour_annee(5,'mai')
#
# ou Eq 1.1
#
n2 = calcul_n(5,5)
print ('n ( 5 mai) = ',n,n2)
del_h = 1.0  # heure avancée
heure = 12.0
minute = 0.0
jour = n
temps_civil = temps(jour,heure,minute)
E = equation_temps(n)
Delta_temps = 4*(lon - Lst)+E
# avec ce delta je peux trouver l'heure solaire. Mais je peux aussi faire appel
# à une fonction qui est dans ma librairie
# calcul de l'heure solaire pour midi le 5 mai
temps_sol = heure_solaire(lon,Lst,del_h,temps_civil)
#
# si on ne met pas d'heure d'horloge, on calcul le temps actulel
#
#
#
# question 1b)
#
# exemple du calcul de l'angle horaire à partir de l'heure solaire
#avec le calcul
omecalc = 15*((temps_sol.heure+temps_sol.minu/60)-12)
#ou avec une fonction

#
# exemple formatté
#
#print ('ome  = ', '%.2f'% ome, '°')
# exemple du calcul de la declinaison solaire 5 mai
delta = 23.45*sind(360*(284+n)/365)
#
# ou appel à la fonction du module solaire
#
deltb  = decl_solaire(n)
# ou appel à la fonction du module solaire en utilisant 1.5 b
#
print ('del = ( 5 mai )','%.2f'% delta,'°')
#
# exemple du calcul du zenith solaire 5 mai à midi heure civile
#
#cosdz = sind(delta)*sind(phi) + cosd(delta)*cosd(phi)*cosd(ome)
#theza = np.arccos(cosdz)*180/pi
#
# ou appel à la fonction du module solaire
#
#thezb = zenith_solaire(phi,delta,ome)
#print ('thez = ','%.2f'% theza, '°')
# exemple du calcul de l'azimuth  solaire 5 mai à midi heure civile
#gams = azimuth_solaire(theza,delta,phi,ome)
#print ('gams = ','%.2f'% gams, '°')
#
#  Exemple 2.1 a)
#
# exemple de l'irradiation extra-terrestr sur le plan horizontal le 5 mai de midi solaire à 1 heure PM solaire
#
Gon =  irradiation_extraterrestre_normale(n)
Gon2 = 1367*(1+0.033*cosd(360*n/365))




print ('Gon ( 5 mai)  = ','%.2f'% Gon, 'W/m2')
om1 = 0
om2 = 15.0
Ioa1 = Gon*12.0*3600/pi*(cosd(phi)*cosd(delta)*(sind(om2)-sind(om1))+pi*(15.0)/180.0*sind(phi)*sind(delta))
#
# example d'applel à la fonction irradiation horaire sur plan horizontal
#
Ioa1b = irradiation_extraterrestre_horaire(n,phi,om1,om2)
Ioa2 = Ioa1/1e6 # transforme en MJ/hr
print('Réponse question a)')
print('\t  I o  (midi- 1 heure  solaire)= ','%.3f' % Ioa2 ,' MJ/hr-m2')
IoakW = Ioa2/3.6 # valeur en kWatts/m2
print('\t  I o  (midi- 1 heure  solaire)= ','%.3f' % IoakW ,' kW/m2')
#
#  Exemple 2.1 b)
#
print('Réponse question b)')
print ('\t heure solaire (midi 5 mai) = ',temps_sol.heure)
print ('\t minute solaire (midi 5 mai)= ',temps_sol.minu)
# exemple de l'irradiation extra-terrestr sur le plan horizontal le 5 mai de midi solaire à 1 heure PM solaire
#
ome = angle_horaire(temps_sol)
om1 = ome
om2 = ome+ 15.0
Iob1 = irradiation_extraterrestre_horaire(n,phi,om1,om2)
Iob2 = Iob1/1e6 # transforme en MJ/hr
print('\t I o  (midi- 1 heure  civile)= ','%.3f' % Iob2 ,' MJ/hr-m2')
Iob4 = irradiation_extraterrestre_horaire(n,phi,-7.5,7.5)
print('\t I o  (max)= ','%.3f' % (Iob4/1e6) ,' MJ/hr-m2')
IobkW = Iob2/3.6 # valeur en kWatts/m2
print('\t I o  (midi- 1 heure  civile)= ','%.3f' % IobkW ,' kW/m2')
#
# Le coucher du soleil est entre 19 et 20 hrs
print('Réponse question b) (mod) entre 19 hrs et 20 heure heure solaire')
om1b = 7*15
om2b = om1b + 15
#
# Calcul de l'irradiation solaire entre 19 hrs et 20 hrs
#
coss = -tand(delta)*tand(phi)
omes = arccosd(coss)
# mauvaise facon
#
Iob1 = Gon*12.0*3600/pi*(cosd(phi)*cosd(delta)*(sind(om2b)-sind(om1b))+pi*(15.0)/180.0*sind(phi)*sind(delta))
#
# bonne facon  ( omega 2 = omes coucher du soleil)
#
Iob2 = Gon*12.0*3600/pi*(cosd(phi)*cosd(delta)*(sind(omes)-sind(om1b))+pi*(omes - om1b)/180.0*sind(phi)*sind(delta))
#
Iob3 = irradiation_extraterrestre_horaire(n,phi,om1b,om2b)
print('\t  I o  (coucher du soleil solaire mauvaise facon))= ','%.3f' % (Iob1/1e6) ,' MJ/hr-m2')
print('\t  I o  (coucher du soleil  solaire)= ','%.3f' % (Iob2/1e6) ,' MJ/hr-m2')
print('\t  I o  (coucher du soleil  solaire)= ','%.3f' % (Iob3/1e6) ,' MJ/hr-m2')

#
# 2.1 c
#
# exemple de l'irradiation extra-terrestre sur le plan horizontal pour le 5 mai
#
#
# calcul de l'angle du coucher du soleil
#
print('Réponse question c)')

print ('\t Angle lever-coucher soleil (5 mai)  = ','%.2f' % omes ,'°')
Ho1 = Gon*24.0*3600/pi*(cosd(phi)*cosd(delta)*sind(omes)+pi*(omes)/180.0*sind(phi)*sind(delta))
Ho = irradiation_extraterrestre_jour(n,phi)
Ho1 = Ho1/1e6
Ho2 = Ho/1e6
print ('\t Ho (5 mai)  = ','%.2f' % Ho1 ,' MJ/jr-m2')
print ('\t Ho (5 mai)  = ','%.2f' % Ho2 ,' MJ/jr-m2')



#
# 2.1 d
print('Réponse question d)')
# exemple de l'irradiation extra-terrestr sur le plan horizontal le mois de mai
mois = 5 # mois de mai
Hom = irradiation_extraterrestre_jour_moyen(mois,phi)
Hom = Hom/1e6    # en MJ
print ('\t Hom ( mois de mai (2.6a))  = ','%.3f' %  Hom,' MJ/jr-m2' )
# exemple de l'irradiation extra-terrestr sur le plan horizontal le mois de mai par le jour type du mois
ntype = j_type[mois-1]
Hom2 = irradiation_extraterrestre_jour(ntype,phi)
delta_type = decl_solaire(ntype)
coss = -tand(delta_type)*tand(phi)
omes_type = arccosd(coss)
Gon_type = 1367*(1+0.033*cosd(360*ntype/365))
Hom3 = Gon_type*24.0*3600/pi*(cosd(phi)*cosd(delta_type)*sind(omes_type)+pi*(omes_type)/180.0*sind(phi)*sind(delta_type))

Hom2 = Hom2/1e6    # en MJ
Hom3 = Hom3/1e6    # en MJ
print ('\t Hom ( mois de mai (2.6b))  = ','%.3f' %  Hom2 ,' MJ/jr-m2' )
print ('\t Hom ( mois de mai (2.6b))  = ','%.3f' %  Hom3 ,' MJ/jr-m2' )
