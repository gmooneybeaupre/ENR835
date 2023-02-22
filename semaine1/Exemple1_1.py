#coding: utf-8
#
# importation des modules utilisés
import numpy as np
from solar_mod import *
pi = np.pi
#
# Vecteurs utilisés dans certains calculs
jm = np.array([17,47,75,105,135,162,198,228,258,288,318,344])      # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])               # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
Gsc = 1367.0
#
#
phi = 45.5  # latitude  Montreal
lon = -73.0 -35.0/60.0  # longitude  Montreal (73 deg 35 min ouest)
Lst = -75        # longitude  méridien HNE -5
#
n = jour_mois_jour_annee(5,'mai')
print ('n= ',n)
del_h = 1.0  # heure avancée
heure = 12.0
minute = 0.0
jour = n
temps_civil = temps(jour,heure,minute)
E = equation_temps(n)
Delta_temps = 4*(lon - Lst)+E
print('Delta T(min) = ',Delta_temps)
# avec ce delta je peux trouver l'heure solaire. Mais je peux aussi faire appel
# à une fonction qui est dans ma librairie
# calcul de l'heure solaire pour midi le 5 mai
temps_sol = heure_solaire(lon,Lst,del_h,temps_civil)
print ('heure solaire (midi 5 mai) = ',temps_sol.heure)
print ('minute solaire (midi 5 mai)= ',temps_sol.minu)
print ('jour  solaire = ',temps_sol.jour)
#
# si on ne met pas d'heure d'horloge, on calcul le temps actulel
#
temps_sol2 = heure_solaire(lon,Lst)
n2 = temps_sol2.jour
print ('heure solaire (aujourd''hui) = ',temps_sol2.heure)
print ('minute solaire  (aujourd''hui)= ',temps_sol2.minu)
print ('jour  solaire  (aujourd''hui) = ',temps_sol2.jour)
#
#
# question 1b)
#
# exemple du calcul de l'angle horaire à partir de l'heure solaire
#avec le calcul
omecalc = 15*((temps_sol.heure+temps_sol.minu/60)-12)
print ('omecalc = ',omecalc)
#ou avec une fonction
ome = angle_horaire(temps_sol)
print ('ome = ',ome)
#
# exemple formatté
#
print ('ome mis en forme = ', '%.2f'% ome, '°')
# exemple du calcul de la declinaison solaire 5 mai
delta = 23.45*sind(360*(284+n)/365)
#
# ou appel à la fonction du module solaire
#
deltb  = decl_solaire(n)
# ou appel à la fonction du module solaire en utilisant 1.5 b
#
deltc  = decl_solaire(n,cas = 2)
print ('del = ','%.2f'% delta,'°')
print ('del = ','%.2f'% deltb,'°')
print ('del = ','%.2f'% deltc,'°')
#
# exemple du calcul du zenith solaire 5 mai à midi heure civile
#
cosdz = sind(delta)*sind(phi) + cosd(delta)*cosd(phi)*cosd(ome)
theza = arccosd(cosdz)
#
# ou appel à la fonction du module solaire
#
thezb = zenith_solaire(phi,delta,ome)
print ('thez = ','%.2f'% theza, '°')
print ('thez = ','%.2f'% thezb, '°')
# exemple du calcul de l'azimuth  solaire 5 mai à midi heure civile
gams = azimuth_solaire(theza,delta,phi,ome)
print ('gams = ','%.2f'% gams, '°')
#
#  exemple 1.1c
#
delt2 = 23.45*sind(360*(284+n2)/365)
cosom = -tand(delt2)*tand(phi)
om = arccosd(cosom)
cosom2 = (cosd(90.833)-sind(delt2)*sind(phi))/(cosd(delt2)*cosd(phi))
om2 = arccosd(cosom2)
lever2,coucher2,duree2 = duree_jour(n2,phi)
lever_civil2 = heure_civile(lon,Lst,0,lever2)
coucher_civil2 = heure_civile(lon,Lst,0,coucher2)
print ('Heure lever  soleil aujourd''hui =', '%.0f'% lever_civil2.heure, 'h')
print ('Minute lever soleil aujourd''hui =', '%.0f'% lever_civil2.minu, 'min')
print ('Heure coucher  soleil  aujourd''hui =', '%.0f'% coucher_civil2.heure, 'h')
print ('Minute coucher soleil  aujourd''hui =', '%.0f'% coucher_civil2.minu, 'min')
print ('Durée du jour  aujourd''hui =', '%.0f'% duree2.heure, 'h', '%.0f'% duree2.minu, 'min')

#
# Exemple 1.1 d
beta = 40.0
gam = 0.0
the1 = normale_solaire(delta,phi,ome,beta,gam) # Eq 1.10
print ('theta  = ',the1)
the2 = normale_solaire2(theza,gams,beta,gam) # Eq 1.11
print ('theta  = ',the2)