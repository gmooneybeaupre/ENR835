#coding: utf-8
#
# exemple solaire 2.3
#
# Calcul de l'ombrage sur un capteur solaire
#
#
from solar_mod import *
import numpy as np
from matplotlib.pyplot import *
from scipy.optimize import newton,brentq


def change_angle(X):
    x1 = 360*(284+1)/365
    x2 = 360*(284+365)/365
    y1 = x1 - 360
    if X < y1:
        Y = X + 720
    else:
        Y = X + 360
    return Y

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
#
#
#
# debut du programme
#
ym = 26.0
zm = 45.0
tps = ym/zm
psi = arctand(tps)

n1 =  jour_mois_jour_annee(17,'oct')
n2 =  jour_mois_jour_annee(17,'aout')
n3 =  jour_mois_jour_annee(17,'dec')
n = n1
delta = decl_solaire(n)
sde = sind(delta)
cde = cosd(delta)
cosom = -tand(delta)*tand(phi)
ome = arccosd(cosom)
thez_midi = zenith_solaire(phi,delta,0)
alps_midi = 90 - thez_midi
print('alpha midi = ',alps_midi)
#
# heure où le soleil devient apparent le 21 septembre
#
costh  = tand(psi)*sind(delta)/(tand(psi)*sind(phi) -cosd(phi))
thez_deb = arccosd(costh)
alps_deb = 90 - thez_deb
print('alpha debut = ',alps_deb)
#
# verification
#
come = (costh - sind(delta)*sind(phi))/(cosd(delta)*cosd(phi))
ome_deb = -arccosd(come)
gams = azimuth_solaire(thez_deb,delta,phi,ome_deb)
talp = tand(alps_deb)/cosd(gams)
alpha_profil = arctand(talp)
print('alph_profil = ',alpha_profil,psi)
debut = heure_angle(ome_deb)
fin = heure_angle(-ome_deb)
print('omega debut',ome_deb)
print('heure debut',debut.heure,debut.minu)
print('heure fin',fin.heure,fin.minu)
#
#
# on recherche quand alpha = psi à midi
#
thezn = 90 - psi
delta = phi - thezn
print('delta = ',delta)
X1 = arcsind(delta/23.45)
#print(X1)
Y1 = change_angle(X1)
n1 = round(365*Y1/360 - 284)
jour1,mois1 = jour_annee_jour_mois(n1)
print('n1 = ',n1)
print(jour1,mois1)
X2 = 180 - X1
Y2 = change_angle(X2)
n2 = round(365*Y2/360 - 284)
jour2,mois2 = jour_annee_jour_mois(n2)
print('n2 = ',n2)
print(jour2,mois2)