#coding: utf-8
#
# exemple solaire 2.5
#
# Calcul de l'ombrage sur un capteur solaire
#
#
from solar_mod import *
import  numpy  as np
from matplotlib.pyplot import *


j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])  # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
Gsc = 1367.0
phi = 45.5  # latitude  Montreal
lon = -73.0 - 35.0/60.0  # longitude  Montreal (73 deg 35 min ouest)
Lst = -75        # longitude  méridien HNE -5
#
#
# a qulle heure ombrage le 21 octobre
#
gamc = 0
L = 2.1
beta = 60
A = L*sind(beta)
B = L*cosd(beta)
nx =jour_mois_jour_annee(21,'dec')
delt  = decl_solaire(nx)
coss = -tand(delt)*tand(phi)
omes = arccosd(coss)
ome1 = 3*15
thez = zenith_solaire(phi,delt,ome1)
gams = azimuth_solaire(thez,delt,phi,ome1)
als = 90-thez
al_midi = 90 - (phi - delt)
X1 = A*cosd(-gams-gamc)/tand(als)
X2 = A*cosd(gams-gamc)/tand(als)
X = max(X1,X2)
Z = L*sind(180 - (beta + al_midi))/sind(al_midi)
Z = X + B
Y = Z - B
psi = arctand(A/Y)
alp = arctand(tand(als)/cosd(gams-gamc))
print('psi = ',psi)
print('alp = ',alp)
omev = np.arange(-ome1,ome1+1)
n = len(omev)
alpv = np.zeros(n)
alsv = np.zeros(n)
psiv = np.zeros(n)
Fv = np.zeros(n)
th1 = 180 - psi - beta
for i in range(0,n):
    ome = omev[i]
    thezn = zenith_solaire(phi,delt,ome)
    gamsn = azimuth_solaire(thezn,delt,phi,ome)
    psiv[i] = psi
    alsv[i] = 90-thezn
    tanp = tand(alsv[i])/cosd(gamsn-gamc)
    alpv[i] = arctand(tanp)
    th2 = psi - alpv[i]
    th3 = 180 - th2 - th1
    th4 = 180 - th3
    d = sind(alpv[i])*Z/sind(th4)
    Fv[i] = max((L - d)/L,0)
figure(1)
plot(omev,alpv,omev,alsv,omev,psiv)
legend(('profil','solaire','psi'))
figure(2)
plot(omev,Fv)
show()
