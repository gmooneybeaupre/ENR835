#coding: latin-1
#
# exemple solaire 2.3
#
# Calcul de l'ombrage sur un capteur solaire
#
#
from solar_mod import *
import  numpy as np
from matplotlib.pyplot import *


j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])  # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
Gsc = 1367.0

def change_angle(X):
    x1 = 360*(284+1)/365
    x2 = 360*(284+365)/365
    y1 = x1 - 360
    if X < y1:
        Y = X + 720
    else:
        Y = X + 360
    return Y

phi = 45.5  # latitude  Montreal
lon = -73.0 -35.0/60.0  # longitude  Montreal (73 deg 35 min ouest)
Lst = -75        # longitude  méridien HNE -5
#
#
def fct(x,delta,gamm):
    y =(cosd(x)*sind(phi) - sind(delta))/(sind(x)*cosd(phi)) - cosd(gamm)
    return y

#
#
ym = 26.0
zm = 45.0
x1 = -8.0
x2 = 52.0
tan1 = ym/zm
psi = arctand(tan1)
gams1 = arctand(x1/zm)
gams2 = arctand(x2/zm)
als1 = arctand(tan1*cosd(gams1))
als2 = arctand(tan1*cosd(gams2))
thez1 = 90 - als1
thez2 = 90 - als2
#
# A lorsque alpA = psi
#

sind1 = cosd(thez1)*sind(phi) - cosd(gams1)*sind(thez1)*cosd(phi)
delta1 = arcsind(sind1)
cosd1 = cosd(delta1)
print('delta point A = ',delta1)
sn = delta1/23.45
X1 = arcsind(sn)
X2 = 180 - X1
Y1 = change_angle(X1)
Y2 = change_angle(X2)
n1 = round(365*Y1/360 - 284)
print('jour de l annee = ',n1)
jour1,mois1 =jour_annee_jour_mois(n1)
print('jour 1 = ',jour1,' du mois  ',mois1)
n2 = round(365*Y2/360 - 284)
print('jour de l annee = ',n2)
jour2,mois2 =jour_annee_jour_mois(n2)
print('jour 2 = ',jour2,' du mois  ',mois2)
#
# calcul de l'heure
#
come = (cosd(thez1) - sind(delta1)*sind(phi))/(cosd(delta1)*cosd(phi))
#come = (cosd(thez1) - sind1*sph)/(cosd1*cph)
ome = np.sign(gams1)*arccosd(come)
dh = ome/15
dm = dh*60
temps1 = heure_angle(ome)
print('temps A  = ',temps1.heure,' heure  ',temps1.minu,' minutes')
print ('omega = ',ome)
#
# point B lorsque alp_A = psi
#

thezn = newton(fct,50,args = (delta1,gams2))
cb = cosd(gams2)*cosd(phi)
sb = sind(phi)
di = -sind(delta1)
dd = np.sqrt(sb**2+cb**2)
x = arctand(sb/cb)
d = (di/dd)
thzB = arcsind(d) + x
print ('thez (PM) = ',thezn,thzB)
alpB = 90 - thzB
talpB = tand(alpB)/cosd(gams2)
alpB = arctand(talpB)
print('alpha profil point B (29 oct) =',alpB)
costhx = cosd(thzB)
comex = (costhx - sind(delta1)*sind(phi))/(cosd(delta1)*cosd(phi))
omex = arccosd(comex)
tempsx = heure_angle(omex)
print('temps X  = ',tempsx.heure,' heure  ',tempsx.minu,' minutes')
print ('omega = ',omex)

#
# B  lorsque alpB = psi
#
sind2 = cosd(thez2)*sind(phi) - cosd(gams2)*sind(thez2)*cosd(phi)
delta2 = arcsind(sind2)
cosd2 = cosd(delta2)
print('delta point B = ',delta2)
sn = delta2/23.45
X1 = arcsind(sn)
X2 = 180 - X1
Y1 = change_angle(X1)
Y2 = change_angle(X2)
n1 = round(365*Y1/360 - 284)
print('jour de l annee = ',n1)
jour1,mois1 =jour_annee_jour_mois(n1)
print('jour 1 = ',jour1,' du mois  ',mois1)
n2 = round(365*Y2/360 - 284)
print('jour de l annee = ',n2)
jour2,mois2 =jour_annee_jour_mois(n2)
print('jour 2 = ',jour2,' du mois  ',mois2)
#
# calcul de l'heure
#
come = (cosd(thez2) - sind(delta2)*sind(phi))/(cosd(delta2)*cosd(phi))
ome = np.sign(gams2)*arccosd(come)
temps2 = heure_angle(ome)
print('temps B  = ',temps2.heure,' heure  ',temps2.minu,' minutes')
print ('omega = ',ome)
#
# point A lorsque alpB = psi
#
thezn = newton(fct,50,args=(delta2,gams1))
cb = cosd(gams1)*cosd(phi)
sb = sind(phi)
di = -sind(delta2)
dd = np.sqrt(sb**2+cb**2)
x = arctand(sb/cb)
d = (di/dd)
thzA = arcsind(d) + x
print ('thez (PM) = ',thezn,thzA)
alpA = 90 - thzA
talpA = tand(alpA)/cosd(gams1)
alpA = arctand(talpA)
print('alpha profil point A (17 oct) =',alpA)
comex = (cosd(thzA) - sind(delta2)*sind(phi))/(cosd(delta2)*cosd(phi))
omex = -arccosd(comex)
tempsx = heure_angle(omex)
print('temps XB  = ',tempsx.heure,' heure  ',tempsx.minu,' minutes')
print ('omega = ',omex)
