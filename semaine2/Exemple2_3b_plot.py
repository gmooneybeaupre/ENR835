#coding: latin-1
from __future__ import unicode_literals
import sys
from solar_mod import *
from numpy import *
from collections import *
from matplotlib.pyplot import *

jm = array([17,47,75,105,135,162,198,228,258,288,318,344])
nm = array([31,28,31,30,31,30,31,31,30,31,30,31])
hrm = array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures ?coul?es apr?s chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
i1 = 0
phi= 45.5
njour = 11
gamso = zeros((7,11))
alpo = zeros((7,11))
nv = zeros(11)
nvv = zeros(11)
nv[0] = jour_mois_jour_annee(21,'dec')
nv[1] = jour_mois_jour_annee(21,'jan')
nv[2] = jour_mois_jour_annee(9,'fev')
nv[3] = jour_mois_jour_annee(23,'fev')
nv[4] = jour_mois_jour_annee(8,'mars')
nv[5] = jour_mois_jour_annee(21,'mars')
nv[6] = jour_mois_jour_annee(3,'avr')
nv[7] = jour_mois_jour_annee(16,'avr')
nv[8] = jour_mois_jour_annee(1,'mai')
nv[9] = jour_mois_jour_annee(21,'mai')
nv[10] = jour_mois_jour_annee(21,'juin')
nvv[0] = jour_mois_jour_annee(21,'dec')
nvv[1]  = jour_mois_jour_annee(22,'nov')
nvv[2]  = jour_mois_jour_annee(3,'nov')
nvv[3]  = jour_mois_jour_annee(17,'oct')
nvv[4]  = jour_mois_jour_annee(6,'oct')
nvv[5]  = jour_mois_jour_annee(21,'sept')
nvv[6]  = jour_mois_jour_annee(10,'sept');
nvv[7]  = jour_mois_jour_annee(28,'aout')
nvv[8]  = jour_mois_jour_annee(12,'aout')
nvv[9]  = jour_mois_jour_annee(24,'juil')
nvv[10]  = jour_mois_jour_annee(21,'juin')
dth = pi/20.0
th = arange(0,2*pi+dth,dth)
nth = len(th)
x = zeros(nth)
y = zeros(nth)
alpc = list()
gamc = list()
#
# plot les cerckles correspondants aux angles de zénith ( 10 à 90 °)
#
fig = figure()
nr = 9
dr = 1.0/nr
t = text(-.25,1.05,'Trajectoire Solaire Montreal', fontsize=15)
for i in range(0,9):
    r = (i+1)*dr
    for j in range(0,nth):
        x[j] = r*cos(th[j])
        y[j] = r*sin(th[j])
    plot(x,y,'k')
    strr = str(90 - 10*(i+1))+'°'
    text(y[0],x[0]+0.02,strr)
#
# plot les lignes correspondants aux angles d'azimuth ( -180 à 180 °)
#
rv = array([0,1])
thv = array([1,1])
x = zeros(2)
y = zeros(2)
for i in range(0,24):
    gamm = -180.0+15.0*i
    th = gamm*pi/180.0*thv
    for j in range(0,2):
        x[j] = rv[j]*cos(th[j])
        y[j] = rv[j]*sin(th[j])
    plot(x,y,'k')
    angle = 15.0*i
    angl = (angle<=180)*angle + (angle>180)*(angle-360.0)
    strr = str(angl)+'°'
    xx =  (sign(y[1])+1)/2.0*(y[1]+0.05)-(sign(y[1])-1)/2.0*(y[1]-0.12)
    text(xx,x[1]+sign(x[1])*0.02,strr)
axis('equal')
#
# plot des positions solaires
#
for i in range(0,11):
    n = nv[i]
    delt  = decl_solaire(n)      # calcul de la declinaison solaire
    com = -tand(phi)*tand(delt)
    oms = arccosd(com)
    om = arange(-oms,oms)
    no = len(om)
    alp = zeros(no)
    gams = zeros(no)
    for j in range(0,no):
        ome = om[j]
        thez = zenith_solaire(phi,delt,ome)     # calcul du zenith solaire
        gams[j] = azimuth_solaire(thez,delt,phi,ome)   # calcul  de l'azimuth solaire
        alp[j] = 90.0-thez
    alpc.append(alp)
    gamc.append(gams)
    i1 = i1+nm[i]
    for id in range(0,7):
        thez = zenith_solaire(phi,delt,-15*(id+1))     # calcul du zenith solaire
        gamso[id,i] = azimuth_solaire(thez,delt,phi,-15*(id+1))   # calcul  de l'azimuth solaire
        alph = 90-thez
        alph = max(alph,0)
        alpo[id,i] = alph
for i in range(0,11):
    rv = 1.0-alpc[i]/90.0
    thv = (gamc[i]-90.0)
    nrv = len(rv)
    x = zeros(nrv)
    y = zeros(nrv)
    for j in range(0,nrv):
        x[j] = rv[j]*cosd(thv[j])
        y[j] = rv[j]*sind(thv[j])
    plot(x,y,'k')
    jour,mois = jour_annee_jour_mois(nv[i])
    strr = str(jour)+' '+mois
    xx = x[0]-.15
    yy = x[len(x)-1]+.04
    if (i>0) & (i< (njour -1)):
        jour,mois = jour_annee_jour_mois(nvv[i])      # modification 19-fev-2020
        strr2 = '%.0f' % jour +' '+ mois
        strr = strr + '-'+ strr2
        xx = x[0]-.35
        yy = x[len(x)-1]+.1
    text(xx,y[0]-0.02,strr)
    text(yy,y[0]-0.02,strr)
A = zeros((11,2))
alph = zeros(11)
gams = zeros(11)
def cherche_del(x,phi,ome):
    thez = zenith_solaire(phi,x,ome)
    return thez - 90

for i in range(0,7):
    A[:,0] = gamso[i,:]
    A[:,1] = alpo[i,:]
    B = A[A[:,0].argsort(),]
    k = 0
    flag = 0
    for j in range(0,11):
        if (B[j,1] != 0):
            alph[k] = B[j,1]
            gams[k] = B[j,0]
            k = k+1
        elif flag == 0:
            flag = 1
            om = -15.0*(i+1)
            dd = brentq(cherche_del,-24.0,24,args=(phi,om))
            alph[k] = 0.0
            gams[k] = azimuth_solaire(90.0,dd,phi,om)
            k = k+1
    rv = 1.0-alph/90.0
    thv = (gams-90.0)*pi/180.0
    if k > 0:
        nrv = k
        x = zeros(nrv)
        y = zeros(nrv)
        for j in range(0,nrv):
            x[j] = rv[j]*cos(thv[j])
            y[j] = rv[j]*sin(thv[j])
        plot(x,y,'r')
        strr = str(i+1)+' '+' PM'
        text(x[0],y[0]+0.05,strr)
        thv = (-gams-90.0)*pi/180.0
        x = zeros(nrv)
        y = zeros(nrv)
        for j in range(0,nrv):
            x[j] = rv[j]*cos(thv[j])
            y[j] = rv[j]*sin(thv[j])
        plot(x,y,'r')
        strr = str(11-i)+' '+' AM'
        text(x[0],y[0]+0.05,strr)
# exemple
xticks([])
yticks([])
axis('equal')


xm = arange(-8,52)
ym = 26.0
zm = 45.0
x1 = -8.0
x2 = 52.0
tan1 = ym/zm
n = len(xm)
gams = zeros(n)
thez = zeros(n)
for i in range(0,n):
    psi = arctand(tan1)
    gams[i] = arctand(xm[i]/zm)
    ta = cosd(gams[i])*tan1
    als = arctand(ta)
    thez[i] = 90.0 - als
rv = thez/90.0
thv = -90.0 - gams
x = zeros(n)
y = zeros(n)
for j in range(0,n):
    x[j] = rv[j]*cosd(thv[j])
    y[j] = rv[j]*sind(thv[j])
def cherche_thez(x,phi,gams,delt):
    y = cosd(gams)*cosd(phi)*sind(x) - cosd(x)*sind(phi) + sind(delt)
    return y
x2 = zeros(n)
y2 = zeros(n)
nj = nv[0]
delt  = decl_solaire(nj)      # calcul de la declinaison solaire
for j in range(0,n):
    thee = brentq(cherche_thez,10,90,args=(phi,gams[j],delt))
    rv2 = thee/90.0
    x2[j] = rv2*cosd(thv[j])
    y2[j] = rv2*sind(thv[j])
x3 = array([x[0],x2[0]])
y3 = array([y[0],y2[0]])
x4 = array([x[n-1],x2[n-1]])
y4 = array([y[n-1],y2[n-1]])
plot(x,y,'b',x2,y2,'b',x3,y3,'b',x4,y4,'b',linewidth=4)

show()
