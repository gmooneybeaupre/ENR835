#coding: utf-8
from __future__ import unicode_literals
import sys
from solar_mod import *
import numpy as np
from matplotlib.pyplot import *

jm = np.array([17,47,75,105,135,162,198,228,258,288,318,344])
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures ?coul?es apr?s chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
i1 = 0
phi= 45.5
njour = 11
gamso = np.zeros((7,11))
alpo = np.zeros((7,11))
nv = np.zeros(11)
nvv = np.zeros(11)
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
th = np.arange(0,2*pi+dth,dth)
nth = len(th)
x = np.zeros(nth)
y = np.zeros(nth)
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
        x[j] = r*np.cos(th[j])
        y[j] = r*np.sin(th[j])
    plot(x,y,'k')
    strr = str(90 - 10*(i+1))+'°'
    text(y[0],x[0]+0.02,strr)
#
# plot les lignes correspondants aux angles d'azimuth ( -180 à 180 °)
#
rv = np.array([0,1])
thv = np.array([1,1])
x = np.zeros(2)
y = np.zeros(2)
for i in range(0,24):
    gamm = -180.0+15.0*i
    th = gamm*pi/180.0*thv
    for j in range(0,2):
        x[j] = rv[j]*np.cos(th[j])
        y[j] = rv[j]*np.sin(th[j])
    plot(x,y,'k')
    angle = 15.0*i
    angl = (angle<=180)*angle + (angle>180)*(angle-360.0)
    strr = str(angl)+'°'
    xx =  (np.sign(y[1])+1)/2.0*(y[1]+0.05)-(np.sign(y[1])-1)/2.0*(y[1]-0.12)
    text(xx,x[1]+np.sign(x[1])*0.02,strr)
axis('equal')

#
# plot des positions solaires
#
for i in range(0,11):
    n = nv[i]
    delt  = decl_solaire(n)      # calcul de la declinaison solaire
    com = -tand(phi)*tand(delt)
    oms = arccosd(com)
    om = np.arange(-oms,oms)
    no = len(om)
    alp = np.zeros(no)
    gams = np.zeros(no)
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
    x = np.zeros(nrv)
    y = np.zeros(nrv)
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
A = np.zeros((11,2))
alph = np.zeros(11)
gams = np.zeros(11)
def cherche_del(x,phi,ome):
    thez = zenith_solaire(phi,x,ome)
    return thez - 90
#
# plot des lignes d'heures
#
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
        x = np.zeros(nrv)
        y = np.zeros(nrv)
        for j in range(0,nrv):
            x[j] = rv[j]*np.cos(thv[j])
            y[j] = rv[j]*np.sin(thv[j])
        plot(x,y,'r')
        strr = str(i+1)+' '+' PM'
        text(x[0],y[0]+0.05,strr)
        thv = (-gams-90.0)*pi/180.0
        x = np.zeros(nrv)
        y = np.zeros(nrv)
        for j in range(0,nrv):
            x[j] = rv[j]*np.cos(thv[j])
            y[j] = rv[j]*np.sin(thv[j])
        plot(x,y,'r')
        strr = str(11-i)+' '+' AM'
        text(x[0],y[0]+0.05,strr)
# exemple
xticks([])
yticks([])
axis('equal')
#figManager = get_current_fig_manager()
#figManager.resize(*figManager.window.maxsize())
#fig.savefig('temp.png', dpi=fig.dpi)

xm = np.arange(-8,52)
ym = 26.0
zm = 45.0
x1 = -8.0
x2 = 52.0
gams = np.arange(-88,88)
tan1 = ym/zm
n = len(gams)
#gams = np.zeros(n)
thez = np.zeros(n)
for i in range(0,n):
    psi = arctand(tan1)
    ta = cosd(gams[i])*tan1
    als = arctand(ta)
    thez[i] = 90.0 - als
rv = thez/90.0
thv = -90.0 - gams
xb = np.zeros(n)
yb = np.zeros(n)
for j in range(0,n):
    xb[j] = rv[j]*cosd(thv[j])
    yb[j] = rv[j]*sind(thv[j])
def cherche_thez(x,phi,gams,delt):
    y = cosd(gams)*cosd(phi)*sind(x) - cosd(x)*sind(phi) + sind(delt)
    return y

nj = nv[0]  # soltice hiver
delt  = decl_solaire(nj)      # calcul de la declinaison solaire
com = -tand(phi)*tand(delt)
oms = arccosd(com)
om = np.arange(-oms,oms)
no = len(om)
x2 = np.zeros(no)
y2 = np.zeros(no)
for j in range(0,no):
    ome = om[j]
    thez = zenith_solaire(phi,delt,ome)     # calcul du zenith solaire
    gam = azimuth_solaire(thez,delt,phi,ome)   # calcul  de l'azimuth solaire
    rv2 = thez/90.0
    thx = -90 - gam
    x2[j] = rv2*cosd(thx)
    y2[j] = rv2*sind(thx)
plot(xb,yb,'b',x2,y2,'b',linewidth=4)
show()
