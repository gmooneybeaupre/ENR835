#coding: utf-8
from solar_mod import *
#from properties_mod import *
import numpy as np
from matplotlib.pyplot import *
from matplotlib.sankey import Sankey

j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])               # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
mois = ['janv','févr','mars','avril','mai','juin','juil','aout','sept','oct','nov','déc']
#
# lecture fichier de sortie
#
N = 2
Ac = 2.874*N
donnees = np.loadtxt("devoir4.txt",skiprows = 1)
temps = donnees[:,0]   # heures
dt = temps[1] - temps[0]    # delta t
It = donnees[:,1]
qu_capteur = donnees[:,2]
qin_res = donnees[:,3]
qou_res = donnees[:,4]
qaux = donnees[:,5]
qpertes  = donnees[:,6]
DeltaEp  = donnees[:,7]      # Variation d'énergie emmagasinée
Eaux = sum(qaux)*dt # kJoules
Epertes = sum(qpertes)*dt  # kJoules
EperteskWh = Epertes/3600
EauxkWh = Eaux/3600
Eout_res = sum(qou_res)*dt
Echarge =  Eout_res + Eaux  # kJoules
EchargekWh = Echarge/3600
Eutile = sum(qu_capteur)*dt # kJoules
EutilekWh = Eutile/3600
Esoleil = sum(It)*dt*Ac # kJoules
Ein_res = -sum(qin_res)*dt # kJoules
Ein_reskWh = Ein_res/3600
DeltaE =  sum(DeltaEp)*dt # kJoules
F = 1  - Eaux/Echarge
print('F = ',F)
rend = Ein_res/Esoleil
print('rend capteurs = ',rend)
Eutile_load = Ein_res - Epertes
rend_sys = Eutile_load/Esoleil
print('rend système = ',rend_sys)
fv = np.zeros(12)
L = np.zeros(12)
i1 = 0
for im in  range(0,12):
    i2 = int(hrr[im]/dt)
    Echarge_res_m = sum(qou_res[i1:i2])*dt
    Eaux_m = sum(qaux[i1:i2])*dt # kJoules
    Echarge_m = Echarge_res_m + Eaux_m# kJoules
    fv[im] = max(1  - Eaux_m/Echarge_m,0)
    L[im] = Echarge_m
    i1 = i2 + 1
    print('Fraction solaire du mois de ','%s' % mois[im], ' = ','%.2f'  % (100*fv[im]),' %' )
F2 = sum(L*fv)/sum(L)
print('Fraction solaire annuelle TRNSYS','%.2f' % (100*F2),' %' )
#
# exemple dpour visualier la fraction soalire

im = np.arange(1,13)

#
bar(im,fv, width = 0.5, color = 'red')
title('Fraction solaire mensuelle')
ax = gca()
ax.set_xticks(im)
ax.set_xticklabels(mois, minor=False, rotation=45)

#
# bilan reservoir
#
Ein = Ein_res
Eout =  Eout_res +  Epertes
Delta_theo = Ein - Eout
print('Variation d''énergie du réservoir par bilan =',Delta_theo, ' KJ')
print('Variation d''énergie du réservoir TRNSYS =',DeltaE, ' KJ')
#
# Flow d'énergie par diagramme de Sankey
#
flag_plot = True
if (flag_plot):
    E1 = Esoleil/1e6  # GJoules
    Eperdue_capteurs = Esoleil - Eutile # pertes capteurs
    E2 = Eperdue_capteurs/1e6 # GJoules
    E4 = Epertes/1e6   # pertes reservoir
    E5 = (Eout_res)/1e6  # Energie solaire utile à la charge
    figure = Sankey(format='%.1f',unit ='GJoules',offset = 4.5,scale = 0.5,margin = 1,gap=0.25, radius=0.1)
    figure.add(flows=[E1, -E2 , -E4, -E5],
           labels=['Solaire', 'Pertes collecteur', 'Pertes reservoir', 'Charge'],
           orientations=[ 0,  1, 1,  0],
           trunklength = 20,
           pathlengths = [0.25, 1,  3.75, 4.25],linewidth = 1)
    figure.finish()
    title("Diagramme de Sankey")
show()


