#coding: utf-8
from solar_mod import *
#from properties_mod import *
import numpy as np
from matplotlib.pyplot import *
from matplotlib.sankey import Sankey

# DONNÉES PRÉLIMINAIRES
j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])               # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
mois = ['janv','févr','mars','avril','mai','juin','juil','aout','sept','oct','nov','déc']
H_b = np.array([1.5817, 2.5884, 3.8525, 4.5808, 5.6846, 6.2114, 6.0228, 4.8676, 4.0267, 2.5372, 1.2235, 1.1742]) # kWh/m^2-jr
H_b = H_b * 3.6e6 # J/m^2-jr
tabb_tan = np.array([0.96, 0.96, 0.96, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.96, 0.96])
Tinf = np.array([-10.3, -8.8, -2.4, 5.7, 12.9, 18, 20.8, 19.4, 14.5, 8.3, 1.6, -6.9])
rhog = np.array([0.7, 0.7, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0,4])
N = 2
phi = 45.5
gamma = 0 
beta = 45
Ac = 2.874*N
ep = 0.9
cp = 4190
mp = 71.32/3600
Ccoll = cp*mp
Cmin = Ccoll
eta0 = 0.717
Fr_tan = eta0
a1 = 14.45076/3.6
b0 = 0.11
b1 = 0.051
FrUL = a1
Tref = 100
Frp_Fr = 1/(1 + ((Ac*FrUL)/Ccoll)*((Ccoll/(ep*Cmin))-1))

# EXTRACTION DES DONNÉES TRNSYS
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

# MÉTHODE F-CHART
print("")
print("Méthode F-chart")
fchart_m = np.zeros(12)
L = L*3600
sum1 = 0
sum2 = 0

for m in range(0,12):
    njt = j_type[m]
    Ho_b = irradiation_extraterrestre_jour(njt,phi)
    Kt_b = H_b[m]/Ho_b
    delt = decl_solaire(njt)
    omesm = angle_sunset(phi,delt)
    Hd_b = Erbs_mois(Kt_b,omesm)
    Hb_b = H_b[m] - Hd_b
    Rb_b = Calcul_Rb_mois(phi,njt,beta,gamma)
    Ht_b,HmtbE,HmtdE,HmtrE = modele_isotropique(H_b[m],Hb_b,Hd_b,beta,Rb_b,rhog[m])

    Xm = FrUL*Frp_Fr*(Tref-Tinf[m])*(3600*24*nm[m])*(Ac/L[m])
    Ym = Fr_tan*tabb_tan[m]*Frp_Fr*nm[m]*(Ac/L[m])*Ht_b
    fchart_m[m] = fchart(Xm,Ym)
    
    sum1 = sum1 + fchart_m[m]*L[m]
    sum2 = sum2 + L[m]
    
    mess = "F-chart pour le mois de " + mois[m] + " : "
    print(mess,'%.2f' % (100*fchart_m[m]), " %")

fchart_annee = sum1/sum2
mess = "F-chart pour l'année : "
print(mess, '%.2f' % (100*fchart_annee), " %")    
print("")

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
flag_plot = False
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


