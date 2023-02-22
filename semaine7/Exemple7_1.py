#coding: utf-8
#
# Exemple 7.1
#
# Calcul de la température de sortie du capteur
#
#
from solar_mod import *
import numpy as np
from matplotlib.pyplot import *
#
# importer les donnnées du fichier
#


heure = np.array([7,8,9,10,11,12])
qpputile = np.array([0,0,300,400,500,200])
charge = np.array([2000,2000,0,0,0,2000])
Ac=5
qutile =  qpputile*Ac  #

Ta=20
mCp=200*4190
UA=12
Tavant1=65
Tavant2=65
dt=3600
n = len(heure)
#initialisation du tableau de sortie
np.insert(heure,0,6)
Tres1 = np.zeros(n+1)
Tres2 = np.zeros(n+1)
Tres1[0] = Tavant1
Tres2[0] = Tavant2
j= 1
ti = heure[0]
heure = np.insert(heure,0,ti-1)
qutile = np.insert(qutile,0,qutile[0])
charge = np.insert(charge,0,charge[0])
for i in range (0,n):
    x = 1 + (dt/mCp)*UA/2
    Tapres1 = Tavant1 + (dt/mCp)*(qutile[j]-charge[j]-UA*(Tavant1-Ta))
    Tres1[j] = Tapres1
    Tavant1 = Tapres1
    Tapres2 = (Tavant2 + (dt/mCp)*(qutile[j]-charge[j]-UA*(Tavant2/2-Ta)))/x
    Tmoy = (Tapres2 + Tavant2)/2
    Tapres= Tavant2 + (dt/mCp)*(qutile[i]-charge[i]-UA*(Tmoy-Ta))
    Tres2[j]=Tapres2
    Tavant2=Tapres2
    j=j+1
print(Tres1)
print(Tres2)

#faire le graphique
fig,ax = subplots()

ax.plot(heure,Tres1,'red',marker='s',label='Temperature 1 du réservoir ')
ax.plot(heure,Tres2,'blue',marker='o',label='Temperature 2 du réservoir ')
ax2=ax.twinx()
ax2.plot(heure,qutile,color='green',marker='+',label='Flux utile')
ax2.plot(heure,charge,color='black',marker='*',label='Flux charge')

#pour légende"
ax.plot(nan,color='green',marker='+',label='Flux utile')
ax.plot(nan,color='black',marker='*',label='Flux charge')


ax.set_xlabel("Heure",fontsize =12)
ax.set_ylabel("Temperature du réservoir",fontsize =12)
ax2.set_ylabel("Energie",fontsize =12)

ax.legend(fontsize =12)
show()