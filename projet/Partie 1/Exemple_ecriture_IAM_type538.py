#coding: utf-8
#
# exemple ecriture fichier
#
import numpy as np
NomFichier = 'Thermomax.txt'
Fichier = open(NomFichier,'w')      # instanciation de l'objet Fichier de la classe file
thetr = np.array([0,10,20,30,40,50,60,90])
thelong = thetr
Kthet  = np.array([1.0,1.01,1.02,1.04,1.04,0.99,0.90,0])
Kthel  = np.array([1.0,1.0,0.99,0.97,0.95,0.91,0.83,0])
n = len(thetr)
for i  in range(0,n):
    Fichier.write(str(thetr[i])+' ')
Fichier.write('\n')
for i  in range(0,n):
    Fichier.write(str(thelong[i])+' ')
Fichier.write('\n')
for i  in range(0,n):
    for j  in range(0,n):
        K = Kthet[i]*Kthel[j]
        Fichier.write(str(K)+'\n')
Fichier.close()