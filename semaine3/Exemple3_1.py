#coding: utf-8
#
# exemple solaire 3.1
#
# Exemple de calcul de la radiation sur plan incliné par les 4 modeles standards
#
from sys import *
from solar_mod import *
import numpy as np

jm = np.array([17,47,75,105,135,162,198,228,258,288,318,344])# jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
Gsc = 1367

phi = 40  # latitude  Madison
#
#
# donnes du problème
#
I = 1.04e6  # Irradiation totale plan horizontal
beta = 60.0 # inclinaison capteur
rhog = 0.6  # albédo
gam = 0.0   # azimuth
#
# Calcul pour le 20 février
#
jour = 20
mois = 'fev'
n = jour_mois_jour_annee(jour,mois);
ome1 = -3.0*15.0    # angle correspondant à 9h00 AM heure solaire
ome2 = ome1+15.0
Io = irradiation_extraterrestre_horaire(n,phi,ome1,ome2)
kt = I/Io
#
# Calcul de I direct et diffus a l'aide des formule de Erbs
#
Idt = Erbs_horaire(kt)
Id = I*Idt
Ib = I - Id
#
# Calcul du coeff Rb
#
omen = (ome1+ome2)/2.0
Rb = Calcul_Rb(phi,n,omen,beta,gam)
#
# modele isotropique
#
Itb1 = Ib*Rb
Itd1 = Id*(1+cosd(beta))/2.0
Itr1 = I*rhog*(1-cosd(beta))/2.0
It1 =  Itb1+Itd1+Itr1
#
# modele HD
#
Ai = Ib/Io
Itb2 = (Ib+Id*Ai)*Rb
Itd2 = Id*(1-Ai)*(1+cosd(beta))/2.0
Itr2 = I*rhog*(1-cosd(beta))/2.0
It2 = Itb2+Itd2+Itr2
# model HDKR
f = np.sqrt(Ib/I)
Itb3 = (Ib+Id*Ai)*Rb
Itd3 = Id*(1-Ai)*(1+cosd(beta))/2.0*(1+f*sind(beta/2.0)**3)
Itr3 = I*rhog*(1-cosd(beta))/2.0
It3 =  Itb3+Itd3+Itr3
print ('modele isotropique', '%.2f'% (It1/1e6),' MJ/hr-m2')
print ('modele Hay davis','%.2f'% (It2/1e6),' MJ/hr-m2')
print ('modele HDKR', '%.2f'% (It3/1e6),' MJ/hr-m2')
#
# 4 modeles par appel de fonctions
#
It1b,Itb1,Itd1,Irt1 = modele_isotropique(I,Ib,Id,beta,Rb,rhog)
It2b,Itb2,Itd2,Itr2 = modele_hay_davis(I,Ib,Id,beta,Rb,rhog,Io)
It3b,Itb3,Itd3,Itr3 = modele_hdkr(I,Ib,Id,beta,Rb,rhog,Io)
print ('modele isotropique', '%.2f'% (It1b/1e6),' MJ/hr-m2')
print ('modele Hay davis','%.2f'% (It2b/1e6),' MJ/hr-m2')
print ('modele HDKR', '%.2f'% (It3b/1e6),' MJ/hr-m2')
# model Perez plus complexe
delt  = decl_solaire(n)
thez = zenith_solaire(phi,delt,omen)
Gon = irradiation_extraterrestre_normale(n)
Ion = Gon*3600.0
the = normale_solaire(delt,phi,omen,beta,gam)
It4,Ibt4,Idt4,Irt4 = modele_perez(I,Ib,Id,beta,Rb,rhog,Io,Ion,thez,the)
print ('modele Perez', '%.2f'% (It4/1e6),' MJ/hr-m2')
