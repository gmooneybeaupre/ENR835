#coding: utf-8
#
# exemple solaire 5.1
#
# Exemple du calcul de la radiation absorbée par un capteur vitre (1 vitre)
#
#
from solar_mod import *
import numpy as np
from matplotlib.pyplot import *

modele = 'iso'
#
# données du problème
#
I = 2.2     # Mj/hr m2
Ib = 1.4
Id = I - Ib  #
rhog = 0.6
phi = 31.0
beta = 60.0      # inclinaison du capteur
gam = 0.0        # azimuth du capteur
n2 = 1.526     # indice de refraction du verre
alpn = 0.93   # coef abs solaire normal
KL = .037      # coef extinction de la vitre
#
# exemple le 11 février entre 11:00 hrs et midi
#
n = jour_mois_jour_annee(11,'fev')
delt  = decl_solaire(n)      # calcul de la declinaison solaire
ome1 = -15  # 11 hrs AM
ome2 = 0    # midi
omen = -7.5
the1 = normale_solaire(delt,phi,omen,beta,gam)   # derection de la radiation directe
thez = zenith_solaire(delt,phi,omen)     # zenith solaire

#
#
Rb = Calcul_Rb(phi,n,omen,beta,gam)
Rbn = cosd(the1)/cosd(thez)
print (' Rb =', Rb,Rbn)
# calcul de la radiation incidente
Itb = Ib*Rb
Itd = Id*(1 + cosd(beta))/2
Itr = I*rhog*(1-cosd(beta))/2
It = Itb + Itd + Itr
It2,Itb2,Itd2,Itr2 = modele_isotropique(I,Ib,Id,beta,Rb,rhog)
#
# Calcul de (Tau*alpha) radiation directe
#
# Calcul au long
#
# calcul de the2 (loi de snell)
#

n1 = 1  # air
the2 = arcsind(n1*sind(the1)/n2)
if the1 == 0:
    rpe = ((n1 - n2)/(n1+n2))**2
    rpa = ((n1 - n2)/(n1+n2))**2
else:

    th12 = the1+the2
    th2m1 = the2-the1
    #
    # Eq 5.4
    #

    rpe = sind(th2m1)**2/sind(th12)**2
    rpa = tand(th2m1)**2/tand(th12)**2
#
# Eq 5.9
#
tau_a = np.exp(-KL/cosd(the2))
rhopa = rpa + rpa*(tau_a*(1-rpa))**2/(1.0 - (rpa*tau_a)**2)
rhope = rpe + rpe*(tau_a*(1.0-rpe))**2/(1.0 - (rpe*tau_a)**2)
rhov = (rhopa + rhope)/2.0
taupa =  tau_a*(1-rpa)**2/(1.0 - (rpa*tau_a)**2)
taupe =  tau_a*(1-rpe)**2/(1.0 - (rpe*tau_a)**2)
tauv = (taupa + taupe)/2.0
alpv = 1.0 - tauv - rhov
#
# Rq 5.11 radiation directe La fonction alp_alpn du module solaire
#
al_aln = alp_alpn(the1)
alpc = alpn*al_aln
tau_alba = tauv*alpc/(1.0-(1.0-alpc)*rhov)
#
# Autre façon utilisant la fonction Calcul_tau_alpha du module solaire fourni
#
tau_alb = Calcul_tau_al(the1,alpn,KL,n2)   # fonction Calcul_tau_al: tau_al = Calcul_tau_al(the1,alpn,KL,n2,n1), si on ne donne pas n1, il prend n1 =1 (air ou vide)

#
#
# Calcul du (tau*alpha) pour la  radiation diffuse
# On ne le fait pas au long,
#
the1d = angle_diffus(beta)  # eq 4.7
tau_ald = Calcul_tau_al(the1d,alpn,KL,n2)
# Calcul du (tau*alpha) pour la
# radiation reflechie
the1r = angle_reflechi(beta)     # eq 4.8
tau_alr = Calcul_tau_al(the1r,alpn,KL,n2)
#
# caclul de la radiation absorbée modele isotrope
#
Sb1 = Itb*tau_alb
Sd1 = Itd*tau_ald
Sr1 = Itr*tau_alr
S1 = Sb1+Sd1+Sr1
# calcul du (tau*alpha) moyen
tau_al_moy = S1/It
print ('tau_al_moy = ',tau_al_moy)
