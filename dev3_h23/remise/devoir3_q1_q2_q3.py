#coding: utf-8

from solar_mod import *
from properties_mod import *
import numpy as np

# QUESTION 1
# DONNÉES PRÉLIMINAIRES
g = 9.8
sig = 5.67e-8
Kd = 0.86
beta = 25.0
gamma = 0
It = 926
Itd = 89
Itr = 39
Itb = It - Itd - Itr
tau_alpha_n = 0.88
e1 = 0.25           # emissivité de la plaque
e2 = 0.88           # emissivité de la vitre
Lair  = 0.030       # epaisseur d'air (metres)
T1c = 52.0           # Temperature de la plaque (Celsius)
T1k = T1c+273        # Temperature de la plaque (Kelvin)
Tinfc = 18
Tinfk = Tinfc + 273
Tskyc = 18
Tskyk = Tskyc + 273
hconve = 12 #W/(m².K)
#
#
# a) Radiation absorbée S
tau_alpha_d = Kd*tau_alpha_n
tau_alpha_r = Kd*tau_alpha_n
Kb =  1 #car angle incidence nul
tau_alpha_b = Kb*tau_alpha_n
S = tau_alpha_b*Itb + tau_alpha_d*Itd + tau_alpha_r*Itr
print ("Question 1a)")
print("La radiation absorbée est de :", S, "W/m²")
#
#
#b) Pertes par unité de surface vers le haut
ok = False
compt = 1
compt_max = 30
delta = 0.01
T2c = 35    #hypothese
while not ok :
    Dt = T1c-T2c
    T2k = T2c + 273.0
    Tairk = (T1k+T2k)/2
    # Calcul du coef de convection interne
    betas = 1/Tairk
    nui = air_prop('nu',Tairk)
    ali = air_prop('al',Tairk)
    ki = air_prop('k',Tairk)
    Ra = g*betas*abs(Dt)*Lair**3/(nui*ali)
    Ra_c = 1708.0/cosd(beta)
    ct = cosd(beta)
    f1 = max(0,1.0-1708.0/(Ra*ct))
    f2 = 1.0-1708*(sind(1.8*beta))**1.6/(Ra*ct)
    f3 = max(0,(Ra*ct/5830)**(1.0/3.0)-1.0)
    Nui = 1.0 + 1.44*f1*f2+f3
    hconvi = Nui*ki/Lair
    # Calcul du coefficient de radiation interne
    hradi = sig*(T1k+T2k)*(T1k**2+T2k**2)/(1.0/e1+1.0/e2-1.0)
    # Résistance équivalente interne
    Ripp = 1.0/(hconvi+hradi)    
    # Calcul des pertes internes
    #qppi = (T1c - T2c)/Ripp
    qconvi = hconvi*(T1c-T2c)
    qradi = hradi*(T1c-T2c)
    qppi2 = qconvi+qradi
    #print ('q" int radiation ',qradi,' W/m2')
    #print ('q" int convection ',qconvi,' W/m2')
    #print ('q" int totales ',qppi,' W/m2'
    #print ('q" int totales ',qppi2,' W/m2')
    #
    # pertes extérieures
    # Coefficient de convection externe donné
    Rcpp = 1/hconve
    # Calcul du coefficient de radiation externe
    T3c = T2c   #car on néglige la résistance de conduction dans le verre
    T3k = T3c + 273
    hrade = e2*sig*(T3k**2+Tskyk**2)*(T3k+Tskyk)
    Rrpp = 1/hrade
    # Calcul des pertes externes
    qconve = hconve*(T3c-Tinfc)
    qrade = hrade*(T3c-Tskyc)
    qppe = qconve+qrade
    #print ('q" ext radiation ',qrade,' W/m2')
    #print ('q" ext convection ',qconve,' W/m2')
    #print ('q" ext totales ',qppe,' W/m2')
    #
    # Bilan à la surface externe de la vitre
    Rintpp = Ripp
    C1 = 1/Rintpp +  1/Rcpp + 1/Rrpp
    T3n = (T1c/Rintpp + Tinfc/Rcpp + Tskyc/Rrpp)/C1
    # vérification des bilans thermiques avec nouvelle valeur de T3
    qconve = hconve*(T3n-Tinfc)
    qrade = hrade*(T3n-Tskyc)
    qppen = qconve+qrade
    qppin = (T1c - T2c)/Ripp
    #
    diff = abs(T3n - T3c)
    if diff < delta:
        ok = True
    else:
        T2c = T3n
        compt = compt + 1
        if compt > compt_max:
            ok = True
            erreur = 1
Uhaut = qppin/(T1c-Tinfc)
#print ('T vitre ',T2c,T3c,T3n)
#print ('U haut ',Uhaut,' W/m2-k')
#print ('q" int correct ',qppin,' W/m2')
print ("Question 1b)")
print("Les pertes par unité de surface vers le haut sont de :", qppen,"W/m²")
#U_Klein =  U_Klein(T1c,Tinfc,beta,hconve,e1,e2,1)
#print ('U Klein ',U_Klein,' W/m2-k')
#
#
# c) Rendement du capteur
eta = (It-qppen)/It
print ("Question 1c)")
print("Le rendement du capteur est de :", eta*100, "%")


# QUESTION 2
# DONNÉES PRÉLIMINAIRES
#
#Y = 2
#H = 0.8
Ac = 0.8*2
W = 0.100          # distance entre les tubes
D = 0.007
Cb = 100 
deltaa = 0.001
ka = 400
#N = int(round(H/W))
N = 8
mp = 0.032
mp1 = mp/N
hf = 1100 
Cp = 4180
#
#
# a) F, F' et Fr
UL = Uhaut   #car on néglige les pertes du bas et des côtés
# Calcul du rendement d'ailette F
m = np.sqrt(UL/(ka*deltaa))
x = m*(W-D)/2
F = np.tanh(x)/x
# Calcul du rendement d'absorbeur F'
Rpjoint = 1.0/Cb
Rpconv = 1.0/(hf*pi*D)
den = W*(1.0/(UL*(D+(W-D)*F))+Rpjoint+Rpconv)
Fp = 1.0/(UL*den)
# Calcul du facteur de récupération Fr
X = UL*Ac*Fp/(mp*Cp)
Fpp = 1/X*(1-np.exp(-X))
Fr = Fp*Fpp
print ("Question 2a)")
print ('F = ',F)
print ("F' = ",Fp)
print ('Fr = ',Fr)
#
#
# b)  Températures d’entrée et de sortie du fluide
Tpc = T1c
qupp = (S - UL*(Tpc-Tinfc))
qu = qupp*Ac
Tfi = -(qupp/Fr-S)/UL + Tinfc
Tfo = Tfi + qu/(mp*Cp)
print ("Question 2b)")
print ("La température d'entrée est : ", Tfi, "°C")
print ("La température de sortie est : ", Tfo, "°C")


# QUESTION 3
# DONNÉES PRÉLIMINAIRES
Gsc = 1367.0
phi = 45.5  # latitude  Montreal
lon = -73.0 -35.0/60.0  # longitude  Montreal
Lst = -75        # longitude  méridien HNE -5
n =jour_mois_jour_annee(6,'juin') # Le 6 juin

# IMPORTATION DU FICHIER MÉTEO
donnees = np.loadtxt("Meteo_horaire_Montreal.txt", skiprows = 1)   # donnes climatiques horaire Montreal  % valeurs en kJ/hr/m2
temps = donnees[:,0]   # heures
Tair = donnees[:,1]
I =  donnees[:,3]/3.6   # valeurs de I total en W/m²
Ib =  donnees[:,4]/3.6   # valeurs de I direct  W/m²
Id = donnees[:,5]/3.6   # valeurs de I diffus en W/m²

print("Question 3")
beta = 25
gam = 0
rhog = 0.3
Kd = 0.86
tau_alpha_n = 0.88
b0 = 0.18
b1 = 0.08
delt = decl_solaire(n)

hr = range((n-1)*24,n*24)
h = 0
eta_glob = 0
tr_ens = 0
for i in hr:
    # Calculs angle horaire et angle d'incidence
    ome1 = 15*(h - 12)
    ome2 = 15*((h+1) - 12)
    omen = (ome1+ome2)/2.0
    thetab = normale_solaire(delt,phi,omen,beta,gam)
    
    # Calcul des radiations
    Rb = Calcul_Rb(phi,n,omen,beta,gam)
    Io = irradiation_extraterrestre_horaire(n,phi,ome1,ome2)
    It,Ibt,Idt,Irt = modele_hdkr(I[i],Ib[i],Id[i],beta,Rb,rhog,Io)
    It = It
    Ibt = Ibt
    Idt = Idt
    Irt = Irt
    
    # Calcul de la radiation absorbée
    Kb = Calcul_Ka(It,Ibt,Idt,Irt,thetab,b0,beta,b1)
    tau_alpha_b = Kb*tau_alpha_n
    S = tau_alpha_b*Ibt+tau_alpha_d*Idt+tau_alpha_r*Irt
    
    # Calcul de la chaleur absorbée
    Tinf = Tair[i]
    qupp = Fr*(S-UL*(Tfi-Tinf))
    
    # Calcul du rendement
    eta = qupp/It
    
    if It > 0.01:
        if eta > 0:
            eta_glob += eta
        else:
            eta_glob += 0
        tr_ens += 1
        print("Le 6 juin de " + str(h) + "h00 à " + str(h+1) + "h00 (Montréal) \t:")
        print("\t Température ambiante\t:\t" + str(Tinf) + "\t °C")
        print("\t Radiation solaire\t:\t" + str(It) + "\t W/m²")
        print("\t Radiation absorbée\t:\t" + str(S) + "\t W/m²")
        print("\t Chaleur récupérée\t:\t" + str(qupp) + "\t W/m²")
        print("\t Rendement\t:\t" + str(eta))
    h+= 1
# Calcul du rendement global pour le 6 juin
eta_glob = eta_glob/tr_ens 
print("RENDEMENT GLOBAL (6 juin)\t:\t" + str(eta_glob))