#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:41:43 2023

@author: gmooneybeaupre
"""
from sys import *
from solar_mod import *
import numpy as np

# DONNÉES PRÉLIMINAIRES
j_type = np.array([17,47,75,105,135,162,198,228,258,288,318,344])  # jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
nom = ['jan','fev','mars','avr','mai','juin','juil','aout','sept','oct','nov','dec']
Gsc = 1367.0
phi = 45.5  # latitude  Montreal
lon = -73.0 -35.0/60.0  # longitude  Montreal
Lst = -75        # longitude  méridien HNE -5
n =jour_mois_jour_annee(6,'juin') # Le 6 juin

# IMPORTATION DU FICHIER MÉTEO
donnees = np.loadtxt("Meteo_horaire_Montreal.txt",skiprows = 1)   # donnes climatiques horaire Montreal  % valeurs en kJ/hr/m2
temps = donnees[:,0]   # heures
Tair = donnees[:,1]
I =  donnees[:,3]/1000   # valeurs de I total en MJoules/m2/ hr
Ib =  donnees[:,4]/1000   # valeurs de I direct  en MJoules/m2/ hr
Id = donnees[:,5]/1000   # valeurs de I diffus en MJoules/m2/ hr

# PARTIE 1 : TROUVER It
print("QUESTION 3")
print("PARTIE 1 : TROUVER It POUR CHAQUE TRANCHES HORAIRES ENSOLEILLÉES LE 6 JUIN")
beta = 25
gam = 0
rhog = 0.3

hr = range((n-1)*24,n*24)
h = 0

for i in hr:
    ome1 = 15*(h - 12)
    ome2 = 15*((h+1) - 12)
    omen = (ome1+ome2)/2.0

    Rb = Calcul_Rb(phi,n,omen,beta,gam)
    Io = irradiation_extraterrestre_horaire(n,phi,ome1,ome2)
    It,Ibt,Idt,Irt = modele_hdkr(I[i],Ib[i],Id[i],beta,Rb,rhog,Io)
    if It > 0.01:
        print("Le 6 juin de " + str(h) + "h00 à " + str(h+1) + "h00 (Montréal) \t:\t" + str(It) + " MJ/hr-m²")
    h+= 1

# PARTIE 2 : AVEC LE CAPTEUR DE LA QUESTION 2
print("PARTIE 2 : AVEC LE CAPTEUR DE LA QUESTION 2")
