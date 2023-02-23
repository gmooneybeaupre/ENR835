#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 18:02:24 2023

@author: gmooneybeaupre
"""

# DONNÉES DU PROBLÈME
modele = 'iso'
g = 9.8
sig = 5.67e-8
Kd = 0.86
beta = 25
gamma = 0
It = 926
Id = 89
Ir = 39
tau_alpha_n = 0.88
Lair = 30               # SOURCE D'ERREUR, À VÉRIF. (en mm)
Tplaquec = 52
Tplaquek = Tplaquec + 273
hext = 12
Tskyc = 18
Tskyk = Tskyc + 273
Tambc = Tskyc
Tambk = Tskyk 
e1 = 0.25
e2 = 0.88

# COMPUTATION DES SOLUTIONS
# a) Trouver la radiation absorbée S
Kb = 
tau_alpha_d = Kd*tau_alpha_n
tau_alpha_r = Kd*tau_alpha_n

