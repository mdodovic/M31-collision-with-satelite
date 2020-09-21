# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 15:05:05 2018

@author: matij
"""
import math as math

kriticnaGustina = 143.5012 #[Msol/kpc^3]
scale_length_halo = 8.18 # [kpc]
delta_c = 429231.966735 # bezdimenzioni parametar - izracunat

def mass(r):
    M = 4*math.pi*delta_c*kriticnaGustina*(scale_length_halo**3)*(math.log((r+scale_length_halo)/scale_length_halo) - (r/(r+scale_length_halo)))
    return M

def virialMass(virialRadius):
    return mass(virialRadius)

if __name__ == '__main__':
    virialRadius = 200.2
    print(virialMass(virialRadius))
