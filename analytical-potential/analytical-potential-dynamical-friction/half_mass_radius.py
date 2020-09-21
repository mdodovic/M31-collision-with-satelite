# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 14:16:11 2018

@author: matij
"""

import numpy as np
import math as math

G = 4.515353945657138e-39 # [kpc^3/(Msol*s^2)]
kriticnaGustina = 277.72/(0.71*0.71) # [Msol/kpc^2]
delta_c = 27e4 # bezdimenzioni parametar
scale_length_halo = 8.18 # [kpc]
scale_length_halo_3 = scale_length_halo**3 # [kpc]


def masaH(r):
    M = 4*math.pi*delta_c*kriticnaGustina*(scale_length_halo_3)*(math.log((r+scale_length_halo)/scale_length_halo) - (r/(r+scale_length_halo)))
    return M


def intenzitet(a,b,c):
    return (a*a + b*b + c*c)**0.5

def halfMassRadiusAnalyticalPotential(Mhalo):
    halfMass = Mhalo/2
    r = 0.
    ukMass = 0
    while r <= 250 :
        if ukMass >= halfMass:
            return r
        
        #print(ukMass)
        #print(sum(massArray[r>radArray]))
        #print(r)
        ukMass = masaH(r)
        print(ukMass)
        r += 0.0001

    #print(halfMass)

def halfMassRadiusN_Body(mass, massArray ,radArray):
    halfMass = mass/2
    #print(halfMass)
    r = 0.
    ukMass = 0
    while r <= 200 :
        if ukMass >= halfMass:
            return r
        
        #print(ukMass)
        #print(sum(massArray[r>radArray]))
        #print(r)
        ukMass = sum(massArray[r>radArray])
        #print(ukMass)
        r += 0.0001

def loadgalaxy(path):
    m,x,y,z,vx,vy,vz = np.loadtxt(path,unpack = True)
    return m,x,y,z,vx,vy,vz

if __name__== "__main__":
    pathFolder = 'C:/Users/matij/Desktop/M31_izolacija/pocetni_uslovi/' 
    pathFileName = 'galaksija.txt'
    m,x,y,z,vx,vy,vz = loadgalaxy(pathFolder + pathFileName)
    radArray = [intenzitet(x[i],y[i],z[i]) for i in range(len(m))]
    radXYArray = [intenzitet(x[i],y[i],0) for i in range(len(m))]
    radArray = np.asarray(radArray)
    m = np.asarray(m)
    z = np.asarray(z)
    radXYArray = np.asarray(radXYArray)
    
    massDisk = m[:108928]
    massBulge = m[108929:108929+96246]
    massHalo = m[108929+96247:]
    
    radDisk = radArray[:108928]
    verticalDisk = z[:108928]
    radialalDisk = radXYArray[:108928]
    radBulge = radArray[108929:108929+96246]
    radHalo = radArray[108929+96247:]

    Mhalo = 88.4591e10 # [Msol]
    print("Half mass radius haloa: " + str(halfMassRadiusAnalyticalPotential(Mhalo)))
    
    #print(max(z[:108928]))
    #print(sum(m[:108928]))
    #print(sum(m[:108928])/2)
    
    #print(verticalDisk == z[:108928])
    
    #print(sum(m[:108928]))#disk
    #print(sum(m[108929:108929+96246])) #disk
    #print(sum(m[108929+96247:])) #halo
    
    #print("Half mass radius diska: " + str(halfMassRadius(sum(massDisk), massDisk, radDisk)))
    #print("Half mass radius bulgea: " + str(halfMassRadius(sum(massBulge), massBulge, radBulge)))
    #print("Half mass radius haloa: " + str(halfMassRadius(sum(massHalo), massHalo, radHalo)))
    #print("Half mass radius diska - vertikalni: " + str(halfMassRadiusN_Body(sum(massDisk), massDisk, verticalDisk)))
    #print("Half mass radius diska - radijalni: " + str(halfMassRadius(sum(massDisk), massDisk, radialalDisk)))
    #print("Half mass radius diska: " + str(halfMassRadius(sum(massDisk), massDisk, radDisk)))
    #r_half_radijalni_deo_diska = 10.681499999988445
    #r_half_ceo_disk = 10.694599999988414
    #x = [5 , 4, 6, 2, 8]
    #a = [10, 9, 8, 7, 6]
    #x = np.asarray(x)
    #a = np.asarray(a)
    #ind = x >= 5
    #b = a[ind]
    #print(b)
    
    #for i in range(len(m)):
    #    if intenzitet(x[i],y[i],z[i]) <=rVir:
    #        M += m[i]*2.325e9
        
    #print(M)
    
    #N=467081      # total number of particles
    #Nd=108929      # disk particles
    #Nb=96247     # bulge particles
    #Nh=261905      # halo particles
    
