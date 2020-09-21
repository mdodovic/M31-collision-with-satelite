# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 00:21:37 2018

@author: matij
"""
import numpy as np
def criticalDensity():
    H = 71.9 #[km/s/Mpc] Hubble constant(Bovin et al. 2017)
    H = H *1000/(3.086e22) # [1/s]
    #G = 6.67e-11 #[m^3/kg/s^2]
    G = 4.515353945657138e-39 # [kpc^3/(Msol*s^2)]
    return 3*H**2/8/np.pi/G

def concentrationParameter(virialRadius,scaleLength):
    return virialRadius/scaleLength

def deltaC(c):
    return 200/3*(c**3)/(np.log(1+c) - (c/(1+c)))

if __name__== "__main__":
    virialRadius = 200.2
    scaleLength = 8.18
    print(concentrationParameter(virialRadius, scaleLength))
    print(deltaC(concentrationParameter(virialRadius,scaleLength)))
    print(criticalDensity())
