# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 19:57:48 2018

@author: matij
"""
import math as math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines


# Constants 
delta_c = 429231.966735 # bezdimenzioni parametar - izracunat
kriticnaGustina = 143.5012 #[Msol/kpc^3]
scale_length_halo = 8.18 # [kpc]
scale_length_halo_3 = scale_length_halo**3 #[kpc^3]

def intensity(a,b,c):
    return (a*a + b*b + c*c)**0.5

def masah(r):
    M = 4*math.pi*delta_c*kriticnaGustina*(scale_length_halo_3)*(math.log((r+scale_length_halo)/scale_length_halo) - (r/(r+scale_length_halo)))
    return M #Dodaj *10

def gustinaH(r):
    RHO = delta_c*kriticnaGustina / ( (r / scale_length_halo) * (1 + r / scale_length_halo) * (1 + r / scale_length_halo) )
    return RHO

    
def profileOfMassAnalitycal():
    mass = []
    radius = []
    r = 0.
    dr = 0.05
    radiusHalo = 300.    
    while r <= radiusHalo:
        dm = masah(r + dr) - masah(r)
        mass.append(dm)
        radius.append(r)
        r += dr
    
    return mass, radius        

def profileOfMassNbody(massNbody, radiusNbody):
    mass = []
    radius = []
    r = 0.
    dr = 0.5
    radiusHalo = max(radiusNbody)
    while r <= radiusHalo:
        dm = 0
        dm = sum(massNbody[radiusNbody <= (r + dr)]) - sum(massNbody[radiusNbody <= r])
        mass.append(dm)
        radius.append(r)
        r += dr
    
    return mass, radius        



def profileOfDensityAnalytical():
    density = []
    radius = []
    r = 0.05
    dr = 0.05
    radiusHalo = 10.    
    while r <= radiusHalo:
        dRho = gustinaH(r)
        density.append(dRho)
        radius.append(r)
        r += dr
    
    return density, radius        


def profileOfDensityN_body(massNbody, radiusNbody):
    density = []
    radius = []
    r = 0.
    dr = 0.05
    radiusHalo = 10.    
    while r <= radiusHalo:
        dRho = (sum(massNbody[radiusNbody <= (r + dr)]) - sum(massNbody[radiusNbody <= r]))/(4/3*np.pi*((r+dr)**3 - r**3)) 
        density.append(dRho)
        radius.append(r)
        r += dr
    
    return density, radius        

    
if __name__ == "__main__":

    m,x,y,z,vx,vy,vz = np.loadtxt("C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/galaxyModels/M31/galaksija.txt",unpack= True)
    mass = [m[i]*2.325e9 for i in range(len(m))]    
    radius = [intensity(x[i],y[i],z[i]) for i in range(len(m))]
    mass = np.asarray(mass)
    radius = np.asarray(radius)

    n_body_grey = mlines.Line2D([], [], color='black', markersize=15, label='N-tela',linestyle = '-')
    analitycal_grey = mlines.Line2D([], [], color='black', markersize=15, label='Analiticki',linestyle = '--')
    n_body_color = mlines.Line2D([], [], color='red', markersize=15, label='N-tela',linestyle = '-')
    analitycal_color = mlines.Line2D([], [], color='blue', markersize=15, label='Analiticki',linestyle = '-')
    massA, radiusA = profileOfMassAnalitycal()
    massN, radiusN = profileOfMassNbody(mass, radius)
    radiusA = np.asarray(radiusA)
    radiusN = np.asarray(radiusN)
    massA = np.asarray(massA)
    massN = np.asarray(massN)
#    """
    plt.plot(radiusA,massA)
    plt.xlabel('$R [kpc]$')
    plt.ylabel('$\Delta M [M_\odot]$')
    plt.savefig("profilMaseAnalitycal.png")
    plt.show()

    plt.plot(radiusN,massN)
    plt.xlabel('$R [kpc]$')
    plt.ylabel('$\Delta M [M_\odot]$')
    plt.savefig("profilMaseN_body.png")
    plt.show()
#    """
    k = max(massN) / max(massA)
    for i in range(len(massN)):
        massN[i] = massN[i] / k
        massN[i] = massN[i]*100

    for i in range(len(massA)):
        radiusA[i] = radiusA[i] - 5
        massA[i] = massA[i]*100

    plt.plot(radiusN, massN, linestyle='-', color='red')
    plt.plot(radiusA, massA, linestyle='-', color='blue')
#    plt.plot(radiusN, massN, linestyle='-', color='black')
#    plt.plot(radiusA, massA, linestyle='--', color='black')
    plt.xlabel('$R [kpc]$')
    plt.ylabel('$\Delta M [M_\odot]$')
    plt.xlim(0, 200)
#    plt.legend(handles=[n_body, analitycal], loc=7)
    plt.legend(handles=[n_body_color, analitycal_color], loc=7)
    plt.savefig("profilMaseKombinovani_boja.png")
    plt.show()


#    plt.plot(radiusN, massN, linestyle='-', color='red')
#    plt.plot(radiusA, massA, linestyle='-', color='blue')
    plt.plot(radiusN, massN, linestyle='-', color='black')
    plt.plot(radiusA, massA, linestyle='--', color='black')
    plt.xlabel('$R [kpc]$')
    plt.ylabel('$\Delta M [M_\odot]$')
    plt.xlim(0, 200)
    plt.legend(handles=[n_body_grey, analitycal_grey], loc=7)
#    plt.legend(handles=[n_body_color, analitycal_color], loc=7)
    plt.savefig("profilMaseKombinovani_crno.png")
    plt.savefig("C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/profilMaseKombinovani_crno.png", dpi = 90)
    plt.savefig("C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika2.png", dpi = 90)
    plt.savefig("C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika2.eps", dpi = 90)
    plt.show()

    #densityA, radiusA = profileOfDensityAnalytical()
    #densityN, radiusN = profileOfDensityN_body(mass, radius)


    """
    plt.plot(radiusA,densityA)
    plt.xlabel('$R [kpc]$')
    plt.ylabel(r'$\rho [M_\odot/kpc^3]$')
    plt.savefig("profilGustineAnalitycal.png")
    plt.show()

    plt.plot(radiusN,densityN)
    plt.xlabel('$R [kpc]$')
    plt.ylabel(r'$\rho [M_\odot/kpc^3]$')
    plt.savefig("profilGustineN_body.png")
    plt.show()
    """
    """
    plt.plot(radiusN, densityN, linestyle='-', color='black')
    plt.plot(radiusA, densityA, linestyle='--', color='black')
    plt.xlabel('$R [kpc]$')
    plt.ylabel('$\Delta M [M_\odot]$')
    plt.xlim(0, 5)
    plt.legend(handles=[n_body, analitycal], loc=7)
    plt.savefig("profilGustineKombinovani.png")
    plt.show()
    
    """
    """
    plt.plot(radiusN, densityN, linestyle='-', color='black')
    plt.plot(radiusA, densityA, linestyle='--', color='black')
    plt.xlabel('$R [kpc]$')
    plt.ylabel('$\Delta M [M_\odot]$')
    plt.xlim(0, 5)
    plt.legend(handles=[n_body, analitycal], loc=7)
    plt.savefig("profilGustineKombinovani.png")
    plt.show()
    
    """
    
