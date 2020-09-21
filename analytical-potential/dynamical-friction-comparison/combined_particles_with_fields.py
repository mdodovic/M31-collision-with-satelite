# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 17:48:55 2018

@author: matij
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

#polja u radovima, jeinice su stepeni
x1=[2.015, 1.745, 1.483, 1.226, 0.969, 0.717, 0.467, 0.219]
y1=[-3.965, -3.525, -3.087, -2.653, -2.264, -1.768, -1.327,-0.886]
x1err=0.33
y1err=0.22

#koriscena polja, jedinice su kpc
x = []
y = []
xerr = []
yerr = []
const = 784 * np.pi / 180
for i in range(len(x1)):
    x.append(const * x1[i])
    y.append(const * y1[i])
    xerr.append(const * x1err)
    yerr.append(const * y1err)


fileNonDinamicalFriction = "C:\\Users\\matij\\Desktop\\Projekat_2018_dinamicko_trenje\\analiticki_potencijal\\JuzniTok\\haloAP_BezTrenja\\"
fileDinamicalFriction = "C:\\Users\\matij\\Desktop\\Projekat_2018_dinamicko_trenje\\analiticki_potencijal\\JuzniTok\\haloAP_SaTrenjem\\"
outputImages = "C:\\Users\\matij\\Desktop\\Projekat_2018_dinamicko_trenje\\analiticki_potencijal\\JuzniTok\\combinedCases\\haloAP\\"
streamField = "inDegrees.txt"
radius = "radijus.txt"
energy = "energija.txt"
angularMomentum = "momentImpulsa.txt"
trajectory = "xyz.txt"


if __name__ == "__main__":

    nonFriction = mlines.Line2D([], [], color='blue', markersize=15, label='bez_trenja', linestyle = '-')
    Friction = mlines.Line2D([], [], color='red', markersize=15, label='sa_trenjem',linestyle = '-')

    xNDF, yNDF, zNDF = np.loadtxt(fileNonDinamicalFriction + trajectory, unpack = True)
    xDF, yDF, zDF = np.loadtxt(fileDinamicalFriction + trajectory, unpack = True)
    plt.plot(xNDF,yNDF, c = 'blue')
    plt.plot(xDF,yDF, c = 'red')
    plt.xlim(-200, 200)
    plt.ylim(-200, 200)
    plt.xlabel("X [kpc]") #, fontsize=15)
    plt.ylabel("Y [kpc]") #, fontsize=15)
    plt.legend(handles=[nonFriction,Friction], loc=7)
    plt.savefig(outputImages + "trajectory.png", dpi = 300)
    plt.show()     

    xNDF, yNDF, zNDF = np.loadtxt(fileNonDinamicalFriction + trajectory, unpack = True)
    xDF, yDF, zDF = np.loadtxt(fileDinamicalFriction + trajectory, unpack = True)
    plt.plot(xNDF,yNDF, c = 'blue')
    plt.plot(xDF,yDF, c = 'red')
    plt.xlim(68.417, -6.8417)
    plt.ylim(-68.417, -6.8417)
    plt.xlabel("X [kpc]") #, fontsize=15)
    plt.ylabel("Y [kpc]") #, fontsize=15)
    plt.legend(handles=[nonFriction,Friction], loc=2)
    plt.errorbar(x, y, yerr,xerr, fmt='', color = 'black',linestyle="None")
    plt.savefig(outputImages + "streamField.png", dpi = 300)
    plt.show() 
    
    tNDF,eNDF = np.loadtxt(fileNonDinamicalFriction + energy, unpack = True)
    tDF,eDF = np.loadtxt(fileDinamicalFriction + energy, unpack = True)
    plt.plot(tNDF,eNDF, c = 'blue')
    plt.plot(tDF,eDF, c = 'red')
    plt.legend(handles=[nonFriction,Friction], loc=4)
    plt.xlabel(r"$T[Gyr]$")
    plt.ylabel(r"$E[M_\odot\ kpc^2 / s]$")    
    plt.savefig(outputImages + "energy.png", dpi = 300)
    plt.show()    
    
    tNDF,lNDF = np.loadtxt(fileNonDinamicalFriction + angularMomentum, unpack = True)
    tDF,lDF = np.loadtxt(fileDinamicalFriction + angularMomentum, unpack = True)
    plt.plot(tNDF, lNDF, c = 'blue')
    plt.plot(tDF, lDF, c = 'red')
    plt.legend(handles=[nonFriction,Friction], loc=4)
    plt.xlabel(r"$T[Gyr]$")
    plt.ylabel(r"$L[M_\odot\ kpc^2/s]$")
    plt.savefig(outputImages + "angularMomentum.png", dpi = 300)
    plt.show()    
        
    tNDF, rNDF = np.loadtxt(fileNonDinamicalFriction + radius, unpack = True)
    tDF, rDF = np.loadtxt(fileDinamicalFriction + radius, unpack = True)
    plt.plot(tNDF, rNDF, c = 'blue')
    plt.plot(tDF, rDF, c = 'red')
    plt.xlabel(r"$T[Gyr]$")
    plt.ylabel(r"$R[kpc]$")
    plt.legend(handles=[nonFriction,Friction], loc=4)
    plt.savefig(outputImages + "radius.png", dpi = 300)
    plt.show()  
