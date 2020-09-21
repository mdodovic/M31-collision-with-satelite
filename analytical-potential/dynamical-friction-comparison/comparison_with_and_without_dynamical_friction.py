# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 17:48:55 2018

@author: matij
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def rotation(x,y,z):
    """
    rotates vector (3 dimension) for angle alpha = 77 and betta = 37 degrees 
    """
    alpha = (np.pi/180)*(77)
    betta  = (np.pi/180)*(37) 
        
    x0 = x
    y0 = y
    z0 = z
    
    x = x0*np.cos(betta) - y0*np.cos(alpha)*np.sin(betta) + z0*np.sin(alpha)*np.sin(betta)
    y = x0*np.sin(betta) + y0*np.cos(alpha)*np.cos(betta) - z0*np.sin(alpha)*np.cos(betta)
    z = y0*np.sin(alpha) + z0*np.cos(alpha)

    return x,y,z

"""
x = -84.41 #kpc
y = 152.47
z = -97.08

x1, y1, z1 = rotation(x, y, z)

print(x1,y1)

zeta0 = (x1/(z1+784))*(180/np.pi)
eta0 = (y1/(z1+784))*(180/np.pi)

print(zeta0, eta0)
"""
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


fileNonDinamicalFriction = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/dwarvesWithDisk/"
fileDinamicalFriction = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/dwarvesWithDisk/"
outputImages = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/combinedCases/testParticleWithDisk/"
#streamField = "inDegrees.txt"
radiusAP = "radijus.txt"
raidusNP = "radius.txt"
energyAP = "energija.txt"
energyNP = "energy_halo.txt"
trajectoryAP = "xyz.txt"
trajectoryNP = "trajectory_kpc.txt"
orjentationAP = "trajectoryDegree.txt"

orjentationNP = "trajectory_degree.txt"

dubinaAP = "dubinaAaliticki.txt"
dubinaNP = "dubinaN_tela.txt"
#angularMomentum = "momentImpulsa.txt"

if __name__ == "__main__":

    nonFriction = mlines.Line2D([], [], color='green', markersize=15, label='bez_trenja', linestyle = '-')
    Friction = mlines.Line2D([], [], color='blue', markersize=15, label='sa_trenjem',linestyle = '-')
    numerical = mlines.Line2D([], [], color='red', markersize=15, label='numericki',linestyle = '-')
    
    #Energy
    tNum,eNum = np.loadtxt(fileNumericalPotential + energyNP, unpack = True)
    tNDF,eNDF = np.loadtxt(fileNonDinamicalFriction + energyAP, unpack = True)
    tDF,eDF = np.loadtxt(fileDinamicalFriction + energyAP, unpack = True)
        
    plt.plot(tNum, eNum, c = 'red')    
    plt.plot(tNDF, eNDF, c = 'green')
    plt.plot(tDF, eDF, c = 'blue')
    plt.legend(handles=[nonFriction,Friction,numerical], loc='lower left', prop={'size': 13})
    plt.xlabel(r"$T$ [Gyr]")
    plt.ylabel(r"$E [M_\odot\ kpc^2 / s^2]$")    
    plt.savefig(outputImages + "energy.png", dpi = 300)
    plt.show()    

    #Radius
    tNum, rNum = np.loadtxt(fileNumericalPotential + raidusNP, unpack = True)
    tNDF, rNDF = np.loadtxt(fileNonDinamicalFriction + radiusAP, unpack = True)
    tDF, rDF = np.loadtxt(fileDinamicalFriction + radiusAP, unpack = True)

    plt.plot(tNum, rNum, c = 'red')
    plt.plot(tNDF, rNDF, c = 'green')
    plt.plot(tDF, rDF, c = 'blue')
    plt.xlabel(r"$T$ [Gyr]")
    plt.ylabel(r"$R$ [kpc]")
    plt.legend(handles=[nonFriction,Friction, numerical], loc='lower left', prop={'size': 13})
    plt.savefig(outputImages + "radius.png", dpi = 300)
    plt.show()  

    #Trajectory in kpc
    xNum, yNum, zNum = np.loadtxt(fileNumericalPotential + trajectoryNP, unpack = True)
    #xNum, yNum, zNum = rotation(xNum, yNum, zNum)
    xNDF, yNDF, zNDF = np.loadtxt(fileNonDinamicalFriction + trajectoryAP, unpack = True)
    #xNDF, yNDF, zNDF = rotation(xNDF, yNDF, zNDF)
    xDF, yDF, zDF = np.loadtxt(fileDinamicalFriction + trajectoryAP, unpack = True)
    #xDF, yDF, zDF = rotation(xDF, yDF, zDF)

    plt.plot(xNum, yNum, c='red')
    plt.plot(xNDF, yNDF, c='green')
    plt.plot(xDF, yDF, c='blue')
    plt.xlim(100, -100)
    #plt.ylim(-70, 70)
    plt.xlabel("X [kpc]")  # , fontsize=15)
    plt.ylabel("Y [kpc]")  # , fontsize=15)
    plt.errorbar(x, y, yerr, xerr, fmt='', color='black', linestyle="None")
    plt.legend(handles=[nonFriction, Friction, numerical], loc=7)
    plt.savefig(outputImages + "trajectory.png", dpi=300)
    plt.show()

    plt.plot(xNum, yNum, c='red')
    plt.plot(xNDF, yNDF, c='green')
    plt.plot(xDF, yDF, c='blue')
    plt.xlim(45, -20)
    plt.ylim(-65, 15)
    plt.xlabel("X [kpc]")  # , fontsize=15)
    plt.ylabel("Y [kpc]")  # , fontsize=15)
    plt.errorbar(x, y, yerr, xerr, fmt='', color='black', linestyle="None")
    plt.legend(handles=[nonFriction, Friction, numerical], loc=7)
    plt.savefig(outputImages + "trajectoryZoom.png", dpi=300)
    plt.show()

    #Trajectory in degree
    etaNum, zetaNum = np.loadtxt(fileNumericalPotential + orjentationNP, unpack = True)
    etaNDF, zetaNDF = np.loadtxt(fileNonDinamicalFriction + orjentationAP, unpack = True)
    etaDF, zetaDF = np.loadtxt(fileDinamicalFriction + orjentationAP, unpack = True)
    
    plt.plot(etaNum, zetaNum, c='red')
    plt.plot(etaNDF, zetaNDF, c='green')
    plt.plot(etaDF, zetaDF, c='blue')
    plt.xlim(100, -100)
    #plt.ylim(-70, 70)
    plt.xlabel("X [kpc]")  # , fontsize=15)
    plt.ylabel("Y [kpc]")  # , fontsize=15)
    plt.errorbar(x, y, yerr, xerr, fmt='', color='black', linestyle="None")
    plt.legend(handles=[nonFriction, Friction, numerical], loc=7)
    plt.savefig(outputImages + "orjentation.png", dpi=300)
    plt.show()

    plt.plot(etaNum, zetaNum, c='red')
    plt.plot(etaNDF, zetaNDF, c='green')
    plt.plot(etaDF, zetaDF, c='blue')
    plt.xlim(45, -20)
    plt.ylim(-65, 15)
    plt.xlabel("X [kpc]")  # , fontsize=15)
    plt.ylabel("Y [kpc]")  # , fontsize=15)
    plt.errorbar(x, y, yerr, xerr, fmt='', color='black', linestyle="None")
    plt.legend(handles=[nonFriction, Friction, numerical], loc=7)
    plt.savefig(outputImages + "orientationZoom.png", dpi=300)
    plt.show()



    #Dubina 
    mNum, dNum = np.loadtxt(fileNumericalPotential + dubinaNP, unpack = True)
    mNDF, dNDF = np.loadtxt(fileNonDinamicalFriction + dubinaAP, unpack = True)
    mDF, dDF = np.loadtxt(fileDinamicalFriction + dubinaAP, unpack = True)

    xDubinaField=[ 1.32, 1.82, 2.38, 2.84, 3.35, 3.86, 4.37]
    yDubinaField=[ 829-785, 836-785, 840-785, 855-785, 860-785, 877-785, 886-785]
    xDubinaField=np.asarray(xDubinaField)
    yDubinaField=np.asarray(yDubinaField)
    #y2=y1-785
    xerrDubina = 0.3
    yerrDubina = 25


    plt.plot(mNDF, dNDF, c='green')
    plt.plot(mDF, dDF, c='blue')
    plt.plot(mNum, dNum, c='red')
    #plt.xlim(45, -20)
    #plt.ylim(-65, 15)
    plt.xlabel("m [\N{DEGREE SIGN}]")  # , fontsize=15)
    plt.ylabel("R [kpc]")  # , fontsize=15)
    plt.errorbar(xDubinaField, yDubinaField, yerrDubina, xerrDubina, fmt='', color = 'black',linestyle="None")
    plt.legend(handles=[nonFriction, Friction, numerical], loc='best')
    plt.savefig(outputImages + "dubina.png", dpi=300)
    plt.show()

    plt.plot(mDF, dDF, c='blue')
    plt.plot(mNDF, dNDF, c='green')
    plt.plot(mNum, dNum, c='red')
    plt.xlim(-1, 7)
    plt.ylim(0, 150)
    plt.xlabel("m [\N{DEGREE SIGN}]")  # , fontsize=15)
    plt.ylabel("R [kpc]")  # , fontsize=15)
    plt.errorbar(xDubinaField, yDubinaField, yerrDubina, xerrDubina, fmt='', color = 'black',linestyle="None")
    plt.legend(handles=[nonFriction, Friction, numerical], loc='best')
    plt.savefig(outputImages + "dubinaZoom.png", dpi=300)
    plt.show()

"""
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

    tNDF, rNDF = np.loadtxt(fileNonDinamicalFriction + radius, unpack = True)
    tDF, rDF = np.loadtxt(fileDinamicalFriction + radius, unpack = True)
    plt.plot(tNDF, rNDF, c = 'blue')
    plt.plot(tDF, rDF, c = 'red')
    plt.xlabel(r"$T[Gyr]$")
    plt.ylabel(r"$R[kpc]$")
    plt.legend(handles=[nonFriction,Friction], loc=4)
    plt.savefig(outputImages + "radius.png", dpi = 300)
    plt.show()  
"""