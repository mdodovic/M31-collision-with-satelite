import numpy as np
import matplotlib.pyplot as plt

# Constatns
G = 4.515353945657138e-39 # Gravitaciona konstanta [kpc^3/(Msol*s^2)]
kriticnaGustina = 143.5012 #[Msol/kpc^3] - izracunat


def intenzitet(a,b,c):
    """
    function calculate intensity of three-component vector

    Input parameters:
    a,b,c - coordinates of vector

    Output parameter
    intensity - intensity of vector
    """
    intensity = (a*a + b*b + c*c)**0.5
    return intensity

def kineticEnergy(pathToFile):
    """
    Function calculate kinetic energy of a particle 
    Formula: Ek = 1/2mv^2 in units 

    Input parameters
    pathToFile - path to snapshot of single particle
                 contain mass [Mass_Of_Sun]; position (x,y,z) [kpc]; velocity (vx,vy,vz) [km/s]

    Output parameter
    ek - kinetic energy [Mass_Of_Sun kpc^2 s^-2]
    """
    m, x, y, z, vx, vy, vz = np.loadtxt(pathToFile, unpack=True)
    vx = vx*3.24078e-17
    vy = vy*3.24078e-17
    vz = vz*3.24078e-17
    ek = 0.5 * m * (intenzitet(vx, vy, vz)**2)
    return ek

def potentialEnergy(pathToFileParticle, pathToFileSystem):
    """
    Function calculate potential energy of a particle with given system of particle
    Formula is: Ep = - G * m1 * m2 / r  in units [Mass_Of_Sun km^2 s^-2]

    Input parameters
    pathToFileParticle - path to snapshot of single particle
                         contain mass [Mass_Of_Sun]; position (x,y,z) [kpc]; velocity (vx,vy,vz) [km/s]

    pathToFIleSystem - path to snapshot of system of particles
                    potential energy of Particle (pathToFIleParticle) is calculated to this system

    Output parameter
    ep - potential energy
    """

    mP, xP, yP, zP, vxP, vyP, vzP = np.loadtxt(pathToFileParticle, unpack=True)
    mS, xS, yS, zS, vxS, vyS, vzS = np.loadtxt(pathToFileSystem, unpack=True)
    n = len(mS)
    ep = 0.
    for i in range(n):
        ep += -G * mP * mS[i] / (intenzitet(xP - xS[i], yP - yS[i], zP - zS[i]))
    return ep




scale_length_halo = 8.18 # [kpc]
delta_c = 429231.966735 # bezdimenzioni parametar - izracunat
potencijalHaloConst = -4*np.pi*G*delta_c*kriticnaGustina*(scale_length_halo**3)

def potencijalH(r):
    FI = potencijalHaloConst/r*np.log((r+scale_length_halo)/scale_length_halo)
    return FI




N = 61
poc = 0
kr = N + 1
#############
folderName = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/numericalPotential/2Case/txtFiles/"
folderOut = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/numericalPotential/2Case/results/"
starsName = "stars_"
diskName = "disk_"
bulgeName = "bulge_"
haloName = "halo_"

energyHalo = []
#energyGalaxy = []
time = []

haloPotentialEnergyConstant = 1.548668791026037
for k in range(61):
    print(k/61*100)
    deo = str(k).zfill(3)
    starsCard = folderName + starsName + deo + ".txt"
    diskCard = folderName + diskName + deo + ".txt"
    bulgeCard = folderName + bulgeName + deo + ".txt"
    haloCard = folderName + haloName + deo + ".txt"

    mS, xS, yS, zS, vxS, vyS, vzS = np.loadtxt(starsCard, unpack=True)
    energyHalo.append(kineticEnergy(starsCard) + potentialEnergy(starsCard, haloCard)* haloPotentialEnergyConstant )
    #energyGalaxy.append(kineticEnergy(starsCard) + potentialEnergy(starsCard, haloCard) + potentialEnergy(starsCard, bulgeCard) + potentialEnergy(starsCard, diskCard))
    time.append(k / 10 * 0.5)
file = open(folderOut + "energy_halo.txt","w")
for i in range(len(time)):
    file.write(str(time[i]) + " " + str(energyHalo[i]) + "\n")
file.close

plt.plot(time, energyHalo)
plt.xlabel("$T$ [Gyr]")
plt.ylabel("$E$ [M$_{sol}$kpc$^2$/s$^2$]")
plt.savefig(folderOut + "energy_halo.png",dpi = 90)
plt.show()
"""
file = open(folderOut + "energy_galaxy.txt","w")
for i in range(len(time)):
    file.write(str(time[i]) + " " + str(energyGalaxy[i]) + "\n")
file.close

plt.plot(time, energyGalaxy)
plt.xlabel("$T$ [Gyr]")
plt.ylabel("$E$ [M$_{sol}$kpc$^2$/s$^2$]")
plt.savefig(folderOut + "energy_galaxy.png",dpi = 90)
plt.show()
"""