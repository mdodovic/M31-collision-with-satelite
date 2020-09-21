import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as math

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

def kiloparsecToDegrees(xInKpcArray, yInKpcArray, zInKpcArray):
    """
    Function converts kiloparsec into degrees
    
    Input parameters:
    xInKpcArray, yInKpcArray, zInKpcArray - arrays of coordinates in kiloparsec

    Output parameters:
    zeta, eta - coordinates in degrees
    """

    zeta = []
    eta = []
    for i in range(len(xInKpcArray)):
        zeta.append((xInKpcArray[i]/(zInKpcArray[i]+784))*(180/np.pi))
        eta.append((yInKpcArray[i]/(zInKpcArray[i]+784))*(180/np.pi))

    return zeta, eta

N=61
poc=0
kr=N+1
#############
folderName = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/numericalPotential/2Case/txtFiles/"
folderOut = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/numericalPotential/2Case/results/"
starsName = "stars_"
xStream_degrees=[2.015, 1.745, 1.483, 1.226, 0.969, 0.717, 0.467, 0.219]
yStream_degrees=[-3.965, -3.525, -3.087, -2.653, -2.264, -1.768, -1.327,-0.886]
xStreamError_degrees=0.33
yStreamError_degrees=0.22

xStream_kpc = []
yStream_kpc = []
xStreamError_kpc = []
yStreamError_kpc = []
const = 784 * np.pi / 180
for i in range(len(xStream_degrees)):
    xStream_kpc.append(const * xStream_degrees[i])
    yStream_kpc.append(const * yStream_degrees[i])
    xStreamError_kpc.append(const * xStreamError_degrees)
    yStreamError_kpc.append(const * yStreamError_degrees)


x_kpc = []
y_kpc = []
z_kpc = []
x_degrees = []
y_degrees = []
radius = []
time = []
for k in range(61):

    if k<=9:
        deo = "00"
    if k<=99 and k>=10:
        deo = "0"
    if k >= 100:
        deo = ""

    starsKart = folderName + starsName + deo + str(k) + ".txt"
    mS, xS, yS, zS, vxS, xyS, vzS = np.loadtxt(starsKart, unpack=True)  # ovde ubacujes bulge
    zeta,eta = kiloparsecToDegrees([xS],[yS],[zS])
    radius.append(intenzitet(xS,yS,zS))
    time.append(k/10*0.5)
    x_kpc.append(xS)
    y_kpc.append(yS)
    z_kpc.append(zS)

    x_degrees.append(zeta[0])
    y_degrees.append(eta[0])

file = open(folderOut + "trajectory_kpc.txt","w")
for i in range(len(x_kpc)):
    file.write(str(x_kpc[i]) + " " + str(y_kpc[i]) + " " + str(z_kpc[i]) + "\n")
file.close

plt.plot(x_kpc, y_kpc)
plt.xlim(50, -120)
plt.ylim(-70,200)
plt.xlabel(r'$X$ [kpc]')
plt.ylabel(r'$Y$ [kpc]')
plt.errorbar(xStream_kpc, yStream_kpc, yStreamError_kpc, xStreamError_kpc, fmt='', color = 'black',linestyle="None")
plt.savefig(folderOut + "trajectory_kpc.png", dpi = 90)
plt.show()

file = open(folderOut + "trajectory_degree.txt","w")
for i in range(len(x_degrees)):
    file.write(str(x_degrees[i]) + " " + str(y_degrees[i]) + "\n")
file.close

plt.plot(x_degrees, y_degrees)
plt.xlim(5, -5)
plt.ylim(-5, 5)
plt.xlabel("x [\N{DEGREE SIGN}]")
plt.ylabel("y [\N{DEGREE SIGN}]")
plt.savefig(folderOut + "trajectorDegrees.png", dpi = 90)
plt.errorbar(xStream_degrees, yStream_degrees, yStreamError_degrees, xStreamError_degrees, fmt='', color = 'black',linestyle="None")
plt.show()

file = open(folderOut + "radius.txt","w")
for i in range(len(time)):
    file.write(str(time[i]) + " " + str(radius[i]) + "\n")
file.close

plt.plot(time,radius)
plt.xlabel(r'$T$ [Gyr]')
plt.ylabel(r'$R$ [kpc]')
plt.savefig(folderOut + "radius.png", dpi = 90)
plt.show()


m = []
for i in range(len(x_degrees)):
    m.append(0.504*(x_degrees[i])-0.864*(y_degrees[i]))
m = np.asarray(m)
radius = np.asarray(radius)
###############################
radius = radius + 10
##############################

file = open(folderOut + "dubinaN_tela.txt","w")
for i in range(len(x_degrees)):
    file.write(str(m[i]) + " " + str(radius[i]) + "\n")
file.close

xDubinaField=[ 1.32, 1.82, 2.38, 2.84, 3.35, 3.86, 4.37]
yDubinaField=[ 829-785, 836-785, 840-785, 855-785, 860-785, 877-785, 886-785]
xDubinaField=np.asarray(xDubinaField)
yDubinaField=np.asarray(yDubinaField)
#y2=y1-785
xerrDubina = 0.3
yerrDubina = 25

#xDubinaField=np.fliplr([xDubinaField])[0]
#yDubinaField=np.fliplr([yDubinaField])[0]


#m = np.fliplr([m])[0]
#radius = np.fliplr([radius])[0]

plt.plot(m,radius)
plt.errorbar(xDubinaField, yDubinaField, yerrDubina, xerrDubina, fmt='', color = 'black',linestyle="None")
#plt.title("Orijentacija toka")
#plt.xlim(0, 5)
#plt.ylim(-30,150)
plt.xlabel("m [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("R [kpc]", fontsize=15)
plt.savefig(folderOut + 'dubinaN_tela', dpi=300 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()


