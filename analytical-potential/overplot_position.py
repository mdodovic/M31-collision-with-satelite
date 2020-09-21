import math as math
import numpy as np
import matplotlib.pyplot as plt

def intenzitet(a,b,c):
    return (a*a + b*b + c*c)**0.5

#Linkovi

fileSveKomponente_in = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/galaxyAp_WithoutFriction/"
fileHalo_in = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/haloAP_BezTrenja/"

t = np.loadtxt(fileHalo_in + "vreme.txt", unpack=True)
t = np.asarray(t)
xS,yS,zS = np.loadtxt(fileSveKomponente_in + "xyz.txt", unpack = True)
xH,yH,zH = np.loadtxt(fileHalo_in + "xyz.txt", unpack = True)

rS = [intenzitet(xS[i],yS[i],zS[i]) for i in range(len(xS))]
rH = [intenzitet(xH[i],yH[i],zH[i]) for i in range(len(xS))]
r = [(rS[i] - rH[i]) / rS[i] for i in range(len(rS))]


print(len(t))
print(len(rS))

plt.plot(t[:3000],np.abs(r),c='black')
#plt.plot(t[:3000],rH)#,c='black')
#plt.plot(xH,yH)#,c='black')
#plt.plot(xH,yH)#,c='black')
plt.xlabel("$T[Gyr]$")
plt.ylabel("$ \delta R $" + "%")
#plt.xlim(70,-70)
#plt.ylim(-70,70)
plt.savefig("odnosRadijusa.png",dpi = 300)
plt.savefig("C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/odnosRadijusa.png", dpi = 90)
plt.savefig("C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika1.png", dpi = 90)
plt.savefig("C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika1.eps", dpi = 90)

plt.show()



"""
x = []
y = []
z = []

print(xS[0])
print(xH[0])
print(yS[0])
print(yH[0])

for i in range(int(len(xH)/2)):
    x.append(((xS[i] - xH[i])))
    y.append(((yS[i] - yH[i])))
    z.append(((zS[i] - zH[i])))


plt.plot(x,y)#,c='black')

#plt.plot(xH,yH)#,c='black')
#plt.plot(xH,yH)#,c='black')

plt.xlabel("$dX[AJ]$")
plt.ylabel("$dY[AJ]$")
#plt.xlim(70,-70)
#plt.ylim(-70,70)
plt.savefig("xy_ravan",dpi = 300)
plt.show()
"""

