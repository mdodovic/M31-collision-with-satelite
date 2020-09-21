import math as math
import numpy as np
import matplotlib.pyplot as plt

#Linkovi

file_in ="C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiSaTrenjemSveKomponente/"

file_out ="C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiSaTrenjemSveKomponente/"

m,x,y,z,vx,vy,vz= np.loadtxt(file_in + "cestica.txt", unpack=True)
t = np.loadtxt(file_in + "vreme.txt", unpack = True)    


def intenzitet(a,b,c):
    return (a*a + b*b + c*c)**0.5

ubrzanje = []
vreme = []

for i in range(1,len(m)):
    v1 = intenzitet(vx[i-1],vy[i-1],vz[i-1])
    v2 = intenzitet(vx[i],vy[i],vz[i])
    t1 = t[i-1]
    t2 = t[i]
    ubrzanje.append((v2-v1)/(t2-t1))
    vreme.append(t1)
    
for i in range(len(ubrzanje)):
    ubrzanje[i] = ubrzanje[i]*1e13
    

file = open(file_out + "ubrzanje.txt","w")
for i in range(len(ubrzanje)):
    file.write(str(vreme[i]) + " " + str(ubrzanje[i]) + "\n")
file.close

plt.plot(vreme,ubrzanje,c='black')
plt.xlabel("$T[Gyr]$")
plt.ylabel(r"$a[10^{-13} kpc/s^2]$")
plt.xlim(min(vreme),max(vreme))
plt.savefig(file_out + "ubrzanje",dpi = 300)
plt.show()
