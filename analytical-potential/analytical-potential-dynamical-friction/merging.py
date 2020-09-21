import math as math
import numpy as np
import matplotlib.pyplot as plt

#Linkovi

file_in ="C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiSaTrenjemSveKomponente/snapshotovi/"

file_out ="C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiSaTrenjemSveKomponente/"


N = 202

M = [];
X = [];
Y = [];
Z = [];
VX = [];
VY = [];
VZ = [];



for i in range(N):
        
    m,x,y,z,vx,vy,vz= np.loadtxt(file_in + "cestica_" + str(i).zfill(4) + ".txt", unpack=True)
    M.append(m)
    X.append(x)
    Y.append(y)
    Z.append(z)
    VX.append(vx)
    VY.append(vy)
    VZ.append(vz)

file = open(file_out + "cestica.txt","w")
for i in range(len(M)):
    file.write(str(M[i]) + " " + str(X[i]) + " " + str(Y[i]) + " " + str(Z[i]) + " " + str(VX[i]) + " " + str(VY[i]) + " " + str(VZ[i]) + "\n")
file.close