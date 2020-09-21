# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:41:24 2018

@author: matij
"""

import math as math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import timeit
import scipy as sp


fileOut = fileIn = "C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/mappingOfMasses/"
fileOut = "C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/mappingOfMasses/final/"
x1=[2.015, 1.745, 1.483, 1.226, 0.969, 0.717, 0.467, 0.219]
y1=[-3.965, -3.525, -3.087, -2.653, -2.264, -1.768, -1.327,-0.886]
xerr=0.33
yerr=0.22

plt.errorbar(x1, y1, yerr,xerr, fmt='', color = 'black',linestyle="None")
plt.xlim(5, -1.5)
plt.ylim(-12, 3)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
for i in range(6,1,-1):
    
    zeta, eta = np.loadtxt(fileIn + "trajectoryMass = " + str(i) + ".txt", unpack = True)    
    plt.plot(zeta,eta, label='$Mass = {i}e10$'.format(i=i))
    #plt.text(2.5, 6, "M = " + textOnImage, fontsize = 15)
    #plt.savefig(file_out + "trajectory - Mass = " + nameOfImage + ".png", dpi = 300)

plt.legend(loc='best')
plt.savefig(fileOut + "trajectory2e10-6e10.png", dpi = 300)
plt.show()



Mass1_black = mlines.Line2D([], [], color='black', label='M = 6e10',marker='.', markersize=5)
Mass2_black = mlines.Line2D([], [], color='black', label='M = 5e10',marker='h', markersize=5)
Mass3_black = mlines.Line2D([], [], color='black', label='M = 4e10',marker='o', markersize=5)
Mass4_black = mlines.Line2D([], [], color='black', label='M = 3e10',marker='^', markersize=5)
Mass5_black = mlines.Line2D([], [], color='black', label='M = 2e10',marker='*', markersize=5)
Mass6_black = mlines.Line2D([], [], color='black', label='M = ',marker='+', markersize=5)
markerArray = ['.', 'h', 'o', '^', '*', '+']

# #############################################################################
plt.errorbar(x1, y1, yerr,xerr, fmt='', color = 'black',linestyle="None")
plt.xlim(5, -1.5)
plt.ylim(-12, 3)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)


for i in range(6,1,-1):
    
    zeta, eta = np.loadtxt(fileIn + "trajectoryMass = " + str(i) + ".txt", unpack = True)    
    plt.plot(zeta[0::50],eta[0::50], label='$Mass = {i}e10$'.format(i=i), color = 'black', marker = markerArray[1-i], markersize = 1)
    #plt.text(2.5, 6, "M = " + textOnImage, fontsize = 15)
    #plt.savefig(file_out + "trajectory - Mass = " + nameOfImage + ".png", dpi = 300)

plt.legend(handles=[Mass1_black,Mass2_black,Mass3_black,Mass4_black,Mass5_black], loc='best')
plt.savefig(fileOut + "trajectory2e10-6e10.png", dpi = 300)
plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/trajectory2e10-6e10_crno.png", dpi = 90)
plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika7Gore.png", dpi = 90)
plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika7Gore.eps", dpi = 90)
plt.show()

########################################################################################
########################################################################################
plt.errorbar(x1, y1, yerr,xerr, fmt='', color = 'black',linestyle="None")
plt.xlim(4, -0.5)
plt.ylim(-9, 3)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)

zeta, eta = np.loadtxt(fileIn + "trajectoryMass = 5" + str(1) + ".txt", unpack = True)
plt.plot(zeta,eta, label='$Mass = 5.{i}e10$'.format(i=1))
for i in range(9,-1,-2):    
    zeta, eta = np.loadtxt(fileIn + "trajectoryMass = 4" + str(i) + ".txt", unpack = True)
    
    plt.plot(zeta,eta, label='$Mass = 4.{i}e10$'.format(i=i))
    #plt.text(2.5, 6, "M = " + textOnImage, fontsize = 15)
    #plt.savefig(file_out + "trajectory - Mass = " + nameOfImage + ".png", dpi = 300)


plt.legend(loc='best')
plt.savefig(fileOut + "trajectory4e10-5e10.png", dpi = 300)
plt.show()

Mass1_black = mlines.Line2D([], [], color='black', label='M = 5.1e10',marker='.', markersize=5)
Mass2_black = mlines.Line2D([], [], color='black', label='M = 4.9e10',marker='h', markersize=5)
Mass3_black = mlines.Line2D([], [], color='black', label='M = 4.7e10',marker='o', markersize=5)
Mass4_black = mlines.Line2D([], [], color='black', label='M = 4.5e10',marker='^', markersize=5)
Mass5_black = mlines.Line2D([], [], color='black', label='M = 4.3e10',marker='*', markersize=5)
Mass6_black = mlines.Line2D([], [], color='black', label='M = 4.1e10',marker='+', markersize=5)
markerArray = ['.', 'h', 'o', '^', '*', '+']

plt.errorbar(x1, y1, yerr,xerr, fmt='', color = 'black',linestyle="None")
plt.xlim(4, -0.5)
plt.ylim(-9, 3)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)

zeta, eta = np.loadtxt(fileIn + "trajectoryMass = 5" + str(1) + ".txt", unpack = True)
plt.plot(zeta[0::50],eta[0::50], label='$Mass = 5.{i}e10$'.format(i=1) , color = 'black', marker = markerArray[0], markersize = 1)
j = 1
for i in range(9,-1,-2):    
    zeta, eta = np.loadtxt(fileIn + "trajectoryMass = 4" + str(i) + ".txt", unpack = True)    
    plt.plot(zeta[0::50],eta[0::50], color = 'black', marker = markerArray[j], markersize = 1)
    #plt.text(2.5, 6, "M = " + textOnImage, fontsize = 15)
    #plt.savefig(file_out + "trajectory - Mass = " + nameOfImage + ".png", dpi = 300)
    j = j+1

plt.legend(handles=[Mass1_black,Mass2_black,Mass3_black,Mass4_black,Mass5_black,Mass6_black], loc='best')

plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/trajectory4e10-5e10_crno.png", dpi = 90)
plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika7Dole.png", dpi = 90)
plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika7Dole.eps", dpi = 90)
plt.savefig(fileOut + "trajectory4e10-5e10.png", dpi = 300)
plt.show()

# #############################################################################
# #############################################################################

plt.errorbar(x1, y1, yerr,xerr, fmt='', color = 'black',linestyle="None")
plt.xlim(3.4, -1.5)
plt.ylim(-7, 3.2)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)


zeta, eta = np.loadtxt(fileIn + "trajectoryMass = 110.txt", unpack = True)
plt.plot(zeta,eta, label='$Mass = 1e11$')

for i in range(90,58,-10):    
    zeta, eta = np.loadtxt(fileIn + "trajectoryMass = " + str(i) + ".txt", unpack = True)
    
    plt.plot(zeta,eta, label='$Mass = {i}e10$'.format(i=int(i/10)))
    #plt.text(2.5, 6, "M = " + textOnImage, fontsize = 15)
    #plt.savefig(file_out + "trajectory - Mass = " + nameOfImage + ".png", dpi = 300)


plt.legend(loc='best')
plt.savefig(fileOut + "trajectory6e10-1e11.png", dpi = 300)
plt.show()


Mass1_black = mlines.Line2D([], [], color='black', label='M = 1e11',marker='.', markersize=5)
Mass2_black = mlines.Line2D([], [], color='black', label='M = 9e10',marker='h', markersize=5)
Mass3_black = mlines.Line2D([], [], color='black', label='M = 8e10',marker='o', markersize=5)
Mass4_black = mlines.Line2D([], [], color='black', label='M = 7e10',marker='^', markersize=5)
Mass5_black = mlines.Line2D([], [], color='black', label='M = 6e10',marker='*', markersize=5)
markerArray = ['.', 'h', 'o', '^', '*']


plt.errorbar(x1, y1, yerr,xerr, fmt='', color = 'black',linestyle="None")
plt.xlim(3.4, -1.5)
plt.ylim(-7, 3.2)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)


zeta, eta = np.loadtxt(fileIn + "trajectoryMass = 110.txt", unpack = True)
plt.plot(zeta[0::10],eta[0::10], color = 'black', marker = markerArray[0], markersize = 1)
j= 1
for i in range(90,58,-10):    
    zeta, eta = np.loadtxt(fileIn + "trajectoryMass = " + str(i) + ".txt", unpack = True)
    plt.plot(zeta[0::10],eta[0::10], color = 'black', marker = markerArray[j], markersize = 1)
    #plt.text(2.5, 6, "M = " + textOnImage, fontsize = 15)
    #plt.savefig(file_out + "trajectory - Mass = " + nameOfImage + ".png", dpi = 300)
    j += 1

plt.savefig(fileOut + "trajectory6e10-1e11.png", dpi = 300)
plt.legend(handles=[Mass1_black,Mass2_black,Mass3_black,Mass4_black,Mass5_black], loc='best')

plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/trajectory6e10-1e11.png", dpi = 90)
plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika7Sredina.png", dpi = 90)
plt.savefig("C:/Users/Matija/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/graficiRad/vektorski/slika7Sredina.eps", dpi = 90)
plt.show()
