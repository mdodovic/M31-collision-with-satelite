# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:44:51 2018

@author: matij
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

if __name__ == "__main__":
    pathFolderWithDF = 'C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiSaTrenjemSveKomponente/'
    pathFolderWithoutDF = 'C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiBezTrenjaSveKomponente/'
    nameFile = 'radijus.txt'
    tNDF, rNDF = np.loadtxt(pathFolderWithoutDF + nameFile, unpack = True)
    tDF, rDF = np.loadtxt(pathFolderWithDF + nameFile, unpack = True)
    plt.plot(tNDF, rNDF)
    plt.plot(tNDF, rDF)
    plt.xlabel('T[Gyr]')
    plt.ylabel('R[kpc]')
    bezTrenja = mlines.Line2D([], [], color='blue', markersize=15, label='bez_trenja', linestyle = '-')
    saTrenjem = mlines.Line2D([], [], color='green', markersize=15, label='sa_trenjem',linestyle = '-')
    plt.legend(handles=[bezTrenja,saTrenjem], loc=4)

    plt.savefig("radiusDF_NDF.png", dpi=90)
    plt.show()
