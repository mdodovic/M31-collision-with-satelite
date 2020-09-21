import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy.stats as sts

# =============================================================================
# Ucitavanje svih podataka - koordinata i IDeva za trenutke 2.0Gyrs - 3.0Gyrs
# =============================================================================

zeta00=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta00.txt", unpack=True)
eta00=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta00.txt", unpack=True)
zeta01=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta01.txt", unpack=True)
eta01=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta01.txt", unpack=True)
zeta02=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta02.txt", unpack=True)
eta02=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta02.txt", unpack=True)
zeta03=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta03.txt", unpack=True)
eta03=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta03.txt", unpack=True)
zeta04=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta04.txt", unpack=True)
eta04=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta04.txt", unpack=True)
zeta05=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta05.txt", unpack=True)
eta05=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta05.txt", unpack=True)
zeta06=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta06.txt", unpack=True)
eta06=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta06.txt", unpack=True)
zeta07=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta07.txt", unpack=True)
eta07=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta07.txt", unpack=True)
zeta08=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta08.txt", unpack=True)
eta08=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta08.txt", unpack=True)
zeta09=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta09.txt", unpack=True)
eta09=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta09.txt", unpack=True)
zeta10=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta10.txt", unpack=True)
eta10=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta10.txt", unpack=True)
zeta11=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta11.txt", unpack=True)
eta11=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta11.txt", unpack=True)
zeta12=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta12.txt", unpack=True)
eta12=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta12.txt", unpack=True)
zeta13=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta13.txt", unpack=True)
eta13=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta13.txt", unpack=True)
zeta14=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta14.txt", unpack=True)
eta14=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta14.txt", unpack=True)
zeta15=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta15.txt", unpack=True)
eta15=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta15.txt", unpack=True)
zeta16=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta16.txt", unpack=True)
eta16=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta16.txt", unpack=True)
zeta17=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta17.txt", unpack=True)
eta17=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta17.txt", unpack=True)
zeta18=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta18.txt", unpack=True)
eta18=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta18.txt", unpack=True)
zeta19=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta19.txt", unpack=True)
eta19=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta19.txt", unpack=True)
zeta20=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta20.txt", unpack=True)
eta20=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta20.txt", unpack=True)
zeta21=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta21.txt", unpack=True)
eta21=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta21.txt", unpack=True)
zeta22=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta22.txt", unpack=True)
eta22=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta22.txt", unpack=True)
zeta23=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta23.txt", unpack=True)
eta23=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta23.txt", unpack=True)
zeta24=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta24.txt", unpack=True)
eta24=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta24.txt", unpack=True)
zeta25=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta25.txt", unpack=True)
eta25=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta25.txt", unpack=True)
zeta26=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta26.txt", unpack=True)
eta26=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta26.txt", unpack=True)
zeta27=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta27.txt", unpack=True)
eta27=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta27.txt", unpack=True)
zeta28=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta28.txt", unpack=True)
eta28=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta28.txt", unpack=True)
zeta29=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta29.txt", unpack=True)
eta29=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta29.txt", unpack=True)
zeta30=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta30.txt", unpack=True)
eta30=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta30.txt", unpack=True)

print("Ucitane zeta i eta...")

id00=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id00.txt", unpack=True)
id01=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id01.txt", unpack=True)
id02=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id02.txt", unpack=True)
id03=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id03.txt", unpack=True)
id04=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id04.txt", unpack=True)
id05=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id05.txt", unpack=True)
id06=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id06.txt", unpack=True)
id07=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id07.txt", unpack=True)
id08=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id08.txt", unpack=True)
id09=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id09.txt", unpack=True)
id10=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id10.txt", unpack=True)
id11=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id11.txt", unpack=True)
id12=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id12.txt", unpack=True)
id13=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id13.txt", unpack=True)
id14=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id14.txt", unpack=True)
id15=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id15.txt", unpack=True)
id16=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id16.txt", unpack=True)
id17=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id17.txt", unpack=True)
id18=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id18.txt", unpack=True)
id19=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id19.txt", unpack=True)
id20=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id20.txt", unpack=True)
id21=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id21.txt", unpack=True)
id22=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id22.txt", unpack=True)
id23=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id23.txt", unpack=True)
id24=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id24.txt", unpack=True)
id25=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id25.txt", unpack=True)
id26=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id26.txt", unpack=True)
id27=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id27.txt", unpack=True)
id28=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id28.txt", unpack=True)
id29=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id29.txt", unpack=True)
id30=np.loadtxt("C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id30.txt", unpack=True)

print("Sve ucitano...")

# =============================================================================
# Definisanje pomsatrackih polja i konstanti
# =============================================================================

Ndisk = 205176 #broj cestica barionske materije Andromede
Nbuldge = 262005 #broj cestica tamne materije Andromede
NAnd = Ndisk + Nbuldge #ukupan broj cestica u Andromedi
Nhalo = 131072 #broj cestica barionske materije patuljka
Nstars = 248809 #broj cestica tamne materije patuljka
mdisk = 3.40944e-005 #masa jedne cestice
mbuldge = 0.000337753 
mhalo = 1.67998e-006
mstars = 1.67980e-005
masapatuljka = Nhalo * mhalo

x1=[2.015, 1.745, 1.483, 1.226, 0.969, 0.717, 0.467, 0.219]
y1=[-3.965, -3.525, -3.087, -2.653, -2.264, -1.768, -1.327,-0.886]
xerr=0.33
yerr=0.22
x1=np.asarray(x1)
y1=np.asarray(y1)

# =============================================================================
# Definisanje nizova i funkcija
# =============================================================================

def skracivanje(zet, et):    #funkcija koja uzima samo delove koji su potrebni 
    zetkon = []               #da bimse posmatracki podaci fitovali
    etkon = []
    for i in range (0, len(zet)):            
        if zet[i] > -0.5 and et[i] < -0.5:
            zetkon.append(zet[i])
            etkon.append(et[i])
        elif zet[i] < 1.3 and et[i] < -4:
            zetkon.append(zet[i])
            etkon.append(et[i])
    etkon = np.asarray(etkon)
    zetkon = np.asarray(zetkon)
    return etkon, zetkon

def centar_mase(ajdi, zeta, eta):
    pomocnizet = 0
    pomocniet = 0
    for i in range(0, Nhalo):
        if eta[int(ajdi[i])-NAnd] > -2:
            pomocnizet = pomocnizet + zeta[int(ajdi[i]-NAnd)] * mhalo
            pomocniet = pomocniet + eta[int(ajdi[i])-NAnd] * mhalo
    pomocnizet = np.asarray(pomocnizet)
    pomocniet = np.asarray(pomocniet)
    izlazzet = 1/masapatuljka * pomocnizet
    izlazet = 1/masapatuljka * pomocniet
    return izlazzet, izlazet

print("Skracivanje ID-eva")

id00 = id00[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id01 = id01[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id02 = id02[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id03 = id03[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id04 = id04[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id05 = id05[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id06 = id06[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id07 = id07[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id08 = id08[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id09 = id09[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id10 = id10[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id11 = id11[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id12 = id12[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id13 = id13[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id14 = id14[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id15 = id15[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id16 = id16[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id17 = id17[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id18 = id18[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id19 = id19[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id20 = id20[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id21 = id21[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id22 = id22[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id23 = id23[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id24 = id24[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id25 = id25[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id26 = id26[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id27 = id27[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id28 = id28[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id29 = id29[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]
id30 = id30[Ndisk+Nbuldge:Ndisk+Nbuldge+Nhalo]

# =============================================================================
# Primena funkcije za centar mase na podatke
# =============================================================================

print("Racanmo centre mase...")

centmaszet00, centmaset00 = centar_mase(id00, zeta00, eta00)
centmaszet01, centmaset01 = centar_mase(id01, zeta01, eta01)
centmaszet02, centmaset02 = centar_mase(id02, zeta02, eta02)
centmaszet03, centmaset03 = centar_mase(id03, zeta03, eta03)
centmaszet04, centmaset04 = centar_mase(id04, zeta04, eta04)
centmaszet05, centmaset05 = centar_mase(id05, zeta05, eta05)
centmaszet06, centmaset06 = centar_mase(id06, zeta06, eta06)
centmaszet07, centmaset07 = centar_mase(id07, zeta07, eta07)
centmaszet08, centmaset08 = centar_mase(id08, zeta08, eta08)
centmaszet09, centmaset09 = centar_mase(id09, zeta09, eta09)
centmaszet10, centmaset10 = centar_mase(id10, zeta10, eta10)
centmaszet11, centmaset11 = centar_mase(id11, zeta11, eta11)
centmaszet12, centmaset12 = centar_mase(id12, zeta12, eta12)
centmaszet13, centmaset13 = centar_mase(id13, zeta13, eta13)
centmaszet14, centmaset14 = centar_mase(id14, zeta14, eta14)
centmaszet15, centmaset15 = centar_mase(id15, zeta15, eta15)
centmaszet16, centmaset16 = centar_mase(id16, zeta16, eta16)
centmaszet17, centmaset17 = centar_mase(id17, zeta17, eta17)
centmaszet18, centmaset18 = centar_mase(id18, zeta18, eta18)
centmaszet19, centmaset19 = centar_mase(id19, zeta19, eta19)
centmaszet20, centmaset20 = centar_mase(id20, zeta20, eta20)
centmaszet21, centmaset21 = centar_mase(id21, zeta21, eta21)
centmaszet22, centmaset22 = centar_mase(id22, zeta22, eta22)
centmaszet23, centmaset23 = centar_mase(id23, zeta23, eta23)
centmaszet24, centmaset24 = centar_mase(id24, zeta24, eta24)
centmaszet25, centmaset25 = centar_mase(id25, zeta25, eta25)
centmaszet26, centmaset26 = centar_mase(id26, zeta26, eta26)
centmaszet27, centmaset27 = centar_mase(id27, zeta27, eta27)
centmaszet28, centmaset28 = centar_mase(id28, zeta28, eta28)
centmaszet29, centmaset29 = centar_mase(id29, zeta29, eta29)
centmaszet30, centmaset30 = centar_mase(id30, zeta30, eta30)

# =============================================================================
# Primena funkcije skracivanja na podatke
# =============================================================================

#etakon20, zetakon20 = skracivanje(zeta20, eta20)
#etakon21, zetakon21 = skracivanje(zeta21, eta21)
#etakon22, zetakon22 = skracivanje(zeta22, eta22)
#etakon23, zetakon23 = skracivanje(zeta23, eta23)
#etakon24, zetakon24 = skracivanje(zeta24, eta24)
#etakon25, zetakon25 = skracivanje(zeta25, eta25)
#etakon26, zetakon26 = skracivanje(zeta26, eta26)
#etakon27, zetakon27 = skracivanje(zeta27, eta27)
#etakon28, zetakon28 = skracivanje(zeta28, eta28)
#etakon29, zetakon29 = skracivanje(zeta29, eta29)
#etakon30, zetakon30 = skracivanje(zeta30, eta30)

# =============================================================================
# Radimo linearne fitove orijentacija streama
# =============================================================================

#slope20, intercept20, r_value, p_value, std_err = sts.linregress(zetakon20, etakon20) 
#slope21, intercept21, r_value, p_value, std_err = sts.linregress(zetakon21, etakon21) 
#slope22, intercept22, r_value, p_value, std_err = sts.linregress(zetakon22, etakon22) 
#slope23, intercept23, r_value, p_value, std_err = sts.linregress(zetakon23, etakon23) 
#slope24, intercept24, r_value, p_value, std_err = sts.linregress(zetakon24, etakon24) 
#slope25, intercept25, r_value, p_value, std_err = sts.linregress(zetakon25, etakon25) 
#slope26, intercept26, r_value, p_value, std_err = sts.linregress(zetakon26, etakon26) 
#slope27, intercept27, r_value, p_value, std_err = sts.linregress(zetakon27, etakon27) 
#slope28, intercept28, r_value, p_value, std_err = sts.linregress(zetakon28, etakon28) 
#slope29, intercept29, r_value, p_value, std_err = sts.linregress(zetakon29, etakon29)
#slope30, intercept30, r_value, p_value, std_err = sts.linregress(zetakon30, etakon30) 

# =============================================================================
# Plotovanje orijentacije
# =============================================================================

print("Pocetak plotovanja...")

plt.scatter(zeta00, eta00, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet00, centmaset00, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon00, intercept00+zetakon00*slope00, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.0Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t00.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta01, eta01, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet01, centmaset01, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon01, intercept01+zetakon01*slope01, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.1Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t01.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta02, eta02, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet02, centmaset02, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon02, intercept02+zetakon02*slope02, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.2Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t02.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta03, eta03, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet03, centmaset03, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon03, intercept03+zetakon03*slope03, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.3Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t03.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta04, eta04, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet04, centmaset04, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon04, intercept04+zetakon04*slope04, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.4Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t04.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta05, eta05, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet05, centmaset05, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon05, intercept05+zetakon05*slope05, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.5Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t05.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta06, eta06, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet06, centmaset06, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon06, intercept06+zetakon06*slope06, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.6Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t06.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta07, eta07, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet07, centmaset07, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon07, intercept07+zetakon07*slope07, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.7Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t07.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta08, eta08, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet08, centmaset08, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon08, intercept08+zetakon08*slope08, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.8Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t08.png', dpi=1090 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta09, eta09, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet09, centmaset09, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon09, intercept09+zetakon09*slope09, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 0.9Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t09.png', dpi=1100 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta10, eta10, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet10, centmaset10, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon10, intercept10+zetakon10*slope10, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.0Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t10.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta11, eta11, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet11, centmaset11, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon11, intercept11+zetakon11*slope11, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.1Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t11.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta12, eta12, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet12, centmaset12, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon12, intercept12+zetakon12*slope12, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.2Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t12.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta13, eta13, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet13, centmaset13, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -6)
plt.ylim(-8, 15)
#plt.plot(zetakon13, intercept13+zetakon13*slope13, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.3Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t13.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta14, eta14, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet14, centmaset14, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon14, intercept14+zetakon14*slope14, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.4Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t14.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta15, eta15, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet15, centmaset15, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon15, intercept15+zetakon15*slope15, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.5Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t15.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta16, eta16, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet16, centmaset16, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon16, intercept16+zetakon16*slope16, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.6Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t16.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta17, eta17 , s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet18, centmaset18, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon17, intercept17+zetakon17*slope17, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.7Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t17.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta18, eta18, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet18, centmaset18, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon18, intercept18+zetakon18*slope18, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.8Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t18.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta19, eta19, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet19, centmaset19, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon19, intercept19+zetakon19*slope19, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 1.9Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t19.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta20, eta20, s=0.1, color='red') #ovo je plotovanje svih cestica
plt.scatter(centmaszet20, centmaset20, color='yellow') #ovo je plotovanje centra mase
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon20, intercept20+zetakon20*slope20, color='red') #ovo je plotovanje srednje vrednosti toka
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.0Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t20.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta21, eta21, s=0.1, color='red')
plt.scatter(centmaszet21, centmaset21, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon21, intercept21+zetakon21*slope21, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.1Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t21.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta22, eta22, s=0.1, color='red')
plt.scatter(centmaszet22, centmaset22, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon22, intercept22+zetakon22*slope22, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.2Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t22.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta23, eta23, s=0.1, color='red')
plt.scatter(centmaszet23, centmaset23, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon23, intercept23+zetakon23*slope23, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.3Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t23.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta24, eta24, s=0.1, color='red')
plt.scatter(centmaszet24, centmaset24, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon24, intercept24+zetakon24*slope24, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.4Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t24.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta25, eta25, s=0.1, color='red')
plt.scatter(centmaszet25, centmaset25, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon25, intercept25+zetakon25*slope25, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.5Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t25.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta26, eta26, s=0.1, color='red')
plt.scatter(centmaszet26, centmaset26, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon26, intercept26+zetakon26*slope26, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.6Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t26.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta27, eta27, s=0.1, color='red')
plt.scatter(centmaszet27, centmaset27, color='yellow')
plt.xlim(4, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon27, intercept27+zetakon27*slope27, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.7Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t27.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta28, eta28, s=0.1, color='red')
plt.scatter(centmaszet28, centmaset28, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon28, intercept28+zetakon28*slope28, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.8Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t28.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta29, eta29, s=0.1, color='red')
plt.scatter(centmaszet29, centmaset29, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon29, intercept29+zetakon29*slope29, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 2.9Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t29.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()

plt.scatter(zeta30, eta30, s=0.1, color='red')
plt.scatter(centmaszet30, centmaset30, color='yellow')
plt.xlim(3, -2)
plt.ylim(-8, 5)
#plt.plot(zetakon30, intercept30+zetakon30*slope30, color='red')
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.title("Orijentacija toka i centar mase")
#plt.xlim(5, -0.5)
#plt.ylim(-5,-0.5)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
disk = mlines.Line2D([], [], color='white', markersize=15, label='t = 3.0Gyrs')
plt.legend(handles=[disk], loc=2, fontsize = 15)
plt.savefig('C:/Users/Pedja/Desktop/Orijentacija plotovi/t30.png', dpi=1200 , edgecolor='w' , facecolor='w' , bbox_inches = 'tight')
plt.show()
