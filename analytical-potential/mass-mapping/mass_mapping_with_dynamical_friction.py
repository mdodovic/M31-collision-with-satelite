import math as math
import numpy as np
import matplotlib.pyplot as plt
import timeit
import scipy as sp


start = timeit.default_timer() #Pocetak programa, sluzi da se odredi duzina trajanja izvrsavanja koda

#Linkovi

#Apsolutna lokacija gde da se cuvaju slike:
file_out = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/mappingOfMasses/"

#Apsolutna lokacija gde da se cuaju snapshotovi (podaci masa, x,y,z,vx,vy,vz za svaki vremenski korak):
#snapshot_out = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/JuzniTok/haloAP_SaTrenjem/snapshotovi/"

#Konstante
G = 4.515353945657138e-39 # Gravitaciona konstanta [kpc^3/(Msol*s^2)]
kriticnaGustina = 143.5012 # Kriticna gustina svemira[Msol/kpc^3]
pi_sqrt = (math.sqrt(math.pi)) #Koren iz pi
dva_sqrt = math.sqrt(2) #Koren iz 2
konstantaDinamickoTrenje = -4*math.pi*G*G #Konstanta za dinamicko trenje, koje se koristi u delu gde se racuna sila dinamickog trenja

# Podaci o analitickom potencijalu
#Mase komponenata
#Mbulge = 3.2e9 # [Msol]
#Mdisk = 3.7e9 # [Msol]
Mhalo = 88e10 # [Msol]
Map = Mhalo #+ Mdisk + Mbulge # [Msol]

rVir = 200.2 #Virijalni radijus [kpc]
massRVir = 964642422540.9172 #masa haloa unutar virijalnog radijusa [Msol] 

#Polumaseni radijusi [kpc]
halfMassRadiusH = 16.37329999997518 #radijus na kojoj je masa haloa pola

#Scale length
scale_length_halo = 8.18 # [kpc]
scale_length_halo_3 = scale_length_halo**3 # scale length na 3[kpc^3] 

#Konstanta za f funkcije, sluze za dinamicko trenje 
fhConst = -1.0 * (G * massRVir)/ (math.log((1+rVir/scale_length_halo) - (rVir / (rVir + scale_length_halo))))

delta_c = 429231.966735 # bezdimenzioni parametar - izracunat

#Konstanta potecijala
potencijalHaloConst = -4*np.pi*G*delta_c*kriticnaGustina*(scale_length_halo**3)



#Pocetni uslovi vremena simulacije

dt = 86400*365*1e6 # korak [s]
t=0 #pocetak simulacije
k = 0 # broj snapshota
tmax=86400*365*1e9 # trajanje simulacije [s] 
n = 3 #broj milijardi godina
tmax = tmax * n # za n milijardi godina



def intenzitet(a,b,c):
    """
    Funkcija za izracunavanje intenziteta vektora sa 3 komponente    
    Input parametri:
    a, b, c - komponente vektora
    Output parametri:
    intenzitet vektra
    """
    return (a*a + b*b + c*c)**0.5

def masah(r):
    """
    Izracunavanje mase haloa unutar radijusa r    
    """
    M = 4*math.pi*delta_c*kriticnaGustina*(scale_length_halo_3)*(math.log((r+scale_length_halo)/scale_length_halo) - (r/(r+scale_length_halo)))
    return M

def masa(r):
    #poziv funkcije masah
    return masah(r)

def potencijalH(r):
    #izracunavanje potecijala haloa na radijusu R
    FI = potencijalHaloConst/r*np.log((r+scale_length_halo)/scale_length_halo)
    return FI

def potencijal(r):
    #Pozvanje funkijce potecijalH
    return potencijalH(r)

# FUNKCIJE GUSTINE
def gustinaH(r):
    #NFW - gustina na radijusu r
    RHO = delta_c*kriticnaGustina / ( (r / scale_length_halo) * (1 + r / scale_length_halo) * (1 + r / scale_length_halo) )
    return RHO

#Dinamicko trenje    
def fh(r):
    #Fh formula iz rada o subhaloima (Penarrubia J. and Benson A. 2005)
    fh = fhConst * (math.log((1+r/scale_length_halo) - (r/(r + scale_length_halo)))) / (r**2) 
    return fh
def sigmaH(donja_granica):
    #Formula za disperziju brzina iz rada (Penarrubia J. and Benson A. 2005)    
    #donja_granica je polozaj cestice 
    gornja_granica = 1000.# beskonacno velika donja granica integrala
    korak = 1. #[kpc]    
    x = []
    y = []
    #print(donja_granica)
        
    r_prim = donja_granica
    #Pravljenje nizova za racunanje integrala
    while r_prim <= gornja_granica:
        x.append(r_prim) #dr'
        y.append(fh(r_prim)*gustinaH(r_prim)) #podintegralna funkcija
        r_prim += korak

    k = -sp.integrate.simps(y, x) #izracunavanje integrala
    y[:] = []
    x[:] = []
    s = (1/gustinaH(r)) * k #Disperzija
    return s
    
def xhi(v, disperzijaBrzine):    
    X = v/(dva_sqrt*disperzijaBrzine)
    return X        

def coulombLog1(r, disperzijaBrzine):
    #Biney i Tremain    
    bmax = r
    rh = halfMassRadiusH
    vtyp = disperzijaBrzine 
    M = Mhalo
    L = math.log(bmax/max(rh,G*M/(vtyp**2)))        
    return L    

def coulombLogH(r,disperzijaBrzine):
    return coulombLog1(r, disperzijaBrzine)

def silaDinamickoTrenjeHalo(m, r, disperzijaBrzine, vx, vy, vz, koordinatnaBrzina):
    v = intenzitet(vx,vy,vz) #intenzitet brzine 
    X = xhi(v, disperzijaBrzine) 
    a = ((konstantaDinamickoTrenje * m * gustinaH(r) * coulombLogH(r, disperzijaBrzine))/(v**3))*(sp.special.erf(X) - (2*X*(math.e**(-X*X)))/pi_sqrt) * koordinatnaBrzina
    return a

def siladinamickoTrenje(m, r, disperzijaBrzine, vx, vy, vz, koordinatnaBrzina):
    #pozivanje sile za halo 
    a = silaDinamickoTrenjeHalo(m, r, disperzijaBrzine, vx, vy, vz, koordinatnaBrzina)
    return a

#Klasa telo
class telo:
    x = 0
    y = 0
    z = 0
    vx = 0
    vxh = 0
    vy = 0
    vyh = 0
    vz = 0
    vzh = 0
    ax = 0
    ay = 0
    az = 0    
    m = 0
    Ek = 0
    Ep = 0
    L = 0
    
    def __init__(self,x,y,z,vx,vy,vz,ax,ay,az,m,Ek,Ep,L):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.ax = ax
        self.ay = ay
        self.az = az
        self.m = m
        self.Ek = Ek
        self.Ep = Ep
        self.L = L

 
tela = []        

#Pocetni uslovi cestice
N=1
pocetno_rast_analiticki_x = 0
pocetno_rast_analiticki_y = 0
pocetno_rast_analiticki_z = 0

#Pocetni uslovi 
x = -84.41 #kpc
y = 152.47
z = -97.08
vx = 123.631622883 * 3.24078e-17  / 100
vy = 1/10 *-49.6091823145 * 3.24078e-17 / 100
vz = 30 * 36.3023759879 * 3.24078e-17 / 100

M = 8.5e10 #Msol
textOnImage = "8.5e10"
nameOfImage = "85"

#import sys
#sys.exit()

#print(M/Map*100) # ako je vece od 20 ne bi valjalo

#Energije i momenti impulsa
Ek = 1/2*M*((intenzitet(vx,vy,vz))**2)
Ep = potencijal(intenzitet(x,y,z))*M
L = M*intenzitet(x,y,z)*intenzitet(vx,vy,vz)


tela.append(telo(x,y,z,vx,vy,vz,0,0,0,M,Ek,Ep,L))

#print(1/2*tela[0].m*intenzitet(tela[0].vx, tela[0].vy, tela[0].vz)**2)###OVO SVUDA
#print(1/2*M*intenzitet(0, vk, 0)**2)

#print(potencijal(intenzitet(tela[0].x,tela[0].y,tela[0].z)))
#print(tela[0].Ek)
#print(tela[0].Ep)
#print(tela[0].Ek + tela[0].Ep)



x = []
y = []
z = []

ukupnaEnergija = []
momentImpulsa = []
radialVelocity = []
vreme = []
radijus = []

pomEk = 0
for i in range(N):
    pomEk += tela[i].Ek + tela[i].Ep
pomL = 0
for i in range(N):
    pomL += tela[i].L


radijus.append(intenzitet(tela[0].x,tela[0].y,tela[0].z))
radialVelocity.append(intenzitet(tela[0].vx,tela[0].vy,0))
ukupnaEnergija.append(pomEk)
momentImpulsa.append(pomL)
vreme.append(t)
k += 1


for i in range(N):
    tela[i].vxh = tela[i].vx + tela[i].ax * dt / 2
    tela[i].vyh = tela[i].vy + tela[i].ay * dt / 2
    tela[i].vzh = tela[i].vz + tela[i].az * dt / 2
    
for i in range(N):    
    tela[i].vx += -(tela[i].ax * dt)
    tela[i].vy += -(tela[i].ay * dt)
    tela[i].vz += -(tela[i].az * dt)

    tela[i].Ek = 1/2*tela[i].m*intenzitet(tela[i].vx, tela[i].vy, tela[i].vz)**2

for i in range(N):
    tela[i].x += tela[i].vxh * dt
    tela[i].y += tela[i].vyh * dt
    tela[i].z += tela[i].vzh * dt
    tela[i].L = tela[i].m * intenzitet(tela[i].vx,tela[i].vy,tela[i].vz) * intenzitet(tela[i].x,tela[i].y,tela[i].z)
    tela[i].Ep = potencijal(intenzitet(tela[i].x,tela[i].y,tela[i].z))*tela[i].m

t+=dt    

#print(potencijal(intenzitet(tela[0].x,tela[0].y,tela[0].z)))
#print(tela[0].Ek)
#print(tela[0].Ep)
#print(tela[0].Ek + tela[0].Ep)



pomEk = 0
for i in range(N):
    pomEk += tela[i].Ek + tela[i].Ep
pomL = 0
for i in range(N):
    pomL += tela[i].L

momentImpulsa.append(pomL)
radijus.append(intenzitet(tela[0].x,tela[0].y,tela[0].z))        
radialVelocity.append(intenzitet(tela[0].vx,tela[0].vy,0))
ukupnaEnergija.append(pomEk)
vreme.append(t)
k += 1




while t<=tmax:
    #print(t/tmax*100)
    for i in range(N):
        Fx_N_body = 0
        Fy_N_body = 0
        Fz_N_body = 0
        for j in range(N):
            if(i != j):
                r = (intenzitet(tela[i].x - tela[j].x,tela[i].y - tela[j].y, tela[i].z - tela[j].z))
                Fx_N_body += - G*tela[i].m * tela[j].m *(tela[i].x - tela[j].x) / (r**3)   
                Fy_N_body += - G*tela[i].m * tela[j].m *(tela[i].y - tela[j].y) / (r**3)
                Fz_N_body += - G*tela[i].m * tela[j].m *(tela[i].z - tela[j].z) / (r**3)
        
        r = (intenzitet(tela[i].x - pocetno_rast_analiticki_x, tela[i].y - pocetno_rast_analiticki_y, tela[i].z - pocetno_rast_analiticki_z))
        Fx_analiticki = - G * tela[i].m * masa(r) * (tela[i].x - pocetno_rast_analiticki_x)/((r)**3)

        Fy_analiticki = - G * tela[i].m * masa(r) * (tela[i].y - pocetno_rast_analiticki_y)/((r)**3)

        Fz_analiticki = - G * tela[i].m * masa(r) * (tela[i].z - pocetno_rast_analiticki_z)/((r)**3) #masa po xyz ili kao radijus
        
        Fx = Fx_N_body + Fx_analiticki
        Fy = Fy_N_body + Fy_analiticki
        Fz = Fz_N_body + Fz_analiticki
        disperzijaBrzine = 425 * 3.2408e-17 #kpc/s        
        #disperzijaBrzine = math.sqrt(1.8970379e-28)#math.sqrt(sigmaH(r))
        #disperzijaBrzine = math.sqrt(sigmaH(r))
        #print(1/3.2408e-17 * disperzijaBrzine *10)
        tela[i].ax = Fx / tela[i].m + siladinamickoTrenje(tela[i].m, r, disperzijaBrzine, tela[i].vx, tela[i].vy, tela[i].vz, tela[i].vx)
        tela[i].ay = Fy / tela[i].m + siladinamickoTrenje(tela[i].m, r, disperzijaBrzine, tela[i].vx, tela[i].vy, tela[i].vz, tela[i].vy)
        tela[i].az = Fz / tela[i].m + siladinamickoTrenje(tela[i].m, r, disperzijaBrzine, tela[i].vx, tela[i].vy, tela[i].vz, tela[i].vz)
        
#ZA POTENCIJALNU ENERGIJU SISTEMA CESTICA
#    ek = 0
#    ep = 0
#    for i in range(N):
#        tela[i].Ek = 1/2*tela[i].m*(intenzitet(tela[i].vx,tela[i].vy,tela[i].vz))**2
#        for j in range(N):
#            if i != j:
#                tela[i].Ep = -G*tela[i].m*tela[j].m/(intenzitet(tela[i].x - tela[j].x,tela[i].y - tela[j].y, tela[i].z - tela[j].z))
#        ek += tela[i].Ek
#        ep += tela[i].Ep
        
#    Kineticka_energija.append(ek)
#    Potencijalna_energija.append(ep)
#    Ukupna_energija.append(ek+ep)

        
    for i in range(N):
        tela[i].vxh += tela[i].ax*dt
        tela[i].vyh += tela[i].ay*dt
        tela[i].vzh += tela[i].az*dt
        
    for i in range(N):
        tela[i].vx = -(tela[i].vxh + tela[i].ax * dt / 2)
        tela[i].vy = -(tela[i].vyh + tela[i].ay * dt / 2)
        tela[i].vz = -(tela[i].vzh + tela[i].az * dt / 2)
        tela[i].Ek = 1/2*tela[i].m*intenzitet(tela[i].vx, tela[i].vy, tela[i].vz)**2

    for i in range(N):
        tela[i].x += tela[i].vxh * dt
        x.append(tela[i].x)
        tela[i].y += tela[i].vyh * dt
        y.append(tela[i].y)
        tela[i].z += tela[i].vzh * dt
        z.append(tela[i].z)
        tela[i].L = tela[i].m * intenzitet(tela[i].vx,tela[i].vy,tela[i].vz) * intenzitet(tela[i].x,tela[i].y,tela[i].z)
        tela[i].Ep = potencijal(intenzitet(tela[i].x,tela[i].y,tela[i].z))*tela[i].m
    #print(tela[0].Ek)
    #print(tela[0].Ep)
    #print(tela[0].Ek + tela[0].Ep)
        

    t+=dt

    pomEk = 0
    for i in range(N):
        pomEk += tela[i].Ek + tela[i].Ep
    pomL = 0
    for i in range(N):
        pomL += tela[i].L
    
    momentImpulsa.append(pomL)
    radijus.append(intenzitet(tela[0].x,tela[0].y,tela[0].z))        
    radialVelocity.append(intenzitet(tela[0].vx,tela[0].vy,0))
    ukupnaEnergija.append(pomEk)
    vreme.append(t)

for i in range(1,N):
    v1 = intenzitet(tela[i-1].vx,tela[i-1].vy,tela[i-1].vz)    
    v2 = intenzitet(tela[i].vx,tela[i].vy,tela[i].vz)    

for i in range(len(vreme)):
    vreme[i] = vreme[i]/(86400*365*1e9)
for i in range(len(vreme)):
    radialVelocity[i] = radialVelocity[i]/(3.24078e-17)


"""
#R(t)

plt.plot(vreme,radijus,c='black')
plt.xlabel("$T[Gyr]$")
plt.ylabel("$R[kpc]$")
plt.xlim(min(vreme),max(vreme))
plt.text(1.3, 170, "M = " + textOnImage, fontsize = 15)
plt.savefig(file_out + "radius - Mass = " + nameOfImage + ".png", dpi = 300)
plt.show()
"""

###########################
#POSMATRACKA POLJA
##########################

x1=[2.015, 1.745, 1.483, 1.226, 0.969, 0.717, 0.467, 0.219]
y1=[-3.965, -3.525, -3.087, -2.653, -2.264, -1.768, -1.327,-0.886]
xerr=0.33
yerr=0.22

zeta = []
eta = []
for i in range(len(x)):
    zeta.append((x[i]/(z[i]+784))*(180/np.pi))
    eta.append((y[i]/(z[i]+784))*(180/np.pi))

plt.plot(zeta,eta)
plt.errorbar(x1,y1,yerr,xerr,fmt='o', color = 'black')
plt.xlim(3, -3)
plt.ylim(-8, 8)
plt.xlabel("x [\N{DEGREE SIGN}]", fontsize=15)
plt.ylabel("y [\N{DEGREE SIGN}]", fontsize=15)
plt.text(2.5, 6, "M = " + textOnImage, fontsize = 15)
plt.savefig(file_out + "trajectory - Mass = " + nameOfImage + ".png", dpi = 300)
plt.show()

file = open(file_out + "trajectoryMass = " + nameOfImage + ".txt","w")
for i in range(len(x)):
    file.write(str(zeta[i]) + " " + str(eta[i]) + "\n")
file.close


#Duzina trajanja programa
stop = timeit.default_timer()
print("Duration of code: " + str(np.round(stop - start,2)) + " seconds") 

file = open("C:/Users/matij/Desktop/djubre/brisi.txt","w")
for i in range(len(x)):
    file.write(str(zeta[i]) + " " + str(eta[i]) + "\n")
file.close
