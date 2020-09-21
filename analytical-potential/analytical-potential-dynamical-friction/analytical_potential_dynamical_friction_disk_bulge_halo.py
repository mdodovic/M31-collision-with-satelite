import math as math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import timeit

start = timeit.default_timer()


#Linkovi

file_out = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiSaTrenjemSveKomponente/"
snapshot_out = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiSaTrenjemSveKomponente/snapshotovi/"

#Konstante
G = 4.515353945657138e-39 # [kpc^3/(Msol*s^2)]
kriticnaGustina = 277.72/(0.71*0.71) # [Msol/kpc^2]
pi_sqrt = (math.sqrt(math.pi))
dva_sqrt = math.sqrt(2)
konstantaDinamickoTrenje = -4*math.pi*G*G
# Podaci o analitickom potencijalu

Mbulge = 3.2e9 # [Msol]
Mdisk = 3.7e9 # [Msol]
Mhalo = 88.46e10 # [Msol]
Map = Mhalo + Mbulge# + Mdisk # [Msol]

rVir = 200.2 #Virijalni radijus
massRVir = 953655483574.0 #masa unutar virijalnog radijusa u masama sunca za celu galaksiju

#Polumaseni radijusi [kpc]
halfMassRadiusH = 16.37329999997518 #18.49558899734247 #radijus na kojoj je masa haloa pola
halfMassRadiusB =  1.6492999999998348#radijus na kojoj je masa bulgea pola

#Skale duzine
scale_length_halo = 8.18 # [kpc]
scale_length_halo_3 = scale_length_halo**3
scale_length_bulge = 1.22 # [kpc]
scale_length_disk_vertikalni = 0.57 # [kpc]
scale_length_disk_radijalni = 6.82 # [kpc]

#Konstanta za fi funkcije, i je h d b
fhConst = -1.0 * (G * massRVir)/ (math.log((1+rVir/scale_length_halo) - (rVir / (rVir + scale_length_halo))))
fbConst = -1.0 * (G * massRVir)/ (math.log((1+rVir/scale_length_bulge) - (rVir / (rVir + scale_length_bulge))))
#fdConst

delta_c = 27e4 # bezdimenzioni parametar
#Sigma_0 = Mdisk/(2*np.pi*scale_length_disk_radijalni**2) # centralna gustina diska

#Pocetni uslovi vremena simulacije

dt = 86400*365*1e6 # milion godina
dt = dt * 50 # za 50 miliona godina 
t=0
k = 0 # broj snapshota
tmax=86400*365*1e9 # milijardu godina
n = 11 #KOLIKO MILIJARDI GODINA, namestiti
tmax = tmax * n # za n milijardi godina

#Pomocne funckije
def intenzitet(a,b,c):
    return (a*a + b*b + c*c)**0.5

#Funkcije mase
def masad(r,z):
    M = Mdisk / (2*scale_length_disk_radijalni) * np.tanh(z/scale_length_disk_vertikalni) * (scale_length_disk_radijalni*(1 - (np.e)**(-r/scale_length_disk_radijalni)) - r*(np.e)**(-r/scale_length_disk_radijalni)) 
    return M

def masab(r):
    M = Mbulge*(r**2)/((scale_length_bulge+r)**2)
    return M
def masah(r):
    M = 4*math.pi*delta_c*kriticnaGustina*(scale_length_halo_3)*(math.log((r+scale_length_halo)/scale_length_halo) - (r/(r+scale_length_halo)))
    return M

def masa(r,x,y,z):
    #r2d = intenzitet(x,y,0)
    return masah(r) #+ masab(r) + masad(r2d,z)

#Funkcije potencijala, koje za sad ne trebaju
#Verovatno nisu tacne
"""
def potencijalH(r):
    FI = -4*np.pi*G*delta_c*kriticnaGustina*(scale_length_halo_3)/r*np.log((r+scale_length_halo)/scale_length_halo)
    return FI
def potencijalD(r):
    FI = -2*np.pi*G*Sigma_0*(scale_length_disk_radijalni**2)*(1-((np.e)**(-r/scale_length_disk_radijalni)))/r
    return FI
def potencijalB(r):
    FI = -G*Mbulge/(r+scale_length_bulge)    
    return FI    
     
def potencijal(r):
    return potencijalH(r) + potencijalD(r) + potencijalB(r)
"""

# FUNKCIJE GUSTINE
def gustinaH(r):
    RHO = delta_c*kriticnaGustina / ( (r / scale_length_halo) * (1 + r / scale_length_halo) * (1 + r / scale_length_halo) )
    return RHO
    
def gustinaB(r):
    RHO = Mbulge * scale_length_bulge/(2*np.pi*r) * 1 / ((r + scale_length_bulge)**3)
    return RHO


#Dinamicko trenje    
def fh(r):
    fh = fhConst * (math.log((1+r/scale_length_halo) - (r/(r + scale_length_halo)))) / (r**2) 
    return fh

def fb(r):
    fb = fbConst * (math.log((1+r/scale_length_bulge) - (r/(r + scale_length_bulge)))) / (r**2) 
    return fb

def sigmaH(donja_granica):
    #donja_granica je polozaj cestice 
    gornja_granica = 10000.# beskonacno velika donja granica integrala
    #r = np.linspace(donja_granica, gornja_granica, gornja_granica, True) #pavljenje niza rastojanja od polozaja cestice do velike vrednosti 1Mpc
    korak = 1. #1kpc    
    x = []
    y = []
    r_prim = donja_granica
    while r_prim <= gornja_granica:
        x.append(r_prim)
        y.append(fh(r_prim)*gustinaH(r_prim))
        r_prim += korak
    k = -sp.integrate.simps(y, x)
    y[:] = []
    x[:] = []
    s = (1/gustinaH(r)) * k
    return s

###DA LI JE HALF MASS RADIUS ZA ANALITICKI POTENCIJAL ILI JE TO IZ MODELA GALAKSIJE???


def sigmaB(donja_granica):
    #donja_granica je polozaj cestice 
    gornja_granica = 10000.# beskonacno velika donja granica integrala
    #r = np.linspace(donja_granica, gornja_granica, gornja_granica, True) #pavljenje niza rastojanja od polozaja cestice do velike vrednosti 1Mpc
    korak = 1. #1kpc    
    x = []
    y = []
    r_prim = donja_granica
    while r_prim <= gornja_granica:
        x.append(r_prim)
        y.append(fb(r_prim)*gustinaB(r_prim))
        r_prim += korak
    k = -sp.integrate.simps(y, x)
    y[:] = []
    x[:] = []
    s = (1/gustinaB(r)) * k
    return s

        
def xhi(r, v, disperzijaBrzine):
    X = v/(dva_sqrt*disperzijaBrzine)
    return X        

def coulombLogH(r, disperzijaBrzine):
    bmax = r
    rh = halfMassRadiusH
    vtyp = disperzijaBrzine 
    M = Mhalo
    L = math.log(bmax/max(rh,G*M/(vtyp**2)))
    return L
    
def coulombLogB(r, disperzijaBrzine):
    bmax = r
    rh = halfMassRadiusB
    vtyp = disperzijaBrzine 
    M = Mbulge
    L = math.log(bmax/max(rh,G*M/(vtyp**2)))
    return L
    
def ubrzanjeDinamickoTrenjeHalo(m, r, vx, vy, vz, koordinatnaBrzina):
    v = intenzitet(vx,vy,vz)
    disperzijaBrzine = sigmaH(r)
    X = xhi(v, r, disperzijaBrzine)
    a = ((konstantaDinamickoTrenje * m * gustinaH(r) * coulombLogH(r, disperzijaBrzine))/(v**3))*(sp.special.erf(X) - (2*X*(math.e**(-X*X)))/pi_sqrt) * koordinatnaBrzina
    return a

def ubrzanjeDinamickoTrenjeBulge(m, r, vx, vy, vz, koordinatnaBrzina):
    v = intenzitet(vx,vy,vz)
    disperzijaBrzine = sigmaB(r)
    X = xhi(v, r, disperzijaBrzine)
    a = ((konstantaDinamickoTrenje * m * gustinaB(r) * coulombLogB(r, disperzijaBrzine))/(v**3))*(sp.special.erf(X) - (2*X*(math.e**(-X*X)))/pi_sqrt) * koordinatnaBrzina
    return a


#def ubrzanjeDinamickoTrenjeDisk(m, r, vx, vy, vz, koordinatnaBrzina):
#    v = intenzitet(vx,vy,vz)
#    disperzijaBrzine = sigmaD(r)
#    X = xhi(v, r, disperzijaBrzine)
#    a = ((konstantaDinamickoTrenje * m * gustinaD(r) * coulombLogD(r, disperzijaBrzine))/(v**3))*(sp.special.erf(X) - (2*X*(math.e**(-X*X)))/pi_sqrt) * koordinatnaBrzina
#    return a

def dinamickoTrenje(m, r, vx, vy, vz, koordinatnaBrzina):
    a = ubrzanjeDinamickoTrenjeHalo(m, r, vx, vy, vz, koordinatnaBrzina) #+ ubrzanjeDinamickoTrenjeBulge(m, r, vx, vy, vz, koordinatnaBrzina)  
    return a

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
R = 200 #kpc

vk = (G*Map/R)**(0.5) #kpc/s
#print(vk)
#M = 4.42e10



#vk = 2.46798*3.24078e-17 #kpc/s
M = 88460300 # paziti na ovo da bude 1/10 od ukupne mase
#M = Map/10
Ek = 1/2*M*vk**2
#Ep = potencijal(R)*M
Ep = 0
L = M*R*vk
pocetno_rast_analiticki_x = 0
pocetno_rast_analiticki_y = 0
pocetno_rast_analiticki_z = 0

tela.append(telo(R,1e-15,1e-15,0,vk,0,0,0,0,M,Ek,Ep,L))

x = []
y = []
z = []

momentImpulsa = []
ukupnaEnergija = []
vreme = []
radijus = []

pomEk = 0
for i in range(N):
    pomEk += tela[i].Ek

pomL = 0
for i in range(N):
    pomL += tela[i].L


radijus.append(R)
ukupnaEnergija.append(pomEk)
momentImpulsa.append(pomL)
vreme.append(t)
file = open(snapshot_out + "cestica_" + str(k).zfill(4) + ".txt","w")
k += 1

for i in range(N):
    file.write(str(tela[i].m) + " " + str(tela[i].x) + " " + str(tela[i].y) + " " + str(tela[i].z) + " " + str(tela[i].vx) + " " + str(tela[i].vy) + " " +str(tela[i].vz) + "\n")
file.close


for i in range(N):
    tela[i].vxh = tela[i].vx + tela[i].ax * dt / 2
    tela[i].vyh = tela[i].vy + tela[i].ay * dt / 2
    tela[i].vzh = tela[i].vz + tela[i].az * dt / 2
    
for i in range(N):    
    tela[i].vx += -(tela[i].ax * dt)
    tela[i].vy += -(tela[i].ay * dt)
    tela[i].vz += -(tela[i].az * dt)
    tela[i].Ek = 1/2*tela[i].m*intenzitet(tela[i].vx, tela[i].vy, tela[i].vz)

for i in range(N):
    tela[i].x += tela[i].vxh * dt
    tela[i].y += tela[i].vyh * dt
    tela[i].z += tela[i].vzh * dt
    tela[i].L = tela[i].m * intenzitet(tela[i].vx,tela[i].vy,tela[i].vz) * intenzitet(tela[i].x,tela[i].y,tela[i].z)

t+=dt    

pomEk = 0
for i in range(N):
    pomEk += tela[i].Ek
pomL = 0
for i in range(N):
    pomL += tela[i].L

momentImpulsa.append(pomL)
radijus.append(intenzitet(tela[i].x,tela[i].y,tela[i].z))        
ukupnaEnergija.append(pomEk)
vreme.append(t)

file = open(snapshot_out + "cestica_" + str(k).zfill(4) + ".txt","w")
k += 1

for i in range(N):
    file.write(str(tela[i].m) + " " + str(tela[i].x) + " " + str(tela[i].y) + " " + str(tela[i].z) + " " + str(tela[i].vx) + " " + str(tela[i].vy) + " " +str(tela[i].vz) + "\n")
file.close



while t<=tmax:
    print(t/tmax*100)
    for i in range(N):
        Fx_N_body = 0
        Fy_N_body = 0
        Fz_N_body = 0
        for j in range(N):
            if(i != j):
                r = (intenzitet(tela[i].x - tela[j].x, tela[i].y - tela[j].y, tela[i].z - tela[j].z))
                Fx_N_body += - G*tela[i].m * tela[j].m *(tela[i].x - tela[j].x) / (r**3)   
                Fy_N_body += - G*tela[i].m * tela[j].m *(tela[i].y - tela[j].y) / (r**3)
                Fz_N_body += - G*tela[i].m * tela[j].m *(tela[i].z - tela[j].z) / (r**3)
        
        r = (intenzitet(tela[i].x - pocetno_rast_analiticki_x, tela[i].y - pocetno_rast_analiticki_y, tela[i].z - pocetno_rast_analiticki_z))
        #masaAP = masa(r, tela[i].x, tela[i].y, tela[i].z)        
        Fx_analiticki = - G * tela[i].m * masa(r, tela[i].x, tela[i].y, tela[i].z) * (tela[i].x - pocetno_rast_analiticki_x)/((r)**3)

        Fy_analiticki = - G * tela[i].m * masa(r, tela[i].x, tela[i].y, tela[i].z) * (tela[i].y - pocetno_rast_analiticki_y)/((r)**3)

        Fz_analiticki = - G * tela[i].m * masa(r, tela[i].x, tela[i].y, tela[i].z) * (tela[i].z - pocetno_rast_analiticki_z)/((r)**3) #masa po xyz ili kao radijus
        
        Fx = Fx_N_body + Fx_analiticki
        
        Fy = Fy_N_body + Fy_analiticki 

        Fz = Fz_N_body + Fz_analiticki 

        tela[i].ax = Fx / tela[i].m + dinamickoTrenje(tela[i].m, r, tela[i].vx, tela[i].vy, tela[i].vz, tela[i].vx)
        tela[i].ay = Fy / tela[i].m + dinamickoTrenje(tela[i].m, r, tela[i].vx, tela[i].vy, tela[i].vz, tela[i].vy)
        tela[i].az = Fz / tela[i].m + dinamickoTrenje(tela[i].m, r, tela[i].vx, tela[i].vy, tela[i].vz, tela[i].vz)

        
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
        tela[i].Ek = 1/2*tela[i].m*intenzitet(tela[i].vx, tela[i].vy, tela[i].vz)

    for i in range(N):
        tela[i].x += tela[i].vxh * dt
        x.append(tela[i].x)
        tela[i].y += tela[i].vyh * dt
        y.append(tela[i].y)
        tela[i].z += tela[i].vzh * dt
        z.append(tela[i].z)
        tela[i].L = tela[i].m * intenzitet(tela[i].vx,tela[i].vy,tela[i].vz) * intenzitet(tela[i].x,tela[i].y,tela[i].z)
        

    t+=dt

    pomEk = 0
    for i in range(N):
        pomEk += tela[i].Ek
        
    pomL = 0
    for i in range(N):
        pomL += tela[i].L
    
    momentImpulsa.append(pomL)
        
    radijus.append(intenzitet(tela[i].x,tela[i].y,tela[i].z))        
    ukupnaEnergija.append(pomEk)
    vreme.append(t)
    
    file = open(snapshot_out + "cestica_" + str(k).zfill(4) + ".txt","w")
    k += 1
    for i in range(N):
        file.write(str(tela[i].m) + " " + str(tela[i].x) + " " + str(tela[i].y) + " " + str(tela[i].z) + " " + str(tela[i].vx) + " " + str(tela[i].vy) + " " +str(tela[i].vz) + "\n")
    file.close

for i in range(1,N):
    v1 = intenzitet(tela[i-1].vx,tela[i-1].vy,tela[i-1].vz)    
    v2 = intenzitet(tela[i].vx,tela[i].vy,tela[i].vz)    


#Ek(t)
for i in range(len(vreme)):
    ukupnaEnergija[i] = ukupnaEnergija[i]*1e4
for i in range(len(vreme)):
    vreme[i] = vreme[i]/(86400*365*1e10)

for i in range(len(vreme)):
    momentImpulsa[i] = momentImpulsa[i]*100


file = open(file_out + "momentImpulsa.txt","w")
for i in range(len(vreme)):
    file.write(str(vreme[i]) + " " + str(momentImpulsa[i]) + "\n")
file.close

plt.plot(vreme,momentImpulsa,c='black')
plt.xlabel("$T[Gyr]$")
plt.ylabel("$L[1/10^2Msolkpc^2/s]$")
plt.xlim(min(vreme),max(vreme))
plt.ylim(min(momentImpulsa[1:len(momentImpulsa)-1])-min(momentImpulsa[1:len(momentImpulsa)-1])/100,  max(momentImpulsa)+max(momentImpulsa)/100)
plt.savefig(file_out + "momentImpulsa(t)",dpi = 300)
plt.show()



file = open(file_out + "energija.txt","w")
for i in range(len(vreme)):
    file.write(str(vreme[i]) + " " + str(ukupnaEnergija[i]) + "\n")
file.close

plt.plot(vreme,ukupnaEnergija,c='black')
plt.xlabel("$T[Gyr]$")
plt.ylabel("$Ek[1/10^4Msolkpc^2/s^2]$")
plt.xlim(min(vreme),max(vreme))
plt.ylim(min(ukupnaEnergija[1:len(ukupnaEnergija)-1])-min(ukupnaEnergija[1:len(ukupnaEnergija)-1])/10,max(ukupnaEnergija)+max(ukupnaEnergija)/10)
plt.savefig(file_out + "energija(t)",dpi = 300)
plt.show()


#R(t)
file = open(file_out + "radijus.txt","w")
for i in range(len(vreme)):
    file.write(str(vreme[i]) + " " + str(radijus[i]) + "\n")
file.close

plt.plot(vreme,radijus,c='black')
plt.xlabel("$T[Gyr]$")
plt.ylabel("$R[kpc]$")
plt.xlim(min(vreme),max(vreme))
plt.ylim(min(radijus[1:len(radijus)-1])-min(radijus[1:len(radijus)-1])/10,max(radijus)+max(radijus)/10)
plt.savefig(file_out + "radijus(t)",dpi = 300)
plt.show()

#Vreme
file = open(file_out + "vreme.txt","w")
for i in range(len(vreme)):
    file.write(str(vreme[i]) + "\n")
file.close

#XYZ koorinate cele simulacije
file = open(file_out + "xyz.txt","w")
for i in range(len(x)):
    file.write(str(x[i]) + " " + str(y[i]) + " " + str(z[i]) + "\n")
file.close

plt.plot(x,y,c='black')
plt.xlabel("$X[kpc]$")
plt.ylabel("$Y[kpc]$")
plt.xlim(-250,250)
plt.ylim(-250,250)
plt.savefig(file_out + "xy_ravan",dpi = 300)
plt.show()

#Duzina trajanja programa
stop = timeit.default_timer()
print(stop - start) 
"""
print(ukupnaEnergija[0])
print(momentImpulsa[0])
print(ukupnaEnergija[1])
print(momentImpulsa[1])
"""