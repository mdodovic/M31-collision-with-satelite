import math as math
import numpy as np
import matplotlib.pyplot as plt
import timeit

start = timeit.default_timer()

#Linkovi

file_out = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiBezTrenjaSveKomponente/"
snapshot_out = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/analitickiBezTrenjaSveKomponente/snapshotovi/"

#Konstante
G = 4.515353945657138e-39 # [kpc^3/(Msol*s^2)]
kriticnaGustina = 277.72/(0.71*0.71) # [Msol/kpc^2]

# Podaci o analitickom potencijalu

Mbulge = 3.2e9 # [Msol]
Mdisk = 3.7e9 # [Msol]
Mhalo = 88e10 # [Msol]
Map = Mdisk + Mbulge + Mhalo # [Msol]




scale_length_halo = 8.18 # [kpc]
scale_length_bulge = 1.22 # [kpc]
scale_length_disk_vertikalni = 0.57 # [kpc]
scale_length_disk_radijalni = 6.82 # [kpc]

delta_c = 27e4 # bezdimenzioni parametar
Sigma_0 = Mdisk/(2*np.pi*scale_length_disk_radijalni**2) # centralna gustina diska


#Pocetni uslovi vremena simulacije

dt = 86400*365*1e6 # milion godina
dt = dt * 50 # za 50 miliona godina 
t=0
k = 0 # broj snapshota
tmax=86400*365*1e9 # milijardu godina
n = 11
tmax = tmax * n # za n milijardi godina




def intenzitet(a,b,c):
    return (a*a + b*b + c*c)**0.5

def masad(r,z):
    M = Mdisk/(2*scale_length_disk_radijalni)*np.tanh(z/scale_length_disk_vertikalni)*(scale_length_disk_radijalni*(1 - (np.e)**(-r/scale_length_disk_radijalni)) - r*(np.e)**(-r/scale_length_disk_radijalni)) 
    return M

def masab(r):
    M = Mbulge*(r**2)/((scale_length_bulge+r)**2)
    return M

def masah(r):
    M = 4*math.pi*delta_c*kriticnaGustina*(scale_length_halo**3)*(math.log((r+scale_length_halo)/scale_length_halo) - (r/(r+scale_length_halo)))
    return M

def masa(r,z):
    return masah(r) + masab(r) + masad(r,z)


def potencijalH(r):
    FI = -4*np.pi*G*delta_c*kriticnaGustina*(scale_length_halo**3)/r*np.log((r+scale_length_halo)/scale_length_halo)
    return FI
def potencijalD(r):
    FI = -2*np.pi*G*Sigma_0*(scale_length_disk_radijalni**2)*(1-((np.e)**(-r/scale_length_disk_radijalni)))/r
    return FI
def potencijalB(r):
    FI = -G*Mbulge/(r+scale_length_bulge)    
    return FI    
    
    
def potencijal(r):
    return potencijalH(r) + potencijalD(r) + potencijalB(r)

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

M = 4.42e10
Ek = 1/2*M*vk**2
Ep = potencijal(R)*M
L = M*R*vk
pocetno_rast_analiticki_x = 0
pocetno_rast_analiticki_y = 0
pocetno_rast_analiticki_z = 0

tela.append(telo(R,0,0,0,vk,0,0,0,0,M,Ek,Ep,L))

x = []
y = []
z = []

ukupnaEnergija = []
momentImpulsa = []
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
    tela[i].vx += (tela[i].ax * dt)
    tela[i].vy += (tela[i].ay * dt)
    tela[i].vz += (tela[i].az * dt)
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
                r = (intenzitet(tela[i].x - tela[j].x,tela[i].y - tela[j].y, tela[i].z - tela[j].z))
                Fx_N_body += - G*tela[i].m * tela[j].m *(tela[i].x - tela[j].x) / (r**3)   
                Fy_N_body += - G*tela[i].m * tela[j].m *(tela[i].y - tela[j].y) / (r**3)
                Fz_N_body += - G*tela[i].m * tela[j].m *(tela[i].z - tela[j].z) / (r**3)
        
        r = (intenzitet(tela[i].x - pocetno_rast_analiticki_x, tela[i].y - pocetno_rast_analiticki_y, tela[i].z - pocetno_rast_analiticki_z))

        Fx_analiticki = - G * tela[i].m * masa(r,tela[i].z) * (tela[i].x - pocetno_rast_analiticki_x)/((r)**3)

        Fy_analiticki = - G * tela[i].m * masa(r,tela[i].z) * (tela[i].y - pocetno_rast_analiticki_y)/((r)**3)

        Fz_analiticki = - G * tela[i].m * masa(r,tela[i].z) * (tela[i].z - pocetno_rast_analiticki_z)/((r)**3) #masa po xyz ili kao radijus
        
        Fx = Fx_N_body + Fx_analiticki
        Fy = Fy_N_body + Fy_analiticki
        Fz = Fz_N_body + Fz_analiticki


        tela[i].ax = Fx / tela[i].m
        tela[i].ay = Fy / tela[i].m
        tela[i].az = Fz / tela[i].m
        
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
        tela[i].vx = (tela[i].vxh + tela[i].ax * dt / 2)
        tela[i].vy = (tela[i].vyh + tela[i].ay * dt / 2)
        tela[i].vz = (tela[i].vzh + tela[i].az * dt / 2)
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

