import numpy as np
import matplotlib.pyplot as plt

def intenzitet(a,b,c):
    return (a*a + b*b + c*c)**0.5

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
    
    def __init__(self,x,y,z,vx,vxh,vy,vyh,vz,vzh,ax,ay,az,m, Ek,Ep):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vxh = vxh
        self.vy = vy
        self.vyh = vyh
        self.vz = vz
        self.vzh = vzh
        self.ax = ax
        self.ay = ay
        self.az = az
        self.m = m
        self.Ek = Ek
        self.Ep = Ep

 
tela = []        
N = 2

tela.append(telo(1.5e11,0,0,0,0,2*3.14159*1.5e11/(365*24*3600),0,0,0,0,0,0,6e24,0,0))
tela.append(telo(0,0,0,0,0,0,0,0,0,0,0,0,2e30,0,0))        
x = []
y = []
z = []
vx = []
vy = []
vz = []

Kineticka_energija = []
Potencijalna_energija = []
Ukupna_energija = []
Vreme = []
G = 6.67e-11 #m^3*kg^-1*s^-2
dt = 86400#s
t=0
tmax=86400*365*10#s


ek = 0
ep = 0
for i in range(N):
    tela[i].Ek = 1/2*tela[i].m*(intenzitet(tela[i].vx,tela[i].vy,tela[i].vz))**2
    epPom = 0
    for j in range(N):
        if i != j:
            epPom += -G*tela[i].m*tela[j].m/(intenzitet(tela[i].x - tela[j].x,tela[i].y - tela[j].y, tela[i].z - tela[j].z))
             
    tela[i].Ep = epPom   
    ek += tela[i].Ek
    ep += tela[i].Ep
    
Kineticka_energija.append(ek)
Potencijalna_energija.append(ep)
Ukupna_energija.append(ek+ep)
Vreme.append(t)


for i in range(N):
    tela[i].vxh = tela[i].vx + tela[i].ax * dt / 2
    tela[i].vyh = tela[i].vy + tela[i].ay * dt / 2
    tela[i].vzh = tela[i].vz + tela[i].az * dt / 2
    
for i in range(N):    
    tela[i].vx += (tela[i].ax * dt)
    vx.append(tela[i].vx)
    tela[i].vy += (tela[i].ay * dt)
    vy.append(tela[i].vy)        
    tela[i].vz += (tela[i].az * dt)
    vz.append(tela[i].vz)

for i in range(N):
    tela[i].x += tela[i].vxh * dt
    x.append(tela[i].x)
    tela[i].y += tela[i].vyh * dt
    y.append(tela[i].y)
    tela[i].z += tela[i].vzh * dt
    z.append(tela[i].z)

while t<=tmax:
    print(t/tmax*100)
    
    for i in range(N):
        Fx = 0
        Fy = 0
        Fz = 0
        for j in range(N):
            if(i != j):
                r = (intenzitet(tela[i].x - tela[j].x,tela[i].y - tela[j].y, tela[i].z - tela[j].z))
                Fx += - G*tela[i].m * tela[j].m *(tela[i].x - tela[j].x) / (r**3)   
                Fy += - G*tela[i].m * tela[j].m *(tela[i].y - tela[j].y) / (r**3)
                Fz += - G*tela[i].m * tela[j].m *(tela[i].z - tela[j].z) / (r**3)
        tela[i].ax = Fx / tela[i].m
        tela[i].ay = Fy / tela[i].m
        tela[i].az = Fz / tela[i].m


    ek = 0
    ep = 0
    for i in range(N):
        tela[i].Ek = 1/2*tela[i].m*(intenzitet(tela[i].vx,tela[i].vy,tela[i].vz))**2
        epPom = 0        
        for j in range(N):
            if i != j:
                epPom += -G*tela[i].m*tela[j].m/(intenzitet(tela[i].x - tela[j].x,tela[i].y - tela[j].y, tela[i].z - tela[j].z))
        
        tela[i].Ep = epPom
        ek += tela[i].Ek
        ep += tela[i].Ep
        
    Kineticka_energija.append(ek)
    Potencijalna_energija.append(ep)
    Ukupna_energija.append(ek+ep)

        
    for i in range(N):
        tela[i].vxh += tela[i].ax*dt
        tela[i].vyh += tela[i].ay*dt
        tela[i].vzh += tela[i].az*dt
        
    for i in range(N):
        tela[i].vx = tela[i].vxh + tela[i].ax * dt / 2
        vx.append(tela[i].vx)
        tela[i].vy = tela[i].vyh + tela[i].ay * dt / 2
        vy.append(tela[i].vy)        
        tela[i].vz = tela[i].vzh + tela[i].az * dt / 2
        vz.append(tela[i].vz)



    for i in range(N):
        tela[i].x += tela[i].vxh * dt
        x.append(tela[i].x)
        tela[i].y += tela[i].vyh * dt
        y.append(tela[i].y)
        tela[i].z += tela[i].vzh * dt
        z.append(tela[i].z)

    t+=dt

    Vreme.append(t)

for i in range(len(x)):
    x[i]=x[i]/1e11
    y[i]=y[i]/1e11

for i in range(len(Ukupna_energija)):
    Ukupna_energija[i] = Ukupna_energija[i]/1e33
    Vreme[i] = Vreme[i] / 86400/365      

putanja_sacuvaj = "C:/Users/matij/Desktop/Projekat_2018_dinamicko_trenje/analiticki_potencijal/leap_frog/"

plt.scatter(x,y,s = 0.1)
plt.xlabel("$X[AJ]$")
plt.ylabel("$Y[AJ]$")
plt.savefig(putanja_sacuvaj + "xy_ravan",dpi=300)
plt.show()

plt.plot(Vreme,Ukupna_energija)
plt.xlabel("$T[Godina]$")
plt.ylabel("$E[1e33 J]$")
plt.ylim(-8.04,-7.9)
plt.savefig(putanja_sacuvaj + "E(t)",dpi=300)
plt.show()

