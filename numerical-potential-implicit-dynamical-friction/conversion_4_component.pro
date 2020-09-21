N = 467082       ; total number of particles
Ndisk = 108929      ; disk particles
Nbulge = 96247     ; bulge particles
Nhalo = 261905      ; halo particles
Nstars = 1       ;"stars"

unit_gal_mass = 2.325e9  
unit_gad_mass = 1.e10

print, N - Ndisk - Nbulge - Nhalo - Nstars

xx=fltarr(7,N)
openu, 2, 'C:\Users\matij\Desktop\Projekat_2018_dinamicko_trenje\numericalPotential\2Case\pocetniUslovi\combinedGalaxies.txt'
readf, 2, xx
close, 2

massconversion= unit_gal_mass/unit_gad_mass
mpd =xx[0,0]
mpb =xx[0,Ndisk]
mph =xx[0,Ndisk + Nbulge]
mps =xx[0,Ndisk + Nbulge + Nhalo]
print, mpd, mpb, mph, mps
 
tempd =fltarr(7,Ndisk)
tempb =fltarr(7,Nbulge)
temph =fltarr(7,Nhalo)
temps =fltarr(7,Nstars)

print, 'galactICS mass in halo =', mph*Nhalo*unit_gal_mass
print, 'galactICS mass in disk =', mpd*Ndisk*unit_gal_mass
print, 'galactICS mass in bulge=', mpb*Nbulge*unit_gal_mass
print, 'galactICS mass in stars=', mps*Nstars*unit_gal_mass


tempd =xx[*, 0                                :  Ndisk - 1]
tempb =xx[*, Ndisk                            :  Ndisk + Nbulge - 1]
temph =xx[*, Ndisk + Nbulge                   :  Ndisk + Nbulge + Nhalo - 1]
temps =xx[*, Ndisk + Nbulge + Nhalo           :  Ndisk + Nbulge + Nhalo + Nstars - 1]

xx[*, 0                               : Nhalo - 1] = temph
xx[*, Nhalo                           : Nhalo + Ndisk - 1] = tempd
xx[*, Nhalo + Ndisk                   : Nhalo + Ndisk + Nbulge - 1] = tempb
xx[*, Nhalo + Ndisk + Nbulge          : Nhalo + Ndisk + Nbulge + Nstars - 1] = temps

;;;;;;;;;;;;;;;;;;;;;;;;;    create halo + disk + bulge 
;;;;;;;;;;;;;;;;;;;;;;;;;    isolated     

    openw, 4, 'C:\Users\matij\Desktop\Projekat_2018_dinamicko_trenje\numericalPotential\2Case\pocetniUslovi\PUIC_5',/f77_unformatted


npart=lonarr(6)
npart[1]=Nhalo
npart[2]=Ndisk
npart[3]=Nbulge
npart[4]=Nstars

;xi=10.1541     ; in case you want to offset the galaxy
;yi=36.6580
;zi=18.2089

pos=fltarr(3,N)
pos[0,*]=xx[1,*];+xi       ; in case you want to offset the galaxy
pos[1,*]=xx[2,*];+yi       ; add xi, yi, zi
pos[2,*]=xx[3,*];+zi

vel=fltarr(3,N)
vel[0,*]=xx[4,*]*100.
vel[1,*]=xx[5,*]*100.
vel[2,*]=xx[6,*]*100.

;     d=sqrt((x0-pos[0,*])^2.+(y0-pos[1,*])^2.+(z0-pos[2,*])^2.)
;     xort=(x0-pos[0,*])/d
;     yort=(y0-pos[1,*])/d
;     zort=(z0-pos[2,*])/d

;     kick=1000.

;     vxkick=kick*xort
;     vykick=kick*yort
;     vzkick=kick*zort

;     vxkick=-4.83452
;     vykick=-62.0415
;     vzkick=-39.0725
;     vel[0,*]=vel[0,*];+vxkick   ; if you want to kick the galaxy
;     vel[1,*]=vel[1,*];+vykick
;     vel[2,*]=vel[2,*];+vzkick

massarr=dblarr(6)

massarr[1]=mph*massconversion
massarr[2]=mpd*massconversion
massarr[3]=mpb*massconversion
massarr[4]=mps*massconversion
print, 'gadget mass in halo=', massarr[1]*Nhalo*unit_gad_mass 
print, 'gadget mass in disk=', massarr[2]*Ndisk*unit_gad_mass 
print, 'gadget mass in bulge=', massarr[3]*Nbulge*unit_gad_mass 
print, 'gadget mass in stars=', massarr[4]*Nstars*unit_gad_mass

id=lindgen(N)

    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)
    npartTotal[1]=Nhalo
    npartTotal[2]=Ndisk
    npartTotal[3]=Nbulge
    npartTotal[4]=Nstars
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)

    writeu,4,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la
    writeu,4,pos
    writeu,4,vel
    writeu,4,id

    close, 4
print,'end'
;;;;;;;;;;;;;;;;;;;;;   create disk + bulge 
;;;;;;;;;;;;;;;;;;;;;   isolated


;;     openw, 5, 'ndb.iso.ics',/f77_unformatted

;;     npart=lonarr(6)
;;     npart[2]=Nd
;;     npart[3]=Nb
;;     massarr=dblarr(6)
;;     massarr[2]=mpd*massconversion
;;     massarr[3]=mpb*massconversion
;;     time=0.0D
;;     redshift=0.0D
;;     flag_sfr=0L
;;     flag_feedback=0L
;;     npartTotal=lonarr(6)
;;     npartTotal[2]=Nd
;;     npartTotal[3]=Nb
;;     bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
;;     la=intarr(bytesleft/2)

;;     writeu,5,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la
;;     N=300000
;;     id=lindgen(N)
;;     writeu,5,pos[*,Nh:Nh+Nd+Nb-1]
;;     writeu,5,vel[*,Nh:Nh+Nd+Nb-1]
;;     writeu,5,id


;;     close, 5


;!P.MULTI=[0,2,0]
;collist=[255*256L*256L+100*256L+100, 255*256L, 255L, 255*256L+255L, 255*256L*256L+255L]
 ;  window,xsize=2000,ysize=2000

;plot, pos[0,0:Nh-1], pos[1,0:Nh-1], psym=3, xrange=[-300,300], yrange=[-300,300], color=collist(1 mod n_elements(collist))
;oplot, pos[0,Nh:Nh+Nd-1], pos[1,Nh:Nh+Nd-1], psym=3, color=collist(2 mod n_elements(collist))
;oplot, pos[0,Nh+Nd:Nh+Nd+Nb-1], pos[1,Nh+Nd:Nh+Nd+Nb-1], psym=3, color=collist(3 mod n_elements(collist))

;plot, pos[0,0:Nh-1], pos[2,0:Nh-1], psym=3, xrange=[-20,20], yrange=[-20,20], color=collist(1 mod n_elements(collist))
;plot, pos[0,Nh:Nh+Nd-1], pos[2,Nh:Nh+Nd-1], psym=3, color=collist(2 mod n_elements(collist))
;oplot, pos[0,Nh+Nd:Nh+Nd+Nb-1], pos[1,Nh+Nd:Nh+Nd+Nb-1], psym=3, color=collist(3 mod n_elements(collist))



;;;DISK ONLY
;plot, pos[0,Nh:Nh+Nd-1], pos[1,Nh:Nh+Nd-1], psym=3, xrange=[-200,200], yrange=[-200,200], color=collist(1 mod n_elements(collist))
;plot, pos[0,Nh:Nh+Nd-1], pos[2,Nh:Nh+Nd-1], psym=3, xrange=[-200,200], yrange=[-200,200], color=collist(1 mod n_elements(collist))


;;;BULGE ONLY
;plot, pos[0,Nh+Nd:Nh+Nd+Nb-1], pos[1,Nh+Nd:Nh+Nd+Nb-1], psym=3, xrange=[-20,20], yrange=[-20,20], color=collist(1 mod n_elements(collist))
;plot, pos[0,Nh+Nd:Nh+Nd+Nb-1], pos[1,Nh+Nd:Nh+Nd+Nb-1], psym=3, xrange=[-20,20], yrange=[-20,20], color=collist(1 mod n_elements(collist))


print,'END'
;breaking
end
