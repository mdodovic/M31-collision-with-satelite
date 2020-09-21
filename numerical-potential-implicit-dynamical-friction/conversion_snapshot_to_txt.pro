;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example showing how the file format of Gadget can    ;;
;; be read-in in IDL                                    ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

frun='C:\Users\matij\Desktop\Projekat_2018_dinamicko_trenje\numericalPotential\2Case\snapshotovi\' 
fout_komp='C:\Users\matij\Desktop\Projekat_2018_dinamicko_trenje\numericalPotential\2Case\txtFilesUnconverted\'
;fout_mass='C:/Users/matij/Desktop/Bliksi_prolaz//'

unit_gal_mass=2.325e9  
unit_gad_mass=1.e10
massconversion=unit_gad_mass  ; kada se pomnozi sa ovim dobija se da je masa u Masama sunca

                             
N = 61
                           
for num=0,N do begin


    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname=frun+"snapshot_"+exts
    fname=strcompress(fname,/remove_all)
    
    fopennamed = fout_komp +"disk_"+exts+".txt"
    fopennameb = fout_komp +"bulge_"+exts+".txt"
    fopennameh = fout_komp +"halo_"+exts+".txt"
    fopennames = fout_komp + "stars_" + exts+".txt"

;    fopennamemd = fout_mass +"mdisk_"+exts+".txt"
;    fopennamemb = fout_mass +"mbulge_"+exts+".txt"
;    fopennamemh = fout_mass +"mhalo_"+exts+".txt"
;    fopennamems = fout_mass +"mstars_"+exts+".txt"
   
    
    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)
    
    openr,1,fname,/f77_unformatted
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la
    
    
    print,"Time= ", time

    N=total(npart)
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
 
    ind=where((npart gt 0) and (massarr eq 0)) 
    if ind(0) ne -1 then begin
        Nwithmass= total(npart(ind))
        mass=fltarr(Nwithmass)
    endif else begin	
        Nwithmass= 0
    endelse

    readu,1,pos
    readu,1,vel
    readu,1,id
    if Nwithmass gt 0 then begin
      readu,1,mass
    endif

    NGas=npart(0)
    NHalo=npart(1)
    Ndisk=npart(2)
    NBulge=npart(3)
    NStars=npart(4)
    NBndry=npart(5)
    
    if Ngas gt 0 then begin
        u=fltarr(Ngas)
        readu,1,u

        rho=fltarr(Ngas)
        readu,1,rho

        if flag_sfr gt 0 then begin
            sfr=fltarr(Ngas)
            mfs=fltarr(Ngas)
            readu,1,sfr
            readu,1,mfs
        endif
    endif
    close,1

    if Ngas gt 0 then begin
        xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas)
        vxgas=fltarr(Ngas) &  vygas=fltarr(Ngas) &  vzgas=fltarr(Ngas)
        mgas=fltarr(Ngas)
        xgas(*)=pos(0,0:Ngas-1)
        ygas(*)=pos(1,0:Ngas-1)
        zgas(*)=pos(2,0:Ngas-1)
        vxgas(*)=vel(0,0:Ngas-1)
        vygas(*)=vel(1,0:Ngas-1)
        vzgas(*)=vel(2,0:Ngas-1)
        if massarr(0) eq 0 then begin
            mgas(*)=mass(0:Ngas-1)	
        endif else begin
            mgas(*)= massarr(0)
	endelse
    endif

    if Nhalo gt 0 then begin
        xhalo=fltarr(Nhalo) &  yhalo=fltarr(Nhalo)  & zhalo=fltarr(Nhalo)
        vxhalo=fltarr(Nhalo) &  vyhalo=fltarr(Nhalo) &  vzhalo=fltarr(Nhalo)
        mhalo=fltarr(Nhalo)
        xhalo(*)=pos(0,0+Ngas:Nhalo+Ngas-1)
        yhalo(*)=pos(1,0+Ngas:Nhalo+Ngas-1)
        zhalo(*)=pos(2,0+Ngas:Nhalo+Ngas-1)
        vxhalo(*)=vel(0,0+Ngas:Nhalo+Ngas-1)
        vyhalo(*)=vel(1,0+Ngas:Nhalo+Ngas-1)
        vzhalo(*)=vel(2,0+Ngas:Nhalo+Ngas-1)
        if massarr(1) eq 0 then begin
	    skip=0L
            for t=0,0 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mhalo(*)=mass(0+skip:Nhalo-1+skip)	
        endif else begin
            mhalo(*)= massarr(1)
	endelse
    endif

    if Ndisk gt 0 then begin
        xdisk=fltarr(Ndisk) &  ydisk=fltarr(Ndisk)  & zdisk=fltarr(Ndisk)
        vxdisk=fltarr(Ndisk) &  vydisk=fltarr(Ndisk) &  vzdisk=fltarr(Ndisk)
        mdisk=fltarr(Ndisk)
        
        xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vxdisk(*)=vel(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vydisk(*)=vel(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vzdisk(*)=vel(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        if massarr(2) eq 0 then begin
	    skip=0L
            for t=0,1 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mdisk(*)=mass(0+skip:Ndisk-1+skip)	
        endif else begin
            mdisk(*)= massarr(2)
	endelse
    endif

    if Nbulge gt 0 then begin
        xbulge=fltarr(Nbulge) &  ybulge=fltarr(Nbulge)  & zbulge=fltarr(Nbulge)
        vxbulge=fltarr(Nbulge) &  vybulge=fltarr(Nbulge) &  vzbulge=fltarr(Nbulge)
        mbulge=fltarr(Nbulge)
        xbulge(*)=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        ybulge(*)=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        zbulge(*)=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        vxbulge(*)=vel(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        vybulge(*)=vel(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        vzbulge(*)=vel(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        if massarr(3) eq 0 then begin
	    skip=0L
            for t=0,2 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mbulge(*)=mass(0+skip:Nbulge-1+skip)	
        endif else begin
            mbulge(*)= massarr(3)
	endelse
    endif


    if Nstars gt 0 then begin
        xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)
        vxstars=fltarr(Nstars) &  vystars=fltarr(Nstars) &  vzstars=fltarr(Nstars)
        mstars=fltarr(Nstars)
        
        xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        vxstars(*)=vel(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        vystars(*)=vel(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        vzstars(*)=vel(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)

        if massarr(4) eq 0 then begin
	    skip=0L
            for t=0,3 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mstars(*)=mass(0+skip:Nstars-1+skip)	
        endif else begin
            mstars(*)= massarr(4)
	endelse
    endif

    if Nbndry gt 0 then begin
      xbndry=fltarr(Nbndry) &  ybndry=fltarr(Nbndry)  & zbndry=fltarr(Nbndry)
      vxbndry=fltarr(Nbndry) &  vybndry=fltarr(Nbndry) &  vzbndry=fltarr(Nbndry)
      mbndry=fltarr(Nbndry)
      
      xbndry(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+Nbndry-1)
      ybndry(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+Nbndry-1)
      zbndry(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+Nbndry-1)
      vxbndry(*)=vel(0,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+Nbndry-1)
      vybndry(*)=vel(1,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+Nbndry-1)
      vzbndry(*)=vel(2,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+Nbndry-1)

      if massarr(5) eq 0 then begin
        skip=0L
        for t=0,4 do begin  
            if (npart(t) gt 0) and (massarr(t) eq 0) then begin
              skip=skip + npart(t)
            endif
        endfor
        mbndry(*)=mass(0+skip:Nbndry-1+skip)  
      endif else begin
        mbndry(*)= massarr(5)
      endelse
    endif
   
    if Ndisk gt 0 then begin
      openw, 2, fopennamed, WIDTH=2500
      nizd = fltarr(7, Ndisk)
      nizd[0,*]= mdisk*massconversion
      nizd[1,*]= xdisk
      nizd[2,*]= ydisk
      nizd[3,*]= zdisk
      nizd[4,*]= vxdisk
      nizd[5,*]= vydisk;
      nizd[6,*]= vzdisk
      printf,2,nizd
      close, 2
    endif

    if Nbulge gt 0 then begin
      openw, 4, fopennameb, WIDTH=2500
      nizb = fltarr(7, Nbulge)
      nizb[0,*]= mbulge*massconversion
      nizb[1,*]= xbulge
      nizb[2,*]= ybulge
      nizb[3,*]= zbulge
      nizb[4,*]= vxbulge;
      nizb[5,*]= vybulge
      nizb[6,*]= vzbulge
      printf, 4, nizb
      close, 4
    endif

    if Nhalo gt 0 then begin
      openw, 6, fopennameh, WIDTH = 2500
      nizh = fltarr(7, Nhalo)
      nizh[0,*]= mhalo*massconversion
      nizh[1,*]= xhalo
      nizh[2,*]= yhalo
      nizh[3,*]= zhalo
      nizh[4,*]= vxhalo
      nizh[5,*]= vyhalo
      nizh[6,*]= vzhalo 
      printf, 6, nizh
      close, 6
    endif
    
    if Nstars gt 0 then begin
      openw, 8, fopennames, WIDTH=2500
      nizs = fltarr(7, Nstars)
      nizs[0,*]= mstars*massconversion
      nizs[1,*]= xstars
      nizs[2,*]= ystars
      nizs[3,*]= zstars
      nizs[4,*]= vxstars
      nizs[5,*]= vystars
      nizs[6,*]= vzstars 
      printf, 8, nizs
      close, 8
    endif

    if Nbndry gt 0 then begin
      openw, 6, fopennamebnd, WIDTH = 2500
      nizbnd = fltarr(7, Nbndry)
      nizbnd[0,*]= mbndry*massconversion
      nizbnd[1,*]= xbndry
      nizbnd[2,*]= ybndry
      nizbnd[3,*]= zbndry
      nizbnd[4,*]= vxbndry
      nizbnd[5,*]= vybndry
      nizbnd[6,*]= vzbndry 
      printf, 6, nizbnd
      close, 6
    endif

  
endfor	

print, "The end of convert"

end

