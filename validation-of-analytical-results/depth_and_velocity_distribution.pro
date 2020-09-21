;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example showing how the file format of Gadget can    ;;
;; be read-in in IDL                                    ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

frun="C:/Users/Pedja/Desktop/Snapshot_ultimate";Pedja/Projekat17-18/Galaksije/Sudari/Snapshot_ultimate"  

window,xsize=600,ysize=333
!P.multi=[0,3,1]

xlen=20


for num=00,00 do begin


    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname=frun+"/snapshot_"+exts 
    fname=strcompress(fname,/remove_all)

    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)


    print,fname
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
    NDisk=npart(2)
    NBulge=npart(3)
    NStars=npart(4)

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
        xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
        xgas(*)=pos(0,0:Ngas-1)
        ygas(*)=pos(1,0:Ngas-1)
        zgas(*)=pos(2,0:Ngas-1)
        if massarr(0) eq 0 then begin
            mgas(*)=mass(0:Ngas-1)	
        endif else begin
            mgas(*)= massarr(0)
	endelse
    endif

    if Nhalo gt 0 then begin
        xhalo=fltarr(NHalo) &  yhalo=fltarr(Nhalo) & zhalo=fltarr(Nhalo) & mhalo=fltarr(Nhalo)
        xhalo(*)=pos(0,0+Ngas:Nhalo+Ngas-1)
        yhalo(*)=pos(1,0+Ngas:Nhalo+Ngas-1)
        zhalo(*)=pos(2,0+Ngas:Nhalo+Ngas-1)
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
        xdisk=fltarr(NDisk) &  ydisk=fltarr(NDisk) &  zdisk=fltarr(NDisk) & mdisk=fltarr(NDisk)
        xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
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
        xbulge=fltarr(NBulge) &  ybulge=fltarr(NBulge) &  zbulge=fltarr(NBulge) & mbulge=fltarr(NBulge)
        xbulge(*)=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        ybulge(*)=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        zbulge(*)=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
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
        xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & mstars=fltarr(Nstars)
        xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        if massarr(4) eq 0 then begin
	    skip=0L
            for t=0,2 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mstars(*)=mass(0+skip:Nstars-1+skip)	
        endif else begin
            mstars(*)= massarr(4)
	endelse
    endif
    
    alpha = (!dpi/180)*(77)
    beta  = (!dpi/180)*(37) 
    
   xdisk0 = xdisk
   ydisk0 = ydisk
   zdisk0 = zdisk

   xdisk = xdisk0*cos(beta) - ydisk0*cos(alpha)*sin(beta) + zdisk0*sin(alpha)*sin(beta)
   ydisk = xdisk0*sin(beta) + ydisk0*cos(alpha)*cos(beta) - zdisk0*sin(alpha)*cos(beta)
   zdisk = ydisk0*sin(alpha) + zdisk0*cos(alpha)
 
    ;xbulge = xbulge*cos(beta) - ybulge*cos(alpha)*sin(beta) + zbulge*sin(alpha)*sin(beta)
    ;ybulge = xbulge*sin(beta) + ybulge*cos(alpha)*cos(beta) - zbulge*sin(alpha)*cos(beta)
    ;zbulge = ybulge*sin(alpha) + zbulge*cos(alpha)
 
   ;xhalo = xhalo - 108.779
   ;yhalo = yhalo + 208.902
   ;zhalo = zhalo + 93.096
    
   xhalo0 = xhalo
   yhalo0 = yhalo
   zhalo0 = zhalo
    
    xhalo = xhalo0*cos(beta) - yhalo0*cos(alpha)*sin(beta) + zhalo0*sin(alpha)*sin(beta)
    yhalo = xhalo0*sin(beta) + yhalo0*cos(alpha)*cos(beta) - zhalo0*sin(alpha)*cos(beta) 
    zhalo = yhalo0*sin(alpha) + zhalo0*cos(alpha)

    ;xstars = xstars*cos(beta) - ystars*cos(alpha)*sin(beta) + zstars*sin(alpha)*sin(beta)
    ;ystars = xstars*sin(beta) + ystars*cos(alpha)*cos(beta) - zstars*sin(alpha)*cos(beta)
    ;zstars = ystars*sin(alpha) + zstars*cos(alpha)

    
   zeta= (xdisk/(zdisk+784))*(180/!dpi)
   eta=(ydisk/(zdisk+784))*(180/!dpi)

   zeta1=(xhalo/(zhalo+784))*(180/!dpi)
   eta1=(yhalo/(zhalo+784))*(180/!dpi)


   x1=[2.015, 1.745, 1.483, 1.226, 0.969, 0.717, 0.467, 0.219]
   y1=[-3.965, -3.525, -3.087, -2.653, -2.264, -1.768, -1.327,-0.886]
   x1err=[0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33]
   y1err=[0.22,0.22,0.22,0.22,0.22,0.22,0.22,0.22]

exts1='00'
    exts1=exts1+strcompress(string(num),/remove_all)
    exts1=strmid(exts1,strlen(exts1)-2,2)
    
;fnamezeta = 'C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdzeta' + exts1 + '.txt'
;fnameeta   = 'C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Orijentacija/lsdeta' + exts1 + '.txt'
;print, fnamezeta
;print, fnameeta 
 
 
;openw, 2, fnamezeta, WIDTH = 16 ;ubacujemo sve pozicije u txt
;    printf, 2, zeta1
;    close, 2
;openw, 3, fnameeta, WIDTH = 16 ;ubacujemo sve pozicije u txt
;    printf, 3, eta1
;    close, 3

;fnameid = 'C:/Users/Pedja/Desktop/Projekat2018/Modeli/Model1/Rezultati/Centar_mase/id' + exts1 + '.txt'
;print, fnameid
;openw, 5, fname1, WIDTH = 16
;printf, 5, id
;close, 5

   !P.MULTI=[0,1,0]
   collist=[255*256L*256L+100*256L+100, 255*256L, 255L, 255*256L+255L, 255*256L*256L+255L]
    ;;;plotujemo Andromedu
    plot, zeta, eta, psym=3, xrange= [xlen,-xlen], yrange= [-xlen,xlen], xstyle=1,ystyle=1, xtitle='x1[stepeni]',ytitle='y1[stepeni]', color=collist(3 mod n_elements(collist))
    
    ;;overplotujemo dwarf
    oplot, zeta1, eta1, psym=3, color=collist(2 mod n_elements(collist)) ;, xrange= [xlen, -xlen], yrange= [-xlen,xlen], xstyle=1,ystyle=1 
   
     ;oplot, x1, y1, psym=3
;     errplot, x1, y1+0.22, y1-0.22;, psym=3, color='red'
         
    wait, 0.5	

endfor	
print, 'gotovo'
end
