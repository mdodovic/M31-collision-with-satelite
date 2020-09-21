;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example showing how the file format of Gadget can    ;;
;; be read-in in IDL                                    ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    
alpha = (!dpi/180)*(77)
beta  = (!dpi/180)*(37) 
;print,alpha
;print,betta   
xdisk0 = -84.41
ydisk0 = 152.47
zdisk0 = -97.08

xdisk = xdisk0*cos(beta) - ydisk0*cos(alpha)*sin(beta) + zdisk0*sin(alpha)*sin(beta)
ydisk = xdisk0*sin(beta) + ydisk0*cos(alpha)*cos(beta) - zdisk0*sin(alpha)*cos(beta)
zdisk = ydisk0*sin(alpha) + zdisk0*cos(alpha)

;print,xdisk
;print,ydisk
;print,zdisk
      
zeta= (xdisk/(zdisk+784))*(180/!dpi)
eta=  (ydisk/(zdisk+784))*(180/!dpi)

print,zeta  
print,eta

end