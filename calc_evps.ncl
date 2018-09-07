;;; NCL script to calculate evaporation rate from hfls & tas

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


;; The following variables should be specified as command-line arguments
;; e.g.: ncl file=\"$file\" script.ncl

fint = addfile(tasfile, "r")
finh = addfile(hflsfile, "r")

system("rm -f "+outfile)
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


; use hfls file as template

; global attributes

att_names = getvaratts(finh) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = finh@$att_names(i)$
end do

history = "Created " + systemfunc("date") + " by "+systemfunc("whoami")+"@"+systemfunc("hostname")+" using NCL script from source files "+hflsfile+" and "+tasfile

fout@history = history

; copy variables

var_names = getfilevarnames (finh) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. "hfls") then   
  fout->$var_names(i)$ = finh->$var_names(i)$
end if
end do




; calculate evaporation


print("Units should be W m-2 and K")
print("Units are: "+finh->hfls@units+", "+fint->tas@units)


;; Because there's not enough room to fit everything in memory at
;; once, we loop on time


;; initialize this way to get correct structure
evps = finh->hfls*0


nt = dimsizes(finh->time)
do i=0, nt-1

  t = fint->tas(i,:,:) - 273.15
  h = finh->hfls(i,:,:)

  ;; calculate latent heat of condensation for water in J/kg:

  ;; Lcw(T) = 6.14342e-5 T^3 + 1.58927e-3 T^2 - 2.36418 T + 2500.79 J/g
  ;; T = temperature in degrees C.  Formula valid for -40C to +40 C.

  a0 = 2500.79
  a1 = -2.36418
  a2 = 1.58927e-3
  a3 = 6.14342e-5

  result = a0 + a1*t
  t = t*t
  result = result + a2 * t
  t = t*t
  result = result + a3 * t
  result = result * 1000    ;; J/g -> J/kg

  ;; divide hfls (W/m^2) by Lcw (J/kg) to get evap rate (kg/m^2/s)
  result = h / result

  evps(i,:,:) = (/result/)

end do


xdim = finh->hfls!2
evps!2 = xdim
evps&$xdim$=finh->$xdim$

ydim = finh->hfls!1
evps!1 = ydim
evps&$ydim$=finh->$ydim$

evps!0 = "time"

evps&time = finh->time

evps@units = "kg m-2 s-1"
evps@long_name = "Surface Evaporation Rate (estimated)"
evps@standard_name = "water_evaporation_flux"
evps@missing_value = 1.e+20
evps@_FillValue = 1.e+20
evps@coordinates = h@coordinates
evps@grid_mapping = h@grid_mapping

fout->evps = evps

exit

;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
