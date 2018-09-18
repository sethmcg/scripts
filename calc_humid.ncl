;;; NCL script to calculate relative humidity and dewpoint from huss, tas, and ps

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Command-line arguments:
;; psfile       required: name of input file with pressure
;; tasfile      required: name of input file with temperature
;; hussfile     required: name of input file with specific humidity
;; hursfile     required: name of output file for relative humidity
;; dewptfile    required: name of output file for dewpoint temperature
;; verbose      boolean option: print progress indicators?


if (.not. isvar("verbose")) then
  verbose = False
end if


finp = addfile(psfile, "r")
fint = addfile(tasfile, "r")
finh = addfile(hussfile, "r")

system("rm -f "+hursfile+" "+dewptfile)
fouth = addfile(hursfile, "c")
foutd = addfile(dewptfile, "c")
filedimdef(fouth,"time",-1,True)  ;; make time dimension unlimited
filedimdef(foutd,"time",-1,True)  ;; make time dimension unlimited

;; use tas as template, since it has a reference height

; global attributes


att_names = getvaratts(fint) 

do i = 0,dimsizes(att_names)-1
  fouth@$att_names(i)$ = fint@$att_names(i)$
  foutd@$att_names(i)$ = fint@$att_names(i)$
end do


;; To track the full history of the file, we need histories of all
;; the ancestors.

BR = tochar(10)

history = "###################################"+BR
history = history + "History of parent file " + hussfile + ":"+BR
history = history + finh@history+BR
history = history + "-------------------------"+BR
history = history + "History of parent file " + psfile + ":"+BR
history = history + finp@history+BR
history = history + "-------------------------"+BR
history = history + "History of parent file " + tasfile + ":"+BR
history = history + fint@history+BR
history = history + "#########################"+BR
history = history + "Created " + systemfunc("date")
history = history + " by "+systemfunc("whoami")+"@"+systemfunc("hostname")
history = history + " using NCL script calc_humid.ncl from parent files listed above."+BR

fouth@history = history
foutd@history = history

;; update unique tracking_ids
fouth@tracking_id = systemfunc("uuidgen")
foutd@tracking_id = systemfunc("uuidgen")


; copy variables

var_names = getfilevarnames (fint) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. "tas") then   
  fouth->$var_names(i)$ = fint->$var_names(i)$
  foutd->$var_names(i)$ = fint->$var_names(i)$
end if
end do



;; check units

if(fint->tas@units .ne. "K") then
  print("tas has units '"+fint->tas@units+"'; needs to be 'K'")
  exit
end if
if(finh->huss@units .ne. "kg kg-1") then
  print("huss has units '"+finh->huss@units+"'; needs to be 'kg kg-1'")
  exit
end if
if(finp->ps@units .ne. "Pa") then
  print("ps has units '"+finp->ps@units+"'; needs to be 'Pa'")
  exit
end if


;; Calculate relative humidity

;; Because there's not enough room to fit everything in memory at
;; once, we loop on time


;; initialize this way to get correct structure

hurs = fint->tas*0
dewpt = hurs


nt = dimsizes(fint->time)
do i=0, nt-1

  t = fint->tas(i,:,:)
  p = finp->ps(i,:,:)
  h = finh->huss(i,:,:)
;; x = mixing ratio
  x = h/(1-h)

  temphurs = relhum(t,x,p)
;; NB: function output is % RH; dewtemp_trh wants % RH, but CF unit is 1
  tempdewpt = dewtemp_trh(t,temphurs)
  temphurs = temphurs/100

  hurs(i,:,:) = (/temphurs/)
  dewpt(i,:,:) = (/tempdewpt/)

   if (verbose .and. (i%100) .eq. 0) then print(""+i) end if 
end do


xdim = fint->tas!2
hurs!2 = xdim
hurs&$xdim$=fint->$xdim$

ydim = fint->tas!1
hurs!1 = ydim
hurs&$ydim$=fint->$ydim$

hurs!0 = "time"
;hurs!1 = "yc"
;hurs!2 = "xc"
;hurs!1 = "lat"
;hurs!2 = "lon"

hurs&time = fint->time
;hurs&xc = fint->xc
;hurs&yc = fint->yc
;hurs&lon = fint->lon
;hurs&lat = fint->lat

hurs@units = "1"
hurs@long_name = "Relative Humidity"
hurs@standard_name = "relative_humidity"
hurs@missing_value = 1.e+20
hurs@_FillValue = 1.e+20
hurs@coordinates = t@coordinates
hurs@grid_mapping = t@grid_mapping

fouth->hurs = hurs

xdim = fint->tas!2
dewpt!2 = xdim
dewpt&$xdim$=fint->$xdim$

ydim = fint->tas!1
dewpt!1 = ydim
dewpt&$ydim$=fint->$ydim$


dewpt!0 = "time"
;dewpt!1 = "yc"
;dewpt!2 = "xc"
;dewpt!1 = "lat"
;dewpt!2 = "lon"

dewpt&time = fint->time
;dewpt&xc = fint->xc
;dewpt&yc = fint->yc
;dewpt&xc = fint->lon
;dewpt&yc = fint->lat

dewpt@units = "K"
dewpt@long_name = "Dew Point Temperature"
dewpt@standard_name = "dew_point_temperature"
dewpt@missing_value = 1.e+20
dewpt@_FillValue = 1.e+20
dewpt@coordinates = t@coordinates
dewpt@grid_mapping = t@grid_mapping

foutd->dewpt = dewpt

exit

;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
