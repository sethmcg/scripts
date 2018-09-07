;;; NCL script to calculate heat index from temp & relative humidity

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "heat_stress.ncl"

;; Command-line arguments:
;; tasfile      required: name of input file with temperature
;; hursfile     required: name of input file with relative humidity
;; hidxfile     required: name of output file for heat index
;; outunits     string option: units of output (degC, K, or degF)

tstrings = (/"degC", "K", "degF"/)

if (.not. isvar("outunits")) then
  outunits = "degC"
end if
ou = ind(outunits .eq. tstrings)
if (ismissing(ou)) then
  print("Invalid output units: "+outunits)
  exit
end if



fint = addfile(tasfile, "r")
finh = addfile(hursfile, "r")

system("rm -f "+hidxfile)
fout = addfile(hidxfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


; global attributes

att_names = getvaratts(fint) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fint@$att_names(i)$
end do

;; To track the full history of the file, we need histories of all
;; the ancestors.

BR = tochar(10)

history = "###################################"+BR
history = history + "History of parent file " + hursfile + ":"+BR
history = history + finh@history+BR
history = history + "-------------------------"+BR
history = history + "History of parent file " + tasfile + ":"+BR
history = history + fint@history+BR
history = history + "#########################"+BR
history = history + "Created " + systemfunc("date")
history = history + " by "+systemfunc("whoami")+"@"+systemfunc("hostname")
history = history + " using NCL script calc_hidx.ncl from parent files listed above."+BR

fout@history = history

; copy variables

var_names = getfilevarnames (fint) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. "tas") then   
  fout->$var_names(i)$ = fint->$var_names(i)$
end if
end do


;; check units

tunits  = fint->tas@units
rhunits = finh->hurs@units

if(.not. any(tunits .eq. tstrings)) then
  print("invalid units for tas: "+tunits)
  exit
end if
if(.not. (rhunits .eq. "1") .or. (rhunits .eq. "%")) then
  print("invalid units for rh: "+rhunits)
  exit
end if


; calculate heat index

;; Note: cannot loop on time, because heat_index_nws_eqns uses
;; "where", which errors on a non-array value

;; initialize this way to get correct structure
hidx = fint->tas*0

iounit = (/ind(tunits .eq. tstrings), ou/)

t = fint->tas
rh = finh->hurs
if(rhunits .eq. "1") then
  rh = rh * 100  ;; convert to %
end if

hidx(:,:,:) = (/heat_index_nws(t, rh, iounit, False)/)


xdim = fint->tas!2
hidx!2 = xdim
hidx&$xdim$=fint->$xdim$

ydim = fint->tas!1
hidx!1 = ydim
hidx&$ydim$=fint->$ydim$

hidx!0 = "time"

hidx&time = fint->time

hidx@units = outunits
hidx@long_name = "NWS Heat Index"
;hidx@standard_name = ""  ;; currently undefined
hidx@missing_value = 1.e+20
hidx@_FillValue = 1.e+20
hidx@coordinates = t@coordinates
hidx@grid_mapping = t@grid_mapping

fout->hidx = hidx

exit

;; Copyright 2016 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
