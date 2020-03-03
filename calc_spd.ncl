;;; NCL script to calculate wind speed from u & v components

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; The following variables should be specified as command-line arguments
;; e.g.: ncl file=\"$file\" script.ncl

;ufile = "uas_CRCM_2001010103.nc"
;vfile = "vas_CRCM_2001010103.nc"

;outfile = "spd_CRCM_2001010103.nc"


finu = addfile(ufile, "r")
finv = addfile(vfile, "r")

system("rm -f "+outfile)
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


; copy/set global attributes

att_names = getvaratts(finu) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = finu@$att_names(i)$
end do

history = "Created " + systemfunc("date") + " by "+systemfunc("whoami")+"@"+systemfunc("hostname")+" using NCL script from source files "+ufile+" and "+vfile

fout@history = history

;; update unique tracking_id
fout@tracking_id = systemfunc("uuidgen")


; copy variables

var_names = getfilevarnames (finu) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. "uas") then   
  fout->$var_names(i)$ = finu->$var_names(i)$
end if
end do

; calculate wind speed
; note: roundabout approach necessary to avoid memory ceiling...

;;spd = finu->uas
;;tmp = finv->vas
;;
;;spd = spd^2
;;tmp = tmp^2
;;spd = spd + tmp
;;spd = sqrt(spd)


spd = finu->uas
spd = spd^2

nt = dimsizes(spd&time)
do t = 0, nt-1
tmp = finv->vas(t,:,:)
tmp = tmp^2
spd(t,:,:) = spd(t,:,:) + tmp
end do

spd = sqrt(spd)



;; get rid of superfluous attributes

delete_VarAtts(spd, -1)

spd@long_name = "Near-Surface Wind Speed"
spd@standard_name = "wind_speed"

varatts = (/"units", "missing_value", "_FillValue", "coordinates", "grid_mapping"/)

do i = 0,dimsizes(varatts)-1      
  spd@$varatts(i)$ = finu->uas@$varatts(i)$
end do


fout->spd = spd

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
