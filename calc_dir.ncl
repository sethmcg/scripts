;;; NCL script to calculate wind direction from u & v components

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; The following variables should be specified as command-line arguments
;; e.g.: ncl file=\"$file\" script.ncl

;ufile = "uas_CRCM_2001010103.nc"
;vfile = "vas_CRCM_2001010103.nc"

;outfile = "dir_CRCM_2001010103.nc"

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


; copy variables

var_names = getfilevarnames (finu) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. "uas") then   
  fout->$var_names(i)$ = finu->$var_names(i)$
end if
end do


; calculate wind direction

rad2deg = 45.0/atan(1.0)
dir = finu->uas 
printVarSummary(dir)
dir = dir * 0.
printVarSummary(dir)
nt = dimsizes(dir&time)
printVarSummary(dir)
do t = 0, nt-1
   tmp  = atan2(finu->uas(t,:,:), finv->vas(t,:,:)) * rad2deg + 180
   dir(t,:,:) = dir(t,:,:) + tmp
end do


varatts = (/"missing_value", "_FillValue", "coordinates", "grid_mapping"/)
do i = 0,dimsizes(varatts)-1
   dir@$varatts(i)$ = finu->uas@$varatts(i)$
end do

dir@units = "degrees"
dir@units_convention = "degrees clockwise from north"
dir@long_name = "10-meter Wind Direction"
dir@standard_name = "wind_from_direction"

fout->dir = dir

exit

;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
