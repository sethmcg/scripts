;;; NCL script to calculate wind speed from u & v components

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


;;; Usage: Invoke using nclwrap or via command line like so: 
;;; ncl -x arg=\"value\" [etc.] script.ncl
;;;
;;; Mandatory command-line arguments:
;;; ufile (netcdf file with data variable 'uas'
;;; vfile (netcdf file with data variable 'vas'
;;; outfile (name of netcdf file to write output to)
;;;    other variables & global attributes are copied from ufile
;;;
;;; optional command-line arguments:
;;; vname (name of output variable, defaults to 'spd')


;; Default name for output variable
if(.not. isvar("vname")) then
  vname = "spd"
end if


finu = addfile(ufile, "r")
finv = addfile(vfile, "r")

system("rm -f "+outfile)
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


; copy/set global attributes

att_names = getvaratts(finu) 

do i = 0,dimsizes(att_names)-1
  if(att_names(i) .eq. "history_of_appended_files") then
    ;; skip
  else if(att_names(i) .eq. "title") then
    title = str_sub_str(finu@title, "Velocity","Speed")
    title = str_sub_str(title, "Eastward","")
    title = str_sub_str(title, "Northward","")
    title = str_sub_str(title, "  "," ")
    fout@title = title
  else 
    fout@$att_names(i)$ = finu@$att_names(i)$
  end if
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
; roundabout looping approach taken to reduce memory usage

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


varatts = (/"units", "missing_value", "_FillValue", "cell_methods", "bias_correction", "remap", "coordinates", "grid_mapping"/)

valid = isatt(finu->uas, varatts)

do i = 0,dimsizes(varatts)-1
  if(valid(i)) then
    spd@$varatts(i)$ = finu->uas@$varatts(i)$
  end if
end do

;fout->spd = spd
fout->$vname$ = spd

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
