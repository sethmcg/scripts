;;; NCL script to vertically interpolate pressure-level data to height

;; Command-line arguments:
;; infile   required: name of input file with data to be interpolated
;; zfile    required: name of input file with zg (geopotential height)
;; outfile  required: name of output file
;; varname  required: data variable name
;; zname    optional: height variable name (defaults to 'zg')
;; height   requires: target height (same units as zfile)

;infile  = "data/raw/ua.hist.MPI-ESM-LR.WRF.mon.NAM-22.raw.nc" 
;zfile   = "data/raw/zg.hist.MPI-ESM-LR.WRF.mon.NAM-22.raw.nc"
;outfile = "data/ua.100m.hist.MPI-ESM-LR.WRF.mon.NAM-22.raw.nc"
;varname = "ua"
;height=100.

if (.not. isvar("zname")) then
  zname = "zg"
end if


fin = addfile(infile, "r")
zin = addfile(zfile, "r")


vin = fin->$varname$
z = zin->zg

;;result = wrf_interp_3d_z(vin, z, height)

opt = True
opt@inc2dlevs=True

result = wrf_user_interp_level(vin, z, height, opt)


;; add coordinate variables

dimvars = getvardims(result)
do i = 0,dimsizes(dimvars)-1
  result&$dimvars(i)$ = vin&$dimvars(i)$
end do

;; copy over variable attributes

varatts = getfilevaratts(fin, varname)
do i = 0,dimsizes(varatts)-1
  result@$varatts(i)$ = fin->$varname$@$varatts(i)$
end do


;;;;;;;;;;;;;;;;;;  output

system("rm -f "+outfile)
fout = addfile(outfile, "c")
filedimdef(fout, "time", -1, True)  ;; make time dimension unlimited



; copy/set global attributes
; (WRF att spam deleted elsewhere)

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  if(att_names(i) .eq. "history_of_appended_files") then
    ;; skip
  else 
    fout@$att_names(i)$ = fin@$att_names(i)$
  end if
end do



;; To track the full history of the file, we need histories of all                    
;; the ancestors.  Following NCO style                                                

BR = str_get_nl()

parents = "Parent files had the following history attributes:"+BR
parents = parents + "input file "+infile+":"+BR
parents = parents + fin@history+BR
parents = parents + "height file "+zfile+":"+BR
parents = parents + zin@history+BR

fout@history_of_parent_files = parents

history = systemfunc("date")+": "+outfile+ " vertically interpolated to "
history = history + height + z@units + " from "
history = history + infile +" and " + zfile + " using NCL function wrf_user_interp_level."

fout@history = history


;; update unique tracking_id
fout@tracking_id = systemfunc("uuidgen")


; copy ancillary & coordinate variables

skipme = venn2_difference(getvardims(result), getvardims(vin))

if(dimsizes(skipme) .ne. 1) then
  print("Error: unhandled case: skipping more than 1 coordinate variable.")
  exit()
end if

var_names = getfilevarnames (fin) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. varname .and. var_names(i) .ne. skipme) then   
  fout->$var_names(i)$ = fin->$var_names(i)$
end if
end do



;;; write result to file

fout->$varname$ = result






exit


;; Copyright 2023 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
