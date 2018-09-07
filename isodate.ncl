;;; NCL script for range of dates in file in ISO-8601 (YYYYMMDD)

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"

;; Pass filename in via command-line.  Example:

;; ncl -x infile=\"table6/prtot_CRCM_climatology.nc\"

format = "%Y%N%D"


;; open file
fin = addfile(infile, "r")

time = fin->time
nt = dimsizes(time)


if(isatt(time,"bounds")) then
  bnds = fin->$time@bounds$
  first = ut_string(bnds(0,0), format)
  last = ut_string(bnds(nt-1,1), format)
else 
  first = ut_string(time(0), format)
  last = ut_string(time(nt-1), format)
end if

print(first+"-"+last)

exit


;; Copyright 2016 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
