;;; NCL script for printing dates of first and last timestep in file

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"

;; Pass filename in via command-line.  Example:

;; ncl -x infile=\"table6/prtot_CRCM_climatology.nc\"

format = "%Y-%N-%D %H:%M:%S"


;; open file
fin = addfile(infile, "r")

time = fin->time
if(time@calendar .eq. "proleptic_gregorian") then time@calendar = "gregorian" end if
nt = dimsizes(time)

infile = systemfunc("basename "+infile)

first = ut_string(time(0), format)
last = ut_string(time(nt-1), format)

print(""+infile+"	"+first+"	"+last)

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
