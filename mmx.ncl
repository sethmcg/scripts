;;; Prints minimum and maximum values of data variable

;; Pass in variable name and input file from command-line: 

;; ncl -x varname=\"prtot\" infile=\"table6/prtot_CRCM_climatology.nc\"

tab  =  "	"

;; open file
fin = addfile(infile, "r")

data = fin->$varname$

print(infile+tab+min(data)+tab+max(data))

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
