;;; Prints minimum and maximum values of data variable, plus
;;; the coordinates where they occur.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Pass in variable name and input file from command-line: 

;; ncl -x varname=\"prtot\" infile=\"table6/prtot_CRCM_climatology.nc\"

tab  =  "	"

;; open file
fin = addfile(infile, "r")

d = dimsizes(fin->$varname$)
data = ndtooned(fin->$varname$)

mini = minind(data)
maxi = maxind(data)

if (ismissing(mini)) then
minloc = d*0-1
mindata = "MISSING"
else
mnl = ind_resolve(mini, d)
minloc = mnl(0,:)
mindata = data(mini)
end if

if (ismissing(mini)) then
maxloc = d*0-1
maxdata = "MISSING"
else
mxl = ind_resolve(maxi, d)
maxloc = mxl(0,:)
maxdata = data(maxi)
end if

;;result = infile + tab + data(mini) + tab + data(maxi) + tab
result = infile + tab + mindata + tab + maxdata + tab

do i = 1,dimsizes(d)
 result = result + minloc(i-1)+","
end do

result = result + tab

do i = 1,dimsizes(d)
 result = result + maxloc(i-1)+","
end do


print(""+result)

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
