;;; Removes the largest or smallest value in the file.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Pass in filename and varname via command-line.  Defaults to
;; removing largest value; to remove smallest value, pass in outfile
;; with value 'neg'.  (This flag is named 'outfile' so this program
;; can be used with nclwrap script without twiddling.)

;;infile = "fixed.nc"
;;varname = "pr"
;;outfile = "neg"

f = addfile(infile, "w")


d = dimsizes(f->$varname$)

if (outfile .eq. "neg") then
mi = minind(ndtooned(f->$varname$))
else
mi = maxind(ndtooned(f->$varname$))
end if

loc = ind_resolve(mi, d)

;;print(loc)
;;print(f->$varname$(loc(0,0), loc(0,1), loc(0,2)))

f->$varname$(loc(0,0), loc(0,1), loc(0,2)) = f->$varname$@missing_value

quit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
