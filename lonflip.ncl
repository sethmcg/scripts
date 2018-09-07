;;; NCL lonflip: reorder longitude from -180:180 -> 0:360

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Specify arguments to script via command-line:
;; ncl lonflip.ncl infile=\"$in\" outfile=\"$out\" varname=\"$var\"

system("rm "+outfile)
fin = addfile(infile, "r")
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited

data = fin->$varname$

odata = lonFlip(data)

printVarSummary(data)
printVarSummary(odata)

fout->$varname$ = odata



; copy global attributes

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fin@$att_names(i)$
end do

;; No need to add history entry - varname@lonFlip attribute records operation

; copy variables

var_names = getfilevarnames (fin) ;

do i = 0,dimsizes(var_names)-1   
  if (var_names(i) .eq. "lat"    .or. \
      var_names(i) .eq. "lon"    .or. \
      var_names(i) .eq. varname       \
      ) then
        print(""+var_names(i)+": skipped")
  else
  print(""+var_names(i))
    fout->$var_names(i)$ = fin->$var_names(i)$
end if
end do


delete(fin)
delete(fout)
delete(var_names)
delete(att_names)



;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
