;;; Removed a degenerate height dimension.
;;; Needs revision - not generic enough

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Note: specify infile and outfile as command-line arguments:
;; ncl -x infile=\"$infile\" outfile=\"$outfile\" remove-degenerate-level.ncl

system("rm "+outfile)
fin = addfile(infile, "r")
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)


; copy global attributes

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fin@$att_names(i)$
end do


; copy variables

var_names = getfilevarnames (fin) ;

do i = 0,dimsizes(var_names)-1

  print(""+var_names(i))

  badnames = (/ "cli", "clw", "hus", "ta", "ua", "va", "wa", "zg"/)


if (any(badnames .eq. var_names(i))  )then 
    fout->$var_names(i)$ = fin->$var_names(i)$(:,0,:,:)
    fout->$var_names(i)$@coordinates = "lon lat level"
  else 

    if  (var_names(i) .eq. "level") then
      level = fin->level(0)
      level!0 = "ncl_scalar"
      fout->level = level
    else 
      fout->$var_names(i)$ = fin->$var_names(i)$
    end if

  end if
end do


delete(fin)
delete(fout)
delete(var_names)
delete(att_names)

;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
