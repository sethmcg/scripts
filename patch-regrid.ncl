load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;infile  = "swe_MM5I_ncep_10days.nc"
;wgtfile = "weights/wgt.MM5I.UDEL.grid.nc"
;outfile = "test.swe.nc"
;varname = "swe"

;; Specify arguments to script via command-line:
;; ncl patch-regrid.ncl infile=\"$in\" wgtfile=\"$wgt\" outfile=\"$out\" varname=\"$var\"

system("rm "+outfile)
fin = addfile(infile, "r")
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited

data = fin->$varname$

opt               = True
opt@CopyVarAtts   = True
opt@CopyVarCoords = True
opt@FixMapMethodForPatch = True  ;; else metadata sez bilinear instead of patch

odata = ESMF_regrid_with_weights (data, wgtfile, opt)

delete(odata@grid_mapping)

fout->$varname$ = odata

fout->lat@standard_name = "latitude"
fout->lon@standard_name = "longitude"


; copy global attributes

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fin@$att_names(i)$
end do



;; append history entry

hstring = systemfunc("date")
hstring = hstring + ": ncl patch-regrid.ncl infile="+infile+" wgtfile="+wgtfile+" outfile="+outfile+" varname="+varname
fout@history = hstring+inttochar(10)+fout@history

;; Copyright 2013 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu 
