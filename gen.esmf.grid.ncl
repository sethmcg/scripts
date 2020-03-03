load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;; create SCRIP gridfile for use in patch regridding

;; Specify arguments to script via command-line:
;; ncl gen.esmf.grid.ncl infile=\"$in\" outfile=\"$out\" maskvar=\"sftlf\"

fin = addfile(infile,"r")
lat = fin->lat
lon = fin->lon

lon = where(lon.lt.0, lon+360, lon)

opt = True
opt@ForceOverwrite = True

if(isvar("maskvar")) then
  ; https://www.ncl.ucar.edu/Document/Functions/ESMF/curvilinear_to_SCRIP.shtml
  landmask = fin->$maskvar$
  opt@GridMask = where(landmask.le.0,0,1)    ; use "Mask2D" in NCL V6.2.1 and earlier
end if


if(dimsizes(dimsizes(lat)) .lt. 2) then  
  rectilinear_to_SCRIP(outfile,lat,lon,opt)
else
  curvilinear_to_SCRIP(outfile,lat,lon,opt)
end if


;; Copyright 2016 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
