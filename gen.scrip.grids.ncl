load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;; create SCRIP gridfiles for each model + UDEL


rcms = (/"CRCM","ECP2","HRM3","MM5I","RCM3","WRFG","UDEL"/)
nm = dimsizes(rcms)

do m=0,nm-1
  model = rcms(m)

  print(""+model)

  infile = "landmask."+model+".nc"
  
  fin = addfile(infile,"r")
  lat = fin->lat
  lon = fin->lon
  land = fin->landmask
  delete(fin)  

  outfile = "scrip.grid."+model+".nc"

  opt = True
  opt@ForceOverwrite = True

  if (model .eq. "UDEL") then
    rectilinear_to_SCRIP(outfile,lat,lon,opt)
  else
    curvilinear_to_SCRIP(outfile,lat,lon,opt)
  end if

  outfile = "scrip.mask."+model+".nc"
  opt@Mask2D = land

  if (model .eq. "UDEL") then
    rectilinear_to_SCRIP(outfile,lat,lon,opt)
  else
    curvilinear_to_SCRIP(outfile,lat,lon,opt)
  end if

  delete([/lat,lon,land,infile,outfile,opt/])

end do


;; Copyright 2013 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
