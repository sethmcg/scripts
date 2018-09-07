;;; NCL script for plotting 1-D netcdf file

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Uncomment the following variables to hard-code them here, or pass
;; values in via command-line, e.g.:

;; ncl varname=\"tas\" infile=\"tas_CRCM_1979010103.nc\" outfile=\"plots/tas_CRCM_1979010103\" title=\"tas_CRCM_1979010103.nc\ \" time-section-plot.ncl


; varname =    "tas"
; infile =     "tas_CRCM_1986010103.nc"
; outfile =    "plot-test"
; title = systemfunc("basename "+infile)+" "

if (.not. isvar("title")) then
  title = systemfunc("basename "+infile)+" "
end if


fin = addfile(infile, "r")

time = fin->time

dims = filevardimsizes(fin,varname)
ny = dims(1)
nx = dims(2)


x0 = nx/2
y0 = ny/2

res = True

res@cnFillOn                = True
res@cnLinesOn               = False
res@cnFillMode              = "RasterFill"

res@lbLabelAutoStride       = True
res@lbBoxLinesOn            = False
res@lbOrientation           = "Vertical"

res@gsnSpreadColors         = True

res@gsnMaximize             = True
res@gsnPaperOrientation     = "landscape"

res@vpWidthF  = 1.0
res@vpHeightF = 0.3

;; create the plot

wks = gsn_open_wks("ps", outfile)

tres = True

tres@gsnMaximize             = True
tres@gsnPaperOrientation     = "landscape"

tres@vpWidthF  = 1.0
tres@vpHeightF = 0.3

tres@tiMainString     = title
tres@xyMarkLineMode   = "Markers"
tres@xyMarker         = 1
tres@xyMarkerSizeF    = 0.01

plot = gsn_csm_xy(wks, time, fin->$varname$(:,y0,x0), tres)

delete(wks)

exit



;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
