;;; NCL script for creating contour plots of 2D time-varying NARCCAP RCM data

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Uncomment the following variables to hard-code them here, or pass
;; values in via command-line.  Note that there are 3 output files,
;; named outfile.tt, outfile.tx, and outfile.ty.  Example:

;; ncl varname=\"tas\" infile=\"tas_CRCM_1979010103.nc\" outfile=\"plots/tas_CRCM_1979010103\" title=\"tas_CRCM_1979010103.nc\ \" time-section-plot.ncl


; varname =    "tas"
; infile =     "tas_CRCM_1986010103.nc"
; outfile =    "plot-test"
; title = systemfunc("basename "+infile)+" "

if (.not. isvar("title")) then
  title = systemfunc("basename "+infile)+" "
end if


;; open file
fin = addfile(infile, "r")

;; read in coordinate variables first

lat = fin->lat
lon = fin->lon
time = fin->time

;;nt = dimsizes(fin->time)
;;nx = dimsizes(fin->xc)
;;ny = dimsizes(fin->yc)
dims = filevardimsizes(fin,varname)
ny = dims(1)
nx = dims(2)

norm = 1.0* max((/nx,ny/))  ;; for plot aspect ratios

if (.not. isvar("x0")) then
x0 = nx/2
end if
if (.not. isvar("y0")) then
y0 = ny/2
end if

;; read in data.

txdata = fin->$varname$(:,y0,:)
tydata = fin->$varname$(:,:,x0)



;; change coordinate variables to play nicely with plotting routines

if(dimsizes(dimsizes(lon)) .gt. 1) then
  txdata!1 = "lon"
  txdata&lon = lon(y0,:)
  
  tydata!1 = "lat"
  tydata&lat = lat(:,x0)
  
  lat00 = lat(y0,x0)
  lon00 = lon(y0,x0)
else
  lat00 = lat(y0)
  lon00 = lon(x0)
end if


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

;; create the plots


wks = gsn_open_wks("ps", outfile+".tx")
gsn_define_colormap(wks,"nrl_sirkes") ;; 21 colors + fg/bg


res@tiMainString     = title+"lat "+lat00

plot = gsn_csm_contour(wks, txdata(lon|:,time|:), res)

delete(wks)



wks = gsn_open_wks("ps", outfile+".ty")
gsn_define_colormap(wks,"nrl_sirkes") ;; 21 colors + fg/bg

res@tiMainString     = title+"lon "+lon00

plot = gsn_csm_contour(wks, tydata(lat|:,time|:), res)

delete(wks)



wks = gsn_open_wks("ps", outfile+".tt")
gsn_define_colormap(wks,"nrl_sirkes") ;; 21 colors + fg/bg

tres = True

tres@gsnMaximize             = True
tres@gsnPaperOrientation     = "landscape"

tres@vpWidthF  = 1.0
tres@vpHeightF = 0.3

tres@tiMainString     = title+"lon "+lon00+", lat "+lat00
tres@xyMarkLineMode   = "Markers"
tres@xyMarker         = 1
;tres@xyMarkerSizeF    = 0.05
tres@xyMarkerSizeF    = 0.01

plot = gsn_csm_xy(wks, time, fin->$varname$(:,y0,x0), tres)

delete(wks)

exit


