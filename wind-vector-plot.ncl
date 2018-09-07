;;; NCL script for creating wind-vector plots of NARCCAP RCM data

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Uncomment the following variables to hard-code them here, or pass values in via command-line.  Example:

;; ncl -x ufile=\"uas_RCM3_cgcm3_2066010103.nc\" vfile=\"vas_RCM3_cgcm3_2066010103.nc\" outfile=\"wind_vec\" timestep=0 wind-vector-plot.ncl

;ufile    = "uas_RCM3_cgcm3_2066010103.nc" 
;vfile    = "vas_RCM3_cgcm3_2066010103.nc" 
;outfile  = "wind_vec" 
;timestep = 0   ;; defaults to last step in file

;; open files
uin = addfile(ufile, "r")
vin = addfile(vfile, "r")

;; read in coordinate variables first

lat = uin->lat
lon = uin->lon
time = uin->time

nx = dimsizes(uin->xc)
ny = dimsizes(uin->yc)

;; default to last timestep in file for QC convenience
if (.not. isvar("timestep")) then
  timestep = dimsizes(time) - 1
end if
;; default title for convenience
if (.not. isvar("title")) then
  date = ut_calendar(time(timestep),0)
  timestamp = ""+date(0,0)+"/"+date(0,1)+"/"+date(0,2)+" "+sprintf("%02.0f", date(0,3))+":"+sprintf("%02.0f", date(0,4))
  title = systemfunc("basename "+ufile)+"~C~"+systemfunc("basename "+vfile)+"~C~"+timestamp
end if


;; read in data.

udata = uin->uas(timestep,:,:)
vdata = vin->vas(timestep,:,:)


;; hook up coordinate variables

udata@lat2d = lat
udata@lon2d = lon
vdata@lat2d = lat
vdata@lon2d = lon


;; read in map projection info
pname = udata@grid_mapping
proj = uin->$pname$


;; set plotting parameters

res = True


;; map projection stuff.

if (lower_case(proj@grid_mapping_name) .eq. "lambert_conformal_conic") then
  res@mpProjection	= "LambertConformal"
  res@mpLambertMeridianF	= proj@longitude_of_central_meridian
  res@mpLambertParallel1F	= proj@standard_parallel(0)
  res@mpLambertParallel2F	= proj@standard_parallel(1)
end if

if (lower_case(proj@grid_mapping_name) .eq. "transverse_mercator") then
  res@mpProjection	= "Mercator"
  res@mpCenterLatF	= proj@latitude_of_projection_origin
  res@mpCenterLonF	= proj@longitude_of_central_meridian
end if

if (lower_case(pname) .eq. "polar_stereographic") then
  res@mpProjection= "Stereographic"
  res@mpCenterLonF= proj@straight_vertical_longitude_from_pole - 360
  res@mpCenterLatF  = proj@latitude_of_projection_origin
end if

if (lower_case(proj@grid_mapping_name) .eq. "rotated_latitude_longitude") then
  res@mpProjection = "Mercator"
  res@mpCenterLonF = proj@grid_north_pole_longitude - 180
  res@mpCenterLatF = 90 - proj@grid_north_pole_latitude
end if


;; /map projection stuff
  

;; map boundaries
;; expand slighltly to be sure we plot things on the borders
res@mpLimitMode             = "Corners"
res@mpLeftCornerLatF        = lat(0,0) -1 ;; SW corner
res@mpLeftCornerLonF        = lon(0,0) -1 ;; SW corner
res@mpRightCornerLatF       = lat(yc|ny-1,xc|nx-1) + 1 ;; NE corner
res@mpRightCornerLonF       = lon(yc|ny-1,xc|nx-1) + 1 ;; NE corner

;; general plotting stuff

res@mpFillOn                = False

res@gsnMaximize             = True

res@vcGlyphStyle            = "CurlyVector"
res@vcRefMagnitudeF         = 10.0
res@vcRefLengthF            = 0.02
res@vcMinDistanceF          = 0.02

res@lbLabelAutoStride       = True
res@lbBoxLinesOn            = False

res@gsnSpreadColors         = True
  
res@tiMainString            = title

;; create the plot

wks = gsn_open_wks("ps", outfile)
gsn_define_colormap(wks,"nrl_sirkes") ;; 21 colors + fg/bg
  
plot = gsn_csm_vector_map(wks, udata, vdata, res)


exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
