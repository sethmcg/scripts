;;; NCL script for creating contour plots of NARCCAP RCM data 

;;; Usage: Invoke using nclwrap or via command line like so: 
;;; ncl -x varname=\"var\" infile=\"filename.nc\" [etc.] plot.ncl
;;;
;;; Other options that can be specified: 
;;; timestep (defaults to last step in file); 
;;; title (defaults to filename, plus timestamp for time-varying);
;;; colormap (defaults to nrl_sirkes); 
;;; cmapfile (cmap from RGB triplets in file, overrides colormap);
;;; creverse (reverse order of cmap)
;;; wkstype (defaults to ps);
;;; fillmode (defaults to raster); 
;;; global (defaults to False, unless data is lat-lon);
;;; drawstates (defaults to False); 
;;; contour (vector of explicit contour levels);
;;; cminmax (vector of contour minval, maxval, spacing); 
;;; zinvert (defaults to False); 
;;; bounds (for global: vector of minlat,maxlat,minlon,maxlon);
;;; projfile (source for map projection info; defaults to infile);
;;;
;;; Auto-detects time-varying vs 2-D data and lat-lon vs projected
;;; data and Does The Right Thing for each case.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; open input file

fin = addfile(infile, "r")

if (.not. isvar("projfile")) then
  pin = fin
else
  pin = addfile(projfile, "r")
end if



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; detect input type and set defaults


;; time-varying or static?

istime = any("time" .eq. getfilevarnames(fin))

if (istime) then
  time = fin->time

  ;; 'julian' is legit (no missing leap year every 100th year), 
  ;; but NCL doesn't know it
  if (isatt(time, "calendar") .and. time@calendar .eq. "julian") then
    time@calendar = "gregorian"
  end if
 
  if (.not. isvar("timestep")) then
    timestep = dimsizes(time) - 1
  end if

  if (.not. isvar("title")) then
    if (isatt(time,"calendar")) then
      date = ut_calendar(time(timestep),0)
      timestamp = ""+date(0,0)+"/"+date(0,1)+"/"+date(0,2)+" "+sprintf("%02.0f", date(0,3))+":"+sprintf("%02.0f", date(0,4))
      title = systemfunc("basename "+infile)+", "+timestamp
    else
      title = systemfunc("basename "+infile)+", "+sprintf("%0.3f",time(timestep))
    end if
  end if

end if


;; overrideable plotting defaults

if (.not. isvar("title")) then
  title = systemfunc("basename "+infile)
end if

;; default colormap is handled after wks is created

if (.not. isvar("outfile")) then
  wkstype = "x11"
  outfile = systemfunc("basename "+infile)
end if
if (.not. isvar("wkstype")) then
  wkstype = "ps"
end if

if (.not. isvar("fillmode")) then
  fillmode = "RasterFill"
end if

if (.not. isvar("drawstates")) then
  drawstates = False
end if




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; read in data


if (istime) then
  data = fin->$varname$(timestep,:,:)
else 
  data = fin->$varname$
end if

if(isvar("zinvert").and.zinvert) then
  data = -data
end if


lat = pin->lat
lon = pin->lon


if (isatt(data, "grid_mapping")) then
  data@lat2d = fin->lat
  data@lon2d = fin->lon
end if


if (.not. isvar("global")) then

  projvarnames = getfilevarnames(pin)
  do pv = 0,dimsizes(projvarnames)-1
    if (isfilevaratt(pin, projvarnames(pv),"grid_mapping"))
      pname = pin->$projvarnames(pv)$@grid_mapping
      proj = pin->$pname$
      nx = dimsizes(pin->xc)
      ny = dimsizes(pin->yc)
      break
    end if
  end do
  
  if (isvar("proj")) then
    global = False
  else
    global = True
  end if  

end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; map projection & bounds

res = True

if(global) then

  res@mpLimitMode             = "LatLon"

  if (.not. isvar("bounds")) then
    res@mpMinLatF               = min(lat)
    res@mpMaxLatF               = max(lat)
    res@mpMinLonF               = min(lon)
    res@mpMaxLonF               = max(lon)
  else
    res@mpMinLatF               = bounds(0)
    res@mpMaxLatF               = bounds(1)
    res@mpMinLonF               = bounds(2)
    res@mpMaxLonF               = bounds(3)
  end if
      
  res@mpCenterLonF            = (res@mpMinLonF + res@mpMaxLonF)/2

else
  
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


  ;; map boundaries
  ;; expand slighltly to be sure we plot things on the borders
  res@mpLimitMode             = "Corners"
  res@mpLeftCornerLatF        = lat(0,0) ;; SW corner
  res@mpLeftCornerLonF        = lon(0,0) ;; SW corner
  res@mpRightCornerLatF       = lat(yc|ny-1,xc|nx-1) ;; NE corner
  res@mpRightCornerLonF       = lon(yc|ny-1,xc|nx-1) ;; NE corner

end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; general plot settings

res@mpFillOn                = False

res@cnFillMode              = fillmode
res@cnFillOn                = True
res@cnLinesOn               = False

res@gsnAddCyclic            = False

res@gsnMaximize             = True
res@gsnSpreadColors         = True

res@lbBoxLinesOn            = False
res@lbLabelAutoStride       = True

if(drawstates) then
res@mpOutlineBoundarySets   = "geophysicalandusstates"
res@mpDataBaseVersion         = "mediumres"             ; select database
res@mpDataSetName             = "Earth..2"
end if

res@tiMainString            = title
if (isvar("leftstring")) then
res@gsnLeftString           = leftstring
end if

if (isvar("cminmax")) then
res@cnLevelSelectionMode    = "ManualLevels"
res@cnMinLevelValF          = cminmax(0)
res@cnMaxLevelValF          = cminmax(1)
res@cnLevelSpacingF         = cminmax(2)
end if

if (isvar("contour")) then
res@cnLevelSelectionMode    = "ExplicitLevels"
res@cnLevels                = contour
end if



if (wkstype .eq. "x11") then
;;  res@gsnFrame = False
end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; actually create the plot

wks = gsn_open_wks(wkstype, outfile)

if (isvar("cmapfile")) then
  colormap = RGBtoCmap(cmapfile)
end if

if (.not. isvar("colormap")) then
  colormap = "nrl_sirkes"
  ;; amwg_blueyellowred ?
end if

gsn_define_colormap(wks,colormap)

if (isvar("creverse") .and. creverse) then
  gsn_reverse_colormap(wks)
end if

plot = gsn_csm_contour_map(wks, data, res)

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
