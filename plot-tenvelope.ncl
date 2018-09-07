;; Creates a temperature envelope plot: a scatter plot of temperatures
;; at particular times with the daily tmin & tmax overlaid for
;; comparison.

;; Inputs: temperature file, tmin file, and tmax file.

;; specify x and y coordinates for a point as array index values.  If
;; not specified, the point nearest the center of the domain is used.

;; specify start time as a date.  If unspecified, the middle of the
;; coverage period is used.

;; specify length of time plotted as number of days.  Default is 60
;; days.


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;****************************************************


fave = addfile(tavefile, "r")
fmin = addfile(tminfile, "r")
fmax = addfile(tmaxfile, "r")


if(.not. isvar("x")) then
  d = getfilevardimsizes(fave,"tas")
  x = d(2)/2
end if
if(.not. isvar("y")) then
  d = getfilevardimsizes(fave,"tas")
  y = d(1)/2
end if
if(.not. isvar("length")) then
  length = 60
end if

time = fave->time

if(.not. isvar("start")) then
  sdate=cd_calendar(time(dimsizes(time)/2),0)
  start = ""+sdate(0,0)+"-"+sdate(0,1)+"-"+sdate(0,2)
end if


ds = stringtointeger(str_split(start,"-"))
cal=0
cal@calendar=time@calendar
t0 = cd_inv_calendar(ds(0),ds(1),ds(2),0,0,0,time@units,cal)
t1 = t0+length

tave = fave->tas({t0:t1},y,x)
tmin = fmin->tasmin({t0:t1},y,x)
tmax = fmax->tasmax({t0:t1},y,x)

tt = fave->time({t0:t1})

dummy = new((/dimsizes(tmin),2/),double)

dummy(:,0) = (/fmin->time_bnds({t0:t1},0)/)
dummy(:,1) = (/fmin->time_bnds({t0:t1},1)/)

tm = ndtooned(dummy)

dummy(:,0) = tmin
dummy(:,1) = tmin

pmin = ndtooned(dummy)

dummy(:,0) = tmax
dummy(:,1) = tmax

pmax = ndtooned(dummy)


  colors = (/"black","blue","red"/)

  wks = gsn_open_wks("ps",outfile) 

  res                  = True 

  res@gsnMaximize      = True
  res@gsnPaperOrientation = "portrait"
  res@gsnDraw          = False
  res@gsnFrame         = False

  res@trYMinF          = min(tmin)    ; You could also just use
  res@trYMaxF          = max(tmax)    ; min/max of u0,u1,u2.


  res@xyMarkLineMode  = "Markers"
  res@xyMarkers       = 1
  res@xyMarkerColor   = colors(0)
  res@xyMarkerSizeF    = 0.02

  plot0 = gsn_csm_xy(wks,tt,tave,res)


  res@xyMarkLineMode  = "Lines"
  res@xyLineThicknessF = 2.0
  res@xyLineColor      = colors(0)

  res@xyLineColor     = colors(1)
  plot1 = gsn_csm_xy(wks,tm,pmin,res)

  res@xyLineColor   = colors(2)

  plot2 = gsn_csm_xy(wks,tm,pmax,res)

  overlay(plot0,plot1)
  overlay(plot0,plot2)

  draw(plot0)
  frame(wks)

;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
