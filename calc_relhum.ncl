;;; NCL script to calculate relative humidity or dewpoint from
;;; specific humidity, temperature, and pressure

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Command-line arguments:
;; pfile     required: name of input file with pressure
;; qfile     required: name of input file with specific humidity
;; tfile     required: name of input file with temperature
;; outfile   required: name of output file for relative humidity
;; verbose   boolean option: print progress indicators?
;; dewpt     boolean option: calculate dewpt instead of relhum
;; pname     optional: name of pressure variable (default: "ps")
;; qname     optional: name of specific humidity variable (default: "huss")
;; tname     optional: name of temperature variable (default: "tas")
;; outname   optional: name of output variable (default: "hurs" / "dewpt")


;; control flags

if (.not. isvar("verbose")) then
  verbose = False
end if

if (.not. isvar("dewpt")) then
  dewpt = False
end if


;; print timestamped message to log progress
procedure message(text) begin
  if(verbose) then
    print(systemfunc("date '+%F %T'") + str_get_tab + text)
  end if
end


;; overrideable variable names

if (.not. isvar("pname")) then
  pname = "ps"
end if

if (.not. isvar("tname")) then
  tname = "tas"
end if

if (.not. isvar("qname")) then
  qname = "huss"
end if

if (.not. isvar("outname")) then
  if(dewpt) then
    outname = "dewpt"
  else 
    outname = "hurs"
  end if
end if


message("starting calc_relhum")

;; add input files

finp = addfile(pfile, "r")
finq = addfile(qfile, "r")
fint = addfile(tfile, "r")

;; check units

tunits = fint->$tname$@units
if(tunits .ne. "K" .and. tunits .ne. "degC") then
  print(tname+" has units '"+tunits+"'; needs to be 'K' or 'degC'")
  exit
end if
qunits = finq->$qname$@units
if(qunits .ne. "kg kg-1" .and. qunits .ne. "1" ) then
  print(qname+" has units '"+qunits+"'; needs to be 'kg kg-1' or '1'")
  exit
end if
punits = finp->$pname$@units
if(punits .ne. "Pa") then
  print(pname+" has units '"+punits+"'; needs to be 'Pa'")
  exit
end if



;; file setup, since we'll write output as we go
message("starting output file setup")
system("rm -f " + outfile)
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


;; use tfile as template, since it has a reference height

;; copy global attributes

att_names = getvaratts(fint) 

do i = 0,dimsizes(att_names)-1
  if(att_names(i) .eq. "history_of_appended_files") then
    ;; skip
  else if(att_names(i) .eq. "title") then
    if (dewpt) then
      fout@title = str_sub_str(fint@title, "2m Temperature", "Dew Point Temperature")
    else
      fout@title = str_sub_str(fint@title, "2m Temperature", "Near-Surface Relative Humidity")
    end if
  else
    fout@$att_names(i)$ = fint@$att_names(i)$
  end if
  end if
end do



;; To track the full history of the file, we need histories of all
;; the ancestors.  Following NCO style

BR = str_get_nl()

parents = "Parent files had the following history attributes:"+BR
parents = parents + "Specific humidity file "+qfile+":"+BR
parents = parents + finq@history+BR
parents = parents + "Temperature file "+tfile+":"+BR
parents = parents + fint@history+BR
parents = parents + "Pressure file "+pfile+":"+BR
parents = parents + finp@history+BR

fout@history_of_parent_files = parents

history = systemfunc("date")+": "+outname+ " calculated from "
history = history + qname +", " + tname + ", and " + pname
history = history + " using NCL function relhum."

fout@history = history

;; update unique tracking_ids
fout@tracking_id = systemfunc("uuidgen")


; copy coordinate and ancillary variables
var_names = getfilevarnames (fint)

do i = 0, dimsizes(var_names)-1      
  if (var_names(i) .ne. tname) then
    fout->$var_names(i)$ = fint->$var_names(i)$
  end if
end do

message("output file setup done")


;; read inputs

message("read q")
q = finq->$qname$
x = q/(1-q)  ;; x = mixing ratio

message("read t")
t = fint->$tname$
if(tunits .eq. "degC") then
  t = t + 273.15
end if

message("read p")
;; allow for static pressure instead of time-varying
if(isfilevardim(finp, pname, "time")) then
    p = finp->$pname$
else 
  p = conform(q, finp->$pname$, (/1,2/))
end if



;; Calculate relative humidity

message("starting relhum calculation")
temphurs = relhum(t,x,p)
message("relhum calculation done")

if(dewpt) then
  message("starting dewtemp calculation")
  result = dewtemp_trh(t,temphurs)
  message("dewtemp calculation done")
else 
  message("clip %RH to (0,100)")
  result = tofloat(temphurs < 100)  ;; clip values above 100%
end if


message("setup ancillary & metadata on result")


;; set missing value

if(isfilevaratt(fint, tname, "missing_value")) then
  missval = fint->$tname$@missing_value
else if(isfilevaratt(fint, tname, "missing_value")) then
  missval = fint->$tname$@_FillValue
else 
  missval = 1e20
end if
end if

result@_FillValue = missval
result@missing_value = missval

xdim = fint->$tname$!2
result!2 = xdim
result&$xdim$=fint->$xdim$

ydim = fint->$tname$!1
result!1 = ydim
result&$ydim$=fint->$ydim$

result!0 = "time"
result&time = fint->time


if(dewpt) then
  result@units = "K"
  result@long_name = "Dew Point Temperature"
  result@standard_name = "dew_point_temperature"
else
  result@units = "%"
  result@long_name = "Near-Surface Relative Humidity"
  result@standard_name = "relative_humidity"
end if


;; copy misc attributes from temperature file (which should have
;; appropriate 2-m height stuff)

varatts = (/"cell_methods", "bias_correction", "remap", "coordinates", "grid_mapping"/)

valid = isatt(fint->$tname$, varatts)

do i = 0,dimsizes(varatts)-1
  if(valid(i)) then
    result@$varatts(i)$ = fint->$tname$@$varatts(i)$
  end if
end do

message("write result")

fout->$outname$ = result

message("calc_relhum done")

exit

;; Copyright 2020 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
