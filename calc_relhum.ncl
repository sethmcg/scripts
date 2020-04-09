;;; NCL script to calculate relative humidity or dewpoint from
;;; specific humidity, temperature, and pressure

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

qfile="/glade/collections/cdg/data/cordex/bcdata2/mbcn-Daymet-ns/NAM-44i/day/CRCM5-UQAM/MPI-ESM-MR/hist/huss/huss.hist.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.mbcn-Daymet-ns.nc"
tfile="/glade/collections/cdg/data/cordex/bcdata2/mbcn-Daymet-ns/NAM-44i/day/CRCM5-UQAM/MPI-ESM-MR/hist/tmean/tmean.hist.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.mbcn-Daymet-ns.nc"
pfile="/glade/p/cisl/risc/narccap/obs/Daymet/ps/sps.Daymet.fixed.NAM-44i.nc"
outfile="~/glade-scratch/test-hurs.nc"
pname="sps"
tname="tmean"
verbose=True



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

;; Changes:

;; use overrideable varnames
;; hursfile -> outfile
;; psfile tasfile hussfile -> pfile tfile qfile
;; flag for dewpt instead of hurs
;; flag for static pressure
;; new nco-style parent histories
;; convert input units if needed


if (.not. isvar("verbose")) then
  verbose = False
end if

if (.not. isvar("dewpt")) then
  dewpt = False
end if

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

system("rm -f " + outfile)
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


;; use tfile as template, since it has a reference height

; global attributes


att_names = getvaratts(fint) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fint@$att_names(i)$
end do



;; To track the full history of the file, we need histories of all
;; the ancestors.  Following NCO style

BR = tochar(10)

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


; copy variables

var_names = getfilevarnames (fint) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. tname) then
  fout->$var_names(i)$ = fint->$var_names(i)$
end if
end do



;; Calculate relative humidity

;; Because there's not enough room to fit everything in memory at
;; once, we loop on time

if(isfilevaratt(fint, tname, "missing_value")) then
  missval = fint->$tname$@missing_value
else if(isfilevaratt(fint, tname, "missing_value")) then
  missval = fint->$tname$@_FillValue
else 
  missval = 1e20
end if
end if


result = new(getfilevardimsizes(fint, tname), float, missval) 
result@missing_value = result@_FillValue


pstatic = .not. isfilevardim(finp, pname, "time")
if(pstatic) then
  p = finp->$pname$
end if

nt = getfilevardimsizes(fint, "time")

do i=0, nt-1

  t = fint->$tname$(i,:,:)
  if(tunits .eq. "degC") then
    t = t + 273.15
  end if

  if(.not.pstatic) then
    p = finp->$pname$(i,:,:)
  end if

  q = finq->$qname$(i,:,:)
  x = q/(1-q)  ;; x = mixing ratio


  ;; NB: function output is % RH; dewtemp_trh wants % RH, but CF unit is 1
  temphurs = relhum(t,x,p)
  if(dewpt) then
    tempdewpt = dewtemp_trh(t,temphurs)
    result(i,:,:) = tofloat((/tempdewpt/))
  else 
    temphurs = temphurs/100
    result(i,:,:) = tofloat((/temphurs/))
  end if

   if (verbose .and. (i%50) .eq. 0) then 
     print(""+i+tochar(9)+systemfunc("date")) 
   end if 
end do


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
  result@units = "1"
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


fout->$outname$ = result


exit

;; Copyright 2020 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
