;;; NCL script to calculate relative humidity or dewpoint from
;;; specific humidity, temperature, and pressure

print(systemfunc("date"))

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

qfile="/glade/scratch/mcginnis/huss-10k.nc"
tfile="/glade/scratch/mcginnis/tmean-10k.nc"
pfile="/glade/scratch/mcginnis/sps.Daymet.fixed.NAM-22i.nc"
outfile="~/glade-scratch/mono-hurs-10k.nc"
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

print(systemfunc("date")+tochar(9)+"addfiles")

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
print(systemfunc("date")+tochar(9)+"create out")

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
print(systemfunc("date")+tochar(9)+"copy vars to out")

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

print(systemfunc("date")+tochar(9)+"allocate result")


result = new(getfilevardimsizes(fint, tname), float, missval) 
result@missing_value = result@_FillValue


print(systemfunc("date")+tochar(9)+"read p")
pstatic = .not. isfilevardim(finp, pname, "time")
if(pstatic) then
  p = conform(result, finp->$pname$, (/1,2/))
else 
  p = finp->$pname$
end if



;;nt = getfilevardimsizes(fint, "time")

;;do i=0, nt-1
print(systemfunc("date")+tochar(9)+"read t")

;;  t = fint->$tname$(i,:,:)
  t = fint->$tname$
  if(tunits .eq. "degC") then
    t = t + 273.15
  end if

;;  if(.not.pstatic) then
;;    p = finp->$pname$(i,:,:)
;;  end if

;;  q = finq->$qname$(i,:,:)
print(systemfunc("date")+tochar(9)+"read q")

  q = finq->$qname$
  x = q/(1-q)  ;; x = mixing ratio

print(systemfunc("date")+tochar(9)+"calc start")

  ;; NB: function output is % RH; dewtemp_trh wants % RH, but CF unit is 1
  temphurs = relhum(t,x,p)
  if(dewpt) then
    tempdewpt = dewtemp_trh(t,temphurs)
;;    result(i,:,:) = tofloat((/tempdewpt/))
  else 
print(systemfunc("date")+tochar(9)+"assign result")

result = tofloat(temphurs/100)
;;    temphurs = temphurs/100
;;    result(i,:,:) = tofloat((/temphurs/))
  end if

;;   if (verbose .and. (i%50) .eq. 0) then 
;;     print(""+i+tochar(9)+systemfunc("date")) 
;;   end if 
;;end do

print(systemfunc("date")+tochar(9)+"calc finish")



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

print(systemfunc("date")+tochar(9)+"write result")

fout->$outname$ = result

print(systemfunc("date")+tochar(9)+"done")

exit

;; Copyright 2020 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
