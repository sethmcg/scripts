;;; Prints value and corresponding date/time for time coordinate

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"

;; Pass filename in via command-line.  Also accepts a flag for
;; printing time_bounds values.  Example:

;; ncl -x infile=\"table6/prtot_CRCM_climatology.nc\" bounds=True datetime.ncl

format = "%Y-%N-%D %H:%M:%S"

if (.not. isvar("bounds")) then
  bounds = False
end if


;; open file
fin = addfile(infile, "r")


time = fin->time
nt = dimsizes(time)

if(time@calendar .eq. "proleptic_gregorian") then
  time@calendar = "gregorian"
  print("WARNING: treating proleptic_gregorian calendar as gregorian")
end if

date = ut_string(time,format)


if(bounds) then
  time_bnds = fin->$time@bounds$
  time_bnds@units = time@units
  time_bnds@calendar = time@calendar

  date_bnds_lo = ut_string(time_bnds(:,0),format)
  date_bnds_hi = ut_string(time_bnds(:,1),format)
end if


tab  = "	"

do i=0,nt-1

  result = ""+i+tab+time(i)+tab+date(i)
  if (bounds) then
    result = result + tab + date_bnds_lo(i) + tab + date_bnds_hi(i)
  end if

  print(""+result)

end do

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
