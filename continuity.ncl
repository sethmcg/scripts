;;; NCL script to check for gaps in time coordinate

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 
function ut_date(time[1])
local date, result
begin
  date = ut_calendar(time,-5)
  result = ""+date(0,0)+"/"+date(0,1)+"/"+date(0,2)+" "+date(0,3)+":"+date(0,4)+":"+date(0,5)
  return(result) 
end


;; Pass in input filename via command-line.  You can also specify timestep length, which defaults to 0.125 (3 hours, 1/8 day)

;; ncl -x varname=\"prtot\" infile=\"table6/prtot_CRCM_climatology.nc\" step=1 continuity.ncl 


 tab = "	"

if (.not. isvar("step")) then
  step = 0.125
end if

;; open file
fin = addfile(infile, "r")

time = fin->time
nt = dimsizes(time)

print(""+infile)
 
 do i=0,nt-2
 
 delta = time(i+1)-time(i)
 
 if (delta .ne. step) then

   print("----")
   print (""+i+tab+time(i)+tab+ut_date(time(i)))
   print (""+(i+1)+tab+time(i+1)+tab+ut_date(time(i+1)))
   print (tab+delta)
 
 end if
 
end do

print("====")
print("")
 
exit



;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
