;;; NCL script to check for gaps in time coordinate

;;; slightly tweaked version: prints output in terse format for
;;; feeding to nctimefill script I wrote

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


;; Uncomment the following variables to hard-code them here, or pass values in via command-line.  Example:

;; ncl -x varname=\"prtot\" infile=\"table6/prtot_CRCM_climatology.nc\" outfile=\"plots/prtot_CRCM_s0\" title=\"CRCM\" timestep=\"0\" narccap-generic-plot.ncl

; infile =     "table6/prtot_CRCM_climatology.nc"

; tab = "	"

if (.not. isvar("step")) then
  step = 0.125
end if

;; open file
fin = addfile(infile, "r")

time = fin->time
nt = dimsizes(time)

;print("file:"+tab+infile)
;print("nsteps:"+tab+(nt-1))
;print("units:"+tab+time@units)
;print("first:"+tab+ut_date(time(0)))
;print("last: "+tab+ut_date(time(nt-1)))
;print("gaps")
;print("----")

;print(""+infile)
 
 do i=0,nt-2
 
 delta = time(i+1)-time(i)
 
 if (delta .ne. step) then

   print (""+i+" "+(i+1)+" "+(delta/step - 1)+" "+step+" "+infile)

;   print("----")
;   print (""+i+tab+time(i)+tab+ut_date(time(i)))
;   print (""+(i+1)+tab+time(i+1)+tab+ut_date(time(i+1)))
;;   print ("    "+delta)
;   print (tab+delta)
;   print("----")
 
 end if
 
end do

;print("====")
;print("")
 
exit



;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
