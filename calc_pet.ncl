;;; NCL script to calculate PET (potential evaptranspiration) from
;;; other variables

load "fao56-crop.ncl"

;;; This script requires nine different input variables, so rather
;;; than specify all of them via the command line, we assume that the
;;; files are all located in . and that they all have the same name
;;; except for the leading variable.

;;; Specify the filename as a command-line argument, e.g.:
;;; ncl fname=\"_RCA4_209012-210011.nc\" calc_pet.ncl

;;; Note: hurs and dewpt can be calculated from huss, ps, and tas
;;; using calc_humid.ncl

;fname = ".nc"

;; read in all the input variables

fin = addfile("sfcWind"+fname, "r")
wind = fin->sfcWind
delete(fin)

fin = addfile("hurs"+fname, "r")
hurs = fin->hurs
delete(fin)

fin = addfile("dewpt"+fname, "r")
dewpt = fin->dewpt
delete(fin)

fin = addfile("ps"+fname, "r")
ps = fin->ps
delete(fin)

fin = addfile("rsds"+fname, "r")
rsds = fin->rsds
delete(fin)

fin = addfile("rlds"+fname, "r")
rlds = fin->rlds
delete(fin)

fin = addfile("rsus"+fname, "r")
rsus = fin->rsus
delete(fin)

fin = addfile("rlus"+fname, "r")
rlus = fin->rlus
delete(fin)

fin = addfile("tas"+fname, "r")
tas = fin->tas



;; PET is calculated using Shuttleworth's 1993 version of the Penman
;;  equation:

;;;  E_mass= (m * R_n + gamma_ * 6.43 * (1+0.536 * U_2 ) * delta_e) /
;;;          (lambda_v * (m + gamma_ ))


;;m = Slope of the saturation vapor pressure curve (kPa K−1)      	 
m = satvpr_slope_fao56 (tas, 1)

;; R_n = Net irradiance (MJ m−2 day−1)
R_n = (rsds + rlds - rsus - rlus) * 86400e-6

;; gamma_ = psychrometric constant = 0.0016286 * P_kPa / lambda_v (kPa K−1)
gamma_ = psychro_fao56(ps, 1)

;; U_2 = wind speed (m s−1)
;; 	 == sfcWind

;; delta = vapor pressure deficit (kPa)
delta_e	= (1 - hurs) * satvpr_tdew_fao56(dewpt, 1)

;; lambda_v = latent heat of vaporization for water (MJ kg−1)
lambda_v = 2.257

evspsblpot = (m * R_n + gamma_ * 6.43 * (1+0.536 * wind ) * delta_e) / (lambda_v * (m + gamma_ ))
evspsblpot = evspsblpot / 86400  ;; convert units from mm/day to kg m-2 s-1



outfile = "evspblpot" + fname


system("rm -f "+outfile)
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


;; copy/set global attributes - note fin is tasfile

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  if(att_names(i) .eq. "history_of_appended_files") then
    ;; skip
  else if(att_names(i) .eq. "title") then
    fout@title = str_sub_str(fin@title, "2m Temperature", "Potential Evapotranspiration")
  else
    fout@$att_names(i)$ = fin@$att_names(i)$
  end if
  end if
end do

history = "Created " + systemfunc("date") + " by "+systemfunc("whoami")+"@"+systemfunc("hostname")+" using NCL script calc_pet.ncl from input files [dewpt hurs ps rlds rlus rsds rsus sfcWind tas]"+fname



fout@history = history

;; update unique tracking_id
fout@tracking_id = systemfunc("uuidgen")


; copy variables

var_names = getfilevarnames (fin) ;

do i = 0,dimsizes(var_names)-1      
if (var_names(i) .ne. "tas") then   
  fout->$var_names(i)$ = fin->$var_names(i)$
end if
end do


;; get rid of superfluous attributes

delete_VarAtts(evspsblpot, -1)

evspsblpot@long_name = "Potential Evapotranspiration"
evspsblpot@standard_name = "water_potential_evaporation_flux"
evspsblpot@units = "kg m-2 s-1"

varatts = (/"missing_value", "_FillValue", "coordinates", "grid_mapping"/)

do i = 0,dimsizes(varatts)-1      
  evspsblpot@$varatts(i)$ = fin->tas@$varatts(i)$
end do


fout->evspsblpot = evspsblpot

exit


;; Copyright 2015 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
