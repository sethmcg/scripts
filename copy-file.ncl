;;; NCL template script; copies infile to outfile

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Note: infile and outfile should be specified as command-line arguments to NCL:

;; ncl -x infile=\"$infile\" outfile=\"$outfile\"  copy-file.ncl


infile = "tas_CRCM_1979010103.nc"
outfile = "tas_fixed.nc"

system("rm "+outfile)
fin = addfile(infile, "r")
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited



; copy global attributes

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fin@$att_names(i)$
end do


; copy variables

var_names = getfilevarnames (fin) ;

do i = 0,dimsizes(var_names)-1      
  fout->$var_names(i)$ = fin->$var_names(i)$
end do


delete(fin)
delete(fout)
delete(var_names)
delete(att_names)
