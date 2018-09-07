;;; NCL coarsening script: regrids high resolution data to lower resolution using hi2lores
;;; %missing allowed is 50%


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Specify arguments to script via command-line:
;; ncl coarsen.ncl datafile=\"$data\" gridfile=\"$grid\" outfile=\"$out\" varname=\"$var\"


system("rm "+outfile)
fin = addfile(datafile, "r")
gin = addfile(gridfile, "r")
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited

data = fin->$varname$
olat = gin->lat
olon = gin->lon

opt=True
opt@critpc = 50
odata = area_hi2lores_Wrap(data&lon, data&lat, data, False, 1, olon, olat, opt)
odata@max_missing_allowed_in_coarsening = "50%"

fout->$varname$ = odata



; copy global attributes

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fin@$att_names(i)$
end do

fout@history = systemfunc("date")+": ncl datafile="+datafile+" gridfile="+gridfile+" outfile="+outfile+" coarsen.ncl"+inttochar(10)+fout@history


; copy variables

var_names = getfilevarnames (fin) ;

do i = 0,dimsizes(var_names)-1   
  if (var_names(i) .eq. "lat"    .or. \
      var_names(i) .eq. "lon"    .or. \
      var_names(i) .eq. varname       \
      ) then
        print(""+var_names(i)+": skipped")
  else
  print(""+var_names(i))
    fout->$var_names(i)$ = fin->$var_names(i)$
end if
end do


delete(fin)
delete(fout)
delete(var_names)
delete(att_names)


