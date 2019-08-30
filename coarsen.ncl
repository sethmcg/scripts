;;; NCL coarsening script: regrids high resolution data to lower
;;; resolution using hi2lores 

;;; % missing allowed is 50%


;; Specify arguments to script via command-line:
;; ncl coarsen.ncl datafile=\"$data\" gridfile=\"$grid\" outfile=\"$out\" varname=\"$var\"

;; optional: critpc = fraction of inputs missing to count output as
;;   missing (defaults to 50%)


system("rm -f "+outfile)
fin = addfile(datafile, "r")
gin = addfile(gridfile, "r")
fout = addfile(outfile, "c")

if(isvar("critpc")) then
  crithist = " critpc="+critpc
else
  critpc = 50
end if



;; time-varying or static?

if(any("time" .eq. getfiledims(fin))) then
  filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited
end if


data = fin->$varname$
olat = gin->lat
olon = gin->lon


opt=True
opt@critpc = critpc
odata = area_hi2lores_Wrap(data&lon, data&lat, data, False, 1, olon, olat, opt)
odata@missing_threshold = critpc+"%"

fout->$varname$ = odata

;; deal with grid bounds

hasbounds = False
latbname = ""
lonbname = ""

if(isatt(olat, "bounds")) then
  hasbounds = True
  latbname = olat@bounds
  fout->$latbname$ = gin->$latbname$
end if

if(isatt(olon, "bounds")) then
  hasbounds = True
  lonbname = olon@bounds
  fout->$lonbname$ = gin->$lonbname$
end if


; copy global attributes

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fin@$att_names(i)$
end do

histring = systemfunc("date")
histring = histring + ": ncl coarsen.ncl"
histring = histring + " datafile="+datafile
histring = histring + " gridfile="+gridfile
histring = histring + " outfile="+outfile
if(isvar("crithist")) then
  histring = histring + crithist
end if

fout@history = histring + inttochar(10) + fout@history

fout@tracking_id = systemfunc("uuidgen")


; copy variables

invars = getfilevarnames (fin) ;

skipme = (/"lat", "lon", latbname, lonbname, varname/)

do i = 0,dimsizes(invars)-1   
  if(any(invars(i) .eq. skipme)) then
;;    print(""+invars(i)+": skipped")
  else
;;  print(""+invars(i))
    fout->$invars(i)$ = fin->$invars(i)$
end if
end do
