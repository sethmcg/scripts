;;; NCL script that adds lat2d & lon2d to a file that only has lat & lon

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;; Note: infile should be specified as command-line arguments to NCL:

;; ncl -x infile=\"$infile\" outfile=\"$outfile\"  copy-file.ncl

print(""+infile)

fin = addfile(infile, "w")

lat = fin->lat
lon = fin->lon

if(any(lon.lt.0)) then
  lon(ind(lon.lt.0)) = lon(ind(lon.lt.0))+360
end if


nx = dimsizes(lon)
ny = dimsizes(lat)

lon2d = conform_dims((/ny,nx/), lon, 1)
lat2d = conform_dims((/ny,nx/), lat, 0)

lon2d!0 = "lat"
lon2d&lat = lat
lon2d!1 = "lon"
lon2d&lon = lon

lat2d!0 = "lat"
lat2d&lat = lat
lat2d!1 = "lon"
lat2d&lon = lon

lon2d@units = "degrees_east"
lon2d@long_name = "longitude"
lon2d@standard_name = "longitude"

lat2d@units = "degrees_north"
lat2d@long_name = "latitude"
lat2d@standard_name = "latitude"

fin->lat2d = lat2d
fin->lon2d = lon2d

;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
