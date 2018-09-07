args<-commandArgs(TRUE)

library(ncdf4)


tlat    = as.numeric(args[1])
tlon    = as.numeric(args[2])
infile  = args[3]
outfile = args[4]

fin<-nc_open(infile)

lat<-ncvar_get(fin,"lat")
lon<-ncvar_get(fin,"lon")

dy = abs(lat - tlat)
dx = abs(lon%%360 - tlon%%360)

#image(dx+dy)

ind = arrayInd(which.min(dx+dy),dim(lat))-1

xind = paste0("-d xc,", ind[1], ",", ind[1])
yind = paste0("-d yc,", ind[2], ",", ind[2])

command = paste("ncks", xind, yind, infile, outfile)

#print(command)

system(command)



