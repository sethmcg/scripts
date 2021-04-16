args<-commandArgs(TRUE)

library(ncdf4)


north   = as.numeric(args[1])
south   = as.numeric(args[2])
east    = as.numeric(args[3])
west    = as.numeric(args[4])
infile  = args[5]

fin<-nc_open(infile, readunlim=FALSE)

lat<-ncvar_get(fin,"lat")
lon<-ncvar_get(fin,"lon")

plot(lon, lat, xlim=c(west, east), ylim=c(south, north))

inbox <- lon <= east & lon >= west & lat <= north & lat >= south


## dim order = lon, lat, [time]

bbox <- apply(which(inbox, TRUE), 2, range)

dnames <-  sapply(fin$var$lat$dim, `[[`, "name")

xsub <- paste0("-d ", dnames[1], ',', bbox[1,1], ',', bbox[2,1])
ysub <- paste0("-d ", dnames[2], ',', bbox[1,2], ',', bbox[2,2])

cat(xsub, ysub, "\n")

