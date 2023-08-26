args<-commandArgs(TRUE)

library(ncdf4)
library(scales)

north   = as.numeric(args[1])
south   = as.numeric(args[2])
east    = as.numeric(args[3])
west    = as.numeric(args[4])
infile  = args[5]


fin<-nc_open(infile, readunlim=FALSE)

lat<-ncvar_get(fin,"lat")
lon<-ncvar_get(fin,"lon")

plot(lon, lat,
     xlim=expand_range(c(west, east), 0.2),
     ylim=expand_range(c(south, north), 0.2))
rect(west, south, east, north, border='red')
inds <- which(lat > -100, TRUE)
text(lon, lat, inds[,1], pos=2)
text(lon, lat, inds[,2], pos=3)

inbox <- lon <= east & lon >= west & lat <= north & lat >= south


## dim order = lon, lat, [time]
     
bbox <- apply(which(inbox, TRUE), 2, range)

dnames <-  sapply(fin$var$lat$dim, `[[`, "name")

## -F because we're counting indexes starting at 1, not 0
     
xsub <- paste0("-F -d ", dnames[1], ',', bbox[1,1], ',', bbox[2,1])
ysub <- paste0("-F -d ", dnames[2], ',', bbox[1,2], ',', bbox[2,2])

cat(xsub, ysub, "\n")

