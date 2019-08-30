args<-commandArgs(TRUE)
library(ncdf4)
library(PCICt)

infile  = args[1]
## For testing:
# infile = "~/DATA/cordex/data/raw/NAM-44/day/WRF/MPI-ESM-LR/rcp85/tasmin/tasmin.rcp85.MPI-ESM-LR.WRF.day.NAM-44.raw.nc"


## Planned for future: 
# Pass further arguments to control output:
# quiet: don't print filename
# range: print ranges instead of enumerating values
# count, date, index: print count/date/index of missing timesteps 

fname <- basename(infile)

varname <- unlist(strsplit(fname, ".", fixed=TRUE))[1]

nc <- nc_open(infile)

bad <- c()

## Looping over timesteps is MUCH faster than reading the whole thing in
for(i in 1:nc$dim$time$len){
  slab <- ncvar_get(nc, varname, start=c(1,1,i), count=c(-1,-1,1))
  if(all(is.na(slab))) { bad <- c(bad,i) }
}


if(length(bad) > 0) {


## Condense vector into runs of sequential values

headtail <- function(x){
  if(length(x) > 1) {c(head(x,1), tail(x,1))} else {x}
}

condense <- function(x){
  runs <- split(x, c(0, cumsum(diff(x) > 1)))
  names(runs) <- NULL
  lapply(runs, headtail)
}


runs <- condense(bad)



## Print count of missing timesteps

cat(fname, ":\t", length(bad)," missing timesteps\n", sep="")



## Print indices of missing timesteps
indexes <- sapply(runs, paste, collapse="-")

cat("index:\t", paste(indexes, collapse=", "), "\n")


## Print dates of missing timesteps

time <- ncvar_get(nc, "time")
runtimes <- rapply(runs, function(x){time[x]}, how="replace")

calendar <- ncatt_get(nc, "time", "calendar")$value
ustring  <- strsplit(ncatt_get(nc, "time", "units")$value, " ")[[1]]
tunits   <- ustring[1]
epoch    <- ustring[3]
if(length(ustring) > 3){ epoch <- paste(epoch, ustring[4]) }
   
timetodate <- function(x){
  conversions <- c(days=24*60*60, hours=60*60, minutes=60, seconds=1)
  if(! tunits %in% names(conversions)){
    stop("Conversion factor from ", tunits, " to seconds for PCICt not known.")
  }
  tdate <- as.PCICt(x*conversions[tunits], cal=calendar, origin=epoch)
  as.character(tdate, format="%Y-%m-%d")
}

dates <- sapply(rapply(runtimes, timetodate, how="replace"), paste, collapse=" to ")

cat("dates:\t", paste(dates, collapse="; "), "\n\n")

}
