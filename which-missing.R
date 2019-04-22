args<-commandArgs(TRUE)
library(ncdf4)

#infile = "prec.hist.GFDL-ESM2M.WRF.day.NAM-44i.raw.nc"

infile  = args[1]

varname <- unlist(strsplit(basename(infile), ".", fixed=TRUE))[1]

fin<-nc_open(infile)

bad <- c()

for(i in 1:fin$dim$time$len){
#  cat(i," ")
  slab <- ncvar_get(fin, varname, start=c(1,1,i), count=c(-1,-1,1))
  if(all(is.na(slab))) { bad <- c(bad,i) }
}

cat(infile, "\t", length(bad),"\n")
cat(bad, "\n")
