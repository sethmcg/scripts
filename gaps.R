library(ncdf4)

space <- "    "

args <- commandArgs(TRUE)

infile <- args[1]
fin    <- nc_open(infile)
time   <- ncvar_get(fin, "time")

dtab <- table(round(diff(time), 3))

gaps <- dtab[-which(names(dtab)=="1")]

if(length(gaps) > 0) {
  result <- paste(names(dtab), dtab, sep=":", collapse=space)
  cat(infile, space, result, "\n")
}

