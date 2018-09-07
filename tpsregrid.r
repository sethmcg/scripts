library(ncdf4, quietly=TRUE)
library(fields, quietly=TRUE)

tpsRegridNetCDF <- function(infile, varname, outfile, gridfile,
                            maskfile="", mask.var="mask", radius=2.5,
                            verbose=FALSE, diagnostics=FALSE,
                            recalc.tps=FALSE, load.tps=FALSE,
                            save.tps=FALSE, tpsfile="",
                            log=FALSE, trace=NULL, fix.lon=TRUE ) {

  ## Regrids netcdf data using the thin-plate-spline algorithm
  ##
  ## Args:
  ##   infile:   name of input netcdf file (string)
  ##   varname:  name of data variable in infile (string)
  ##   outfile:  name of output netcdf file (string)
  ##   gridfile: name of netcdf file w/ lat/lon to project to (string)
  ##   maskfile: name of netcdf file with location mask (string)
  ##   radius:   radius for TPS fit.
     
  ##   verbose:          print progress messages? (bool)
  ##   diagnostics:      interactive sanity-check as we go? (bool)
     
  ##   recalc.tps:       recalculate fit object each timestep? (bool)
  ##   load.tps:         load fit object from file? (bool)
  ##   save.tps:         save fit object to file?   (bool)
  ##   tpsfile:  name of file to load/save fit object from/to (string)

  ##   log:      perform log transformation (bool)

  ##   trace:    values below this threshold set to zero (numeric)

  ##   fix.lon:  normalize longitudes to 0:360
  
  ## Returns: nothing

  
  ## Infile and gridfile are presumed to have horizontal coordinate
  ## dimensions xc & yc with accompanying 2D auxiliary coordinate
  ## variables lat & lon.

  ## Note that gridfile MUST contain a data variable of some kind; the
  ## ncdf4 library can't handle netcdf files that contain only
  ## coordinate variables.
  
  ## Maskfile should have exactly the same coordinates as gridfile; set
  ## maskfile to "" for no masking.
  
  ## A radius value of 2.5 (degrees) is a good default for half-degree /
  ## 50-km data.

  ## Non-essential auxiliary information like map projection, time
  ## bounds, and scalar level coordinate does not get copied into
  ## outfile, so outfile will not be cf-compliant.
  
  ## If outfile already exists, it will be clobbered with no warning.
  
  ## The TPS fit depends only on the geometry of the gridpoints, so if
  ## the input coordinates are identical, it can be reused instead of
  ## recalculated, which saves a lot of time.  The only time you
  ## wouldn't want to do this is if there are missing values in the
  ## input data that change from timestep to timestep.
  
  ## If you want to take it one step further, you can reuse a saved fit
  ## object with different input files, as long as they all have
  ## exactly the same grid locations.  Note that we only load/save a
  ## single fit object to file, so this is incompatible with the
  ## recalc.tps option.

  ## A log transform can be used to ensure that the interpolated values
  ## are always positive (e.g., for precipitation).  The interpolation
  ## is performed on ln(data) and exp(result) is returned.

  ## Trace is the threshold of negligibility.  Values below this
  ## threshold are set to trace/10 before interpolation, and to 0 after
  ## interpolation.  Used in conjunction with the log option, it means
  ## you don't have to worry about zero or negative values.  A good
  ## trace value for precipitation is 0.01.  Whether mm or inches, over
  ## a period of days or longer, a trace threshold of 0.01 will be
  ## negligible.  Meteorological data will often have precip floored to
  ## zero below 0.01 inches.
   
  if (recalc.tps && (load.tps || save.tps)) {
    print("tpsRegridNetCDF error: recalc.tps option incompatible with load.tps and save.tps.  Aborting.")
    return()
  }
  
  if (verbose) print("Reading input data...")
    
  if (diagnostics) par(ask=TRUE)


  vname <- list()
  
  ## Get coordinates from input data file
  
  fin <- nc_open(infile)

  for (v in fin$dim){
    va <- ncatt_get(fin, v$name)
    if (identical(va$standard_name, "time")) {
      time <- ncvar_get(fin,v$name)
      vname[["time"]] <- v$name
    }
    if (identical(va$standard_name, "latitude")) {
      iy <- ncvar_get(fin,v$name)
      vname[["iy"]] <- v$name
    }
    if (identical(va$standard_name, "longitude")) {
      ix <- ncvar_get(fin,v$name)
      vname[["ix"]] <- v$name
    }
  }
  
  static <- is.null(vname[["time"]])
  if (static) {
    nt <- 1
  } else {
    nt <- length(time)
  }


  if(is.null(vname[["ix"]]) || is.null(vname[["iy"]])) {
    for (v in fin$var){
      va <- ncatt_get(fin,v$name)
      if (identical(va$standard_name, "latitude")) {
        ilat <- ncvar_get(fin,v$name)
        ix <- v$dim[1][[1]]$vals
        iy <- v$dim[2][[1]]$vals
        vname[["ilat"]] <- v$name
        vname[["ix"]] <- v$dim[[1]]$name
        vname[["iy"]] <- v$dim[[2]]$name
      }
      if (identical(va$standard_name, "longitude")) {
        ilon <- ncvar_get(fin,v$name)
        vname[["ilon"]] <- v$name
      }
    }
  } else {
    ilat <- t(matrix(iy,length(iy),length(ix)))
    ilon <-   matrix(ix,length(ix),length(iy))    
  }
  

  ## Read in data
  ## this slight convolution avoids dropping of degenerate dimensions
  var <-  fin$var[[varname]]
  data <- ncvar_get(fin,var)

  dim(data) <- var$varsize

  if(static){
    dim(data) <- c(dim(data),1)
  }

  
  ## Get target coordinates from grid file
  
  gin <- nc_open(gridfile)

  for (v in gin$dim){
    va <- ncatt_get(gin, v$name)
    if (identical(va$standard_name, "latitude")) {
      yc <- ncvar_get(gin,v$name)
      vname[["yc"]] <- v$name
    }
    if (identical(va$standard_name, "longitude")) {
      xc <- ncvar_get(gin,v$name)
      vname[["xc"]] <- v$name
    }
  }

  twod <- is.null(vname[["xc"]]) || is.null(vname[["yc"]]) 

  if(twod){
    for (v in gin$var){
      va <- ncatt_get(gin,v$name)
      if (identical(va$standard_name, "latitude")) {
        olat <- ncvar_get(gin,v$name)
        xc <- v$dim[1][[1]]$vals
        yc <- v$dim[2][[1]]$vals
        vname[["olat"]] <- v$name
        vname[["xc"]] <- v$dim[[1]]$name
        vname[["yc"]] <- v$dim[[2]]$name
      }
      if (identical(va$standard_name, "longitude")) {
        olon <- ncvar_get(gin,v$name)
        vname[["olon"]] <- v$name
      }
    }
  } else {
    olat <- t(matrix(yc,length(yc),length(xc)))
    olon <-   matrix(xc,length(xc),length(yc))    
  }

  if(fix.lon){
    ilon = ilon %% 360
    olon = olon %% 360
  }
  
  ## If you wanted to dump the data to R format, do something like this:
  ##  precip <- list( x=lon2d, y=lat2d, z=data, t=time)
  ##  save(precip, file="precip.rda")
  
  if (diagnostics) {
    ## check that data looks sane
    image.plot(ilon,ilat,data[,,1])
    title("raw data check")
  }
  
  ## reform the 2-D lat/lon values into N-by-2 arrays
  iloc <- cbind( c(ilon), c(ilat))
  oloc <- cbind( c(olon), c(olat))
  
  ## mask output locations
  if (maskfile != ""){
    min <- nc_open(maskfile)
    mask <- ncvar_get(min, mask.var)
    mask <- ifelse(mask == 1, TRUE, FALSE)
    oloc <- oloc[c(mask),]
  } else {
    mask <- rep(TRUE,length(olat))
  }
  
  
  ## create output array
  nt <- length(time)
  nx <- length(xc)
  ny <- length(yc)
  result <- array(dim=c(nx,ny,nt))
  
  
  ## load fit object
  if(load.tps) {
    if(verbose) print(paste("loading fit object from file ",tpsfile))
    load(tpsfile)
  }

  ## data transforms
  if(!is.null(trace)){
    data[which(data < trace)] = trace/10.0
  }

  if(log) data <- log(data)


  ## main interpolation loop
  if(verbose) print("Calculating...")

  for (t in 1:nt){
  
    ## missing values make tps unhappy, so subset inputs to non-missing
    ind <- c(!is.na(data[,,t]))
    x <- iloc[ind,]
    y <- c(data[,,t])[ind]
    ## Note: works b/c loc & data unroll in same order when flattened w/ c().
    
    if(verbose) print(paste("time = ",t,": ",(length(ind)-sum(ind))," missing",sep=""))
    
    ## visual double-check:
    if(diagnostics){
      quilt.plot(iloc, c(data[,,t]))
      title("iloc vs data slice")
      quilt.plot(x,y)
      title("x vs y")
    }
    
    if(recalc.tps || !exists("fit")) {
    
      if(diagnostics){
        ## check radius.  Ideally you want at least 20 neighbors per point
        look <- nearest.dist(x,x, delta=radius)
        NN <- diff(look@rowpointers) ## nonzero distances each row = num neighbors
        stats(NN)
        hist(NN)
        
        lookdim <- dim(data) 
        look <- matrix(NA, lookdim[1], lookdim[2])
        look[ind] <- NN
        look[look > 20] <- 20
        image.plot(ilon, ilat, look, zlim=c(1,20))
        title("neighbors")
      }
      if(verbose) print("calculating fit")
      fit <- fastTps(x, y, theta=radius)
    }
      
    
      ## interpolate fitted data to new locs.
      if(verbose) print("interpolating data")
      field <- array(dim=c(nx,ny))
      field[mask] <- predict(fit, xnew=oloc, ynew=y)
      result[,,t] <- field
  }

  if(log) result <- exp(result)

  if(!is.null(trace)){
    result[which(result < trace)] = 0
  }
  
  
  if(verbose) print("Finished interpolating; writing output...")
  
  ####################################
  ## And now we do all the output
  
  
  #####
  
  if (save.tps){
    if(verbose) print(paste("saving fit object to file ",tpsfile))
    save(fit,file=tpsfile)
  }
  
  #####
  
  if(verbose) print("creating netcdf output")

  ## convenience function
  copy.attributes <- function(fromfile, fromvar, tofile, tovar, skip=''){
    attlist <- ncatt_get(fromfile, fromvar)
    for (attname in names(attlist)) {
      if(!any(attname == skip)){
        ncatt_put(tofile, tovar, attname, attlist[[attname]])
      }
    }
  }

  
  ## define dimensions  
  ndim <- list()

  ndim[[vname[["xc"]]]] <- ncdim_def(vname[["xc"]], "dummy", xc)
  ndim[[vname[["yc"]]]] <- ncdim_def(vname[["yc"]], "dummy", yc)
  if(!static) {
    ndim[[vname[["time"]]]] <- ncdim_def(vname[["time"]], "dummy", time, unlim=TRUE)
  }
  ndim[[varname]] <- ndim
  if(twod){
    ndim[["lat"]] <- ndim[["lon"]] <- ndim[c(1,2)]
  }

  
  ## define variables
  nvar <- list()
  #!# Note: get missval from file...
  if(twod){
    nvar[["lat"]] <- ncvar_def(name=vname[["olat"]], units="dummy", dim=ndim[["lat"]], missval=1e20)
    nvar[["lon"]] <- ncvar_def(name=vname[["olon"]], units="dummy", dim=ndim[["lon"]], missval=1e20)
  }
  nvar[[varname]] <- ncvar_def(name=varname, units="dummy", dim=ndim[[varname]], missval=1e20)
 
  
  if(verbose) print(paste("writing netcdf file ",outfile))
  
  ## create file
  fout <- nc_create(outfile, nvar)

    
  ## write variables to file

  if(twod){
    ncvar_put(fout,vname[["olat"]],olat)
    ncvar_put(fout,vname[["olon"]],olon)
  }
  if(static){
    ncvar_put(fout,varname,result[,,1])
  } else {
    ncvar_put(fout,varname,result)
  }
 
 
  ## copy attributes

  copy.attributes(fin,varname,fout,varname,skip="grid_mapping")
  copy.attributes(gin,vname[["xc"]],fout,vname[["xc"]])
  copy.attributes(gin,vname[["yc"]],fout,vname[["yc"]])
  if(!static){
    copy.attributes(fin,vname[["time"]],fout,vname[["time"]])
  }
  if(twod){
    copy.attributes(gin,vname[["olat"]],fout,vname[["olat"]])
    copy.attributes(gin,vname[["olon"]],fout,vname[["olon"]])
  }

  copy.attributes(fin,0,fout,0)

  
  ## add history entry
  
  history <- ncatt_get(fin, 0, attname="history")$value
  history <- paste(history, "\n",date(),":", sep="")
  history <- paste(history, " regridded data from ",basename(infile), sep="")
  history <- paste(history, " to grid from ",basename(gridfile), sep="")
  history <- paste(history, " using R function fastTps() in library 'fields'", sep="")
  if(load.tps){
    history <- paste(history," using fit from file ",basename(tpsfile), sep="")  
  }else{
    history <- paste(history," with theta=",radius, sep="")  
  }
  
  ncatt_put(fout, 0, "history", history)
  
  
  ## done
  nc_close(fout)
  if(verbose) print("done!")
  
  return()
}

# Copyright 2010-2012 Univ. Corp. for Atmos. Research
# Author: Seth McGinnis, mcginnis@ucar.edu
