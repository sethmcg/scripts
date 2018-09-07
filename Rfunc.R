if(! exists("tim.colors")){
  library(fields)
}
  
rescale <- function(x, y){
    rx = range(x)
    ry = range(y)
    z = ((x - rx[1])/(rx[2]-rx[1]))*(ry[2]-ry[1])+ry[1]
    return(z)
}

colorize <- function(z,ncol=101,zlim=range(z)){
    zcolor = pmax( zlim[1], pmin( z, zlim[2]))
    colors = tim.colors(ncol)
    zcolor = colors[rescale(zcolor,c(1,ncol))]
    return(zcolor)
}

colorbar <- function(xmin,xmax,ymin,ymax,z,ncol=100,nlab=11,zlim=range(z),...){    
    zbar = rbind(rescale(1:ncol,zlim))
    x = c(xmin,xmax)
    y = seq(ymin,ymax,length.out=ncol+1)
    image(x,y,zbar,col=tim.colors(ncol),add=TRUE)
    labels = signif(seq(min(zlim),max(zlim),len=nlab),digits=2)
    text(xmax,seq(ymin,ymax,len=nlab),labels,pos=4,...)
}

rv <- function(){
  par(bg="black", fg="white", col.axis="white", col.lab="white", col.main="white", col.sub="white")
}

overlays <- function(){
  contour(glon, glat, is.na(landmask), lev=0.5, add=TRUE, label='', col="black")
  contour(glon, glat, is.na(station.density), lev=0.5, add=TRUE, col="white", lty=3, label='', lwd=2)
  contour(glon, glat, station.density, lev=50, add=TRUE, col="white", lwd=2, label='')
}

scatterplot<-function(x,y,z,zlim=range(z),...){
  zcol = colorize(z,zlim=zlim)
  plot(x, y, col=zcol,...)
  colorbar(330,335,10,75,station$data,zlim=zlim)
}
