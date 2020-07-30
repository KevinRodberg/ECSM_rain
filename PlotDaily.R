#--
#   package management:
#     provide automated means for first time use of script to automatically 
#	  install any new packages required for this code, with " calls 
#	  wrapped in a for loop.
#--
pkgChecker <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}

list.of.pkgs <-  c("reshape2","readr","dplyr", "data.table","readxl",
                   "rgeos","sp","dismo","lattice","gridExtra","rasterVis","maptools",
                   "raster","fields","automap","gstat","rgdal",
                   "future","listenv","geosphere")

suppressWarnings(pkgChecker(list.of.pkgs))
# 
# library(reshape2)
# library(raster)
# library(dplyr)
# library(rasterVis)
# library(sp)
# library()

geomSeries <- function(base, max) {
  (base^(0:floor(log(max, base)))-1)/100
}
#-------------------------------------------------
# Set up county boundry shapefile for overlay 
# on raster maps
#-------------------------------------------------
gClip <- function(shp, bb) {
  if (class(bb) == "matrix")
    b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else
    b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

fixDecimals <- function(DF,decPlaces){
  is.num <-sapply(DF,is.numeric)
  DF[is.num] <- lapply(DF[is.num], round,decPlaces)
  return(DF)
}


basePath = '//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/ECSM/Data'
# Clear raster Stack
rasStack <-stack()

# Calculate raster extents
# res =1000
# xmin = 672000.00
# ymin = 143000.00
# nrows = 1059
# ncols = 313


res = 5280*2
xmin = 664767
ymin = 292038
nrows = 88
ncols = 30

xmax = xmin+ (ncols*res)
ymax = ymin + (nrows*res)
UTM17m = CRS("+init=epsg:26917")
HARNSP17ft  = CRS("+init=epsg:2881")

WMDbnd.Path <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)

ECSMras <- raster(resolution=res, xmn=xmin, xmx=xmax,ymn=ymin, ymx=ymax,crs=HARNSP17ft)
proj4string(ECSMras )=CRS("+init=epsg:2881")

rasExt <- extent(ECSMras)
clpBnds2 <- gClip(WMDbnd, rasExt)

# year = 2012
yearList <- seq(1987,2018, by=1)
setwd(basePath)
for (year in yearList){
  yearPath = paste('./rasterPlots',year,sep='/') 
  cat(yearPath,sep='\n')
  dir.create(yearPath, recursive=TRUE, showWarnings=TRUE)
}
for (year in yearList){
  # Clear raster Stack
  rasStack <-stack()
  yearPath= paste(basePath,'rasterPlots',year,sep='/')
  fileDaily <- paste(basePath,"/wmmRain",year,".csv",sep="")
  #wmdRain <-read.table(fileDaily, header=TRUE, sep="\t",  na.strings="0.00")
  wmdRain <-read.csv(fileDaily, header=TRUE, na.strings="0.00")
  # wmdRain$X <-(xmin+(res/2))+ (wmdRain$Column_*(res))
  # wmdRain$Y <-(ymin+(res/2))+ ((nrows+1-wmdRain$Row)*(res))
  coordinates(wmdRain)=~Xcoord+Ycoord
  #	coordinates(wmdRain)=~Longitude+Latitude
  proj4string(wmdRain)=HARNSP17ft 
  # 
  # 	write.csv(cbind(wmdRain[wmdRain$Annual==0.0,]$variable,
  # 	                wmdRain[wmdRain$Annual==0.0,]@coords),'H:/badWMDcells.csv')
  dailyColNames = names(wmdRain[,-c(1,2,length(names(wmdRain)))])
  bom = seq(as.Date(paste0(year,"/1/1")), by = "month", length.out = 12)-as.Date(paste0(year,"/1/1"))+1
  eom =seq(as.Date(paste0(year,"/2/1")), by = "month", length.out = 12)-as.Date(paste0(year,"/1/1"))
  # create points from each column and convert 
  # them to rasters and add each to the raster stack
  for (daily in dailyColNames ){
    ECSMras_1 <-rasterize(wmdRain,ECSMras ,daily,fun=mean)
    rasStack <- stack(rasStack, ECSMras_1)
  }
  
  names(rasStack)<-dailyColNames 
  # Eliminate Scientific Notation
  options(scipen=10000)
  myTheme=rasterTheme(region=brewer.pal('Blues', n=9))

  #for (i in 1:nlayers(rasStack[[1:5]])){
  for (i in 1:nlayers(rasStack)){
    filename=paste(yearPath,"/WMD_",year,names(rasStack)[i],".png",sep="")
    rasPlt <- rasStack[[i]]
    rasPltName <-names(rasStack)[i]
    # pntsPlt <-pointList[[i]]
    NRDplot=( levelplot(rasPlt, par.settings=myTheme, 
                        main=paste0(rasPltName, " ",year," WMD"), 
                        at=geomSeries(base=2, max=500), xlab = NULL, margin=F,
                        layout=c(1,1),contour=FALSE) +
                latticeExtra::layer(sp.polygons(clpBnds2, col='darkgray')))
    
    trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
    # print(Biasplot, split    = c(1,1,4,1),more=TRUE)
    print(NRDplot)
    # print(NRDplot, split     = c(2,1,1,1),more=TRUE)
    # print(NRDBiasplot, split = c(3,1,4,1),more=TRUE)
    # print(diffRasPlt, split  = c(4,1,4,1))
    dev.off()
  }
}

