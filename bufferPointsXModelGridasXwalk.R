#--------------------------------------------------------------------------
#  Program: bufferPointsXModelGrid.R
#  Author:  Kevin A. Rodberg, Science Supervisor, Water Supply Bureau
#           5/14/2018; 7/29/2020; 
#           9/16/2024: results were coming out striped other than ET
#                      see NRDsub change
#  Purpose: Crosswalk/Lookup table for NexRad Pixel to Model Row and Column
#           Using Nearest Neighbor searching algorithm (FNN::get.knn)
#--------------------------------------------------------------------------

list.of.pkgs <-  c("data.table","terra","sp","rgdal","sf",
                   "FNN","raster","future","tcltk2")
#-------------------------------------------
# Set model name you are working with today:
#-------------------------------------------
model = 'ECSM'

pkgChecker <- function(x){
	#--
	#   package management:
	#     provide automated means for first time use of script to automatically 
	#	  install any new packages required for this code, with " calls 
	#	  wrapped in a for loop.
	#--
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}
suppressWarnings(pkgChecker(list.of.pkgs))
rgdal::set_thin_PROJ6_warnings(TRUE)
readPoints <- function(csvFile){
  #-------------------------------------------------
  #  Read NexRad Pixels coordinates into
  #  data frame
  #-------------------------------------------------
  NRD <- utils::read.csv(csvFile)
 # NRDsub <- NRD[!NRD[,length(NRD)]==0,]
  NRDsub <- NRD  # Previous code resulted in striped RECH, RUNOFF and RAIN
  head(NRDsub)
  names(NRDsub)<- c("Pixel_id", "X", "Y" )
	NRDsub$PIX_INDX<-seq.int(nrow(NRDsub))
  #PixelCoords <- NRDsub[c("Pixel_id", "X", "Y")]
  PixelCoords <- NRDsub[c("Pixel_id", "X", "Y","PIX_INDX")]
  return(PixelCoords)
}

readgridPoints<- function(Modelgrd.Path,Model.Shape){
  #-------------------------------------------------
  #  -Read arcGIS shape file of model mesh
  #  -Convert mesh polygons to points 
  #     referencing the center of the model cells
  #  -Change spatial reference if necessary
  #  -Convert point coordinates into a data frame
  #-------------------------------------------------
  #Modelgrd %<-% rgdal::readOGR(Modelgrd.Path,Model.Shape)
  #modCol <- grep("Column_",colnames(Modelgrd@data),ignore.case=T)
  #modRow <- grep("^Row$",colnames(Modelgrd@data),ignore.case=T)
  Modelgrd <-st_read(paste0(Modelgrd.Path,'/',Model.Shape,'.shp'))
  modCol <- grep("Column_",colnames(Modelgrd),ignore.case=T)
  modRow <- grep("^Row$",colnames(Modelgrd),ignore.case=T)
  gridCentroids <- st_centroid(Modelgrd,byid=TRUE)
  print(st_crs(Modelgrd))
  if (!raster::compareCRS(HARNSP17ft,st_crs(Modelgrd))) {
    gridCentroids <-st_transform(gridCentroids,HARNSP17ft)
    # gridCentroids <- sp::spTransform(gridCentroids,HARNSP17ft)
  }
  
  # Convert model cell points to a dataframe
  asSFGC <- sf::st_as_sf(gridCentroids)
  ModelGridCoords <- do.call(base::rbind,sf::st_geometry(asSFGC))
  ModelGridCoords <-base::cbind(st_drop_geometry(Modelgrd[,c(modRow,modCol)]),
                                st_drop_geometry(ModelGridCoords))
  
  return(ModelGridCoords)
}

#-------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#-------------------------------------------------
suppressWarnings(HARNSP17ft  <- sp::CRS("+init=epsg:2881"))
suppressWarnings(HARNUTM17Nm  <- sp::CRS("+init=epsg:3747"))
suppressWarnings(latlongs <- sp::CRS("+proj=longlat +datum=WGS84"))

#-------------------------------------------------
# Provide default basepath with option to change
#-------------------------------------------------
#basePath <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS/"
basePath <-"//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/ECSM/Data"
basePath <-  tk_choose.dir(default = basePath, 
                           caption = "Select input directory for biasNRD files")

#-------------------------------------------------
#  Prepare environment to support multiprocessing
#-------------------------------------------------
cat(paste('Establishing mutliprocessor Plan','\n'))
future::plan(multisession,.skip=TRUE)

#-------------------------------------------------
#  Select input and Read one year for Pixels
#  using a background process
#-------------------------------------------------

#csvFile <- paste0(basePath, paste0("biasNRD2003.csv"))
csvFile <- utils::choose.files(default=paste(basePath,'*.csv',sep='/'))
cat(paste('"+" indicates Progress Reading',csvFile,'\n'))
#f <- future({readPoints(csvFile)})
f<- readPoints(csvFile )
#-------------------------------------------------
# Read Model grid shapefile as polys 
# and convert to points and convert to data frame
# using a background process  
# - readgridPoints() runs async to readPoints()
#-------------------------------------------------

if (model =='ECFTX'){
  Modelgrd.Path <- 
    "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/CFWI/Data/From_SW_SJ"
  Model.Shape <-"ECFTX_GRID_V3"
} else if (model =='ECSM'){
  Modelgrd.Path <- 
    "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DataLib/ModelData/ECSM"
  Model.Shape <-"ECSM_Model_Grid"
} else {
  Modelgrd.Path <- 
    "//ad.sfwmd.gov/dfsroot/data/wsd/PLN/Felipe/NEXRAD/Weekly_Exe/LWC/LWC_NRD_Data/Model_shapefiles"
  Model.Shape <-"LWCSIM_Model_Grid"
  
}
  setwd(Modelgrd.Path)

cat(paste('":" indicates Progress Reading Model shapefile',Model.Shape,'\n'))
if(!exists("ModelGridCoords")) {g <- future({readgridPoints(Modelgrd.Path,Model.Shape)})}
#g <- readgridPoints(Modelgrd.Path,Model.Shape)
#-------------------------------------------------
# Wait for values from futures
#-------------------------------------------------
cat(paste('Waiting for background processing to complete','\n'))

while (!resolved(g)){
  if (!resolved(f)){
  cat("+")
  }
  cat(":")
}
cat("\n")

#PixelCoords <-future::value(f)
gPrime <-future::value(g) 
ModelGridCoords <-future::value(g)
PixelCoords <-f
# ModelGridCoords <-g

cat(paste('Background processing is complete','\n'))

#-------------------------------------------------
# Find closest points using nearest neighbor
#-------------------------------------------------
NearNeighbor <- FNN::get.knnx(PixelCoords[,2:3],ModelGridCoords[,3:4],1)
ClosestPixels <- NearNeighbor[["nn.index"]]
ClosestDistance <- NearNeighbor[["nn.dist"]]

CpixDist <- as.data.frame(cbind(as.character(PixelCoords[ClosestPixels,1]),ClosestDistance))
CpixDist <- cbind(CpixDist,PixelCoords[ClosestPixels,4])
#CpixDist <- as.data.frame(cbind(as.character(PixelCoords[ClosestPixels,c(1,4)]),ClosestDistance))
ModelGridCoords <- cbind(CpixDist,ModelGridCoords)

#-------------------------------------------------
# Output data
#-------------------------------------------------
ModelGCdf <- as.data.frame(ModelGridCoords)
names(ModelGCdf) <- c("PixelID","distance","PIX_INDX","row","col","X","Y")
#names(ModelGCdf) <- c("PIX_INDX","PixelID","row","col")
ModelGCdf$row_col <-paste0(ModelGCdf$row,'_',ModelGCdf$col)
head(ModelGCdf)
ModelGCdf<-ModelGCdf[,c(3,1,4,5,8)]
head(ModelGCdf)
ModelGCdf$row<-as.numeric(ModelGCdf$row)
ModelGCdf$col<-as.numeric(ModelGCdf$col)

#ModelCellfile <- paste(Modelgrd.Path,'/ModelCells',,'.csv',sep='')
ModelCellfile <- gsub('.csv','NN.csv',csvFile)
utils::write.csv(as.data.frame(ModelGCdf),ModelCellfile,row.names =FALSE)
cat(paste('Crosswalk/Lookup table for NexRad Pixel to Model ROw and Column','\n'))
cat(paste('saved to: ',ModelCellfile,'\n'))

gc()

#-------------------------------------------------
# FNN::get.knnx simple example
#
# x1 is like Pixels
# x2 is like modelgrid
#-------------------------------------------------
#   x1 = cbind(runif(5),runif(5))
#   x2 = cbind(runif(10),runif(10))
#   nn = get.knnx(x1,x2,1)
#-------------------------------------------------
# returns index and distance of x1 closest to x2
#-------------------------------------------------
#   nn[["nn.index"]]
#   nn[["nn.dist"]]
