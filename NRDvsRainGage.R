
#==================================================================================
# \\ad.sfwmd.gov\dfsroot\data\wsd\sup\devel\source\R\ECSM_rain\NRDvsRainGage.R
#==================================================================================
# Modified for ECSM:  Kevin A. Rodberg -  April 2020
#
# Modified for LWC:   LWCPrepDataUsingMonthlyBiasMP.R
# Programmer:         Felipe Amorano - 10/02/2018 
#
# Original code:      PrepDataUsingBiasMP.R
# Programmer:         Kevin A. Rodberg  - 05/11/2018
#==================================================================================
# Execution with Multiprocessors significantly reduces execution time
#==================================================================================
# -Calculates monthly bias multipliers from RainGage vs NexRad Pixels.
# -Bias multipliers at RainGages are interpolated using Ordinary Kriging
#     producing an annual 'bias' raster.  
# -Daily NexRad Pixels are rasterized and multiplied by the 'bias' raster 
#     producing 'bias adjusted' Daily NexRad rasters
# -'bias adjusted' Daily NexRad rasters are converted back to pixels 
#     by extracting raster values from the rasters at pixel pointslocations.  
# -Daily pixels values by row are combined as columns and exported to csv 
#     with an Monthly total column added.
# Finally Raster plots of Monthly data are also produced showing:
#     Interpolated Bias, Uncorrected NexRad, Bias Corrected NexRad, 
#       Difference in Original vs Corrected
#==================================================================================

#==================================================================================
# Package Management
#
#   pkgChecker function provides automated means for first time use of script to  
#	  automatically install any new packages required for this code, with " calls 
#	  wrapped in a for loop.
#==================================================================================

pkgChecker <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}

list.of.packages <-  c("reshape2","readr","dplyr","tidyr", "data.table","readxl","rgeos","sp",
                       "dismo","lattice","rasterVis","maptools","raster","fields","automap",
                       "gstat","future","listenv","ggplot2","RANN","geosphere")

suppressWarnings(pkgChecker(list.of.packages))

#---------------------------------------------------
#  "Not In" Function
#---------------------------------------------------
'%!in%' <- function(x,y)!('%in%'(x,y))

#---------------------------------------------------
#  Set number of decimal places for columns in a dataframe
#---------------------------------------------------
fixDecimals <- function(DF,decPlaces){
  is.num <-sapply(DF,is.numeric)
  DF[is.num] <- lapply(DF[is.num], round,decPlaces)
  return(DF)
}

#---------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#---------------------------------------------------
HARNSP17ft  = CRS("+init=epsg:2881")
latlongs = CRS("+proj=longlat +datum=WGS84")

#------------------------------------------------------------
# Set up county boundry shapefile for overlay on raster maps
#------------------------------------------------------------
gClip <- function(shp, bb) {
  if (class(bb) == "matrix")
    b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else
    b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

#------------------------------------------------------------
# Used for clipping results to model bounary
#------------------------------------------------------------
ECSM.Path <- "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/ECSM"
ECSM.Shape <- "ECSM_bnd.shp"
setwd(ECSM.Path)
ECSM <- readShapePoly(ECSM.Shape, proj4string = HARNSP17ft)

WMDbnd.Path <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)

#------------------------------------------------------------
# SQL used to query data to csv: ECSM_RainGage_1985to2019.sql
#------------------------------------------------------------
rainGageData<-read.csv("//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/ECSM/Data/ECSM_RainGageV2.csv",stringsAsFactors = FALSE)
rainGageData$YEAR = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y")
rainGageData$MONTH = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%m")
rainGageData$DAILY_DATE = as.Date(format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y-%m-%d"))

rainGages<- unique(rainGageData[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
rownames(rainGages) <- NULL

#------------------------------------------------------------
# Function  processes one year of data at a time to
# provide async processing of function via "future" package
#------------------------------------------------------------
calcAnnualRainStats <- function(yr,i) {
  cat(yr,sep='/n')
  #------------------------------------------------------------
  # Monthly NRD data files created with:  MonthlyNRDS.R
  #   Daily produced using NRD_ECSMviaDBHydro.R 
  #        (ECSM specific edits required for path, years, 
  #           &  data type (NRD vs ECSM))
  #     which calls queryTotalNRDbyParams.R
  #------------------------------------------------------------
  oneYr<-read.csv(paste0('G:/ECSM/Data/NRDrainMonthly',yr,'.csv'))
  closest<- nn2(oneYr[,3:4],rainGages[,3:4],1)
  index= closest[[1]]
  distance=closest[[2]]
  nearestPixel<-as.data.frame(cbind(rainGages$DBKEY,oneYr[index,]$Pixel_id,  distance))
  names(nearestPixel)<-c('DBKEY', 'Pixel_id', 'Distance')
  nearestPixel$Pixel_id <- as.numeric(as.character(nearestPixel$Pixel_id))
  nearestPixel$DBKEY <- as.character(nearestPixel$DBKEY)
  nearestPixel$Distance <- as.numeric(as.character(nearestPixel$Distance))
  
  ECSM_NRD<- read_csv(paste0("G:/ECSM/Data/ECSM_NRD_",yr,".csv"))
  
  meltNRD<-melt(ECSM_NRD,id=c("ROWnum","Pixel_id","X","Y"))
  meltNRD$daily_date = as.character(meltNRD$variable)
  meltNRD$variable = as.character(meltNRD$variable)
  LongWdates<-cbind(meltNRD[meltNRD$variable < "A",] %>% separate(variable, c("year", "month", "day")))
  LongWdates$daily_date <- as.Date(LongWdates$daily_date, format="%Y-%m-%d")
  #
  # Subset Gage data for the year and filter out specific codes
  #
  RainG<- merge(x=rainGageData[rainGageData$YEAR==yr & rainGageData$CODE %!in% c('X','M','N','PT'),],
                y=nearestPixel, by.x= 'DBKEY', by.y='DBKEY')
  
  Rain<-merge(x=RainG,  y=LongWdates[,c('Pixel_id','daily_date','value')],
              by.x=c('Pixel_id','DAILY_DATE'),by.y=c('Pixel_id','daily_date'))
  
  names(Rain)<-c("Pixel_id","DAILY_DATE", "DBKEY","STATION","AGENCY","XCOORD","YCOORD","Gage","CODE",
                 "YEAR","MONTH","Distance","NexRad" )
  
  Rain$STATION <- as.factor(Rain$STATION)
  
  #-------------------------------------------------
  #  Remove extreme >20 inches and
  #  those really small values that are inconsistent
  #-------------------------------------------------
  Rain<-filter(Rain, Gage >= 0.001 & Gage < 20 & 
                        NexRad >= 0.00 & NexRad < 20 |
                 round(Gage,2) == round(NexRad,2)) 
  #-------------------------------------------------
  # Rain<-filter(Rain, (round(Gage,2) >= 0.01 & Gage < 20 & 
  #                       round(NexRad,2) >= 0.01 & NexRad < 20) )
  
  Rain <- na.omit(Rain)
  
  #-------------------------------------------------
  # Calculate Stats for Gage and NexRad bias
  # Filtered by a number of observations and higher correlation
  #-------------------------------------------------
  Rain.monthly.sum<- Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize_at(vars(Gage, NexRad), sum) %>% 
    rename(Gage.sum = Gage, NexRad.sum = NexRad)%>%
    subset(NexRad.sum > 0)
  
  Rain.monthly.obs<- Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize_at(vars(Gage), ~sum(. >-.001)) %>% 
    rename(Gage.obs = Gage) %>%
    subset(Gage.obs >10)
  
  Corr_GageNexRad<-Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize(correlation = cor(Gage, NexRad, method = "kendall")) %>%
    drop_na() %>% 
    subset(correlation >= .75 & correlation <= 1.0)
  
  RainStats<-  merge(rainGages,merge(Corr_GageNexRad,
                                     merge(x=Rain.monthly.sum,y=Rain.monthly.obs),all.y=TRUE),all.x =TRUE)
  RainStats[is.na(RainStats$correlation),]$correlation = 0.0
  RainStats[is.na(RainStats$Gage.sum),]$Gage.sum = 0.0
  RainStats[is.na(RainStats$NexRad.sum),]$NexRad.sum = 0.0
  RainStats <- na.omit(RainStats)
  
  write.csv(RainStats,paste0('h:/rainstats',yr,'.csv'))
  

  RainStats$bias = RainStats$Gage.sum/RainStats$NexRad.sum
  if( sum(round(RainStats$Gage.sum,2) == round(RainStats$NexRad.sum,2)) >0 ){
    RainStats[round(RainStats$Gage.sum,2)==round(RainStats$NexRad.sum,2),]$bias=1.0
  }
  RainStats<-fixDecimals(RainStats,3)
  
  #-------------------------------------------------
  # Prevent over-correction by limiting 
  # range for bias adjustment factor
  #-------------------------------------------------
  if( sum(is.nan(RainStats$bias)) >0 ){
    RainStats[is.nan(RainStats$bias),]$bias <-1
  }
  if( sum(is.na(RainStats$bias)) >0 ){
    RainStats[is.na(RainStats$bias),]$bias <-1
  }
  if (dim(RainStats[RainStats$bias<.5,])[1] >0){
    RainStats[RainStats$bias<.5,]$bias <- .5
  }
  if (dim(RainStats[RainStats$bias>1.5,])[1] >0){
    RainStats[RainStats$bias>1.5,]$bias <- 1.5
  }
  RainStats$YEAR<-yr
  RainStats<-merge(nearestPixel, merge(rainGages,RainStats))
  
  Rain<-merge(Rain,RainStats[,c("DBKEY","Pixel_id","STATION","MONTH","bias")])
  Rain$NRD<-Rain$NexRad*Rain$bias
  rStats<- RainStats[,c("STATION","YEAR","MONTH","XCOORD","YCOORD","Gage.obs","Gage.sum","NexRad.sum","bias")]
  
  ToPlot<- Rain[, !names(Rain) %in% c('STATION',"AGENCY","XCOORD","YCOORD","CODE","YEAR","MONTH","Distance","bias")]
  testPlot=melt(ToPlot,id=c('Pixel_id','DBKEY','DAILY_DATE','Gage'))
  
  # Swap Pixel_id and Station for plot grouping changes
    # ToPlot<- Rain[, !names(Rain) %in% c('Pixel_id',"AGENCY","XCOORD","YCOORD","CODE","YEAR","MONTH","Distance","bias")]
    # testPlot=melt(ToPlot,id=c('STATION','DBKEY','DAILY_DATE','Gage'))
  #--
  # Plot Rain Gage vs NexRad correlation
  #--
  # for (stn in unique(testPlot$Pixel_id)){
  #   cat (stn,sep='\n')
  # 
  #   filename = paste("G:/ECSM/Data/graphNRDcorrection/",stn,"_",yr,".png", sep = "" )
  #   png(  file = filename, width = 3000,height = 3000,units = "px",  res=300)
  # 
  #   p <- ggplot(testPlot[testPlot$Pixel_id==stn,],aes(x=value, y=Gage,
  #                                            color = paste(DBKEY,variable,sep='_'),shape=DBKEY)) +
  #     labs(title =stn, color = 'NexRad') +
  #     geom_point() +geom_smooth(method = "lm")
  #   print(p)
  #   dev.off()
  # }
  return(rStats)
}

plan(multiprocess)
processed= listenv(NULL)

i=0

for (yr in seq(2000,2018)){
    i = i + 1
    processed[[i]] <- future({calcAnnualRainStats(yr,i)})
}
x=i
my_data <- list()

#-------------------------------------------------
# value function waits for results to become available
# for each process
#-------------------------------------------------
for (i in seq(1:x)){
  my_data[[i]] <-value(processed[[i]])
}

allStats<-do.call(rbind,my_data)
names(allStats)<- c('RainGage','year','month','XCOORD','YCOORD','count','sum_Rainfall','sum_NRD','adjF')

basePath <-  "G:/ECSM/Data/"
ECSM_NRD<- read_csv(paste0("G:/ECSM/Data/ECSM_NRD_",yr,".csv"))
meltNRD<- na.omit(melt(ECSM_NRD,id=c("ROWnum","Pixel_id","X","Y")))
PixelCoords <- meltNRD[c("Pixel_id", "X", "Y")]

#-------------------------------------------------
# Add additional melted NRD data sets 
# to calculate bias based on more than one year
#-------------------------------------------------
NRD <- do.call("rbind", list(meltNRD))

#-------------------------------------------------
# Calculate raster extents
#-------------------------------------------------
xmin = floor(min(PixelCoords[c('X')]))
xmax = ceiling(max(PixelCoords[c('X')]))
ymin = floor(min(PixelCoords[c('Y')]))
ymax = ceiling(max(PixelCoords[c('Y')]))

#-------------------------------------------------
# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
#-------------------------------------------------
rasRows <- (ymax - ymin) / 6561.679
rasCols <- (xmax - xmin) / 6561.679

#-------------------------------------------------
# define raster and map extents using NRD pixel data extents
#-------------------------------------------------
ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,
              ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(ECSM, rasExt)

#-------------------------------------------------
# FUNCTION: dayBiasFn
# Multiplies Daily NRD rasters by bias  
#   calculating adjusted daily NexRad Raster
# and returns 
#-------------------------------------------------
dayBiasFn <- function(DailyNRD,biasRas){
  DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("X", "Y")],
                                         data = DailyNRD,proj4string = HARNSP17ft)
  NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean) * biasRas
  NRDBiasPnts <- raster::extract(NRDras,DailyNRD.pnts,fun=mean,df=TRUE)
  NRDBiasPnts$Pixel_id <- DailyNRD.pnts$Pixel_id
  NRDBiasPnts$X<-DailyNRD.pnts$X
  NRDBiasPnts$Y<-DailyNRD.pnts$Y
  return(NRDBiasPnts[c(3,4,5,2)])
}

#-------------------------------------------------
# FUNCTION: biasByYearMon
#   [Works well with future function for multiprocessing]
# Processes RainVsGage data by year
# Creating CSV files with NRDstat, biasNRD and 
# updates NRD with annual totals
#-------------------------------------------------
biasByYearMon <-function(yearStr,monStr,x){
  ECSM_NRDbyYr<- read_csv(paste0(basePath,"ECSM_NRD_",yearStr,".csv"))
  
  dateList<-names(ECSM_NRDbyYr)[-c(1,2,3,4,length(names(ECSM_NRDbyYr)))]
  NRDbyYr <- melt(ECSM_NRDbyYr, id = c("ROWnum","Pixel_id", "X", "Y"))
  NRDbyYr$ROWnum<-NULL
  
  #-------------------------------------------------
  # read and organize daily NexRad data
  #-------------------------------------------------
  for (d in dateList[1]) {
    DailyNRD <- NRDbyYr[NRDbyYr$variable == d, ]
    DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("X", "Y")],
                                           data = DailyNRD, proj4string = HARNSP17ft)
    
    NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean) 
    #    NRDBiasPnts <- data.frame(extract(NRDras, DailyNRD.pnts))
    DailyNRD.grid <- as(NRDras, "SpatialGridDataFrame")
  }
  NRDPixel<- data.frame(Pixel_id=DailyNRD.pnts$Pixel_id, 
                        X=DailyNRD.pnts$X, Y=DailyNRD.pnts$Y)
  
  mon <- as.numeric(monStr)
  RGdata <- na.omit(allStats[allStats$year == yearStr & 
                       allStats$month == monStr
                             & allStats$adjF <400,])
  
  #   If has no monthly RGdata values ... Setting some dummy values
  
  if (nrow(RGdata)==0) {
    cat (paste('Trying to skip',yearStr, monStr,' Due to lack of data\n'))
    RGdata <- na.omit(allStats[allStats$year == yearStr & 
                                 allStats$month =='01'
                               & allStats$adjF <400,])
    RGdata$month = '04'; RGdata$sum_NRD = 0.0; RGdata$sum_Rainfall = 0.0; RGdata$adjF=1.0
    rainGage.pnts <- SpatialPointsDataFrame(coords = RGdata[, c("XCOORD", "YCOORD")],
                                            data = RGdata,proj4string = HARNSP17ft)
    # biasStuff <-list("B_ras"=NRDras,
    #                  "R_pnts"=DailyNRD.pnts,
    #                  "C_pnts"=DailyNRD.pnts,
    #                  "MonthlyRas"=NRDras,
    #                  "MonNRDRas"=NRDras,
    #                  "year"=as.numeric(yearStr),
    #                  "month"=as.numeric(monStr),
    #                  "monthlyNRDbias"=as.data.frame(NRDPixel[,c(1,2,3,ncol(NRDPixel))]),
    #                  "dailyNRDbias" = as.data.frame(NRDPixel[,-c(ncol(NRDPixel))]))
    # biasStuff <-list("B_ras"=biasRas,
    #                  "R_pnts"=rainGage.pnts,
    #                  "C_pnts"=rainClustData,
    #                  "MonthlyRas"=MonRas,
    #                  "MonNRDRas"=NRDMonRas,
    #                  "year"=as.numeric(yearStr),
    #                  "month"=as.numeric(monStr),
    #                  "monthlyNRDbias"=as.data.frame(NRDbiasPixels[,c(1,2,3,ncol(NRDbiasPixels))]),
    #                  "dailyNRDbias" = as.data.frame(NRDbiasPixels[,-c(ncol(NRDbiasPixels))])
    #)
  }
    
    #-------------------------------------------------
    # Interpolate Bias from RainGages to NexRad pixels
    # Make NRD data correction using Bias
    #-------------------------------------------------
    rainGage.pnts <- SpatialPointsDataFrame(coords = RGdata[, c("XCOORD", "YCOORD")],
                                            data = RGdata,proj4string = HARNSP17ft)
    latlongPnts <- spTransform(rainGage.pnts,latlongs)
    distMatrix <-distm(latlongPnts)
    hc <- hclust(as.dist(distMatrix), method="complete")
    
    # 8329.369=distHaversine(c(-80.94646,25.2545),c(-80.88909,25.2006))
    # 5000 is about 5 kilometer
    #latlongPnts$clust <-cutree(hc,h=2000)
    latlongPnts$clust <-cutree(hc,h=5000)
    tempdf <- as.data.frame(latlongPnts)
    
    biasVals2 <- aggregate(cbind(XCOORD, YCOORD, adjF)~clust, tempdf, mean,na.rm=TRUE) 
    gageClustVals <- aggregate(cbind(XCOORD, YCOORD, sum_Rainfall)~clust, tempdf, mean,na.rm=TRUE) 
    
    biasVals<-biasVals2[, c("XCOORD", "YCOORD", "adjF")]
    rainVals<-gageClustVals[, c("XCOORD", "YCOORD", "sum_Rainfall")]
    
    biasVals <- biasVals[biasVals$XCOORD> NRDras@extent[1] 
                         & biasVals$XCOORD< NRDras@extent[2]
                         & biasVals$YCOORD> NRDras@extent[3] 
                         & biasVals$YCOORD< NRDras@extent[4],]
    
    rainVals <- rainVals[rainVals$XCOORD> NRDras@extent[1] 
                         & rainVals$XCOORD< NRDras@extent[2]
                         & rainVals$YCOORD> NRDras@extent[3] 
                         & rainVals$YCOORD< NRDras@extent[4],]
    
    biasData<-biasVals
    coordinates(biasData) =  ~ XCOORD + YCOORD
    proj4string(biasData) = HARNSP17ft
    biasData$XHARN <- coordinates(biasData)[, 1]
    biasData$YHARN <- coordinates(biasData)[, 2]
    biasData <- spTransform(biasData,HARNSP17ft)
    
    rainClustData<-rainVals
    coordinates(rainClustData) =  ~ XCOORD + YCOORD
    proj4string(rainClustData) = HARNSP17ft
    rainClustData$XHARN <- coordinates(rainClustData)[, 1]
    rainClustData$YHARN <- coordinates(rainClustData)[, 2]
    rainClustData <- spTransform(rainClustData,HARNSP17ft)
    
    #NRDras <- rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean)
    
    #-------------------------------------------------
    #  Theisen Polygon and raster code:
    # theisPoly <- voronoi(rainGage.pnts)
    # TheisRas <- rasterize(theisPoly, NRDras, theisPoly$bias, fun = mean)
    
    #-------------------------------------------------
    #  autoKrige implemented for Ordinary kriging
    #-------------------------------------------------
    if (min(biasData@data$adjF)==max(biasData@data$adjF)){
    #### create empty or single value raster) ####
    #  IDW raster code:
      gs <- gstat(formula=adjF~1, locations=biasData)
      biasRas <- interpolate(ras, gs)
    } else {
      surf <- autoKrige(formula=adjF ~ 1, input_data=biasData, new_data = DailyNRD.grid)
      biasRas <- raster(surf$krige_output)
    }
    # surf <- autoKrige(formula=adjF ~ 1, input_data=biasData, new_data = DailyNRD.grid)
    # biasRas <- raster(surf$krige_output)
    rainGage.pnts$bias <- raster::extract(biasRas,rainGage.pnts,fun=mean,df=TRUE)[,2]
    rainClustData$bias <- raster::extract(biasRas,rainClustData,fun=mean,df=TRUE)[,2]

    
    NRDbiasPixels <-NRDPixel
    names(NRDbiasPixels)<-c( "Pixel_id" ,"X" ,"Y")
    #-------------------------------------------------
    #  Process NexRad for each day calling "dayBiasFn"
    #  NRD with bias results are merged into a single table 
    #-------------------------------------------------
    monthFilter =paste0(yearStr,'-',monStr)
    monList <- dateList[startsWith(dateList,monthFilter)]
    
    for (d in monList) {
      table = dayBiasFn(NRDbyYr[NRDbyYr$variable == d, ],biasRas)
      names(table) <- c( "Pixel_id" ,"X" ,"Y",d  )
      NRDbiasPixels <- merge(NRDbiasPixels, table, by =c("Pixel_id","X","Y"))
    }
    
    #-------------------------------------------------
    # Add final column for monthly total NRD with Bias correction
    # and export to csv
    #-------------------------------------------------
    
    NRDbiasPixels$Monthly<-rowSums(NRDbiasPixels[,-c(1,2,3)],na.rm=TRUE)
    names(NRDbiasPixels)[length(names(NRDbiasPixels))]<- sprintf('Mon%02d',mon)
    NRDbiasPixels <- fixDecimals(NRDbiasPixels,4)
    #csvFile <- paste0(basePath, sprintf("biasNRD%s%02d.csv",yearStr,mon))
    #fwrite(NRDbiasPixels, csvFile) 
    
    #-------------------------------------------------
    # Create Raster for bias corrected Monthly NexRAD rain
    #-------------------------------------------------
    Monthly <-NRDbiasPixels[,c(1,2,3,ncol(NRDbiasPixels))]
    xy <- Monthly[,c(2,3)]
    Monthlypdf <-SpatialPointsDataFrame(coords=xy,data=Monthly,proj4string=HARNSP17ft)
    MonRas <- rasterize(Monthlypdf,ras,Monthlypdf@data[,4],fun=mean)
    rainClustData$NRD <- raster::extract(MonRas,rainClustData,fun=mean,df=TRUE)[,2]
    
    #-------------------------------------------------
    # Add final column for monthly total NRD
    # and export to csv
    #-------------------------------------------------
    col.num <- which(colnames(ECSM_NRDbyYr) %in% monList)
    ECSM_NRDbyYr$Monthly<-rowSums(ECSM_NRDbyYr[,col.num],na.rm=TRUE)
    #csvFile <- paste0(basePath, sprintf("NRD%s%02d.csv",yearStr,mon))
    ECSM_NRDbyYr <- fixDecimals(ECSM_NRDbyYr,4)
    #fwrite(NRDwCoords[,c(1,6,unlist(col.num),ncol(NRDwCoords))], csvFile) 
    
    #-------------------------------------------------
    # Create Raster for Monthly NEXRad rain
    #-------------------------------------------------
    MonNRDs<-ECSM_NRDbyYr[,c(2,3,4,ncol(ECSM_NRDbyYr))]
    xy <- MonNRDs[,c(2,3)]
    MonNRDsspdf <-SpatialPointsDataFrame(coords=xy,data=MonNRDs,proj4string=HARNSP17ft)
    NRDMonRas <- rasterize(MonNRDsspdf,ras,MonNRDsspdf$Monthly,fun=mean)
    rainClustData$Nexrad <- raster::extract(NRDMonRas,rainClustData,fun=mean,df=TRUE)[,2]
    
    
    #-------------------------------------------------
    #  Return a list of 4 spatial objects 
    #     (3 rasters & 1 points)
    # MonRas is corrected NRDMonRas is uncorrected
    #     and some data frames
    #   rainClustData are also points
    #-------------------------------------------------
    biasStuff <-list("B_ras"=biasRas,
                     "R_pnts"=rainGage.pnts,
                     "C_pnts"=rainClustData,
                     "MonthlyRas"=MonRas,
                     "MonNRDRas"=NRDMonRas,
                     "year"=as.numeric(yearStr),
                     "month"=as.numeric(monStr),
                     "monthlyNRDbias"=as.data.frame(NRDbiasPixels[,c(1,2,3,ncol(NRDbiasPixels))]),
                     "dailyNRDbias" = as.data.frame(NRDbiasPixels[,-c(ncol(NRDbiasPixels))])
                     )
  

  return(biasStuff)
}

#-------------------------------------------------
# Set up for multiprocessing function calls
#-------------------------------------------------
processed= listenv(NULL)
yrList=list()

yearStr <- as.character(2014)
#-------------------------------------------------
# Define range of years to process
#-------------------------------------------------
processYears <- seq(2000, 2018)
x=0
yr = processYears[1]
for (yr in processYears) {
  for (mon in seq(1,12)){
    yearStr <- as.character(yr)
    monStr <-sprintf("%02d",mon)
    x=x+1
    cat (paste(yearStr,':',monStr,"\n"))
    #-------------------------------------------------
    # Call FUNCTION "biasByYear" with futures multiprocessing
    # wrapper function
    #-------------------------------------------------
    processed[[x]] <- future({biasByYearMon(yearStr,monStr,x)})
    #processed[[x]] <- biasByYearMon(yearStr,monStr,x)
  }
}

mpList<-list()
rList<-list()
stackList <-list()

#-------------------------------------------------
# value function waits for results to become available
# for each process
#-------------------------------------------------
for (i in seq(1:x)){
  mpList <-value(processed[[i]])
  rList[[i]]<-mpList
}
#-------------------------------------------------
#  Uncomment for not parallel processing
#-------------------------------------------------
#for (i in seq(1:x)){
#  mpList <-processed[[i]]
#  rList[[i]]<-mpList
#}

#-------------------------------------------------
# unlist results returned from FUNCTION "biasByYearMon"
#-------------------------------------------------
stackList = unlist(lapply(rList,"[[",1))
pointList = unlist(lapply(rList,"[[",2))
ClustPntList = unlist(lapply(rList,"[[",3))
AnnRainList = unlist(lapply(rList,"[[",4))
AnnNRDList = unlist(lapply(rList,"[[",5))
yrList = unlist(lapply(rList,"[[",6))
MonthList = unlist(lapply(rList,"[[",7))
filePath <- basePath
AnnualTotals <- list()
iy = 0
for (iyr in processYears){
  iy = iy + 1
  NRDbiasMonthly <-rList[[1]]$monthlyNRDbias[,c(1,2,3)]
  NRDbiasDaily <- rList[[1]]$dailyNRDbias[,c(1,2,3)]
  for (i in seq(from=1, to=x)){
    if (rList[[i]]$year == iyr){
      NRDbiasMonthly <- merge(NRDbiasMonthly, rList[[i]]$monthlyNRDbias[,-c(2,3)], by =c("Pixel_id"))
      NRDbiasDaily <- merge(NRDbiasDaily, rList[[i]]$dailyNRDbias[,-c(2,3)], by =c("Pixel_id"))
    }
  }
  NRDbiasMonthly$Annual<-rowSums(NRDbiasMonthly[,-c(1,2,3)],na.rm=TRUE)
  cat(paste("exporting",sprintf("MonthlybiasNRD%04d.csv",iyr),'\n'))
  csvFile <- paste0(filePath, sprintf("MonthlybiasNRD%04d.csv",iyr))
  cat(csvFile)
  cat('\n')
  fwrite(NRDbiasMonthly, csvFile) 
  
  NRDbiasDaily$Annual<-rowSums(NRDbiasDaily[,-c(1,2,3)],na.rm=TRUE)
  cat(paste("exporting",sprintf("DailybiasNRD%04d.csv",iyr),'\n'))
  csvFile <- paste0(filePath, sprintf("DailybiasNRD%04d.csv",iyr))
  cat(csvFile)
  cat('\n')
  fwrite(NRDbiasDaily, csvFile) 
  AnnualTotals[[iy]]<- NRDbiasDaily$Annual
}
Annuals<-as.data.frame(t(do.call("rbind",AnnualTotals)))
names(Annuals)<- processYears
Annuals<-cbind(NRDbiasMonthly[,1:3],Annuals)
csvFile <- paste0(filePath, "ECSM_biasNRD_Totals.csv")
fwrite(Annuals, csvFile) 

rasStack <-stack()
rasStack <-stack(stackList)
MonRasStack <-stack()
MonRasStack <-stack(AnnRainList)  # corrected
MonNRDStack <-stack()
MonNRDStack <-stack(AnnNRDList)   # uncorrected

yearList <-as.character(processYears)
StackNames = list()
x = 0
yr = processYears[1]
for (yr in processYears) {
  for (mon in seq(1, 12)) {
    x = x + 1
    StackNames[x] <- sprintf("%4d%02d",yr, mon)
  }
}
names(rasStack)<- StackNames
names(MonRasStack)<- StackNames
names(MonNRDStack)<- StackNames

#myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))

myTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
divTheme = rasterTheme(region = brewer.pal('Spectral', n = 11))
adjTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
difTheme = rasterTheme(region = brewer.pal('Spectral', n = 11))

#-------------------------------------------------
# Create plot files for each raster type by month
#-------------------------------------------------
plotOneRas <- function(filename, rasPlt, rasPltName, points, pltTheme, atVals,clpBnds2,WMDbnd){
  myplot=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                     at=atVals,layout=c(1,1),contour=FALSE, margin=F) +
             latticeExtra::layer(sp.polygons(clpBnds2)) +                
             latticeExtra::layer(sp.polygons(WMDbnd)) +
             latticeExtra::layer(sp.text(coordinates(points),txt=points$RainGage,pos=1,cex=.5 )) +
             latticeExtra::layer(sp.points(points, col = "red"))
  )
  return(myplot)
}

MonNRDLayer=MonNRDStack[[1]];points=pointList[[1]];ClustPoints=ClustPntList[[1]]
rasLayer=rasStack[[1]];MonrasLayer=MonRasStack[[1]]
pltStackItems<-function(MonNRDLayer,points,ClustPoints,rasLayer,MonrasLayer){
    
    # Plot Monthly NexRad
    filename=paste(basePath,"BiasPlots/points",names(MonNRDLayer),".png",sep="")
    rasPlt <- MonNRDLayer
    rasPltName <-paste('Uncorrected NexRad\n',names(MonNRDLayer))
    if( sum(is.na(ClustPoints@data$bias)) >0 ){
      ClustPoints@data[is.na(ClustPoints@data$bias),c('bias', 'NRD', 'Nexrad')]<-0
    }
    pltTheme <- adjTheme
    atVals <-c(seq(0,28,length=29),35)
    NRDplot=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                        at=atVals, xlab = NULL, margin=F,
                        layout=c(1,1),contour=FALSE) +
                latticeExtra::layer(sp.polygons(clpBnds2)) +
                latticeExtra::layer(sp.polygons(WMDbnd)) +
                latticeExtra::layer(sp.text(coordinates(ClustPoints),
                                            txt=as.character(round((ClustPoints$sum_Rainfall-
                                                                      ClustPoints$Nexrad),2)),pos=1,cex=.5 )) +
                latticeExtra::layer(sp.points(ClustPoints,pch=1, 
                                              cex=round((ClustPoints$sum_Rainfall-ClustPoints$Nexrad)*.65,2),
                                              col = "black")) +             
                latticeExtra::layer(sp.points(ClustPoints,pch=1, 
                                              cex=round(((ClustPoints$sum_Rainfall-ClustPoints$Nexrad)*-1)*.65,2) ,
                                              col = "red"))
    )
    #  trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
    # print(NRDplot)
    # dev.off()
    
    # Plot Monthly Bias Raster 
    rasPlt <- rasLayer
    rasPltName <-paste("Bias Adjustment\n",names(rasLayer))
    min(points$adjF)
    max(points$adjF)
    pltTheme <- divTheme
    #atVals <-3^((-10:10)/7.8)
    #atVals <-3^((-10:10)/10)
    atVals <-c(0,2.5^((-10:10)/10))
    Biasplot=plotOneRas(filename, rasPlt, rasPltName,points, pltTheme, atVals,clpBnds2,WMDbnd)
    # trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
    # print(Biasplot)
    # dev.off()
    
    #   Plot Monthly Bias Adjusted NexRad
    rasPlt <- MonrasLayer
    rasPltName <-paste("Bias Adjusted NexRad\n",names(MonrasLayer))
    pltTheme <- adjTheme
    atVals <-c(seq(0,28,length=29),35)
    NRDBiasplot=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                            at=atVals, xlab = NULL, margin=F,
                            layout=c(1,1),contour=FALSE) +
                    latticeExtra::layer(sp.polygons(clpBnds2)) +
                    latticeExtra::layer(sp.polygons(WMDbnd)) +
                    latticeExtra::layer(sp.text(coordinates(ClustPoints),
                                                txt=as.character(round(ClustPoints$sum_Rainfall-
                                                                         (ClustPoints$NRD),2)),
                                                pos=1,cex=.5 )) +
                    latticeExtra::layer(sp.points(ClustPoints,pch=1, 
                                                  cex=round(((ClustPoints$sum_Rainfall-
                                                                (ClustPoints$NRD))*.65),2),
                                                  col = "black")) +             
                    latticeExtra::layer(sp.points(ClustPoints,pch=1, 
                                                  cex=round((((ClustPoints$sum_Rainfall-
                                                                 (ClustPoints$NRD))*-1)*.65),2),
                                                  col = "red"))
    )
    # trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
    # print(NRDBiasplot)
    # dev.off()
    
    #   Plot Monthly Bias Adjusted NexRad - NRD
    rasPlt <-MonrasLayer-MonNRDLayer
    rasPltName <-paste("Adj NexRad - Uncorrected\n",names(MonrasLayer))
    pltTheme <- difTheme
    atVals <-c(-6,-5,-4,-3,-2,-1,-.5,0,.5,1,2,4,8,12,16,20,24)
    #atVals <-c(-6,seq(-3,3,length=24),6)
    if( sum(is.na(points@data$adjF)) >0 ){
      points@data[is.na(points@data$adjF),]$adjF = 1
    }
    diffRasPlt=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                           at=atVals, xlab = NULL, margin=F,
                           layout=c(1,1),contour=FALSE) +
                   latticeExtra::layer(sp.polygons(clpBnds2)) +                
                   latticeExtra::layer(sp.polygons(WMDbnd)) +
                   latticeExtra::layer(sp.text(coordinates(ClustPoints),
                                               txt=as.character(round(ClustPoints$sum_Rainfall-
                                                                        (ClustPoints$NRD),2)),
                                               pos=1,cex=.5 )) +                 
                   latticeExtra::layer(sp.points(ClustPoints,pch=1, 
                                                 cex=round(((ClustPoints$sum_Rainfall-
                                                               (ClustPoints$NRD))*.65),2),
                                                 col = "black")) +             
                   latticeExtra::layer(sp.points(ClustPoints,pch=1, 
                                                 cex=round((((ClustPoints$sum_Rainfall-
                                                                (ClustPoints$NRD))*-1)*.65),2),
                                                 col = "red"))
    )
    
    # trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
    # print(diffRasPlt)
    # dev.off()
    trellis.device(device="png", filename=filename, width=4800,height=2400,units="px",res=300)
    print(Biasplot, split    = c(1,1,4,1),more=TRUE)
    print(NRDplot, split     = c(2,1,4,1),more=TRUE)
    print(NRDBiasplot, split = c(3,1,4,1),more=TRUE)
    print(diffRasPlt, split  = c(4,1,4,1))
    dev.off()

}

processed= listenv(NULL)
# for (i in 1:5){
  
for (i in 1:nlayers(rasStack)){
  #pltStackItems(MonNRDStack[[i]],pointList[[i]],ClustPntList[[i]],rasStack[[i]],MonRasStack[[i]])
  # processed[[i]] <- future({pltStackItems(MonNRDStack[[i]],pointList[[i]],
  #                                         ClustPntList[[i]],rasStack[[i]],MonRasStack[[i]])})
  MonNRDLayer=MonNRDStack[[i]];points=pointList[[i]];ClustPoints=ClustPntList[[i]]
  rasLayer=rasStack[[i]];MonrasLayer=MonRasStack[[i]]
  if( sum(is.na(ClustPoints@data$bias)) >0 ){
    ClustPoints@data[is.na(ClustPoints@data$bias),c('bias', 'NRD', 'Nexrad')]<-0
  }
  processed[[i]] <- future({pltStackItems(MonNRDLayer,points,ClustPoints,rasLayer,MonrasLayer)})
  cat(paste(names(MonNRDStack[[i]]),'\n'))
}

for (x in seq(1:i)){
  done <-value(processed[[x]])
  cat('.')
}
cat('\n')
