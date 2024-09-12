#===============================================================================
#  Program: Xtra_NRDvsRainGage.R 
#           \\ad.sfwmd.gov\dfsroot\data\wsd\sup\devel\source\R\ECSM_rain\
#===============================================================================
# Code History:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Original:       PrepDataUsingBiasMP.R            Kevin A. Rodberg - 05/11/2018
# update  LWC:    LWCPrepDataUsingMonthlyBiasMP.R  Felipe Zamorano -  10/02/2018
# update  ECSM:   NRDvsRainGage.R                  Kevin A. Rodberg - April 2020
# update  ECSM:   Xtra_WMMvsRainGage.R             Kevin A. Rodberg -  July 2020
#                   Extend WMM pixels South of grid  
#                   more functions setup to use multiprocessing 
# update  ECSM:   Xtra_NRDvsRainGage.R             Kevin A. Rodberg -  July 2020
#                   Extend NexRAD pixels South of grid  
#                   setup remaining functions to use multiprocessing 
# update  ECSM:   Xtra_NRDvsRainGage2024.R         Kevin A. Rodberg -  Sept 2024
#                   update code removing deprecated packages
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# >>> Execution with Multiprocessors significantly reduces execution time
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#  General Description:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# -Calculates monthly bias multipliers from RainGage vs NexRad Pixels.
# -Bias multipliers at Rain Gages are interpolated using Ordinary Kriging
#     producing a monthly 'bias' raster.   
# -Daily NexRad Pixels are rasterized and multiplied by the 'bias' raster 
#     producing 'bias adjusted' Daily NexRad rasters
# -'bias adjusted' Daily NexRad rasters are converted back to NRD pixels 
#     by extracting raster values from the rasters at pixel points locations.  
# -Daily pixels values by row are combined as columns and exported to csv 
#     with an Monthly total column added.
# -Finally, Raster plots of Monthly data are also produced showing:
#     Interpolated Bias, Uncorrected NexRad, Bias Corrected NexRad, 
#       Difference in https://cran.r-project.org/bin/windows/Rtools/Original vs Corrected
#===============================================================================

list.of.packages <-  c("reshape2","readr","dplyr","tidyr", "data.table",
                       "readxl","terra","sp", "dismo", "lattice","rasterVis",
                       "maptools","raster","fields","automap", "gstat",
                       "future","listenv","ggplot2","RANN","geosphere",
                       "tcltk2", "RColorBrewer")

pkgChecker <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}
suppressWarnings(pkgChecker(list.of.packages))

yrSeq<-seq(2000,2018)

'%!in%' <-  function(x,y)!('%in%'(x,y))

fixDecimals <- function(DF,decPlaces){
  is.num <-sapply(DF,is.numeric)
  colNames = names(DF)
  #DF[is.num] <- lapply(DF[is.num], round,decPlaces)
  for (col in colNames[is.num]){
    DF[[col]]=round(DF[[col]],3)
  }
  return(DF)
}

#------------------------------------------------------------
# Set up county boundry shapefile for overlay on raster maps
#------------------------------------------------------------
gClip <- function(shp, bb) {
  if (class(bb) == "matrix")
    b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else
    b_poly <- as(extent(bb), "SpatialPolygons")
  #gIntersection(shp, b_poly, byid = T)
  terra::crop(shp,b_poly)
}

#---------------------------------------------------
# Define GIS variables:
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
# and WGS84 Lat/Long Spatial Reference variables
#---------------------------------------------------
HARNSP17ft  = CRS("+init=epsg:2881")
latlongs = CRS("+proj=longlat +datum=WGS84")

WMDbnd.Path <-"//ad.sfwmd.gov/dfsroot/data/hpcc_shared/krodberg/NexRadTS"

WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)
setwd("H:/")
calcRainStats<-function(basePath,yr,wideXtra,RG){
  # FUNCTION: calcRainStats
  #   [Works well with future function for multiprocessing]
  # Calculate Series of Monthly Rain Stats for 1 year from gage data
  #---------------------------------------------------
  rainGages<- unique(RG[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rownames(rainGages) <- NULL
  
  ECSM_NRD<- read.csv(paste0(basePath,"ECSM_NRD_",yr,".csv"))
  ECSM_NRD$ROWnum = NULL
  ECSM_NRD[,c(length(ECSM_NRD))] = NULL
  
  # Drop column for day Feb 29 in non-Leap years and get rid of Total column
  col.num=0
  tempWideXtra<-wideXtra
  if( yr/4 - (yr%/%4) != 0) {
    col.num <- which( colnames(tempWideXtra) %in%'X2000.02.29' )
    tempWideXtra[,col.num]<- NULL
  }
  tempWideXtra[,length(tempWideXtra)]<- NULL
  names(tempWideXtra)<-names(ECSM_NRD)
  ECSM_NRD<- rbind(ECSM_NRD,tempWideXtra)
  
  meltNRD<-reshape2::melt(ECSM_NRD,id=c("Pixel_id","X","Y"))
  meltNRD<-meltNRD[meltNRD$variable != 'Annual',]
  meltNRD[is.na(meltNRD$value),]$value = 0.0
  meltNRD <- na.omit(meltNRD)
 
  meltNRD$variable = gsub('X','',as.character(meltNRD$variable))
  meltNRD$daily_date = as.Date(gsub('\\.','-',as.character(meltNRD$variable)))
  meltNRD$variable = as.character(meltNRD$daily_date)
  
  LongWdates<-cbind(meltNRD[meltNRD$variable < "A",] %>% separate(variable, 
                     c("year", "month", "day")))
  LongWdates$daily_date <- as.Date(LongWdates$daily_date, format="%Y-%m-%d")
  LongWdates$value <- as.numeric(LongWdates$value)
  
  #
  # Subset Gage data for the year and filter out specific codes
  #
  RainG<- merge(x=RG,   y=nnRowCo, by.x= 'DBKEY', by.y='DBKEY')
  Rain<-merge(x=RainG,  y=LongWdates[,c('Pixel_id','daily_date','value')],
              by.x=c('Pixel_id','DAILY_DATE'),by.y=c('Pixel_id','daily_date'))
  names(Rain)<-c("Pixel_id","DAILY_DATE", "DBKEY","STATION","AGENCY","XCOORD",
                 "YCOORD","Gage","CODE","YEAR","MONTH","dist","NRDRain" )
  #Rain$STATION <- as.factor(Rain$STATION)
  Rain<-filter(Rain, (Gage >= 0.01 & Gage < 20 & NRDRain >= 0.01 & NRDRain < 20) | 
                 (Gage <= 0.01 & NRDRain <= 0.01))
  Rain <- na.omit(Rain)
  
  #
  # Calculate Stats for Gage and NRDRain bias
  #
  Rain.monthly.sum<- Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize_at(vars(Gage, NRDRain), sum) %>% 
    rename(Gage.sum = Gage, NRDRain.sum = NRDRain)
  
  Rain.monthly.obs<- Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize_at(vars(Gage), ~sum(. >-.001)) %>% 
    rename(Gage.obs = Gage)  %>% 
    subset(Gage.obs >= 14)  
  
  noNAsdByGage<-Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize(sdGage = sd(Gage),.groups='keep') %>%
    drop_na() %>%
    subset(sdGage > 0)
  x1 = as.data.table(noNAsdByGage)
  
  noNAsdByNexRad<-Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize(sdNexRad = sd(NRDRain),.groups='keep') %>%
    drop_na() %>%
    subset(sdNexRad > 0) 
  x2 = as.data.table(noNAsdByNexRad)
  
  Rain <- left_join(Rain, x1, by = c('STATION','DBKEY','MONTH'))
  Rain <- left_join(Rain, x2, by = c('STATION','DBKEY','MONTH'))
  Rain<- na.omit(Rain)

  Rain.monthly.corr<-Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize(correlation = cor(Gage, NRDRain, method = "kendall"),
              .groups='keep') %>%
    drop_na() %>% 
    subset(correlation >= .75 & correlation <= 1.0)
  
  RainStats<-  merge(rainGages,merge(Rain.monthly.corr,
                                     merge(x=Rain.monthly.sum,y=Rain.monthly.obs),all.y=TRUE),all.x =TRUE)
  if (nrow(RainStats[is.na(RainStats$correlation),])> 0) {
    RainStats[is.na(RainStats$correlation),]$correlation = 1.0
  }
  if (nrow(RainStats[is.na(RainStats$Gage.sum),])> 0) {
    RainStats[is.na(RainStats$Gage.sum),]$Gage.sum = 0.0
  }
  if (nrow(RainStats[is.na(RainStats$NRDRain.sum),])> 0) {
    RainStats[is.na(RainStats$NRDRain.sum),]$NRDRain.sum = 0.0
  }
  RainStats<-na.omit(RainStats)
  RainStats[RainStats$YCOORD <= SouthMost,]$NRDRain.sum<-RainStats[RainStats$YCOORD <= SouthMost,]$Gage.sum
  rainGages<- unique(RG[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rainGages<-na.omit(rainGages)
  
  #
  # Rainstats maybe useful for optimizing bias correction
  #
  write.csv(RainStats[!is.na(RainStats$MONTH),],paste0('h:/NRDrainstats',yr,'.csv'))
  
  #
  # Filter RainStats
  #
  RainStats$bias = RainStats$Gage.sum/RainStats$NRDRain.sum
  RainStats[is.nan(RainStats$bias),]$bias = 1.0
  RainStats$bias<-round(RainStats$bias,3)
  
  #
  # Prevent over-correction
  #
  if( sum(is.nan(RainStats$bias)) >0 ){
    RainStats[is.nan(RainStats$bias),]$bias <-1
  }
  if( sum(is.na(RainStats$bias)) >0 ){
    RainStats[is.na(RainStats$bias),]$bias <-1
  }
  if (dim(RainStats[RainStats$bias<.75,])[1] >0){
    RainStats[RainStats$bias<.75,]$bias <- .75
  }
  # if (dim(RainStats[RainStats$bias>1.5,])[1] >0){
  #   RainStats[RainStats$bias>3,]$bias <- 1.5
  if (dim(RainStats[RainStats$bias>4.0,])[1] >0){
    RainStats[RainStats$bias>4.0,]$bias <- 4.0
  }
  # RainStats[RainStats$Gage.obs < SignificantDays,]$bias <- 1
  RainStats$YEAR<-yr
  #   RainG<- merge(x=RG,   y=nnRowCo, by.x= 'DBKEY', by.y='DBKEY')
  x3 = merge(rainGages,RainStats)
  x4<-merge(nnRowCo, x3)
  RainStats = x4
  #RainStats<-merge(nnRowCo, merge(rainGages,RainStats))
  #Rain<-merge(Rain,RainStats[,c("DBKEY","Pixel_id","STATION","MONTH","bias")])
  x5<-merge(x=Rain,y=RainStats[,c("DBKEY","Pixel_id","STATION","MONTH","bias")],
            by.x=c("DBKEY","Pixel_id","STATION","MONTH"), 
            by.y=c("DBKEY","Pixel_id","STATION","MONTH") )
  Rain = x5
  Rain$NRD<-Rain$NRDRain*Rain$bias
  
  keepCols = c("STATION","YEAR","MONTH","XCOORD","YCOORD","Gage.obs",
               "Gage.sum","NRDRain.sum","bias")
  #--------------------------------------------------------------
  #  Return results as list 
  #  as appropriate for multiprocessing "futures" function calls
  #--------------------------------------------------------------
  rainStatStuff=list("Rain"=Rain,
                     "RainStats"=RainStats[,keepCols])
  return(rainStatStuff)
  
  #--------------------------------------------------------------
  #   The following code can be copied into function to 
  #   produce XY charts although it may need modifications 
  #   as variables may have evolved since intial coding
  #--------------------------------------------------------------
  
  # ToPlot<- Rain[, !names(Rain) %in% c('STATION',"AGENCY","XCOORD","YCOORD",
  #                                     "CODE","YEAR","MONTH","dist","bias")]
  # ToPlot<- Rain[, !names(Rain) %in% c('Pixel_id',"AGENCY","XCOORD","YCOORD",
  #                                     "CODE","YEAR","MONTH","dist","bias")]
  # testPlot=melt(ToPlot,id=c('STATION','DBKEY','DAILY_DATE','Gage'))
  # testPlot=melt(ToPlot,id=c('Pixel_id','DBKEY','DAILY_DATE','Gage'))
  
  # for (stn in unique(testPlot$Pixel_id)){
  #   cat (stn,sep='\n')
  # 
  #   fileName = paste0("G:/ECSM/Data/graphNRDcorrection/",stn,"_",yr,".png")
  #   png(  file = fileName, width = 3000,height = 3000,units = "px", res=300)
  # 
  #   p <- ggplot(testPlot[testPlot$Pixel_id==stn,],aes(x=value, y=Gage,
  #                                            color = paste(DBKEY,variable,
  #                                            sep='_'),shape=DBKEY)) +
  #     labs(title =stn, color = 'NRDRain') + 
  #     geom_point() +
  #     geom_smooth(method = "lm")
  #   print(p)
  #   dev.off()
  # }
}

dayBiasFn <- function(DailyNRD,biasRas){
  #-------------------------------------------------
  # FUNCTION: dayBiasFn
  # Multiplies Daily NRD rasters by bias  
  #   calculating adjusted daily NRDRain Raster
  # and returns 
  #-------------------------------------------------  
  DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("X", "Y")],
                    data = DailyNRD,proj4string = HARNSP17ft)
  NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, 
                     fun = mean) * biasRas
  # gs <- gstat(formula=adjF~1, locations=biasData)
  # idw <- interpolate(NRDras, gs)
  NRDBias.pnts <- raster::extract(NRDras,DailyNRD.pnts,fun=mean,df=TRUE)
  NRDBias.pnts$Pixel_id <- DailyNRD.pnts$Pixel_id
  NRDBias.pnts$X<-DailyNRD.pnts$X
  NRDBias.pnts$Y<-DailyNRD.pnts$Y
  return(NRDBias.pnts[c(3,4,5,2)])
}

biasByYearMon <-function(basePath,allStats,Rain,southRain,yearStr,monStr,x){
  yr = as.numeric(yearStr)
  #-------------------------------------------------
  # FUNCTION: biasByYearMon
  #   [Works well with future function for multiprocessing]
  # Processes RainVsGage data by year
  # Creating CSV files with NRDstat, biasNRD and 
  # updates NRD with annual totals
  #-------------------------------------------------  
  # read and organize daily NRDRain data
  #-------------------------------------------------
  NRDrainByYr<- read.csv(paste0(basePath,"ECSM_NRD_",yr,".csv"),
                         stringsAsFactors=F)
  
  NRDbyYr <- reshape2::melt(NRDrainByYr, id = c("ROWnum","Pixel_id", "X", "Y"))
  NRDbyYr$ROWnum<-NULL
  coordinates(NRDrainByYr) =  ~ X + Y
  proj4string(NRDrainByYr) = HARNSP17ft
  NRDrainByYr$XHARN <- coordinates(NRDrainByYr)[, 1]
  NRDrainByYr$YHARN <- coordinates(NRDrainByYr)[, 2]
  NRDrainByYr <- spTransform(NRDrainByYr,HARNSP17ft)
  
  dateList <-as.character(unique(NRDbyYr$variable))
  dayList <- paste0('day',seq(1,length(dateList)-1))
  dateFormList = as.Date(paste0(as.character(yr),'-01-01') ) + 
     as.numeric(gsub('day','',as.character(dayList)))-1
  dayLUp<-as.data.frame(cbind(dayList,as.character(dateFormList)),
                        stringsAsFactors=F)
  dayLUp$V2<-as.character(dayLUp$V2)
  dayLUp <- na.omit(dayLUp)
  
  for (d in dayList[1]) {
    #-------------------------------------------------
    #  Theisen Polygon and raster code to fill in southern area:
    firstDay<-unique(southRain$DAILY_DATE)[1]
    theisPoly <- 
      voronoi(southRain[southRain$DAILY_DATE==firstDay,])
    TheisRas <- rasterize(theisPoly, ras, theisPoly$VALUE, fun = mean)
    
    DailyNRD <- NRDbyYr[NRDbyYr$variable == dateList[1], ]
    DailyNRD.pnts <-SpatialPointsDataFrame(coords = DailyNRD[, c("X", "Y")],
                      data = DailyNRD, proj4string = HARNSP17ft)
    
    NRDras <-rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean) 
    # Do raster overlay Here:
    s<-stack(NRDras,TheisRas)
    mergeRas<-calc(s,fun=function(x) ifelse(is.na(x[1]), 
                                            ifelse(is.na(x[2]),0,x[2]),x[1]))
    NRDBias.pnts <- data.frame(raster::extract(mergeRas, oneYr))
    NRDBias.pnts$Pixel_id <- oneYr$Pixel_id
    NRDBias.pnts$X<-oneYr$Xcoord
    NRDBias.pnts$Y<-oneYr$Ycoord    
    
    DailyNRD.grid <- as(mergeRas, "SpatialGridDataFrame")
  }
  
  NRDPixel<-data.frame(DailyNRD.pnts$Pixel_id,DailyNRD.pnts$X,DailyNRD.pnts$Y)
  NRDxPixel<-data.frame(oneYr$Pixel_id, oneYr$XHARN, oneYr$YHARN)
  names(NRDxPixel)<-c("Pixel_id","X","Y")
  names(NRDPixel)<-c("Pixel_id","X","Y")
  NRDxPixel<- NRDxPixel[NRDxPixel$Y <= SouthMost,]
  mon <- as.numeric(monStr)
  RGdata <- na.omit(allStats[allStats$year == yearStr & 
                               allStats$month == mon &
                              allStats$adjF <400,])
  #   April 2001 has no RGdata values 
  if (nrow(RGdata)==0) {
    cat (paste("No Rain Gage data for",yearStr, monStr,'\n'))
    biasStuff <-
      list("Bias.ras"=NRDras,
           "gageR.pnts"=DailyNRD.pnts,
           "MonthlyRas"=NRDras,
           "MonNRDRas"=NRDras,
           "year"=as.numeric(yearStr),
           "month"=as.numeric(monStr),
           "monthlyNRDbias"=as.data.frame(NRDPixel[,c(1,2,3,ncol(NRDPixel))]),
           "dailyNRDbias" = as.data.frame(NRDPixel[,-c(ncol(NRDPixel))])
           )
  } else {
    #-------------------------------------------------
    # Interpolate Bias from RainGages to NRDRain pixels
    # Make NRD data correction using Bias
    #-------------------------------------------------
    rainGage.pnts<-SpatialPointsDataFrame(coords=RGdata[,c("XCOORD","YCOORD")],
                       data = RGdata,proj4string = HARNSP17ft)
    latlongPnts <- spTransform(rainGage.pnts,latlongs)
    distMatrix <-distm(latlongPnts)
    hc <- hclust(as.dist(distMatrix), method="complete")
    latlongPnts$clust <-cutree(hc,h=2000)
    latlongPnts$clust <-cutree(hc,h=4500)
    tempdf <- as.data.frame(latlongPnts)
    biasVals2 <- aggregate(cbind(XCOORD, YCOORD, adjF)~clust,tempdf,
                           mean,na.rm=TRUE) 
    biasVals<-biasVals2[, c("XCOORD", "YCOORD", "adjF")]
    biasVals <- biasVals[biasVals$XCOORD> NRDras@extent[1] 
                         & biasVals$XCOORD< NRDras@extent[2]
                         & biasVals$YCOORD> NRDras@extent[3] 
                         & biasVals$YCOORD< NRDras@extent[4],]
    biasData<-biasVals
    coordinates(biasData) =  ~ XCOORD + YCOORD
    proj4string(biasData) = HARNSP17ft
    RGdata<-fixDecimals(RGdata,3)
    biasData$XHARN <- coordinates(biasData)[, 1]
    biasData$YHARN <- coordinates(biasData)[, 2]
    biasData <- spTransform(biasData,HARNSP17ft)
    
    #-------------------------------------------------
    #  autoKrige implemented for Ordinary kriging
    #-------------------------------------------------
    surf <- autoKrige(formula=adjF ~ 1, input_data=biasData, 
                      new_data = DailyNRD.grid)
    biasRas <- raster(surf$krige_output)
    rainGage.pnts$bias <- raster::extract(biasRas,rainGage.pnts,
                                          fun=mean,df=TRUE)[,2]
    
    #-------------------------------------------------
    #  IDW raster code:
    # gs <- gstat(formula=adjF~1, locations=biasData)
    # idw <- interpolate(NRDras, gs)
    
    NRDbiasPixels <-NRDPixel
    names(NRDbiasPixels)<-c( "Pixel_id" ,"X" ,"Y")
    NRDbiasxPixels <-NRDxPixel
    names(NRDbiasxPixels)<-c( "Pixel_id" ,"X" ,"Y")
    
    #-------------------------------------------------
    #  Process NRDRain for each day calling "dayBiasFn"
    #  NRD with bias results are merged into a single table 
    #-------------------------------------------------
    monthFilter =paste0(yearStr,'-',monStr)
    dateFormList = as.Date(paste0(as.character(yr),'-01-01') ) + 
      as.numeric(gsub('day','',as.character(dayList)))-1
    dayLUp<-as.data.frame(cbind(dayList,
                                as.character(dateFormList)),
                          stringsAsFactors=F)
    dayLUp$V2<-as.character(dayLUp$V2)
    monList <- dayLUp[startsWith(dayLUp$V2 ,monthFilter),]$dayList
    iday = 0
    for (d in monList) {
      iday = iday + 1

      NRDras <- rasterize(DailyNRD.pnts, ras, DailyNRD.pnts$value, fun = mean)
      
      NRDxPnts<-as.data.frame(oneYr[oneYr$YHARN<=SouthMost,]$Pixel_id)
      NRDxPnts$X<-oneYr[oneYr$YHARN<=SouthMost,]$XHARN
      NRDxPnts$Y<-oneYr[oneYr$YHARN<=SouthMost,]$YHARN
      NRDxPnts$ID <- NULL
      names(NRDxPnts) <-c("Pixel_id","X","Y")
      OneDayRain<-Rain[!is.na(Rain$Gage) & 
                         Rain$DAILY_DATE == dayLUp[dayLUp$dayList==d,]$V2,c(1,7,8,9)]
      names(OneDayRain)<- c("Pixel_id","X","Y","value")
      table = dayBiasFn(NRDbyYr[NRDbyYr$variable == 
                                  dateList[as.numeric(rownames(dayLUp[dayLUp$dayList==d,]))], ],biasRas)
      names(table)<- c( "Pixel_id" ,"X" ,"Y",d)
      #-------------------------------------------------
      #  Theisen Polygon and raster code:
      if (nrow(southRain[southRain$DAILY_DATE==
                         dayLUp[dayLUp$dayList==d,]$V2,])>0){
        theisPoly <- 
          voronoi(southRain[southRain$DAILY_DATE==dayLUp[dayLUp$dayList==d,]$V2,])
        TheisRas <- rasterize(theisPoly, ras, theisPoly$VALUE, fun = mean)
        
        #-------------------------------------------------
        #   Extract Southern area pixels from TheisRas
        #-------------------------------------------------
        #NRDxPnts <- raster::extract(TheisRas,southRain,fun=mean,df=TRUE)
        NRDxPnts <-raster::extract(TheisRas,oneYr[oneYr$YHARN<=SouthMost,],
                                   fun=mean,df=TRUE)
        if (nrow(NRDxPnts[is.na(NRDxPnts$layer),])>0){
          NRDxPnts[is.na(NRDxPnts$layer),]$layer = 0.0
        }
     #   cat(max(NRDxPnts$layer),sep='\n')
        names(NRDxPnts)[names(NRDxPnts) =="layer"]<-c(d)
        
      } else {
        NRDxPnts$val<-0.0
        names(NRDxPnts)[names(NRDxPnts) =="val"]<-c(d)
      }
      NRDbiasxPixels<-cbind(NRDbiasxPixels,NRDxPnts[,c(d)])
      NRDbiasPixels<-merge(NRDbiasPixels, table, by =c("Pixel_id","X","Y"))
    }
    names(NRDbiasxPixels)<-c("Pixel_id","X","Y",monList)
    #-------------------------------------------------
    # Add final column for monthly total NRD with Bias correction
    # for export to csv
    #-------------------------------------------------
    
    NRDbiasPixels$Monthly<-rowSums(NRDbiasPixels[,-c(1,2,3)],na.rm=TRUE)
    NRDbiasxPixels$Monthly<-rowSums(NRDbiasxPixels[,-c(1,2,3)],na.rm=TRUE)
    names(NRDbiasPixels)[length(names(NRDbiasPixels))]<- sprintf('Mon%02d',mon)
    names(NRDbiasxPixels)[length(names(NRDbiasxPixels))]<-sprintf('Mon%02d',mon)
    NRDbiasPixels <- fixDecimals(NRDbiasPixels,4)
    NRDbiasxPixels <- fixDecimals(NRDbiasxPixels,4)
    #csvFile <- paste0(basePath, sprintf("biasNRD%s%02d.csv",yearStr,mon))
    #fwrite(NRDbiasPixels, csvFile) 
    NRDbiasPixels<-rbind(NRDbiasPixels,NRDbiasxPixels)
    
    #-------------------------------------------------
    # Add final column for monthly total NRD
    # for export to csv
    #-------------------------------------------------
    which(colnames(NRDrainByYr)%in%monList)
    #col.num <- which(colnames(NRDrainByYr@data) %in% dateList)
    oneMonth<-colnames(NRDrainByYr@data)[c(FALSE,FALSE,dayLUp$dayList %in% monList,FALSE,FALSE,FALSE)]
    NRDrainByYr$Monthly<-rowSums(NRDrainByYr@data[,oneMonth],na.rm=TRUE)
    NRDrainByYr <- fixDecimals(NRDrainByYr@data,4)
    lastCol=ncol(NRDrainByYr)
    MonNRDrain<-NRDrainByYr[,c("Pixel_id","XHARN","YHARN","Monthly")] 

    lastCol=ncol(NRDbiasxPixels)
    MonXtraNRDrain<-NRDbiasxPixels[,c(1,2,3,lastCol)]
    names(MonXtraNRDrain)<- names(MonNRDrain)
    MonNRDrain<- rbind(MonNRDrain,MonXtraNRDrain)
    coordinates(MonNRDrain) = ~ XHARN + YHARN
    proj4string(MonNRDrain) = HARNSP17ft
    MonNRDrain <- spTransform(MonNRDrain,HARNSP17ft) 
    
    coordinates(NRDbiasPixels) =  ~ X + Y
    proj4string(NRDbiasPixels) = HARNSP17ft
    
    NRDbiasPixels <- spTransform(NRDbiasPixels,HARNSP17ft)    
    #-------------------------------------------------
    # Create Raster for Monthly NRDRain rain
    #-------------------------------------------------
    MonNRDs<-NRDbiasPixels[,c("Pixel_id",sprintf('Mon%02d',mon))]
    monthlyVals <- as.data.frame(NRDbiasPixels[,sprintf('Mon%02d',mon)])
    NRDMonRas <- rasterize(MonNRDrain,ras, MonNRDrain$Monthly)
    
    #---------------------------------------------------------
    # Create Raster for bias corrected Monthly NRDRain rain
    #---------------------------------------------------------
    Monthly <-NRDbiasPixels[,c("Pixel_id",sprintf('Mon%02d',mon))]
    MonRas <- rasterize(Monthly,ras,Monthly@data[,2],fun=mean)    
    Ncols=ncol(NRDbiasPixels)
    monthlyNRDbias=as.data.frame(NRDbiasPixels[,c(1,Ncols)])
    dailyNRDbias=as.data.frame(NRDbiasPixels[,-c(Ncols)])
    if ('coords.x1' %in%  names(monthlyNRDbias)){
      setnames(monthlyNRDbias, 
               old =c('coords.x1', 'coords.x2'), 
               new = c('X','Y'),skip_absent=TRUE)
      setnames(monthlyNRDbias, 
               old =c('coords.x1', 'coords.x2'), 
               new = c('X','Y'),skip_absent=TRUE)
    }
    if ('coords.x1' %in%  names(dailyNRDbias)){
      setnames(dailyNRDbias, 
               old =c('coords.x1', 'coords.x2'), 
               new = c('X','Y'),skip_absent=TRUE)
      setnames(dailyNRDbias, 
               old =c('coords.x1', 'coords.x2'), 
               new = c('X','Y'),skip_absent=TRUE)
    }
    #---------------------------------------------------------
    #  Return a list of 4 spatial objects 
    #     (3 rasters & 1 points),
    #   Year and month
    #    and 2 dataframes    
    Ncols=ncol(NRDbiasPixels)
    biasStuff <-
      list("Bias.ras"=biasRas,
           "gageR.pnts"=rainGage.pnts,
           "MonthlyRas"=MonRas,
           "MonNRDRas"=NRDMonRas,
           "year"=as.numeric(yearStr),
           "month"=as.numeric(monStr),
           "monthlyNRDbias"=monthlyNRDbias,
           "dailyNRDbias" = dailyNRDbias
    )
  }
  return(biasStuff)
}

aggPnts<- function(gageR.pnts){
  #---------------------------------------------------
  # FUNCTION: aggPnts
  # Aggegates data using bias point location clusters
  #   calculating values to be used for bubble plots
  #--------------------------------------------------- 
  latlongPnts <- spTransform(gageR.pnts,latlongs)
  distMatrix <-distm(latlongPnts)
  hc <- hclust(as.dist(distMatrix), method="complete")
  #latlongPnts$clust <-cutree(hc,h=2000)
  latlongPnts$clust <-cutree(hc,h=4500)
  tempdf <- as.data.frame(latlongPnts)
  biasPnts <- aggregate(cbind(XCOORD, YCOORD, adjF, bias, sum_Rain,
                              sum_NRD)~clust, tempdf, mean,na.rm=TRUE) 
  biasPnts$Corrected <- biasPnts$sum_NRD*biasPnts$bias
  biasPnts$ObsVsAdj <- biasPnts$sum_Rain-biasPnts$Corrected
  biasPnts$ObsPcnt <- 100*(biasPnts$sum_Rain-biasPnts$sum_NRD)/biasPnts$sum_Rain
  # Fix division by 0
  if (nrow(biasPnts[!is.finite(biasPnts$ObsPcnt),])>0){
    biasPnts[!is.finite(biasPnts$ObsPcnt),]$ObsPcnt<-0
  }
  # limit values divided by very small values
  if (nrow(biasPnts[biasPnts$ObsPcnt>500,])>0){
    biasPnts[biasPnts$ObsPcnt>500,]$ObsPcnt<-0
  }  
  # ignore values divided by very small values
  if (nrow(biasPnts[biasPnts$ObsPcnt< -500,])>0){
    biasPnts[biasPnts$ObsPcnt< -500,]$ObsPcnt<- 0
  }
  # Fix division by 0
  biasPnts$AdjPcnt <- 100*(biasPnts$ObsVsAdj)/biasPnts$sum_Rain
  if (nrow(biasPnts[!is.finite(biasPnts$AdjPcnt),])>0){
    biasPnts[!is.finite(biasPnts$AdjPcnt),]$AdjPcnt<-0
  }  
  # limit values divided by very small values
  if (nrow(biasPnts[biasPnts$AdjPcnt>500,])>0){
    biasPnts[biasPnts$AdjPcnt>500,]$AdjPcnt<-0
  }  
  # limit values divided by very small values
  if (nrow(biasPnts[biasPnts$AdjPcnt< -500,])>0){
    biasPnts[biasPnts$AdjPcnt< -500,]$AdjPcnt<- -0
  }
  coordinates(biasPnts) =  ~ XCOORD + YCOORD
  proj4string(biasPnts) = HARNSP17ft
  return(biasPnts)
}  

plt1ras <- 
  function(fileName,bias.ras,bias.rasName,gageR.pnts,pltTheme,atVals,clpBnds2){
  #---------------------------------------
  # Create plot panel for a single raster 
  #---------------------------------------  
  panel1 = paste0(bias.rasName,' Bias Multiplier ', 
                  '\nRed \'+\' = Gage Location')  
  myplot=( levelplot(bias.ras, par.settings=pltTheme, main=panel1, 
                     colorkey=list(space="left"),
                     # ylab=list(label="\n",cex=.75),
                     scales = list(x=list(rot=45),y=list(rot=45),cex=.5),
                     at=atVals,layout=c(1,1),contour=FALSE, margin=F) +
             latticeExtra::layer(sp.polygons(clpBnds2)) +
             latticeExtra::layer(sp.text(coordinates(gageR.pnts),
                           txt=gageR.pnts$RainGage,pos=1,cex=.4 )) +
             latticeExtra::layer(sp.points(gageR.pnts, col = "red"))
  )
  return(myplot)
}

plt4ras <-
  function(outPath,pltName,nexR.ras,bias.ras,adjR.ras,diffR.ras,gageR.pnts,aggR.pnts){
  #-------------------------------------------------------------
  # Create plot files showing 4 raster comparison types by month
  #-------------------------------------------------------------
  fileName=paste(outPath,"BiasPlotsNRD/points",pltName,"X",".png",sep="")
  titleStr = paste(month.abb[as.numeric(substr(gsub("X","",pltName),5,6))],
                  substr(gsub("X","",pltName),1,4))
  #diffR.ras <-adjR.ras -nexR.ras
  aggR.pnts <-aggPnts(gageR.pnts)
  myTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
  divTheme = rasterTheme(region = brewer.pal('PiYG', n = 11))
  adjTheme1 = rasterTheme(region = brewer.pal('GnBu', n = 9))
  adjTheme2 = rasterTheme(region = brewer.pal('GnBu', n = 9))
  difTheme = rasterTheme(region = brewer.pal('RdBu', n = 11))
  # ------------------------
  # Monthly NRD Rain Raster 
  # ------------------------  
  pltTheme <- adjTheme1
  atVals <-c(seq(0,5,length=11),seq(6,15,length=4),seq(20,35,length=2))
  panel2 = paste0('Monthly NexRAD','\nBlack O = % Underestimate')

  NRDplot=( levelplot(nexR.ras, 
	                    par.settings=pltTheme, 
	                    main=panel2,at=atVals, 
											xlab = NULL, margin=F,
                      colorkey=FALSE,
                      scales = list(y=list(at=NULL),x=list(rot=45),cex=.5),
                      ylab=list(label="\n\n\n\n",cex=.75),
                      contour=FALSE) +
              latticeExtra::layer(sp.polygons(clpBnds2)) +
              latticeExtra::layer(sp.text(coordinates(aggR.pnts),
                            txt=as.character(round(aggR.pnts$ObsPcnt,0)),
                            pos=1,cex=.5 )) +
              latticeExtra::layer(sp.points(aggR.pnts,pch=1, 
                            cex=round(aggR.pnts$ObsPcnt/-100.,2),
                            col = "black")) + 
              latticeExtra::layer(sp.points(aggR.pnts,pch=1, 
                            cex=round(aggR.pnts$ObsPcnt/100.,2),
                            col = "red"))
            )
  # ------------------------
  # Monthly Bias Raster 
  # ------------------------
  pltTheme <- divTheme
  atVals <-c(.50,.55,.60,.65,.70,.75,.80,.85,.90,.95,
             1.0,1.1,1.2,1.3,1.4,1.5,2.,2.5, 3.,4.)  
  # uses "plt1ras" function which may be adapted to each of the panels
  Biasplot<-
    plt1ras(fileName, bias.ras, titleStr,gageR.pnts, pltTheme, atVals,clpBnds2)
  
  # -------------------------------------  
  #   Monthly Bias Adjusted NRDRain
  # -------------------------------------  
  pltTheme <- adjTheme2
  pal = brewer.pal(n=9,"RdGy")
  atVals <-c(seq(0,5,length=11),seq(6,15,length=4),seq(20,35,length=2))
  panel3 = paste0('Adjusted Monthly NexRAD','\nRed O = % Overestimate')
	
  NRDBiasplot=( levelplot(adjR.ras, par.settings=pltTheme, main=panel3, 
                          colorkey=list(space="left"),
                          at=atVals, xlab = NULL, 
													margin=FALSE,
                          scales = list(y=list(at=NULL),
                                        x=list(rot=45),cex=.5),
                          ylab=list(label="\n",cex=0.75),
                          contour=FALSE) +
                  latticeExtra::layer(sp.polygons(clpBnds2)) +
                  latticeExtra::layer(sp.text(coordinates(aggR.pnts),
                                txt=as.character(round(aggR.pnts$AdjPcnt,0)),
                                pos=1,cex=.5 )) +
                  latticeExtra::layer(sp.points(aggR.pnts,pch=1, 
                                cex=2*round(aggR.pnts$AdjPcnt/-100.,2),
                                col="black")) +             
                  latticeExtra::layer(sp.points(aggR.pnts,pch=1,
                                cex=2*round(aggR.pnts$AdjPcnt/100.,2),
                                col="red"))
                )
  
  # ------------------------------------------------------  
  # Monthly Difference [Bias Adjusted NRDrain - NRD]
  # ------------------------------------------------------  
  pltTheme <- difTheme
  atVals <-c(-20,-10,seq(-5,5,length=21),10,20)

  # Note = 'Black circles = pixels < observed'
  # Note = 'Red circles = pixels > observed'
  panel4 = paste0('Inches of Change','\nValues= +/-"@Gage Clusters')

  diffPlt=(levelplot(diffR.ras,par.settings=pltTheme,main=panel4,at=atVals,
                     xlab = NULL,
                     colorkey=list(space="right"),margin=FALSE,
                     scales = list(y=list(at=NULL),
                                   x=list(rot=45),cex=.5),
                     contour=FALSE) +
             latticeExtra::layer(sp.polygons(clpBnds2)) +
             latticeExtra::layer(sp.text(coordinates(aggR.pnts),
                           txt=as.character(round(aggR.pnts$ObsVsAdj,2)), 
                           pos=1,cex=.5 )) +
             latticeExtra::layer(sp.points(aggR.pnts,pch=1,
                           cex=round(aggR.pnts$ObsVsAdj*-1,2),
													 col = "red")) +             
             latticeExtra::layer(sp.points(aggR.pnts,pch=1,
                           cex=round(aggR.pnts$ObsVsAdj*1,2), 
													 col = "blue"))
           )
  # ------------------------------------------  
  # combine the plot panels to print
  # ------------------------------------------    
  trellis.device(device="png", filename=fileName, 
                 width=4000,height=2400,units="px",res=300)
  
  par.main.text=trellis.par.get('par.main.text')
  par.main.text$just = "left"
  par.main.text$cex = .5
  par.main.text$font = 1
  par.main.text$x = grid::unit(.4, "in")
  trellis.par.set('par.main.text',par.main.text)
  print(Biasplot, split    = c(1,1,4,1),more=TRUE)

  par.main.text$just = "right"
  par.main.text$x = grid::unit(3.0, "in")
  trellis.par.set('par.main.text',par.main.text)  
  print(NRDplot, split     = c(2,1,4,1),more=TRUE)

  par.main.text$just = "left"
  par.main.text$x = grid::unit(.75, "in")
  trellis.par.set('par.main.text',par.main.text)
  print(NRDBiasplot, split = c(3,1,4,1),more=TRUE)

  par.main.text$just = "right"
  par.main.text$x = grid::unit(2.75, "in")
  trellis.par.set('par.main.text',par.main.text)
  print(diffPlt, split  = c(4,1,4,1))
  dev.off()
  return(pltName)
}

########################
##    PROGRAM BODY    ##
########################

#------------------------------------------------------------
# Initialize variables for multiprocessing function calls
#------------------------------------------------------------
processed= listenv(NULL)
yrList=list()


GISPath <- "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/"
ProjPath <- "ECSM/Data/"
basePath <- paste0(GISPath,ProjPath)
outPath <-  paste0(basePath,"BiasPltNRD_0.75-4.0/")

#---------------------------------------------------
#  Read Rain Gage Data if not already in memory
#---------------------------------------------------
if(exists("gageR.df") && object.size(gageR.df) > 350000000 ){
  cat(paste('Skipped Reading data: ECSM_RainGageV2.csv','\n'))  
  } else {
  cat(paste('Reading data: ECSM_RainGageV2.csv','\n'))
  gageR.df = fread(paste0(basePath,"/ECSM_RainGageV2.csv"),sep=',')
  gageR.df$DAILY_DATE = as.Date(gageR.df$DAILY_DATE, format="%m/%d/%Y")
  gageR.df$YEAR = lubridate::year(gageR.df$DAILY_DATE)
  gageR.df$MONTH = lubridate::month(gageR.df$DAILY_DATE)  
  # gageR.df<-
  #   read.csv(paste0(basePath,"ECSM_RainGageV2.csv"),stringsAsFactors=FALSE)
  # gageR.df$YEAR<-
  #   format(as.Date(gageR.df$DAILY_DATE, format="%m/%d/%Y"),"%Y")
  # gageR.df$MONTH<-
  #   format(as.Date(gageR.df$DAILY_DATE, format="%m/%d/%Y"),"%m")
  # gageR.df$DAILY_DATE<-
  #   as.Date(format(as.Date(gageR.df$DAILY_DATE, 
  #                          format="%m/%d/%Y"),"%Y-%m-%d"))
  rainGages<- 
    unique(gageR.df[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rownames(rainGages) <- NULL
  rainGages<-na.omit(rainGages)
}
 
#---------------------------------------------------------
# Define raster mapping extents 
#   and
# Create dataSet: "oneYr" as framework for ECSM pixels
# based upon NRD cell spacing plus rain data South of WWM
#---------------------------------------------------------
cat(paste('Creating Rain pixel structure','\n'))
#---------------------------------------------------------
# Random Leap Year
#---------------------------------------------------------
yr =2000
oneYr<-read.csv(paste0(basePath,'ECSM_NRD_',yr,'.csv'))
i=0

#---------------------------------------------------------
# Define pixel Size
#   5280 * 2 = feet in 2 mile NRD grid
#   6561.679 = feet in 2 Kilometer NRD pixel spacing
#---------------------------------------------------------
halfPixel = 5280
halfPixel = 3280.84
pixSz = halfPixel*2

#---------------------------------------------------------
# calculate raster extents & number of rows and columns
#---------------------------------------------------------
xmin = floor(min(oneYr[c('X')])-halfPixel)
xmax = ceiling(max(oneYr[c('X')])+halfPixel)
ymin = floor(min(oneYr[c('Y')])-halfPixel - (24*(pixSz)))
ymax = ceiling(max(oneYr[c('Y')])+halfPixel)

rasRows <- floor((ymax - ymin) /(pixSz))
rasCols <- floor((xmax - xmin) /(pixSz))

ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,
              ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)

SouthMost<-min(oneYr$Y)-(pixSz) 
oneYr$Annual<-NULL
oneYr$ROWnum<-NULL
dateList<-names(oneYr)[-c(1,2,3)]
iletter =0
wideXtra<-NULL

for (l in LETTERS[24:1]){
  y= SouthMost-(iletter*pixSz)
  iletter = iletter + 1
  for (cols in  seq(1:rasCols)){
    x= xmin + ((cols-1)*(pixSz))
    emptyRow<-c(paste0(l,cols),x,y,replicate(length(dateList),NA))
    OneRowDF <- data.frame(t(emptyRow),stringsAsFactors=F)
    wideXtra<-rbind(wideXtra,OneRowDF)
  }
}

names(wideXtra) <-c("Pixel_id","X","Y",dateList)
oneYr<- rbind(oneYr,wideXtra)  
oneYr$X <- as.numeric(oneYr$X)
oneYr$Y <- as.numeric(oneYr$Y)

closest<- nn2(oneYr[,2:3],rainGages[,3:4],1)
index= closest[[1]]
dist=closest[[2]]

nnRowCo<-as.data.frame(cbind(rainGages$DBKEY,
                                     oneYr[index,]$Pixel_id,dist))
names(nnRowCo)<-c('DBKEY', 'Pixel_id', 'dist')
nnRowCo$DBKEY <- as.character(nnRowCo$DBKEY)
nnRowCo$dist <- as.numeric(as.character(nnRowCo$dist))

coordinates(oneYr) =  ~ X + Y
proj4string(oneYr) = HARNSP17ft
oneYr$XHARN <- coordinates(oneYr)[, 1]
oneYr$YHARN <- coordinates(oneYr)[, 2]
oneYr <- spTransform(oneYr,HARNSP17ft)

cat(paste('Establishing mutliprocessor Plan','\n'))
options(future.rng.onMisue = "ignore")

plan(multisession(workers = 8))
calcRainMP <- listenv()

#---------------------------------------------
# Process calcRainStats:
#   Setup status bar printing calcRainStats
#---------------------------------------------
cat(paste('Calculating Rain statistics by year','\n'))
cat(paste('\nWaiting for ',length(yrSeq), 
          'years of NRD rain to process','\n'))
for (dig in seq(1,4)){
  for (yr in seq(1,length(yrSeq))){
    cat(substring(as.character(yrSeq)[[yr]],dig:dig,dig:dig))
    cat(' ')
  }
  cat('\n')
}

iy=0
codeFilter=c('X','M','N','PT','?')
for (yr in yrSeq){
  cat('] ')
  iy = iy + 1
  yr = yrSeq[iy]
  RG<-gageR.df[gageR.df$YEAR==yr & 
               gageR.df$CODE %!in% codeFilter &
               gageR.df$VALUE >= 0.0,]
  # KAR 2024 note:  Recent changes seem to limit the use of 
  #         future to 5 workers
  #         for this portion of the bias correction due to memory
  #         limitiations at least for ECSM:
  # calcRainMP[[iy]] <- calcRainStats(basePath,yr,wideXtra,RG)
  # gc()
  calcRainMP[[iy]] <- future({calcRainStats(basePath,yr,wideXtra,RG)},.GlobalEnv=FALSE)
} 
cat('\n')

nyr = iy
rainStatData<-list()

for (i in seq(1:nyr)){
  rainStatData[[i]] <-future::value(calcRainMP[[i]])
  cat('^ ')
}
cat('\n')

gc()
# process calcRainStats function return results when not using futures
# for (i in seq(1:iy)){
#  mpList <-calcRainMP[[i]]
#  rainStatData[[i]] <-mpList
# }
cat('Saving Rain Stats image',sp='\n')
save.image("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ECSM_rain/rain75v40.RData")
#####
#  Read rain.RData
#####
# load("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ECSM_rain/rain75v40.RData")
# Remember to rerun the package checker
#-------------------------------------------------------
# unlist results returned from FUNCTION "calcRainStats"
#-------------------------------------------------------
#Rain <-rainStatData[[1]]$Rain
RainStats <- rainStatData[[1]]$RainStats
if (nyr > 1){
  for (iy in seq(2,nyr)){
    #Rain <-rbind(Rain,rainStatData[[iy]]$Rain)
    RainStats <- rbind(RainStats,rainStatData[[iy]]$RainStats)
  }
}
allStats<-RainStats 
names(allStats)<- c('RainGage','year','month','XCOORD','YCOORD',
                    'count','sum_Rain','sum_NRD','adjF')

yearStr = '2018'
monStr = '12'
year = 2018
mon = 12
southRain<-gageR.df[gageR.df$YEAR==year & 
                    gageR.df$MONTH==mon & 
                    gageR.df$CODE %!in% codeFilter &
                    gageR.df$YCOORD <= SouthMost,]

#Prepare 1 year of RaingGage Data South of NRD for voronoi process
coordinates(southRain) =  ~ XCOORD + YCOORD
proj4string(southRain) = HARNSP17ft
southRain$XHARN <- coordinates(southRain)[, 1]
southRain$YHARN <- coordinates(southRain)[, 2]
southRain <- spTransform(southRain,HARNSP17ft)

cat(paste('Establishing New mutliprocessor Plan','\n'))
#plan(multisession(workers = 10))
#calcRainMP <- listenv()

#-------------------------------------------------
# Define range of years to process
#-------------------------------------------------
 
processYears <- yrSeq
cat(paste('\nProcessing bias correction for each month of',
          length(yrSeq), 'years\n'))

options(future.globals.maxSize= 1073741824/2 ) # 1/2 Gig

for (processYears in yrSeq[1:19]) {
  gc()
  cat(paste('Establishing New mutliprocessor Plan','\n'))
  plan(multisession(workers = 9))

  processed <- listenv()
  cat(paste('year',processYears),sep='\n')
  x=0
  yr = processYears[1]
  mon = 1
  for (yr in processYears) {
    yearStr <- as.character(yr)
    cat (paste(yearStr,':'))
    Rain <-data.frame()
    Rain <-rbind(Rain,rainStatData[[match(yr, yrSeq)]]$Rain)
    
    for (mon in seq(1,12)){
      monStr <-sprintf("%02d",mon)
      x=x+1
      cat (paste(monStr," "))
      # Prepare single year of Rain Gage data south of NRD pixel area
      #   for Theissen polygon processing
      southRain<-gageR.df[gageR.df$YEAR==yearStr & 
                          gageR.df$MONTH==mon &
                          gageR.df$YCOORD < SouthMost,]
      coordinates(southRain) =  ~ XCOORD + YCOORD
      proj4string(southRain) = HARNSP17ft
      southRain$XHARN <- coordinates(southRain)[, 1]
      southRain$YHARN <- coordinates(southRain)[, 2]
      southRain <- spTransform(southRain,HARNSP17ft)
      
      
      #-----------------------------------------------------------    
      # multiprocessor function call    
      #-----------------------------------------------------------    
      processed[[x]] <- future({biasByYearMon(basePath,allStats,as.data.frame(Rain[Rain$MONTH== mon,]),
                                              southRain,yearStr,monStr,x)})
      #-----------------------------------------------------------    
      # single processor function call    
      #-----------------------------------------------------------    
      # processed[[x]] <- biasByYearMon(basePath,allStats,Rain,
      #                                 southRain,yearStr,monStr,x)
    }
    cat('\n')
  }
  
  #-------------------------------------------------------
  # value function waits for processed results to become 
  # available for each process
  #-------------------------------------------------------
  cat(paste('Waiting for bias correction processes to finalize','\n'))
  mpList<-list()
  rList<-list()
  # retrieve multiprocessor function call result  values 
  for (i in seq(1:x)){
    mpList <-future::value(processed[[i]])
    rList[[i]]<-mpList
  }
  
  # process function call reurn results when not using futures
  # for (i in seq(1:x)){
  #  mpList <-processed[[i]]
  #  rList[[i]]<-mpList
  # }
  
  #-------------------------------------------------
  # unlist results returned from FUNCTION "biasByYearMon"
  #-------------------------------------------------
  bias.ras.List = unlist(lapply(rList,"[[",1))
  pointList = unlist(lapply(rList,"[[",2))
  AnnRainList = unlist(lapply(rList,"[[",3))
  AnnNRDList = unlist(lapply(rList,"[[",4))
  yrList = unlist(lapply(rList,"[[",5))
  MonthList = unlist(lapply(rList,"[[",6))
  
  filePath<- outPath
  cat(paste('Exporting csv files','\n'))
  for (iyr in processYears){
    NRDbiasMonthly <-rList[[1]]$monthlyNRDbias[,c("Pixel_id","X","Y")]
    NRDbiasDaily <- rList[[1]]$dailyNRDbias[,c("Pixel_id","X","Y")]
    
  
    for (i in seq(from=1, to=x)){
      if (rList[[i]]$year == iyr){
        NRDbiasMonthly <- merge(NRDbiasMonthly, 
   rList[[i]]$monthlyNRDbias[,!names(rList[[i]]$monthlyNRDbias) 
                             %in% c("X","Y")],by =c("Pixel_id"))
        NRDbiasDaily <- merge(NRDbiasDaily, 
   rList[[i]]$dailyNRDbias[,!names(rList[[i]]$dailyNRDbias) 
                           %in% c("X","Y")],by =c("Pixel_id"))
      }
    }
    if(length(NRDbiasMonthly)>4){
      NRDbiasMonthly$Annual<-
        rowSums(NRDbiasMonthly[,!names(NRDbiasMonthly) 
                               %in% c("Pixel_id","X","Y")],na.rm=TRUE)
    }
    csvFile <- paste0(filePath, sprintf("MonthlybiasXNRD%04d.csv",iyr))
    cat(sprintf("MonthlybiasXNRD%04d.csv",iyr))
    cat(' ')
    fwrite(NRDbiasMonthly, csvFile) 
    NRDbiasDaily$Annual<-
      rowSums(NRDbiasDaily[,!names(NRDbiasDaily) 
                           %in% c("Pixel_id","X","Y")],na.rm=TRUE)
    csvFile <- paste0(filePath, sprintf("DailybiasXNRD%04d.csv",iyr))
    cat(sprintf("DailybiasXNRD%04d.csv",iyr))
    cat('\n')
    fwrite(NRDbiasDaily, csvFile) 
  }
  
  # Do not save this image.  Currently trying to output more than 2.9GB
  #save.image("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ECSM_rain/bias.RData")
  
  cat(paste('producing maps','\n'))
  bias.ras.Stack <-stack()
  bias.ras.Stack <-stack(bias.ras.List)
  adjNex.ras.Stack <-stack()
  adjNex.ras.Stack <-stack(AnnRainList)
  nexR.ras.Stack <-stack()
  nexR.ras.Stack <-stack(AnnNRDList)
  
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
  names(bias.ras.Stack)<- StackNames
  names(adjNex.ras.Stack)<- StackNames
  names(nexR.ras.Stack)<- StackNames
  
  ps<-listenv()
  cat(paste('Plotting', nlayers(bias.ras.Stack), 'rain comparison maps \n'))
  
  for (i in 1:nlayers(bias.ras.Stack)){
    gageR.pnts<- pointList[[i]]
    aggR.pnts <-aggPnts(gageR.pnts)
    # Plot Monthly NRDRain
    pltName <-names(nexR.ras.Stack)[i]
    nexR.ras <- nexR.ras.Stack[[i]]
    bias.ras <- bias.ras.Stack[[i]]
    adjR.ras <- adjNex.ras.Stack[[i]]
    diffR.ras <- adjR.ras - nexR.ras
  
    #------------------------------   
    # multiprocessor function call    
    #------------------------------  
    ps[[i]] %<-% plt4ras(outPath,pltName,
                         nexR.ras,bias.ras,adjR.ras,diffR.ras,
                         gageR.pnts,aggR.pnts)
    #------------------------------   
    # single processor function call    
    #------------------------------   
    # ps[[i]] <- plt4ras(outPath,pltName,
    #                    nexR.ras,bias.ras,adjR.ras,diffR.ras,
    #                    gageR.pnts,aggR.pnts)  
    #------------------------------   
    txt<-sprintf('%03d',i)
    if (i%/%10 == (i/4)){cat(txt)} else { cat('.') }
    if (i%/%60 == (i/12)){cat('\n')}
  }
  cat(paste0(txt,'\n'))
}
