pkgChecker <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}

list.of.packages <-  c("reshape2","readr","dplyr","tidyr", "data.table",
                       "readxl","rgeos","sp", "dismo","lattice","rasterVis",
                       "maptools","raster","fields","automap", "gstat",
                       "future","listenv","ggplot2","RANN","geosphere")

suppressWarnings(pkgChecker(list.of.packages))
yrSeq<-seq(1985,2018)
#yrSeq<-seq(1985,1986)

'%!in%' <- function(x,y)!('%in%'(x,y))

fixDecimals <- function(DF,decPlaces){
  is.num <-sapply(DF,is.numeric)
  DF[is.num] <- lapply(DF[is.num], round,decPlaces)
  return(DF)
}

gClip <- function(shp, bb) {
  #------------------------------------------------------------
  # Set up county boundry shapefile for overlay on raster maps
  #------------------------------------------------------------
  if (class(bb) == "matrix")
    b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else
    b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

calcRainStats<-function(basePath,yr,wideXtra,RG){
  # FUNCTION: calcRainStats
  #   [Works well with future function for multiprocessing]
  # Calculate Series of Monthly Rain Stats for 1 year from gage data
  #---------------------------------------------------
  rainGages<- unique(RG[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rownames(rainGages) <- NULL
  #ECSM_WMM<- read_csv(paste0("G:/ECSM/Data/wmmRain",yr,".csv"))
  ECSM_WMM<- read.csv(paste0(basePath,"wmmRain",yr,".csv"))
  ECSM_WMM<-rename(ECSM_WMM,ROWnum=X)
  ECSM_WMM<-rename(ECSM_WMM,wmmRoCo=variable)
  ECSM_WMM$ROWnum = NULL
  ECSM_WMM$Annual = NULL
  col.num=0
  # Drop column day 366 from non-Leap years
  if( yr/4 - (yr%/%4) != 0) {
    col.num <- which(colnames(wideXtra) %!in% colnames(ECSM_WMM))
    wideXtra[,col.num]<- NULL
  }
  ECSM_WMM<- rbind(ECSM_WMM,wideXtra)
  
  meltWMM<-melt(ECSM_WMM,id=c("wmmRoCo","Xcoord","Ycoord"))
  meltWMM<-meltWMM[meltWMM$variable != 'Annual',]
  meltWMM[is.na(meltWMM$value),]$value = 0.0
  meltWMM <- na.omit(meltWMM)
  meltWMM$daily_date = as.Date(paste0(as.character(yr),'-01-01') ) + 
    as.numeric(gsub('day','',as.character(meltWMM$variable)))-1
  meltWMM$variable = as.character(meltWMM$daily_date)
  
  LongWdates<-cbind(meltWMM[meltWMM$variable < "A",] %>% separate(variable, c("year", "month", "day")))
  LongWdates$daily_date <- as.Date(LongWdates$daily_date, format="%Y-%m-%d")
  LongWdates$value <- as.numeric(LongWdates$value)
  
  #
  # Subset Gage data for the year and filter out specific codes
  #
  RainG<- merge(x=RG,   y=nnRowCo, by.x= 'DBKEY', by.y='DBKEY')
  Rain<-merge(x=RainG,  y=LongWdates[,c('wmmRoCo','daily_date','value')],
              by.x=c('wmmRoCo','DAILY_DATE'),by.y=c('wmmRoCo','daily_date'))
  names(Rain)<-c("wmmRoCo","DAILY_DATE", "DBKEY","STATION","AGENCY","XCOORD",
                 "YCOORD","Gage","CODE","YEAR","MONTH","dist","wmmRain" )
  Rain$STATION <- as.factor(Rain$STATION)
  Rain<-filter(Rain, (Gage >= 0.01 & Gage < 20 & wmmRain >= 0.01 & wmmRain < 20) | 
                 (Gage <= 0.01 & wmmRain <= 0.01))
  Rain <- na.omit(Rain)
  
  #
  # Calculate Stats for Gage and wmmRain bias
  #
  Rain.monthly.sum<- Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize_at(vars(Gage, wmmRain), sum) %>% 
    rename(Gage.sum = Gage, wmmRain.sum = wmmRain)
  
  Rain.monthly.obs<- Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize_at(vars(Gage), ~sum(. >-.001)) %>% 
    rename(Gage.obs = Gage)  %>% 
    subset(Gage.obs >= 14)  
  
  Corr_GagewmmRain<-Rain %>%
    group_by(STATION,DBKEY,MONTH) %>%
    summarize(correlation = cor(Gage, wmmRain, method = "kendall")) %>%
    drop_na() %>% 
    subset(correlation >= .75 & correlation <= 1.0)
  
  RainStats<-  merge(rainGages,merge(Corr_GagewmmRain,
                                     merge(x=Rain.monthly.sum,y=Rain.monthly.obs),all.y=TRUE),all.x =TRUE)
  if (nrow(RainStats[is.na(RainStats$correlation),])> 0) {
    RainStats[is.na(RainStats$correlation),]$correlation = 1.0
  }
  if (nrow(RainStats[is.na(RainStats$Gage.sum),])> 0) {
    RainStats[is.na(RainStats$Gage.sum),]$Gage.sum = 0.0
  }
  if (nrow(RainStats[is.na(RainStats$wmmRain.sum),])> 0) {
    RainStats[is.na(RainStats$wmmRain.sum),]$wmmRain.sum = 0.0
  }
  RainStats<-na.omit(RainStats)
  RainStats[RainStats$YCOORD <= SouthMost,]$wmmRain.sum<-RainStats[RainStats$YCOORD <= SouthMost,]$Gage.sum
  rainGages<- unique(RG[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rainGages<-na.omit(rainGages)
  
  #
  # Rainstats maybe useful for optimizing bias correction
  #
  write.csv(RainStats[!is.na(RainStats$MONTH),],paste0('h:/WMMrainstats',yr,'.csv'))
  
  #
  # Filter RainStats
  #
  RainStats$bias = RainStats$Gage.sum/RainStats$wmmRain.sum
  RainStats[is.nan(RainStats$bias),]$bias = 1.0
  RainStats<-fixDecimals(RainStats,3)
  
  #
  # Prevent over-correction
  #
  if( sum(is.nan(RainStats$bias)) >0 ){
    RainStats[is.nan(RainStats$bias),]$bias <-1
  }
  if( sum(is.na(RainStats$bias)) >0 ){
    RainStats[is.na(RainStats$bias),]$bias <-1
  }
  if (dim(RainStats[RainStats$bias<.5,])[1] >0){
    RainStats[RainStats$bias<.5,]$bias <- .5
  }
  if (dim(RainStats[RainStats$bias>3.,])[1] >0){
    RainStats[RainStats$bias>3,]$bias <- 3.0
  }
  # RainStats[RainStats$Gage.obs < SignificantDays,]$bias <- 1
  RainStats$YEAR<-yr
  RainStats<-merge(nnRowCo, merge(rainGages,RainStats))
  Rain<-merge(Rain,RainStats[,c("DBKEY","wmmRoCo","STATION","MONTH","bias")])
  Rain$WMM<-Rain$wmmRain*Rain$bias
  
  keepCols = c("STATION","YEAR","MONTH","XCOORD","YCOORD","Gage.obs",
               "Gage.sum","wmmRain.sum","bias")
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
  
  # ToPlot<- Rain[, !names(Rain) %in% c('STATION',"AGENCY","XCOORD","YCOORD","CODE","YEAR","MONTH","dist","bias")]
  # #ToPlot<- Rain[, !names(Rain) %in% c('wmmRoCo',"AGENCY","XCOORD","YCOORD","CODE","YEAR","MONTH","dist","bias")]
  # #testPlot=melt(ToPlot,id=c('STATION','DBKEY','DAILY_DATE','Gage'))
  # testPlot=melt(ToPlot,id=c('wmmRoCo','DBKEY','DAILY_DATE','Gage'))
  
  # for (stn in unique(testPlot$wmmRoCo)){
  #   cat (stn,sep='\n')
  # 
  #   filename = paste("G:/ECSM/Data/graphWMMcorrection/",stn,"_",yr,".png", sep = "" )
  #   png(  file = filename, width = 3000,height = 3000,units = "px",  res=300)
  # 
  #   p <- ggplot(testPlot[testPlot$wmmRoCo==stn,],aes(x=value, y=Gage,
  #                                            color = paste(DBKEY,variable,sep='_'),shape=DBKEY)) +
  #     labs(title =stn, color = 'wmmRain') +
  #     geom_point() +geom_smooth(method = "lm")
  #   print(p)
  #   dev.off()
  # }
}

dayBiasFn <- function(DailyWMM,biasRas){
  #-------------------------------------------------
  # FUNCTION: dayBiasFn
  # Multiplies Daily WMM rasters by bias  
  #   calculating adjusted daily wmmRain Raster
  # and returns 
  #-------------------------------------------------  
  DailyWMM.pnts <-SpatialPointsDataFrame(coords = DailyWMM[, c("X", "Y")],
                                         data = DailyWMM,proj4string = HARNSP17ft)
  WMMras <-rasterize(DailyWMM.pnts, ras, DailyWMM.pnts$value, fun = mean) * biasRas
  # gs <- gstat(formula=adjF~1, locations=biasData)
  # idw <- interpolate(WMMras, gs)
  WMMBiasPnts <- raster::extract(WMMras,DailyWMM.pnts,fun=mean,df=TRUE)
  WMMBiasPnts$wmmRoCo <- DailyWMM.pnts$wmmRoCo
  WMMBiasPnts$X<-DailyWMM.pnts$X
  WMMBiasPnts$Y<-DailyWMM.pnts$Y
  return(WMMBiasPnts[c(3,4,5,2)])
}

biasByYearMon <-function(basePath,allStats,Rain,southernRain,
                         yearStr,monStr,x){
  yr = as.numeric(yearStr)
  #-------------------------------------------------
  # FUNCTION: biasByYearMon
  #   [Works well with future function for multiprocessing]
  # Processes RainVsGage data by year
  # Creating CSV files with WMMstat, biasWMM and 
  # updates WMM with annual totals
  #-------------------------------------------------  
  # read and organize daily wmmRain data
  #-------------------------------------------------
  WMMrainByYr<- read.csv(paste0(basePath,"wmmRain",yr,".csv"),
                         stringsAsFactors=F)
  
  WMMrainByYr<-rename(WMMrainByYr,ROWnum=X)
  WMMrainByYr<-rename(WMMrainByYr,wmmRoCo=variable)
  WMMrainByYr<-rename(WMMrainByYr,X=Xcoord)
  WMMrainByYr<-rename(WMMrainByYr,Y=Ycoord)
  
  WMMbyYr <- melt(WMMrainByYr, id = c("ROWnum","wmmRoCo", "X", "Y"))
  WMMbyYr$ROWnum<-NULL
  # SouthMost<-min(WMMrainByYr$Y)-(2*5280)
  coordinates(WMMrainByYr) =  ~ X + Y
  proj4string(WMMrainByYr) = HARNSP17ft
  WMMrainByYr$XHARN <- coordinates(WMMrainByYr)[, 1]
  WMMrainByYr$YHARN <- coordinates(WMMrainByYr)[, 2]
  WMMrainByYr <- spTransform(WMMrainByYr,HARNSP17ft)
  
  dateFormList = as.Date(paste0(as.character(yr),'-01-01') ) + 
    as.numeric(gsub('day','',as.character(dateList)))-1
  dayLUp<-as.data.frame(cbind(dateList,
                              as.character(dateFormList)),
                        stringsAsFactors=F)
  dayLUp$V2<-as.character(dayLUp$V2)
  
  for (d in dateList[1]) {
    #-------------------------------------------------
    #  Theisen Polygon and raster code to fill in southern area:
    firstDay<-unique(southernRain$DAILY_DATE)[1]
    theisPoly <- 
      voronoi(southernRain[southernRain$DAILY_DATE==firstDay,])
    TheisRas <- rasterize(theisPoly, ras, theisPoly$VALUE, fun = mean)
    
    DailyWMM <- WMMbyYr[WMMbyYr$variable == d, ]
    DailyWMM.pnts <-SpatialPointsDataFrame(coords = DailyWMM[, c("X", "Y")],
                                           data = DailyWMM, proj4string = HARNSP17ft)
    
    WMMras <-rasterize(DailyWMM.pnts, ras, DailyWMM.pnts$value, fun = mean) 
    # Do raster overlay Here:
    s<-stack(WMMras,TheisRas)
    mergeRas<-calc(s,fun=function(x) ifelse(is.na(x[1]), 
                                            ifelse(is.na(x[2]),0,x[2]),x[1]))
    WMMBiasPnts <- data.frame(raster::extract(mergeRas, oneYr))
    WMMBiasPnts$wmmRoCo <- oneYr$wmmRoCo
    WMMBiasPnts$X<-oneYr$Xcoord
    WMMBiasPnts$Y<-oneYr$Ycoord    
    
    DailyWMM.grid <- as(mergeRas, "SpatialGridDataFrame")
  }
  WMMPixel<- data.frame(DailyWMM.pnts$wmmRoCo, DailyWMM.pnts$X, DailyWMM.pnts$Y)
  WMMxPixel<- data.frame(oneYr$wmmRoCo, oneYr$XHARN, oneYr$YHARN)
  names(WMMxPixel)<-c("wmmRoCo","X","Y")
  names(WMMPixel)<-c("wmmRoCo","X","Y")
  WMMxPixel<- WMMxPixel[WMMxPixel$Y <= SouthMost,]
  mon <- as.numeric(monStr)
  RGdata <- na.omit(allStats[allStats$year == yearStr & 
                               allStats$month == monStr
                             & allStats$adjF <400,])
  #   April 2001 has no RGdata values 
  
  if (nrow(RGdata)==0) {
    cat (paste("No Rain Gage data for",yearStr, monStr,'\n'))
    biasStuff <-list("B_ras"=WMMras,
                     "R_pnts"=DailyWMM.pnts,
                     "MonthlyRas"=WMMras,
                     "MonWMMRas"=WMMras,
                     "year"=as.numeric(yearStr),
                     "month"=as.numeric(monStr),
                     "monthlyWMMbias"=as.data.frame(WMMPixel[,c(1,2,3,ncol(WMMPixel))]),
                     "dailyWMMbias" = as.data.frame(WMMPixel[,-c(ncol(WMMPixel))])
    )
  } else
  {
    #-------------------------------------------------
    # Interpolate Bias from RainGages to wmmRain pixels
    # Make WMM data correction using Bias
    #-------------------------------------------------
    rainGage.pnts <- SpatialPointsDataFrame(coords = RGdata[,c("XCOORD", "YCOORD")],
                                            data = RGdata,proj4string = HARNSP17ft)
    #<<<< Delete this if it works>>>>
    # rainGage.pnts <- SpatialPointsDataFrame(coords = RGdata[RGdata$YCOORD > SouthMost, c("XCOORD", "YCOORD")],
    #                                         data = RGdata,proj4string = HARNSP17ft)    
    latlongPnts <- spTransform(rainGage.pnts,latlongs)
    distMatrix <-distm(latlongPnts)
    hc <- hclust(as.dist(distMatrix), method="complete")
    #latlongPnts$clust <-cutree(hc,h=2000)
    latlongPnts$clust <-cutree(hc,h=4500)
    tempdf <- as.data.frame(latlongPnts)
    
    biasVals2 <- aggregate(cbind(XCOORD, YCOORD, adjF)~clust, tempdf, mean,na.rm=TRUE) 
    biasVals<-biasVals2[, c("XCOORD", "YCOORD", "adjF")]
    
    biasVals <- biasVals[biasVals$XCOORD> WMMras@extent[1] 
                         & biasVals$XCOORD< WMMras@extent[2]
                         & biasVals$YCOORD> WMMras@extent[3] 
                         & biasVals$YCOORD< WMMras@extent[4],]
    
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
    surf <- autoKrige(formula=adjF ~ 1, input_data=biasData, new_data = DailyWMM.grid)
    biasRas <- raster(surf$krige_output)
    rainGage.pnts$bias <- raster::extract(biasRas,rainGage.pnts,fun=mean,df=TRUE)[,2]
    
    #-------------------------------------------------
    #  IDW raster code:
    # gs <- gstat(formula=adjF~1, locations=biasData)
    # idw <- interpolate(WMMras, gs)
    
    WMMbiasPixels <-WMMPixel
    names(WMMbiasPixels)<-c( "wmmRoCo" ,"X" ,"Y")
    WMMbiasxPixels <-WMMxPixel
    names(WMMbiasxPixels)<-c( "wmmRoCo" ,"X" ,"Y")
    
    #-------------------------------------------------
    #  Process wmmRain for each day calling "dayBiasFn"
    #  WMM with bias results are merged into a single table 
    #-------------------------------------------------
    monthFilter =paste0(yearStr,'-',monStr)
    dateFormList = as.Date(paste0(as.character(yr),'-01-01') ) + 
      as.numeric(gsub('day','',as.character(dateList)))-1
    dayLUp<-as.data.frame(cbind(dateList,
                                as.character(dateFormList)),
                          stringsAsFactors=F)
    dayLUp$V2<-as.character(dayLUp$V2)
    monList <- dayLUp[startsWith(dayLUp$V2 ,monthFilter),]$dateList
    iday = 0
    for (d in monList) {
      iday = iday + 1
      #-------------------------------------------------
      #  Theisen Polygon and raster code:
      theisPoly <- 
        voronoi(southernRain[southernRain$DAILY_DATE==
                               dayLUp[dayLUp$dateList==d,]$V2,])
      TheisRas <- rasterize(theisPoly, ras, theisPoly$VALUE, fun = mean)
      
      #-------------------------------------------------
      #   Extract Southern area pixels from TheisRas
      #-------------------------------------------------
      #WMMxPnts <- raster::extract(TheisRas,southernRain,fun=mean,df=TRUE)
      WMMxPnts <- raster::extract(TheisRas,oneYr[oneYr$YHARN<=SouthMost,],fun=mean,df=TRUE)
      WMMras <- rasterize(DailyWMM.pnts, ras, DailyWMM.pnts$value, fun = mean)
      WMMxPnts$wmmRoCo<-oneYr[oneYr$YHARN<=SouthMost,]$wmmRoCo
      WMMxPnts$X<-oneYr[oneYr$YHARN<=SouthMost,]$XHARN
      WMMxPnts$Y<-oneYr[oneYr$YHARN<=SouthMost,]$YHARN
      WMMxPnts$ID <- NULL
      names(WMMxPnts) <-c(d,"wmmRoCo","X","Y")
      OneDayRain<-Rain[!is.na(Rain$Gage) & 
                         Rain$DAILY_DATE == dateFormList[iday],c(1,7,8,9)]
      names(OneDayRain)<- c("wmmRoCo","X","Y","value")
      table = dayBiasFn(WMMbyYr[WMMbyYr$variable == d, ],biasRas)
      names(table) <- c( "wmmRoCo" ,"X" ,"Y",d  )
      WMMbiasxPixels <- merge(WMMbiasxPixels, WMMxPnts, by =c("wmmRoCo","X","Y"))
      WMMbiasPixels <- merge(WMMbiasPixels, table, by =c("wmmRoCo","X","Y"))
    }
    #-------------------------------------------------
    # Add final column for monthly total WMM with Bias correction
    # for export to csv
    #-------------------------------------------------
    
    WMMbiasPixels$Monthly<-rowSums(WMMbiasPixels[,-c(1,2,3)],na.rm=TRUE)
    WMMbiasxPixels$Monthly<-rowSums(WMMbiasxPixels[,-c(1,2,3)],na.rm=TRUE)
    names(WMMbiasPixels)[length(names(WMMbiasPixels))]<- sprintf('Mon%02d',mon)
    names(WMMbiasxPixels)[length(names(WMMbiasxPixels))]<- sprintf('Mon%02d',mon)
    WMMbiasPixels <- fixDecimals(WMMbiasPixels,4)
    WMMbiasxPixels <- fixDecimals(WMMbiasxPixels,4)
    #csvFile <- paste0(basePath, sprintf("biasWMM%s%02d.csv",yearStr,mon))
    #fwrite(WMMbiasPixels, csvFile) 
    WMMbiasPixels<-rbind(WMMbiasPixels,WMMbiasxPixels)
    
    #-------------------------------------------------
    # Add final column for monthly total WMM
    # for export to csv
    #-------------------------------------------------
    col.num <- which(colnames(WMMrainByYr@data) %in% monList)
    WMMrainByYr$Monthly<-rowSums(WMMrainByYr@data[,col.num],na.rm=TRUE)
    WMMrainByYr <- fixDecimals(WMMrainByYr@data,4)
    lastCol=ncol(WMMrainByYr)
    MonWMMrain<-WMMrainByYr[,c(2,lastCol-2,lastCol-1,lastCol)] 
    lastCol=ncol(WMMbiasxPixels)
    MonXtraWMMrain<-WMMbiasxPixels[,c(1,2,3,lastCol)]
    names(MonXtraWMMrain)<- names(MonWMMrain)
    MonWMMrain<- rbind(MonWMMrain,MonXtraWMMrain)
    coordinates(MonWMMrain) = ~ XHARN + YHARN
    proj4string(MonWMMrain) = HARNSP17ft
    MonWMMrain <- spTransform(MonWMMrain,HARNSP17ft) 
    
    coordinates(WMMbiasPixels) =  ~ X + Y
    proj4string(WMMbiasPixels) = HARNSP17ft
    
    WMMbiasPixels <- spTransform(WMMbiasPixels,HARNSP17ft)    
    #-------------------------------------------------
    # Create Raster for Monthly wmmRain rain
    #-------------------------------------------------
    MonWMMs<-WMMbiasPixels[,c("wmmRoCo",sprintf('Mon%02d',mon))]
    monthlyVals <- as.data.frame(WMMbiasPixels[,sprintf('Mon%02d',mon)])
    WMMMonRas <- rasterize(MonWMMrain,ras, MonWMMrain$Monthly)
    
    #---------------------------------------------------------
    # Create Raster for bias corrected Monthly wmmRain rain
    #---------------------------------------------------------
    Monthly <-WMMbiasPixels[,c(1,ncol(WMMbiasPixels))]
    MonRas <- rasterize(Monthly,ras,Monthly@data[,2],fun=mean)    
    #---------------------------------------------------------
    #  Return a list of 4 spatial objects 
    #     (3 rasters & 1 points),
    #   Year and month
    #    and 2 dataframes    
    biasStuff <-list("B_ras"=biasRas,
                     "R_pnts"=rainGage.pnts,
                     "MonthlyRas"=MonRas,
                     "MonWMMRas"=WMMMonRas,
                     "year"=as.numeric(yearStr),
                     "month"=as.numeric(monStr),
                     "monthlyWMMbias"=as.data.frame(WMMbiasPixels[,c(1,ncol(WMMbiasPixels))]),
                     "dailyWMMbias" = as.data.frame(WMMbiasPixels[,-c(ncol(WMMbiasPixels))])
    )
  }
  
  #cat(paste(names(as.data.frame(WMMbiasPixels[,-c(ncol(WMMbiasPixels))])),'\n'))
  return(biasStuff)
}

aggPnts<- function(pntsPlt){
  #-------------------------------------------------
  # FUNCTION: aggPnts
  # Aggegates data using bias point location clusters
  #   calculating values to be used for bubble plots
  #------------------------------------------------- 
  latlongPnts <- spTransform(pntsPlt,latlongs)
  distMatrix <-distm(latlongPnts)
  hc <- hclust(as.dist(distMatrix), method="complete")
  #latlongPnts$clust <-cutree(hc,h=2000)
  latlongPnts$clust <-cutree(hc,h=4500)
  tempdf <- as.data.frame(latlongPnts)
  biasPnts <- aggregate(cbind(XCOORD, YCOORD, adjF, bias, sum_Rain, sum_WMM)~clust, tempdf, mean,na.rm=TRUE) 
  biasPnts$Corrected <- biasPnts$sum_WMM*biasPnts$bias
  biasPnts$ObsVsAdj <- biasPnts$sum_Rain-biasPnts$Corrected
  coordinates(biasPnts) =  ~ XCOORD + YCOORD
  proj4string(biasPnts) = HARNSP17ft
  return(biasPnts)
}  

plotOneRas <- function(filename, rasPlt, rasPltName, pntsPlt, 
                       pltTheme, atVals,clpBnds2){
  #-------------------------------------------------
  # Create plot files for each raster type by month
  #-------------------------------------------------  
  panel1 = paste('Bias Multiplier',rasPltName, '\nRed Crosses = Rain Gages')
  
  myplot=( levelplot(rasPlt, par.settings=pltTheme, main=panel1, 
                     at=atVals,layout=c(1,1),contour=FALSE, margin=F) +
             latticeExtra::layer(sp.polygons(clpBnds2)) +
             latticeExtra::layer(sp.text(coordinates(pntsPlt),txt=pntsPlt$RainGage,pos=1,cex=.5 )) +
             latticeExtra::layer(sp.points(pntsPlt, col = "red"))
  )
  return(myplot)
}


plt4rasters <-function(outPath,PltName,
                       MonWMMPlt,rasPlt,AdjPlt,diffPlt,pntsPlt){

  filename=paste(outPath,"BiasPlotsWMM/points",PltName,"X",".png",sep="")
  #diffPlt <-AdjPlt -MonWMMPlt
  aggPts <-aggPnts(pntsPlt)
  myTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
  divTheme = rasterTheme(region = brewer.pal('PiYG', n = 11))
  adjTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
  difTheme = rasterTheme(region = brewer.pal('RdBu', n = 11))
  
  pltTheme <- adjTheme
  atVals <-c(0,seq(5,20,length=16),30)
  panel2 = paste('Uncorrected ',PltName,'\n','Black O = Underestimate')
  WMMplot=( levelplot(MonWMMPlt, par.settings=pltTheme, main=panel2, 
                      at=atVals, xlab = NULL, margin=F,
                      layout=c(1,1),contour=FALSE) +
              latticeExtra::layer(sp.polygons(clpBnds2)) +
              latticeExtra::layer(sp.text(coordinates(aggPts),
                                          txt=as.character(round((aggPts$sum_Rain-aggPts$sum_WMM),2)),
                                          pos=1,cex=.5 )) +
              latticeExtra::layer(sp.points(aggPts,pch=1, 
                                            cex=round((aggPts$sum_Rain-aggPts$sum_WMM)*.6,2),
                                            col = "black")) +             
              latticeExtra::layer(sp.points(aggPts,pch=1, 
                                            cex=round((aggPts$sum_Rain-aggPts$sum_WMM)*-.6,2),
                                            col = "red"))
  )
  
  # Plot Monthly Bias Raster 
  pltTheme <- divTheme
  #atVals <-seq(0,4,length=21)
  atVals <-c(.50,.55,.60,.65,.70,.75,.80,.85,.90,.95,
             1.0,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.0)
  
  Biasplot=plotOneRas(filename, rasPlt, PltName,pntsPlt, pltTheme, atVals,clpBnds2)
  
  #   Plot Monthly Bias Adjusted wmmRain
  
  pltTheme <- adjTheme
  atVals <-c(0,seq(5,20,length=16),30)
  panel3 = paste('Adjusted ',PltName,'\n','Red O = Overestimate')
  WMMBiasplot=( levelplot(AdjPlt, par.settings=pltTheme, main=panel3, 
                          at=atVals, xlab = NULL, margin=F,
                          layout=c(1,1),contour=FALSE) +
                  latticeExtra::layer(sp.polygons(clpBnds2)) +
                  latticeExtra::layer(sp.text(coordinates(aggPts),
                                              txt=as.character(round(aggPts$ObsVsAdj,2)),
                                              pos=1,cex=.5 )) +
                  latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                                cex=round((aggPts$ObsVsAdj*.6),2),col="black")) +             
                  latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                                cex=round((aggPts$ObsVsAdj*-.6),2),col="red"))
  )
  #   Plot Monthly Bias Adjusted wmmRain - WMM
  pltTheme <- difTheme
  atVals <-c(-15,-10,-6,seq(-3,3,length=20),6,10,15)
  
  if( sum(is.na(pntsPlt@data$bias)) >0 ){
    pntsPlt@data[is.na(pntsPlt@data$bias),]$bias = 1
  }
  # Note = 'Black circles = pixels < observed'
  # Note = 'Red circles = pixels > observed'
  panel4 = paste('Inches of Change\n',PltName)
  diffPlt=(levelplot(diffPlt,par.settings=pltTheme,main=panel4,at=atVals,
                     xlab = NULL,margin=F,layout=c(1,1),contour=FALSE) +
             latticeExtra::layer(sp.polygons(clpBnds2)) +
             latticeExtra::layer(sp.text(coordinates(aggPts),
                                         txt=as.character(round(aggPts$ObsVsAdj,2)), 
                                         pos=1,cex=.5 )) +
             latticeExtra::layer(sp.points(pntsPlt,pch=1,
                                           cex=aggPts$ObsVsAdj*.6, col = "black")) +             
             latticeExtra::layer(sp.points(pntsPlt,pch=1,
                                           cex=aggPts$ObsVsAdj*-.6, col = "red"))
  )
  trellis.device(device="png", filename=filename, 
                 width=4800,height=2400,units="px",res=300)
  print(Biasplot, split    = c(1,1,4,1),more=TRUE)
  print(WMMplot, split     = c(2,1,4,1),more=TRUE)
  print(WMMBiasplot, split = c(3,1,4,1),more=TRUE)
  print(diffPlt, split  = c(4,1,4,1))
  dev.off()
  return(PltName)
}

#---------------------------------------------------
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
#---------------------------------------------------
HARNSP17ft  = CRS("+init=epsg:2881")
latlongs = CRS("+proj=longlat +datum=WGS84")
WMDbnd.Path <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)
GISPath <- "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/"
ProjPath <- "ECSM/Data/"
basePath <- paste0(GISPath,ProjPath)
#---------------------------------------------------
#  Read Rain Gage Data if not already in memory
#---------------------------------------------------
if(exists("rainGageData") && object.size(rainGageData) > 350000000 ){
  cat(paste('Skipped Reading data: ECSM_RainGageV2.csv','\n'))  
} else 
{
  cat(paste('Reading data: ECSM_RainGageV2.csv','\n'))
  rainGageData<-read.csv(paste0(basePath,"ECSM_RainGageV2.csv"),stringsAsFactors = FALSE)
  rainGageData$YEAR = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y")
  rainGageData$MONTH = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%m")
  rainGageData$DAILY_DATE = as.Date(format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y-%m-%d"))
  rainGages<- unique(rainGageData[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rownames(rainGages) <- NULL
  rainGages<-na.omit(rainGages)
}
 

#---------------------------------------------------------
# Define raster mapping extents 
#   and
# Create dataSet: "oneYr" as framework for ECSM pixels
# based upon WMM cell spacing plus rain data South of WWM
#---------------------------------------------------------
cat(paste('Creating Rain pixel structure','\n'))

i=0
#---------------------------------------------------------
# Random Leap Year
#---------------------------------------------------------
yr =2000
oneYr<-read.csv(paste0(basePath,'wmmRain',yr,'.csv'))

#---------------------------------------------------------
# Define pixel Size
#   5280 * 2 = feet in 2 mile WMM grid
# 6561.679 = feet in 2 Kilometer NRD pixel spacing
#---------------------------------------------------------
halfPixel = 5280
pixSz = halfPixel*2

#---------------------------------------------------------
# calculate raster extents & number of rows and columns
#---------------------------------------------------------
xmin = floor(min(oneYr[c('Xcoord')])-halfPixel)
xmax = ceiling(max(oneYr[c('Xcoord')])+halfPixel)
ymin = floor(min(oneYr[c('Ycoord')])-halfPixel - (15*(pixSz)))
ymax = ceiling(max(oneYr[c('Ycoord')])+halfPixel)

rasRows <- floor((ymax - ymin) /(pixSz))
rasCols <- floor((xmax - xmin) /(pixSz))

ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,
              ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)

SouthMost<-min(oneYr$Ycoord)-(pixSz) 
oneYr<-rename(oneYr,wmmRoCo=variable)
oneYr$Annual<-NULL
oneYr$X<-NULL
dateList<-names(oneYr)[-c(1,2,3)]
iletter =0
wideXtra<-NULL
for (l in LETTERS[15:1]){
  y= SouthMost-(iletter*pixSz)
  iletter = iletter + 1
  for (cols in  seq(1:rasCols)){
    x= xmin + ((cols-1)*(pixSz))
    emptyRow<-c(paste0(l,cols),x,y,replicate(length(dateList),NA))
    OneRowDF <- data.frame(t(emptyRow),stringsAsFactors=F)
    wideXtra<-rbind(wideXtra,OneRowDF)
  }
}
names(wideXtra) <-c("wmmRoCo","Xcoord","Ycoord",dateList)
oneYr<- rbind(oneYr,wideXtra)  
oneYr$Xcoord <- as.numeric(oneYr$Xcoord)
oneYr$Ycoord <- as.numeric(oneYr$Ycoord)

closest<- nn2(oneYr[,2:3],rainGages[,3:4],1)

index= closest[[1]]
dist=closest[[2]]

nnRowCo<-as.data.frame(cbind(rainGages$DBKEY,
                                     oneYr[index,]$wmmRoCo,dist))
names(nnRowCo)<-c('DBKEY', 'wmmRoCo', 'dist')
nnRowCo$DBKEY <- as.character(nnRowCo$DBKEY)
nnRowCo$dist <- as.numeric(as.character(nnRowCo$dist))

coordinates(oneYr) =  ~ Xcoord + Ycoord
proj4string(oneYr) = HARNSP17ft
oneYr$XHARN <- coordinates(oneYr)[, 1]
oneYr$YHARN <- coordinates(oneYr)[, 2]
oneYr <- spTransform(oneYr,HARNSP17ft)

cat(paste('Establishing mutliprocessor Plan','\n'))
plan(multiprocess,.skip=TRUE)
calcRainMP <- listenv()

#---------------------------------------------
# Process calcRainStats:
#   Setup status bar printing calcRainStats
#---------------------------------------------
cat(paste('Calculating Rain statistics by year','\n'))
cat(paste('\nWaiting for ',nyr, 
          'years of WMM rain bias correction to process','\n'))
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
  RG<-rainGageData[rainGageData$YEAR==yr & 
                     rainGageData$CODE %!in% codeFilter &
                     rainGageData$VALUE >= 0.0,]
  #calcRainMP[[iy]] <- calcRainStats(basePath,yr,wideXtra,RG)
  calcRainMP[[iy]] <- future({calcRainStats(basePath,yr,wideXtra,RG)})
} 

nyr = iy
rainStat_data<-list()

for (i in seq(1:nyr)){
  rainStat_data[[i]] <-value(calcRainMP[[i]])
}
#-------------------------------------------------------
# unlist results returned from FUNCTION "calcRainStats"
#-------------------------------------------------------
Rain <-rainStat_data[[1]]$Rain
RainStats <- rainStat_data[[1]]$RainStats
if (nyr > 1){
  for (iy in seq(2,nyr)){
    Rain <-rbind(Rain,rainStat_data[[iy]]$Rain)
    RainStats <- rbind(RainStats,rainStat_data[[iy]]$RainStats)
  }
}
allStats<-RainStats 
names(allStats)<- c('RainGage','year','month','XCOORD','YCOORD',
                    'count','sum_Rain','sum_WMM','adjF')

myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))

yearStr = '2017'
monStr = '09'
southernRain<-rainGageData[rainGageData$YEAR==yearStr & 
                         rainGageData$MONTH==monStr &
                         rainGageData$YCOORD <= SouthMost,]

#Prepare 1 year of RaingGage Data South of WMM for voronoi process
coordinates(southernRain) =  ~ XCOORD + YCOORD
proj4string(southernRain) = HARNSP17ft
southernRain$XHARN <- coordinates(southernRain)[, 1]
southernRain$YHARN <- coordinates(southernRain)[, 2]
southernRain <- spTransform(southernRain,HARNSP17ft)

#------------------------------------------------------------
# Initialize variables for multiprocessing function calls
#------------------------------------------------------------

processed= listenv(NULL)
yrList=list()

yearStr <- as.character(2017)
#-------------------------------------------------
# Define range of years to process
#-------------------------------------------------
processYears <- yrSeq

cat(paste('Processing bias correction for each month of',
          length(yrSeq), 'years\n'))
x=0
yr = processYears[1]
for (yr in processYears) {
  yearStr <- as.character(yr)
  cat (paste(yearStr,':'))
  for (mon in seq(1,12)){
    monStr <-sprintf("%02d",mon)
    x=x+1
    cat (paste(monStr," "))
    # Prepare single year of RaingGage Data South of WMM 
    #   for Theisen polygon processing
    southernRain<-rainGageData[rainGageData$YEAR==yearStr & 
                                 rainGageData$MONTH==monStr &
                                 rainGageData$YCOORD < SouthMost,]
    coordinates(southernRain) =  ~ XCOORD + YCOORD
    proj4string(southernRain) = HARNSP17ft
    southernRain$XHARN <- coordinates(southernRain)[, 1]
    southernRain$YHARN <- coordinates(southernRain)[, 2]
    southernRain <- spTransform(southernRain,HARNSP17ft)
    
    #-----------------------------------------------------------
    # Call FUNCTION "biasByYearMon" 
    #   with futures multiprocessing wrapper function
    #-----------------------------------------------------------    
    # multiprocessor function call    
    processed[[x]] <- future({biasByYearMon(basePath,allStats,Rain,
                                            southernRain,yearStr,monStr,x)})
    # single processor function call    
    # processed[[x]] <- biasByYearMon(basePath,allStats,Rain,
    #                                 southernRain,yearStr,monStr,x)
  }
  cat('\n')
}

#-------------------------------------------------
# value function waits for results to become available
# for each process
#-------------------------------------------------
cat(paste('Waiting for bias correction processes to return','\n'))

mpList<-list()
rList<-list()

# retrieve multiprocessor function call result  values 
for (i in seq(1:x)){
  mpList <-value(processed[[i]])
  rList[[i]]<-mpList
}

# process function call reurn results
# for (i in seq(1:x)){
#  mpList <-processed[[i]]
#  rList[[i]]<-mpList
# }

#-------------------------------------------------
# unlist results returned from FUNCTION "biasByYearMon"
#-------------------------------------------------

stackList = unlist(lapply(rList,"[[",1))
pointList = unlist(lapply(rList,"[[",2))
AnnRainList = unlist(lapply(rList,"[[",3))
AnnWMMList = unlist(lapply(rList,"[[",4))
yrList = unlist(lapply(rList,"[[",5))
MonthList = unlist(lapply(rList,"[[",6))

outPath <-  paste0(basePath,"BiasPlotsWMMb0.5-3.0/")
filePath<- outPath
cat(paste('Exporting csv files','\n'))
for (iyr in processYears){
  WMMbiasMonthly <-rList[[1]]$monthlyWMMbias[,c("wmmRoCo","X","Y")]
  WMMbiasDaily <- rList[[1]]$dailyWMMbias[,c("wmmRoCo","X","Y")]

  for (i in seq(from=1, to=x)){
    if (rList[[i]]$year == iyr){
      WMMbiasMonthly <- merge(WMMbiasMonthly, 
 rList[[i]]$monthlyWMMbias[,!names(rList[[i]]$monthlyWMMbias) 
                           %in% c("X","Y")],by =c("wmmRoCo"))
      WMMbiasDaily <- merge(WMMbiasDaily, 
 rList[[i]]$dailyWMMbias[,!names(rList[[i]]$dailyWMMbias) 
                         %in% c("X","Y")],by =c("wmmRoCo"))
    }
  }
  if(length(WMMbiasMonthly)>4){
    WMMbiasMonthly$Annual<-
      rowSums(WMMbiasMonthly[,!names(WMMbiasMonthly) 
                             %in% c("wmmRoCo","X","Y")],na.rm=TRUE)
  }
  csvFile <- paste0(filePath, sprintf("MonthlybiasXWMM%04d.csv",iyr))
  cat(sprintf("MonthlybiasXWMM%04d.csv",iyr))
  cat(' ')
  fwrite(WMMbiasMonthly, csvFile) 
  WMMbiasDaily$Annual<-
    rowSums(WMMbiasDaily[,!names(WMMbiasDaily) 
                         %in% c("wmmRoCo","X","Y")],na.rm=TRUE)
  csvFile <- paste0(filePath, sprintf("DailybiasXWMM%04d.csv",iyr))
  cat(sprintf("DailybiasXWMM%04d.csv",iyr))
  cat('\n')
  fwrite(WMMbiasDaily, csvFile) 
}

cat(paste('producing maps','\n'))
rasStack <-stack()
rasStack <-stack(stackList)
MonRasStack <-stack()
MonRasStack <-stack(AnnRainList)
MonWMMStack <-stack()
MonWMMStack <-stack(AnnWMMList)

#yearList <-as.character(processYears)
#names(rasStack)<- yearList
#names(AnnrasStack)<- yearList
#names(AnnWMMStack)<- yearList

yearList <-as.character(processYears)
StackNames = list()
x = 0
yr = processYears[1]
for (yr in processYears) {
  for (mon in seq(1, 12)) {
    #for (mon in seq(9, 12)) {
    x = x + 1
    StackNames[x] <- sprintf("%4d%02d",yr, mon)
  }
}
names(rasStack)<- StackNames
names(MonRasStack)<- StackNames
names(MonWMMStack)<- StackNames


ps<-listenv()
cat(paste('Plotting', nlayers(rasStack), 'rain comparison maps \n'))
for (i in 1:nlayers(rasStack)){
  pntsPlt<- pointList[[i]]
  # Plot Monthly wmmRain
  PltName <-names(MonWMMStack)[i]
  MonWMMPlt <- MonWMMStack[[i]]
  rasPlt <- rasStack[[i]]
  AdjPlt <- MonRasStack[[i]]
  diffPlt <- AdjPlt - MonWMMPlt
  ps[[i]] <- plt4rasters(outPath,PltName,
                           MonWMMPlt,rasPlt,AdjPlt,diffPlt,pntsPlt)
  if (i%/%10 == (i/10)){cat(i)} else { cat('.') }
}
cat('\n')
