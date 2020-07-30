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



'%!in%' <- function(x,y)!('%in%'(x,y))

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

WMDbnd.Path <- "//whqhpc01p/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)

rainGageData<-read.csv("//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/ECSM/Data/ECSM_RainGage.csv",stringsAsFactors = FALSE)
rainGageData$YEAR = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y")
rainGageData$MONTH = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%m")
rainGageData$DAILY_DATE = as.Date(format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y-%m-%d"))

rainGages<- unique(rainGageData[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
rownames(rainGages) <- NULL
my_data <- list()
i=0
#for (yr in seq(2000,2018)){
for (yr in seq(1985,2018)){
#for (yr in seq(1985,1986)){
    i = i + 1
  cat(yr,sep='/n')
  oneYr<-read.csv(paste0('G:/ECSM/Data/wmmRain',yr,'.csv'))
  oneYr<-rename(oneYr,wmmRoCo=variable)
  closest<- nn2(oneYr[,3:4],rainGages[,3:4],1)
  index= closest[[1]]
  distance=closest[[2]]
  nearestwWMMRoCo<-as.data.frame(cbind(rainGages$DBKEY,oneYr[index,]$wmmRoCo,  distance))
  names(nearestwWMMRoCo)<-c('DBKEY', 'wmmRoCo', 'Distance')
  nearestwWMMRoCo$wmmRoCo <- as.numeric(as.character(nearestwWMMRoCo$wmmRoCo))
  nearestwWMMRoCo$DBKEY <- as.character(nearestwWMMRoCo$DBKEY)
  nearestwWMMRoCo$Distance <- as.numeric(as.character(nearestwWMMRoCo$Distance))
  
  ECSM_WMM<- read_csv(paste0("G:/ECSM/Data/wmmRain",yr,".csv"))
  ECSM_WMM<-rename(ECSM_WMM,ROWnum=X1)
  ECSM_WMM<-rename(ECSM_WMM,wmmRoCo=variable)
  ECSM_WMM<-rename(ECSM_WMM,X=Xcoord)
  ECSM_WMM<-rename(ECSM_WMM,Y=Ycoord)
  
  meltWMM<-melt(ECSM_WMM,id=c("ROWnum","wmmRoCo","X","Y"))
  meltWMM<-meltWMM[meltWMM$variable != 'Annual',]
  meltWMM$daily_date = as.Date(paste0(as.character(yr),'-01-01') ) + 
    as.numeric(gsub('day','',as.character(meltWMM$variable)))-1
  meltWMM$variable = as.character(meltWMM$daily_date)
  LongWdates<-cbind(meltWMM[meltWMM$variable < "A",] %>% separate(variable, c("year", "month", "day")))
  LongWdates$daily_date <- as.Date(LongWdates$daily_date, format="%Y-%m-%d")
  #
  # Subset Gage data for the year and filter out specific codes
  #
  RainG<- merge(x=rainGageData[rainGageData$YEAR==yr & rainGageData$CODE %!in% c('X','M','N','PT'),],
                y=nearestwWMMRoCo, by.x= 'DBKEY', by.y='DBKEY')
  
  Rain<-merge(x=RainG,  y=LongWdates[,c('wmmRoCo','daily_date','value')],
              by.x=c('wmmRoCo','DAILY_DATE'),by.y=c('wmmRoCo','daily_date'))
  
  names(Rain)<-c("wmmRoCo","DAILY_DATE", "DBKEY","STATION","AGENCY","XCOORD","YCOORD","Gage","CODE",
                 "YEAR","MONTH","Distance","wmmRain" )
  
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
  RainStats[is.na(RainStats$correlation),]$correlation = 1.0
  RainStats[is.na(RainStats$Gage.sum),]$Gage.sum = 0.0
  RainStats[is.na(RainStats$wmmRain.sum),]$wmmRain.sum = 0.0
  write.csv(RainStats,paste0('h:/WMMrainstats',yr,'.csv'))
  
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
  if (dim(RainStats[RainStats$bias>4.,])[1] >0){
    RainStats[RainStats$bias>4,]$bias <- 4.0
  }
 # RainStats[RainStats$Gage.obs < SignificantDays,]$bias <- 1
  RainStats$YEAR<-yr
  RainStats<-merge(nearestwWMMRoCo, merge(rainGages,RainStats))
  
  Rain<-merge(Rain,RainStats[,c("DBKEY","wmmRoCo","STATION","MONTH","bias")])
  Rain$WMM<-Rain$wmmRain*Rain$bias
  my_data[[i]] <- RainStats[,c("STATION","YEAR","MONTH","XCOORD","YCOORD","Gage.obs","Gage.sum","wmmRain.sum","bias")]
  
  ToPlot<- Rain[, !names(Rain) %in% c('STATION',"AGENCY","XCOORD","YCOORD","CODE","YEAR","MONTH","Distance","bias")]
  #ToPlot<- Rain[, !names(Rain) %in% c('wmmRoCo',"AGENCY","XCOORD","YCOORD","CODE","YEAR","MONTH","Distance","bias")]
  #testPlot=melt(ToPlot,id=c('STATION','DBKEY','DAILY_DATE','Gage'))
  testPlot=melt(ToPlot,id=c('wmmRoCo','DBKEY','DAILY_DATE','Gage'))
  
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

allStats<-do.call(rbind,my_data)
names(allStats)<- c('RainGage','year','month','XCOORD','YCOORD','count','sum_Rainfall','sum_WMM','adjF')

myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))

basePath <-  "G:/ECSM/Data/"

meltWMM<- na.omit(melt(ECSM_WMM,id=c("ROWnum","wmmRoCo","X","Y")))
PixelCoords <- meltWMM[c("wmmRoCo", "X", "Y")]

#-------------------------------------------------
# Add additional melted WMM data sets 
# to calculate bias based on more than one year
#-------------------------------------------------
WMM <- do.call("rbind", list(meltWMM))

#-------------------------------------------------
# Calculate raster extents
#-------------------------------------------------
xmin = floor(min(PixelCoords[c('X')])-5280)
xmax = ceiling(max(PixelCoords[c('X')])+5280)
ymin = floor(min(PixelCoords[c('Y')])-5280)
ymax = ceiling(max(PixelCoords[c('Y')])+5280)

#-------------------------------------------------
# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
#-------------------------------------------------
# rasRows <- (ymax - ymin) / 6561.679
# rasCols <- (xmax - xmin) / 6561.679
rasRows <- (ymax - ymin) /(5280*2)
rasCols <- (xmax - xmin) /(5280*2)

#-------------------------------------------------
# define raster and map extents using WMM pixel data extents
#-------------------------------------------------
ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)



#-------------------------------------------------
# FUNCTION: dayBiasFn
# Multiplies Daily WMM rasters by bias  
#   calculating adjusted daily wmmRain Raster
# and returns 
#-------------------------------------------------
dayBiasFn <- function(DailyWMM,biasRas){
  DailyWMM.pnts <-SpatialPointsDataFrame(coords = DailyWMM[, c("X", "Y")],
                                         data = DailyWMM,proj4string = HARNSP17ft)
  WMMras <-rasterize(DailyWMM.pnts, ras, DailyWMM.pnts$value, fun = mean) * biasRas
  WMMBiasPnts <- raster::extract(WMMras,DailyWMM.pnts,fun=mean,df=TRUE)
  WMMBiasPnts$wmmRoCo <- DailyWMM.pnts$wmmRoCo
  WMMBiasPnts$X<-DailyWMM.pnts$X
  WMMBiasPnts$Y<-DailyWMM.pnts$Y
  return(WMMBiasPnts[c(3,4,5,2)])
}

#-------------------------------------------------
# FUNCTION: biasByYearMon
#   [Works well with future function for multiprocessing]
# Processes RainVsGage data by year
# Creating CSV files with WMMstat, biasWMM and 
# updates WMM with annual totals
#-------------------------------------------------
yearStr = '2001'
monStr = '05'

biasByYearMon <-function(yearStr,monStr,x){
  #ECSM_WMMbyYr<- read_csv(paste0(basePath,"ECSM_WMM_",yearStr,".csv"))
  ECSM_WMMbyYr<- read_csv(paste0("G:/ECSM/Data/wmmRain",yr,".csv"))
  ECSM_WMMbyYr<-rename(ECSM_WMMbyYr,ROWnum=X1)
  ECSM_WMMbyYr<-rename(ECSM_WMMbyYr,wmmRoCo=variable)
  ECSM_WMMbyYr<-rename(ECSM_WMMbyYr,X=Xcoord)
  ECSM_WMMbyYr<-rename(ECSM_WMMbyYr,Y=Ycoord)
  
  dateList<-names(ECSM_WMMbyYr)[-c(1,2,3,4,length(names(ECSM_WMMbyYr)))]
  WMMbyYr <- melt(ECSM_WMMbyYr, id = c("ROWnum","wmmRoCo", "X", "Y"))
  WMMbyYr$ROWnum<-NULL
  
  #-------------------------------------------------
  # read and organize daily wmmRain data
  #-------------------------------------------------
  for (d in dateList[1]) {
    DailyWMM <- WMMbyYr[WMMbyYr$variable == d, ]
    DailyWMM.pnts <-SpatialPointsDataFrame(coords = DailyWMM[, c("X", "Y")],
                                           data = DailyWMM, proj4string = HARNSP17ft)
    
    WMMras <-rasterize(DailyWMM.pnts, ras, DailyWMM.pnts$value, fun = mean) 
    #    WMMBiasPnts <- data.frame(extract(WMMras, DailyWMM.pnts))
    DailyWMM.grid <- as(WMMras, "SpatialGridDataFrame")
  }
  WMMPixel<- data.frame(DailyWMM.pnts$wmmRoCo, DailyWMM.pnts$X, DailyWMM.pnts$Y)
  
  mon <- as.numeric(monStr)
  RGdata <- na.omit(allStats[allStats$year == yearStr & 
                       allStats$month == monStr
                             & allStats$adjF <400,])
#
  #
    #   April 2001 has no RGdata values 
  
  if (nrow(RGdata)==0) {
    cat (paste(yearStr, monStr,'\n'))
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
    # biasData <- na.omit(RGdata[, c("XCOORD", "YCOORD", "adjF")])
    
    #  coordinates(biasData) =  ~ XCOORD + YCOORD
    #  proj4string(biasData) = HARNSP17ft
    
    #  biasData$X <- coordinates(biasData)[, 1]
    #  biasData$Y <- coordinates(biasData)[, 2]
    
    rainGage.pnts <- SpatialPointsDataFrame(coords = RGdata[, c("XCOORD", "YCOORD")],
                                            data = RGdata,proj4string = HARNSP17ft)
    latlongPnts <- spTransform(rainGage.pnts,latlongs)
    distMatrix <-distm(latlongPnts)
    hc <- hclust(as.dist(distMatrix), method="complete")
    latlongPnts$clust <-cutree(hc,h=2000)
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
    
    #WMMras <- rasterize(DailyWMM.pnts, ras, DailyWMM.pnts$value, fun = mean)
    
    #-------------------------------------------------
    #  Theisen Polygon and raster code:
    # theisPoly <- voronoi(rainGage.pnts)
    # TheisRas <- rasterize(theisPoly, WMMras, theisPoly$bias, fun = mean)
    
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
    #-------------------------------------------------
    #  Process wmmRain for each day calling "dayBiasFn"
    #  WMM with bias results are merged into a single table 
    #-------------------------------------------------
    monthFilter =paste0(yearStr,'-',monStr)
    #cat(paste0('X',yearStr,'.',monStr,'\n'))
    #cat(paste(monthFilter,'\n'))
    dateFormList = as.Date(paste0(as.character(yr),'-01-01') ) + 
      as.numeric(gsub('day','',as.character(dateList)))-1
    dayLU<-as.data.frame(cbind(dateList,as.character(dateFormList)))
    dayLU$V2<-as.character(dayLU$V2)
    monList <- dayLU[startsWith(dayLU$V2 ,monthFilter),]$dateList
    #cat(monList)
    for (d in monList) {
      #   cat(paste("."))
      #cat(paste0(d,'\n'))
      table = dayBiasFn(WMMbyYr[WMMbyYr$variable == d, ],biasRas)
      names(table) <- c( "wmmRoCo" ,"X" ,"Y",d  )
      WMMbiasPixels <- merge(WMMbiasPixels, table, by =c("wmmRoCo","X","Y"))
    }
    #-------------------------------------------------
    # Add final column for monthly total WMM with Bias correction
    # and export to csv
    #-------------------------------------------------
    
    WMMbiasPixels$Monthly<-rowSums(WMMbiasPixels[,-c(1,2,3)],na.rm=TRUE)
    names(WMMbiasPixels)[length(names(WMMbiasPixels))]<- sprintf('Mon%02d',mon)
    WMMbiasPixels <- fixDecimals(WMMbiasPixels,4)
    #csvFile <- paste0(basePath, sprintf("biasWMM%s%02d.csv",yearStr,mon))
    #fwrite(WMMbiasPixels, csvFile) 
    
    #-------------------------------------------------
    # Create Raster for bias corrected Monthly wmmRain rain
    #-------------------------------------------------
    Monthly <-WMMbiasPixels[,c(1,2,3,ncol(WMMbiasPixels))]
    xy <- Monthly[,c(2,3)]
    Monthlypdf <-SpatialPointsDataFrame(coords=xy,data=Monthly,proj4string=HARNSP17ft)
    MonRas <- rasterize(Monthlypdf,ras,Monthlypdf@data[,4],fun=mean)
    
    #-------------------------------------------------
    # Add final column for monthly total WMM
    # and export to csv
    #-------------------------------------------------
    col.num <- which(colnames(ECSM_WMMbyYr) %in% monList)
    #names(ECSM_WMMbyYr)
    ECSM_WMMbyYr$Monthly<-rowSums(ECSM_WMMbyYr[,col.num],na.rm=TRUE)
    #csvFile <- paste0(basePath, sprintf("WMM%s%02d.csv",yearStr,mon))
    ECSM_WMMbyYr <- fixDecimals(ECSM_WMMbyYr,4)
    #fwrite(WMMwCoords[,c(1,6,unlist(col.num),ncol(WMMwCoords))], csvFile) 
    
    #-------------------------------------------------
    # Create Raster for Monthly wmmRain rain
    #-------------------------------------------------
    MonWMMs<-ECSM_WMMbyYr[,c(2,3,4,ncol(ECSM_WMMbyYr))]
    xy <- MonWMMs[,c(2,3)]
    # MonWMMsspdf <-SpatialPointsDataFrame(coords=xy,data=MonWMMs,proj4string=latlongs)
    MonWMMsspdf <-SpatialPointsDataFrame(coords=xy,data=MonWMMs,proj4string=HARNSP17ft)
    #  MonWMMsspdf <-spTransform(MonWMMsspdf,HARNSP17ft)
    WMMMonRas <- rasterize(MonWMMsspdf,ras,MonWMMsspdf$Monthly,fun=mean)
    #-------------------------------------------------
    biasStuff <-list("B_ras"=biasRas,
                     "R_pnts"=rainGage.pnts,
                     "MonthlyRas"=MonRas,
                     "MonWMMRas"=WMMMonRas,
                     "year"=as.numeric(yearStr),
                     "month"=as.numeric(monStr),
                     "monthlyWMMbias"=as.data.frame(WMMbiasPixels[,c(1,2,3,ncol(WMMbiasPixels))]),
                     "dailyWMMbias" = as.data.frame(WMMbiasPixels[,-c(ncol(WMMbiasPixels))])
    )
  }
  #-------------------------------------------------
  #  Return a list of 4 spatial objects 
  #     (3 rasters & 1 points)


  #cat(paste(names(as.data.frame(WMMbiasPixels[,-c(ncol(WMMbiasPixels))])),'\n'))
  return(biasStuff)
}

#-------------------------------------------------
# Set up for multiprocessing function calls
#-------------------------------------------------
plan(multiprocess)
processed= listenv(NULL)
yrList=list()

yearStr <- as.character(1985)
#-------------------------------------------------
# Define range of years to process
#-------------------------------------------------
processYears <- seq(1985, 2018)
#processYears <- seq(1985, 1986)
#processYears <- seq(2003, 2004)
#processYears <- seq(2014, 2014)
x=0
yr = processYears[1]
for (yr in processYears) {
  for (mon in seq(1,12)){
    #for (mon in seq(9,12)){
    yearStr <- as.character(yr)
    #monStr <- as.character(mon)
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


filePath <- basePath


for (iyr in processYears){
  WMMbiasMonthly <-rList[[1]]$monthlyWMMbias[,c(1,2,3)]
  WMMbiasDaily <- rList[[1]]$dailyWMMbias[,c(1,2,3)]
  # names(WMMbiasMonthly)<-c( "wmmRoCo" ,"X" ,"Y")
  # names(WMMbiasDaily)<-c( "wmmRoCo" ,"X" ,"Y")
  for (i in seq(from=1, to=x)){
    if (rList[[i]]$year == iyr){
      WMMbiasMonthly <- merge(WMMbiasMonthly, rList[[i]]$monthlyWMMbias[,-c(2,3)], by =c("wmmRoCo"))
      WMMbiasDaily <- merge(WMMbiasDaily, rList[[i]]$dailyWMMbias[,-c(2,3)], by =c("wmmRoCo"))
    }
  }
  WMMbiasMonthly$Annual<-rowSums(WMMbiasMonthly[,-c(1,2,3)],na.rm=TRUE)
  cat(paste("exporting",sprintf("MonthlybiasWMM%04d.csv",iyr),'\n'))
  csvFile <- paste0(filePath, sprintf("MonthlybiasWMM%04d.csv",iyr))
  cat(csvFile)
  cat('\n')
  fwrite(WMMbiasMonthly, csvFile) 
  
  WMMbiasDaily$Annual<-rowSums(WMMbiasDaily[,-c(1,2,3)],na.rm=TRUE)
  cat(paste("exporting",sprintf("DailybiasWMM%04d.csv",iyr),'\n'))
  csvFile <- paste0(filePath, sprintf("DailybiasWMM%04d.csv",iyr))
  cat(csvFile)
  cat('\n')
  fwrite(WMMbiasDaily, csvFile) 
}

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
myTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
divTheme = rasterTheme(region = brewer.pal('Spectral', n = 11))
adjTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
difTheme = rasterTheme(region = brewer.pal('Spectral', n = 11))

#-------------------------------------------------
# Create plot files for each raster type by month
#-------------------------------------------------
plotOneRas <- function(filename, rasPlt, rasPltName, pntsPlt, pltTheme, atVals,clpBnds2){
  myplot=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                     at=atVals,layout=c(1,1),contour=FALSE, margin=F) +
             latticeExtra::layer(sp.polygons(clpBnds2)) +
             latticeExtra::layer(sp.text(coordinates(pntsPlt),txt=pntsPlt$RainGage,pos=1,cex=.5 )) +
             latticeExtra::layer(sp.points(pntsPlt, col = "red"))
  )
  return(myplot)
}

for (i in 1:nlayers(rasStack)){
  
  # Plot Monthly wmmRain
  #  filename=paste(basePath,"/rasterPlotsLowerRange/MonWMM",names(MonWMMStack)[i],".png",sep="")
  filename=paste(basePath,"BiasPlotsWMM/points",names(MonWMMStack)[i],"d2000",".png",sep="")
  rasPlt <- MonWMMStack[[i]]
  rasPltName <-names(MonWMMStack)[i]
  pntsPlt <-pointList[[i]]
  pltTheme <- adjTheme
  atVals <-c(seq(0,28,length=29),30)
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  WMMplot=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                      at=atVals, xlab = NULL, margin=F,
                      layout=c(1,1),contour=FALSE) +
              latticeExtra::layer(sp.polygons(clpBnds2)) +
              latticeExtra::layer(sp.text(coordinates(pntsPlt),txt=as.character(round((pntsPlt$sum_Rainfall-pntsPlt$sum_WMM),2)),
                                          pos=1,cex=.6 )) +
              latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                            cex=round((pntsPlt$sum_Rainfall-pntsPlt$sum_WMM)*.65,2),
                                            col = "black")) +             
              latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                            cex=round(((pntsPlt$sum_Rainfall-pntsPlt$sum_WMM)*-1)*.65,2) ,
                                            col = "red"))
  )
  #  trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
  # print(WMMplot)
  # dev.off()
  
  # Plot Monthly Bias Raster 
  #  filename=paste(basePath,"/rasterPlotsLowerRange/bias",names(rasStack)[i],".png",sep="")
  rasPlt <- rasStack[[i]]
  rasPltName <-names(rasStack)[i]
  pntsPlt <-pointList[[i]]
  min(pntsPlt$adjF)
  max(pntsPlt$adjF)
  pltTheme <- divTheme
  #atVals <-seq(0,4,length=21)
  atVals <-3^((-10:10)/7.8)
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  Biasplot=plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2)
  # trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
  # print(WMMplot)
  # dev.off()
  
  #   Plot Monthly Bias Adjusted wmmRain
  #  filename=paste(basePath,"/rasterPlotsLowerRange/MonWMMwBias",names(MonRasStack)[i],".png",sep="")
  rasPlt <- MonRasStack[[i]]
  rasPltName <-names(MonRasStack)[i]
  pntsPlt <-pointList[[i]]
  pltTheme <- adjTheme
  #atVals <-c(seq(0,15,length=16),30)
  atVals <-c(seq(0,28,length=29),30)
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  WMMBiasplot=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                          at=atVals, xlab = NULL, margin=F,
                          layout=c(1,1),contour=FALSE) +
                  latticeExtra::layer(sp.polygons(clpBnds2)) +
                  latticeExtra::layer(sp.text(coordinates(pntsPlt),
                                              txt=as.character(round(pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias),2)),
                                              pos=1,cex=.6 )) +
                  latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                                cex=round(((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*.65),2),
                                                col = "black")) +             
                  latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                                cex=round((((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*-1)*.65),2),
                                                col = "red"))
  )
  #   Plot Monthly Bias Adjusted wmmRain - WMM
  rasPlt <-MonRasStack[[i]]-MonWMMStack[[i]] 
  rasPltName <-names(MonRasStack)[i]
  pntsPlt <-pointList[[i]]
  pltTheme <- difTheme
  atVals <-c(seq(-20,20,length=28))
  #atVals <-c(-6,-3,-2,-1,-.5,0,.5,1,2,3,6)
  #atVals <-c(-6,seq(-3,3,length=24),6)
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  
  if( sum(is.na(pntsPlt@data$bias)) >0 ){
    pntsPlt@data[is.na(pntsPlt@data$bias),]$bias = 1
  }
  diffRasPlt=( levelplot(rasPlt, par.settings=pltTheme, main=rasPltName, 
                         at=atVals, xlab = NULL, margin=F,
                         layout=c(1,1),contour=FALSE) +
                 latticeExtra::layer(sp.polygons(clpBnds2)) +
                 latticeExtra::layer(sp.text(coordinates(pntsPlt),
                                             txt=as.character(round(pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias),2)),
                                             pos=1,cex=.5 )) +
                 latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                               cex=round(((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*.65),2),
                                               col = "black")) +             
                 latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                               cex=round((((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*-1)*.65),2),
                                               col = "red"))
  )
  
  # trellis.device(device="png", filename=filename, width=2400,height=2400,units="px",res=300)
  # print(WMMBiasplot)
  # dev.off()
  trellis.device(device="png", filename=filename, width=4800,height=2400,units="px",res=300)
  print(Biasplot, split    = c(1,1,4,1),more=TRUE)
  print(WMMplot, split     = c(2,1,4,1),more=TRUE)
  print(WMMBiasplot, split = c(3,1,4,1),more=TRUE)
  print(diffRasPlt, split  = c(4,1,4,1))
  dev.off()
  
  # print(grid.arrange(Biasplot,WMMplot,WMMBiasplot, ncol=3))
  # for (i in seq(1,20)){
  #   dev.off()
  # }
}
