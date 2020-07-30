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
cat(paste('Reading data: ECSM_RainGageV2.csv','\n'))
rainGageData<-read.csv("//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/ECSM/Data/ECSM_RainGageV2.csv",stringsAsFactors = FALSE)
rainGageData$YEAR = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y")
rainGageData$MONTH = format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%m")
rainGageData$DAILY_DATE = as.Date(format(as.Date(rainGageData$DAILY_DATE, format="%m/%d/%Y"),"%Y-%m-%d"))

rainGages<- unique(rainGageData[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
rownames(rainGages) <- NULL
wmm_data <- list()
calcRainStat<-function(yr){
  
    oneYr<-read.csv(paste0('G:/ECSM/Data/wmmRain',yr,'.csv'))
    
    xmin = floor(min(oneYr[c('Xcoord')])-5280)
    xmax = ceiling(max(oneYr[c('Xcoord')])+5280)
    ymin = floor(min(oneYr[c('Ycoord')])-5280 - (15*(5280*2)))
    ymax = ceiling(max(oneYr[c('Ycoord')])+5280)
    rasRows <- (ymax - ymin) /(5280*2)
    rasCols <- (xmax - xmin) /(5280*2)
    SouthMost<-min(oneYr$Ycoord)-(2*5280) 
    oneYr<-rename(oneYr,wmmRoCo=variable)
    oneYr$Annual<-NULL
    oneYr$X<-NULL
    dateList<-names(oneYr)[-c(1,2,3)]
    lnum =0
    wideXtra<-NULL
    for (l in LETTERS[15:1]){
      y= SouthMost-(lnum*5280*2)
      lnum = lnum + 1
      for (cols in  seq(1:rasCols)){
        x= xmin + ((cols-1)*(5280*2))
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
    
    # Exclude rain gages South of the WMM grid cells
    #  closest<- nn2(oneYr[,3:4],rainGages[rainGages$YCOORD > SouthMost,3:4],1)
    index= closest[[1]]
    distance=closest[[2]]
    nearestwWMMRoCo<-as.data.frame(cbind(rainGages$DBKEY,oneYr[index,]$wmmRoCo,distance))
    names(nearestwWMMRoCo)<-c('DBKEY', 'wmmRoCo', 'Distance')
    nearestwWMMRoCo$DBKEY <- as.character(nearestwWMMRoCo$DBKEY)
    nearestwWMMRoCo$Distance <- as.numeric(as.character(nearestwWMMRoCo$Distance))
    
    coordinates(oneYr) =  ~ Xcoord + Ycoord
    proj4string(oneYr) = HARNSP17ft
    oneYr$XHARN <- coordinates(oneYr)[, 1]
    oneYr$YHARN <- coordinates(oneYr)[, 2]
    oneYr <- spTransform(oneYr,HARNSP17ft)
    
    ECSM_WMM<- read_csv(paste0("G:/ECSM/Data/wmmRain",yr,".csv"))
    ECSM_WMM<-rename(ECSM_WMM,ROWnum=X1)
    ECSM_WMM<-rename(ECSM_WMM,wmmRoCo=variable)
    ECSM_WMM$ROWnum = NULL
    ECSM_WMM$Annual = NULL
    
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
    RainStats[RainStats$YCOORD <= SouthMost,]$wmmRain.sum<-RainStats[RainStats$YCOORD <= SouthMost,]$Gage.sum
    # from 51
    rainGages<- unique(rainGageData[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
    
    write.csv(RainStats[!is.na(RainStats$MONTH),],paste0('h:/WMMrainstats',yr,'.csv'))
    
    #
    # Filter RainStats
    #
    RainStats$bias = RainStats$Gage.sum/RainStats$wmmRain.sum
    RainStats[is.nan(RainStats$bias),]$bias = 1.0
    RainStats<-fixDecimals(RainStats,3)
    RainStats<-na.omit(RainStats)
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
    RainStats<-merge(nearestwWMMRoCo, merge(rainGages,RainStats))
    
    Rain<-merge(Rain,RainStats[,c("DBKEY","wmmRoCo","STATION","MONTH","bias")])
    Rain$WMM<-Rain$wmmRain*Rain$bias
    wmm_data[[i]] <- RainStats[,c("STATION","YEAR","MONTH","XCOORD","YCOORD","Gage.obs","Gage.sum","wmmRain.sum","bias")]
    
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
  
}
i=0
yr =2000
cat(paste('Calculating Rain statistics by year','\n'))
#for (yr in seq(2000,2018)){
for (yr in seq(1985,2018)){
#for (yr in seq(2017,2018)){
#for (yr in seq(2000,2000)){
    #for (yr in seq(1985,1986)){
  i = i + 1
  cat(yr,sep='/n')
  
  calcRainStats(yr)
  oneYr<-read.csv(paste0('G:/ECSM/Data/wmmRain',yr,'.csv'))
  
  xmin = floor(min(oneYr[c('Xcoord')])-5280)
  xmax = ceiling(max(oneYr[c('Xcoord')])+5280)
  ymin = floor(min(oneYr[c('Ycoord')])-5280 - (15*(5280*2)))
  ymax = ceiling(max(oneYr[c('Ycoord')])+5280)
  rasRows <- (ymax - ymin) /(5280*2)
  rasCols <- (xmax - xmin) /(5280*2)
  SouthMost<-min(oneYr$Ycoord)-(2*5280) 
  oneYr<-rename(oneYr,wmmRoCo=variable)
  oneYr$Annual<-NULL
  oneYr$X<-NULL
  dateList<-names(oneYr)[-c(1,2,3)]
  lnum =0
  wideXtra<-NULL
  for (l in LETTERS[15:1]){
    y= SouthMost-(lnum*5280*2)
    lnum = lnum + 1
    for (cols in  seq(1:rasCols)){
      x= xmin + ((cols-1)*(5280*2))
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
  
# Exclude rain gages South of the WMM grid cells
#  closest<- nn2(oneYr[,3:4],rainGages[rainGages$YCOORD > SouthMost,3:4],1)
  index= closest[[1]]
  distance=closest[[2]]
  nearestwWMMRoCo<-as.data.frame(cbind(rainGages$DBKEY,oneYr[index,]$wmmRoCo,distance))
  names(nearestwWMMRoCo)<-c('DBKEY', 'wmmRoCo', 'Distance')
  nearestwWMMRoCo$DBKEY <- as.character(nearestwWMMRoCo$DBKEY)
  nearestwWMMRoCo$Distance <- as.numeric(as.character(nearestwWMMRoCo$Distance))
  
  coordinates(oneYr) =  ~ Xcoord + Ycoord
  proj4string(oneYr) = HARNSP17ft
  oneYr$XHARN <- coordinates(oneYr)[, 1]
  oneYr$YHARN <- coordinates(oneYr)[, 2]
  oneYr <- spTransform(oneYr,HARNSP17ft)
  
  ECSM_WMM<- read_csv(paste0("G:/ECSM/Data/wmmRain",yr,".csv"))
  ECSM_WMM<-rename(ECSM_WMM,ROWnum=X1)
  ECSM_WMM<-rename(ECSM_WMM,wmmRoCo=variable)
  ECSM_WMM$ROWnum = NULL
  ECSM_WMM$Annual = NULL

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
  RainStats[RainStats$YCOORD <= SouthMost,]$wmmRain.sum<-RainStats[RainStats$YCOORD <= SouthMost,]$Gage.sum
  # from 51
  rainGages<- unique(rainGageData[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  
  write.csv(RainStats[!is.na(RainStats$MONTH),],paste0('h:/WMMrainstats',yr,'.csv'))
  
  #
  # Filter RainStats
  #
  RainStats$bias = RainStats$Gage.sum/RainStats$wmmRain.sum
  RainStats[is.nan(RainStats$bias),]$bias = 1.0
  RainStats<-fixDecimals(RainStats,3)
  RainStats<-na.omit(RainStats)
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
  RainStats<-merge(nearestwWMMRoCo, merge(rainGages,RainStats))
  
  Rain<-merge(Rain,RainStats[,c("DBKEY","wmmRoCo","STATION","MONTH","bias")])
  Rain$WMM<-Rain$wmmRain*Rain$bias
  wmm_data[[i]] <- RainStats[,c("STATION","YEAR","MONTH","XCOORD","YCOORD","Gage.obs","Gage.sum","wmmRain.sum","bias")]
  
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

# curerntly only one dataframe in list however
allStats<-do.call(rbind,wmm_data)
names(allStats)<- c('RainGage','year','month','XCOORD','YCOORD','count','sum_Rainfall','sum_WMM','adjF')

myTheme = rasterTheme(region = brewer.pal('Blues', n = 9))

basePath <-  "G:/ECSM/Data/"

meltWMM<- na.omit(melt(ECSM_WMM,id=c("wmmRoCo","Xcoord","Ycoord")))
PixelCoords <- meltWMM[c("wmmRoCo", "Xcoord", "Ycoord")]
PixelCoords$Xcoord <- as.numeric(PixelCoords$Xcoord)
PixelCoords$Ycoord <- as.numeric(PixelCoords$Ycoord)
#-------------------------------------------------
# Add additional melted WMM data sets 
# to calculate bias based on more than one year
#-------------------------------------------------
# curerntly only one dataframe in list however
WMM <- do.call("rbind", list(meltWMM))

#-------------------------------------------------
# Calculate raster extents  including an extra 15 rows to the South
# for rain gages covering the Keys  outside the WMM grid
#-------------------------------------------------
 xmin = floor(min(PixelCoords[c('Xcoord')])-5280)
 xmax = ceiling(max(PixelCoords[c('Xcoord')])+5280)
 ymin = floor(min(PixelCoords[c('Ycoord')])-5280 - (15*(5280*2)))
 ymax = ceiling(max(PixelCoords[c('Ycoord')])+5280)
# xmin = 672000,xmax = 985000, ymin = 143000, ymax = 1202000

#-------------------------------------------------
# calculate number of rows and columns
# 6561.679 = feet in 2 Kilometer pixel spacing
#-------------------------------------------------
# rasRows <- (ymax - ymin) / 6561.679
# rasCols <- (xmax - xmin) / 6561.679
rasRows <- floor((ymax - ymin) /(5280*2))
rasCols <- floor((xmax - xmin) /(5280*2))

#-------------------------------------------------
# define raster and map extents using WMM pixel data extents
#-------------------------------------------------
ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)

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

yearStr = '2017'
monStr = '09'
southernRain<-rainGageData[rainGageData$YEAR==yearStr & 
                         rainGageData$MONTH==monStr &
                         rainGageData$YCOORD < SouthMost,]
#Prepare 1 year of RaingGage Data South of WMM for voronoi process
coordinates(southernRain) =  ~ XCOORD + YCOORD
proj4string(southernRain) = HARNSP17ft
southernRain$XHARN <- coordinates(southernRain)[, 1]
southernRain$YHARN <- coordinates(southernRain)[, 2]
southernRain <- spTransform(southernRain,HARNSP17ft)

biasByYearMon <-function(southernRain,yearStr,monStr,x){
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
  WMMrainByYr<- read.csv(paste0("G:/ECSM/Data/wmmRain",yr,".csv"),
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
    WMMBiasPnts <- data.frame(extract(mergeRas, oneYr))
    WMMBiasPnts$wmmRoCo <- oneYr$wmmRoCo
    WMMBiasPnts$X<-oneYr$Xcoord
    WMMBiasPnts$Y<-oneYr$Ycoord    
   
    DailyWMM.grid <- as(mergeRas, "SpatialGridDataFrame")
  }
  WMMPixel<- data.frame(DailyWMM.pnts$wmmRoCo, DailyWMM.pnts$X, DailyWMM.pnts$Y)
  WMMxPixel<- data.frame(oneYr$wmmRoCo, oneYr$XHARN, oneYr$YHARN)
  names(WMMxPixel)<-c("wmmRoCo","X","Y")
  WMMxPixel<- WMMxPixel[WMMxPixel$Y <= SouthMost,]
  mon <- as.numeric(monStr)
  RGdata <- na.omit(allStats[allStats$year == yearStr & 
                       allStats$month == monStr
                             & allStats$adjF <400,])
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
        voronoi(southernRain[southernRain$DAILY_DATE==dayLUp[dayLUp$dateList==d,]$V2,])
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
    # and export to csv
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
    # and export to csv
    #-------------------------------------------------
    col.num <- which(colnames(WMMrainByYr@data) %in% monList)
    WMMrainByYr$Monthly<-rowSums(WMMrainByYr@data[,col.num],na.rm=TRUE)
    WMMrainByYr <- fixDecimals(WMMrainByYr@data,4)
    # csvFile <- paste0(basePath, sprintf("WMM%s%02d.csv",yearStr,mon))
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
    # MonWMMs<-WMMrainByYr[,c(2,3,4,ncol(WMMrainByYr))]
    # WMMMonRas <- rasterize(MonWMMs,ras,MonWMMs$Monthly,fun=mean)
    MonWMMs<-WMMbiasPixels[,c("wmmRoCo",sprintf('Mon%02d',mon))]
    monthlyVals <- as.data.frame(WMMbiasPixels[,sprintf('Mon%02d',mon)])
    
    WMMMonRas <- rasterize(MonWMMrain,ras, MonWMMrain$Monthly)
    #WMMMonRas <- rasterize(MonWMMs,ras, monthlyVals[,c(sprintf('Mon%02d',mon))])
    #-------------------------------------------------
    # Create Raster for bias corrected Monthly wmmRain rain
    #-------------------------------------------------
    Monthly <-WMMbiasPixels[,c(1,ncol(WMMbiasPixels))]
    # Monthly <-WMMbiasPixels[,c(1,2,3,ncol(WMMbiasPixels))]
    # xy <- Monthly[,c(2,3)]
    # Monthlypdf <-SpatialPointsDataFrame(coords=xy,data=Monthly,proj4string=HARNSP17ft)
    # MonRas <- rasterize(Monthlypdf,ras,Monthlypdf@data[,4],fun=mean)    
    MonRas <- rasterize(Monthly,ras,Monthly@data[,2],fun=mean)    
    #-------------------------------------------------
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
  #-------------------------------------------------
  #  Return a list of 4 spatial objects 
  #     (3 rasters & 1 points)

  #cat(paste(names(as.data.frame(WMMbiasPixels[,-c(ncol(WMMbiasPixels))])),'\n'))
  return(biasStuff)
}

#-------------------------------------------------
# Set up for multiprocessing function calls
#-------------------------------------------------
cat(paste('Establishing mutliprocessor Plan','\n'))
plan(multiprocess,.skip=TRUE)
processed= listenv(NULL)
yrList=list()

yearStr <- as.character(2017)
#-------------------------------------------------
# Define range of years to process
#-------------------------------------------------
processYears <- seq(1985, 2018)
#processYears <- seq(2017, 2018)
#processYears <- seq(2000, 2000)

cat(paste('Processing bias correction by year','\n'))
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
    southernRain<-rainGageData[rainGageData$YEAR==yearStr & 
                             rainGageData$MONTH==monStr &
                             rainGageData$YCOORD < SouthMost,]
    #Prepare 1 year of RaingGage Data South of WMM for voronoi process
    coordinates(southernRain) =  ~ XCOORD + YCOORD
    proj4string(southernRain) = HARNSP17ft
    southernRain$XHARN <- coordinates(southernRain)[, 1]
    southernRain$YHARN <- coordinates(southernRain)[, 2]
    southernRain <- spTransform(southernRain,HARNSP17ft)
    
    processed[[x]] <- future({biasByYearMon(southernRain,yearStr,monStr,x)})
    #processed[[x]] <- biasByYearMon(southernRain,yearStr,monStr,x)
  }
}

mpList<-list()
rList<-list()
stackList <-list()

#-------------------------------------------------
# value function waits for results to become available
# for each process
#-------------------------------------------------
cat(paste('Waiting for bias correction processes to return','\n'))
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

basePath <-  "G:/ECSM/Data/test/"
filePath<- basePath
cat(paste('Exporting csv files','\n'))
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
  length(WMMbiasMonthly)
  if(length(WMMbiasMonthly)>4){
    WMMbiasMonthly$Annual<-rowSums(WMMbiasMonthly[,-c(1,2,3)],na.rm=TRUE)
  }
  cat(paste("exporting",sprintf("MonthlybiasXWMM%04d.csv",iyr),'\n'))
  csvFile <- paste0(filePath, sprintf("MonthlybiasXWMM%04d.csv",iyr))
  cat(csvFile)
  cat('\n')
  fwrite(WMMbiasMonthly, csvFile) 
  
  WMMbiasDaily$Annual<-rowSums(WMMbiasDaily[,-c(1,2,3)],na.rm=TRUE)
  cat(paste("exporting",sprintf("DailybiasXWMM%04d.csv",iyr),'\n'))
  csvFile <- paste0(filePath, sprintf("DailybiasXWMM%04d.csv",iyr))
  cat(csvFile)
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
myTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
divTheme = rasterTheme(region = brewer.pal('Spectral', n = 11))
divTheme = rasterTheme(region = brewer.pal('PiYG', n = 11))
adjTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
difTheme = rasterTheme(region = brewer.pal('Spectral', n = 11))
difTheme = rasterTheme(region = brewer.pal('RdBu', n = 11))

#-------------------------------------------------
# Create plot files for each raster type by month
#-------------------------------------------------
plotOneRas <- function(filename, rasPlt, rasPltName, pntsPlt, pltTheme, atVals,clpBnds2){
  panel1 = paste('Bias Multiplier',rasPltName, '\nRed Crosses = Rain Gages')
  
  myplot=( levelplot(rasPlt, par.settings=pltTheme, main=panel1, 
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
  filename=paste(basePath,"BiasPlotsWMM/points",names(MonWMMStack)[i],"X",".png",sep="")
  rasPlt <- MonWMMStack[[i]]
  rasPltName <-names(MonWMMStack)[i]
  pntsPlt <-pointList[[i]]
  pltTheme <- adjTheme
  #atVals <-c(seq(0,28,length=29),30)
  atVals <-c(0,seq(5,20,length=16),30)
  panel2 = paste('Uncorrected ',rasPltName,'\n','Black O = Pixels < Observed')
  
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  WMMplot=( levelplot(rasPlt, par.settings=pltTheme, main=panel2, 
                      at=atVals, xlab = NULL, margin=F,
                      layout=c(1,1),contour=FALSE) +
              latticeExtra::layer(sp.polygons(clpBnds2)) +
              latticeExtra::layer(sp.text(coordinates(pntsPlt),txt=as.character(round((pntsPlt$sum_Rainfall-pntsPlt$sum_WMM),2)),
                                          pos=1,cex=.6 )) +
              latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                            cex=round((pntsPlt$sum_Rainfall-pntsPlt$sum_WMM)*.45,2),
                                            col = "black")) +             
              latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                            cex=round(((pntsPlt$sum_Rainfall-pntsPlt$sum_WMM)*-1)*.45,2) ,
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
  #atVals <-c(seq(0,28,length=29),30)
  atVals <-c(0,seq(5,20,length=16),30)
  #  print(plotOneRas(filename, rasPlt, rasPltName,pntsPlt, pltTheme, atVals,clpBnds2))
  # Note = 'Red circles = pixels > observed'
  panel3 = paste('Adjusted ',rasPltName,'\n','Red O = Pixels > Observed')
  WMMBiasplot=( levelplot(rasPlt, par.settings=pltTheme, main=panel3, 
                          at=atVals, xlab = NULL, margin=F,
                          layout=c(1,1),contour=FALSE) +
                  latticeExtra::layer(sp.polygons(clpBnds2)) +
                  latticeExtra::layer(sp.text(coordinates(pntsPlt),
                                              txt=as.character(round(pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias),2)),
                                              pos=1,cex=.6 )) +
                  latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                                cex=round(((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*.45),2),
                                                col = "black")) +             
                  latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                                cex=round((((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*-1)*.45),2),
                                                col = "red"))
  )
  #   Plot Monthly Bias Adjusted wmmRain - WMM
  rasPlt <-MonRasStack[[i]]-MonWMMStack[[i]] 
  rasPltName <-names(MonRasStack)[i]
  pntsPlt <-pointList[[i]]
  pltTheme <- difTheme
  #atVals <-c(seq(-10,10,length=28))
  #atVals <-c(-6,-3,-2,-1,-.5,0,.5,1,2,3,6)
  atVals <-c(-10,-6,seq(-3,3,length=24),6,10)

  if( sum(is.na(pntsPlt@data$bias)) >0 ){
    pntsPlt@data[is.na(pntsPlt@data$bias),]$bias = 1
  }
  # Note = 'Black circles = pixels < observed'
  # Note = 'Red circles = pixels > observed'
  panel4 = paste('Inches of Change\n',rasPltName)
  diffRasPlt=( levelplot(rasPlt, par.settings=pltTheme, main=panel4, 
                         at=atVals, xlab = NULL, margin=F,
                         layout=c(1,1),contour=FALSE) +
                 latticeExtra::layer(sp.polygons(clpBnds2)) +
                 latticeExtra::layer(sp.text(coordinates(pntsPlt),
                                             txt=as.character(round(pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias),2)),
                                             pos=1,cex=.5 )) +
                 latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                               cex=round(((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*.45),2),
                                               col = "black")) +             
                 latticeExtra::layer(sp.points(pntsPlt,pch=1, 
                                               cex=round((((pntsPlt$sum_Rainfall-(pntsPlt$sum_WMM*pntsPlt$bias))*-1)*.45),2),
                                               col = "red") )
               
  )
  
  trellis.device(device="png", filename=filename, width=4800,height=2400,units="px",res=300)
  print(Biasplot, split    = c(1,1,4,1),more=TRUE)
  print(WMMplot, split     = c(2,1,4,1),more=TRUE)
  print(WMMBiasplot, split = c(3,1,4,1),more=TRUE)
  print(diffRasPlt, split  = c(4,1,4,1))
  dev.off()
  
label=c("",
        'Black circles = pixels < observed',
        'Red circles = pixels > observed')
}
