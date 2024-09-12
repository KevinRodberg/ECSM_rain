#===============================================================================
#  Program: PlotDailyRasByMonth.R 
#           \\ad.sfwmd.gov\dfsroot\data\wsd\sup\devel\source\R\ECSM_rain\
#===============================================================================
# Code History:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Original:       PlotDailyRasByMonth.R            Kevin A. Rodberg - 10/02/2020
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# >>> Execution with Multiprocessors significantly reduces execution time. 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#  General Description:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# -Daily pixels values are read from  Annual csv's before or after bias correction
#     and filtered for a specified month
# -A single Monthly Raster plot of data is produced showing:
#     Single Monthly Map with Krigged Rain gage map of monthly data 
#     with up to 31 panels of daily WMM.  Monthly map highlights rain differences
# -Indiviual DayByDay maps are also produced for more detailed view.
#===============================================================================

yrSeq<-seq(1985,2000)
monthList <-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

#---------------------------------------------------
# Package management
#---------------------------------------------------
list.of.packages <-  c("reshape2","readr","dplyr","tidyr", "data.table",
                       "readxl",
			#	"rgeos",
				"sp", "dismo", "lattice","rasterVis",
                       "maptools","raster","fields","automap", "gstat",
                       "future","listenv","ggplot2","RANN","geosphere",
                       "tcltk2", "RColorBrewer")
pkgChecker <- function(x){
  for( i in x ){    if( ! require( i , character.only = TRUE ) ){
    install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )    }  }
  }
suppressWarnings(pkgChecker(list.of.packages))

#---------------------------------------------------
# Misc function definitions:
#---------------------------------------------------
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")

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

#---------------------------------------------------
# Define GIS variables:
# NAD83 HARN StatePlane Florida East FIPS 0901 Feet
# and WGS84 Lat/Long Spatial Reference variables
#---------------------------------------------------
HARNSP17ft  = sp::CRS("+init=epsg:2881")
latlongs = sp::CRS("+proj=longlat +datum=WGS84")
WMDbnd.Path <- "//ad.sfwmd.gov/dfsroot/data/hpcc_shared/krodberg/NexRadTS"
WMDbnd.Shape <- "CntyBnds.shp"
setwd(WMDbnd.Path)
WMDbnd <- readShapePoly(WMDbnd.Shape, proj4string = HARNSP17ft)
setwd("H:/")
#---------------------------------------------------
# Define Primary data processing functions
#---------------------------------------------------
readPoints <- function(NRD,year){
  #-------------------------------------------------
  #  Read NexRad Pixels coordinates into  data frame 
  #  [minor variations in original vs bias csv file formats]
  #-------------------------------------------------
  date_range <-as.character(seq(as.Date(paste(year,'01','01',sep='-')),as.Date(paste(year,'12','31',sep='-')) ,by=1))
  if(names(NRD)[1]=="X"){
    NRDsub <- NRD[!NRD[,length(NRD)]==0,-c(1)]
  }   else      {  
    NRDsub <- NRD[!NRD[,length(NRD)]==0,]
  }
  names(NRDsub)<- c("Pixel_id", "X", "Y",date_range,"annual" )
  coordinates(NRDsub)<- ~X+Y
  proj4string(NRDsub)<-HARNSP17ft
  return(NRDsub)
}
calcRainStats<-function(basePath,yr,wideXtra,RG){
  # FUNCTION: calcRainStats
  #   [Works well with future function for multiprocessing]
  # Calculate Series of Monthly Rain Stats for 1 year from gage data
  #---------------------------------------------------
  rainGages<- unique(RG[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rownames(rainGages) <- NULL
  #ECSM_WMM<- read_csv(paste0("G:/ECSM/Data/wmmRain",yr,".csv"))
  ECSM_WMM<- read_csv(paste0(basePath,"/wmmRain",yr,".csv"))
  ECSM_WMM<-rename(ECSM_WMM,ROWnum=X1)
  ECSM_WMM<-rename(ECSM_WMM,wmmRoCo=variable)
  ECSM_WMM$ROWnum = NULL
  ECSM_WMM$Annual = NULL
  cols=names(ECSM_WMM)
  
  # Drop column day 366 from non-Leap years
  col.num=0
  if( yr/4 - (yr%/%4) != 0) {
    col.num <- which(colnames(wideXtra) %!in% colnames(ECSM_WMM))
    wideXtra[,col.num]<- NULL
  }
  names(wideXtra)<-cols
  ECSM_WMM<- rbind(ECSM_WMM,wideXtra)
  ECSMdata <- read.csv(csvFile)
  ECSM_WMM<-ECSM_WMM
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
  LongWdates$wmmRoCo <-as.factor(LongWdates$wmmRoCo)
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
  if (dim(RainStats[RainStats$bias>4.,])[1] >0){
    RainStats[RainStats$bias>4.,]$bias <- 4.0
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
}

radioChoose <- function(radioList) {
  win2 <- tktoplevel()
  frame1 <-  tk2frame(win2,borderwidth = 5,relief = "sunken",padding = 10)
  frame3 <-  tk2frame(win2,borderwidth = 3,relief = "sunken",padding = 10)
  lbl.MonthSelect <- tk2label(win2, text = "Month to Create \n   Daily Rasters: ",font = fontHeading)
  tkpack(lbl.MonthSelect,  side = "top",  expand = TRUE,  ipadx = 5,  ipady = 5,  fill = "x")
  tkpack(frame3,side = "bottom",expand = FALSE,fill = "both")
  tkpack(frame1,side = "left",expand = FALSE,fill = "both")
  rBtnVal1 = tclVar(trimws(radioList[[1]]))
  
  # create matrix arrangement of radioList radioButtons
  btns.f1 = vector()
  ncols = 6
  for (num in seq(1, length(radioList))) {
    r= floor((num-1)/ncols)
    c= num-(ncols*r)-1
    btn <- tk2radiobutton(frame1, width=7)
    tkconfigure(btn, variable = rBtnVal1, value = num)
    tkgrid(tk2label(frame1, text = trimws(radioList[[num]])),btn,row=r,column=c,padx=2,pady=15)
    btns.f1 = append(btns.f1, btn)
  }
  
  tkgrid(tk2button(frame3,text ="Cancel",width = -6,command = fnCncl),
         tk2button(frame3,text ="OK",    width = -6,command = fnOK  ), padx = 10,pady = c(5, 15))
  tkbind(win2, "<Return>", fnOK)
  
  tkraise(win2)
  tkwait.variable(done)
  tkdestroy(win2)
  if (tclvalue(done) != 1) {
    exit("User canceled Month Selection")
  }
  n1  <- as.integer((tclvalue(rBtnVal1)))
  return(list(n1=n1))
  
}

#pnts<-MonthlyRain.pnts

plt1ras <- function(rainRas, rainRasName, Plot.pnts, rainTheme, 
                    atVals,clpBnds2,circles=FALSE){
    #---------------------------------------
    # Create plot panel for a single raster 
    #---------------------------------------  
    Diff.pnts<- Plot.pnts[abs(Plot.pnts$sum_Rain-Plot.pnts$sum_WMM)>1,]
    panel1 = paste0(rainRasName,' Rain')
    if (circles){
        myplot=  levelplot(rainRas, par.settings=rainTheme,
                  main=panel1,
                  colorkey=list(space="left"),
                  scales = list(x=list(rot=45),y=list(rot=45),cex=.5),
                  at=atVals,layout=c(1,1),contour=FALSE, margin=F,
                  panel=function(...){
          panel.levelplot.raster(...)
          sp.polygons(clpBnds2)
          sp.text(coordinates(Diff.pnts),txt=Diff.pnts$RainGage,pos=1,cex=.5 )
          sp.points(Diff.pnts,pch=1,
                    cex=round((Diff.pnts$sum_Rain - Diff.pnts$sum_WMM) / -2,2),col = "black")
          sp.points(Diff.pnts,pch=1, 
                    cex=round((Diff.pnts$sum_Rain-Diff.pnts$sum_WMM)/2,2),col = "red")
        })
      }
    else    {
      myplot= levelplot(rainRas, par.settings=rainTheme, main=panel1, 
                         colorkey=list(space="left"),
                         scales = list(x=list(rot=45),y=list(rot=45),cex=.5),
                         at=atVals,layout=c(1,1),contour=FALSE, margin=F,
                        panel=function(...){
                          panel.levelplot.raster(...)
                 sp.polygons(clpBnds2)   
                 sp.points(Plot.pnts, pch=20,cex=.25, col="blue")
                 sp.points(Diff.pnts,pch=1,
                           cex=round((Diff.pnts$sum_Rain-Diff.pnts$sum_WMM)/-2,2),col = "black") 
                 sp.points(Diff.pnts,pch=1,
                           cex=round((Diff.pnts$sum_Rain-Diff.pnts$sum_WMM)/2,2),col = "red")

      })
    }
   
     return(myplot)
  }

########################
##    PROGRAM BODY    ##
########################
GISPath <- "//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/"
ProjPath <- "ECSM/Data/"
basePath <- paste0(GISPath,ProjPath)
cat(paste("Attention to Pop-up Windows:\n",
          "Select input directory for daily Rain data files:\n"))
basePath <-  tk_choose.dir(default = basePath, 
                           caption = "Select input directory for daily Rain data files")
cat(gsub("/",'\\',basePath,fixed=TRUE),sep='\n')

#-------------------------------------------------
#  Select input and Read one year for Pixels
#  using a background process
#-------------------------------------------------
cat("Select daily rain csvfile\n")
csvFile <- utils::choose.files(default=paste(basePath,'*.csv',sep='/'))
cat(gsub("/",'\\',csvFile,fixed=TRUE),sep='\n')

cat("Select Month using Radio BUttons\n")
options <- radioChoose(monthList)
monthNum=options$n1
month <- monthList[[monthNum]]
cat(paste('\t',month,'\n'))

cat("Select output directory for daily Rain plots\n")
outPath <- tk_choose.dir(default = basePath, 
                         caption = "Select output directory for daily Rain plots")
cat(gsub("/",'\\',outPath,fixed=TRUE),sep='\n')

NRD <- utils::read.csv(csvFile)
if (names(NRD)[1]=='wmmRoCo'){
  NRD<-rename(NRD,variable=wmmRoCo)
  NRD<-rename(NRD,Xcoord=X)
  NRD<-rename(NRD,Ycoord=Y)
}
#---------------------------------------------------
# Extract year for the name of the csv file
#---------------------------------------------------
end=gregexpr(pattern='.csv',csvFile)[[1]]-1
start = end - 3
yr=substr(csvFile,start,end)


cat(paste('Establishing mutliprocessor Plan','\n'))
# if error presents regarding invalid connection try rerunning this:
#   plan(multiprocess)

plan(multiprocess,.skip=TRUE)
calcRainMP <- listenv()

#---------------------------------------------------
# Convert data read from csv into points
#---------------------------------------------------
f <- future({readPoints(NRD,yr)})

yr <- as.numeric(yr)
#---------------------------------------------------
#  Read Rain Gage Data if not already in memory
#---------------------------------------------------
if(exists("gageR.df") && object.size(gageR.df) > 350000000 ){
  cat(paste('Skipped Reading data: ECSM_RainGageV2.csv','\n'))  
} else 
{
  cat(paste('Reading Rain Gage data: ECSM_RainGageV2.csv','\n', 'Approximately 4.4 million records may take 2-3 minutes. \n'))
  gageR.df<-read.csv(paste0(basePath,"/ECSM_RainGageV2.csv"),stringsAsFactors = FALSE)
  altPath = '//ad.sfwmd.gov/dfsroot/data/wsd/MOD/YA/Kevin/ECSM_RainGageV2_UP102720_CSV.csv'
  gageR.df<-read.csv(paste0(basePath,"/ECSM_RainGageV2.csv"),stringsAsFactors = FALSE)
  gageR2.df<-read.csv(altPath,stringsAsFactors = FALSE,header=FALSE)
  
  
  names(gageR2.df)<-names(gageR.df)
  gageR2.df$YEAR = format(as.Date(gageR2.df$DAILY_DATE, format="%m/%d/%Y"),"%Y")
  gageR2.df$MONTH = format(as.Date(gageR2.df$DAILY_DATE, format="%m/%d/%Y"),"%m")
  gageR2.df$DAILY_DATE = as.Date(format(as.Date(gageR2.df$DAILY_DATE, 
                                                format="%m/%d/%Y"),"%Y-%m-%d"))
  gageR2.df[gageR2.df$VALUE=='NULL',]$VALUE<-NA
  gageR2.df$VALUE<- as.numeric(gageR2.df$VALUE)
  gageR2.df[gageR2.df$CODE=='NULL',]$CODE<-''

  gageR.df$YEAR = format(as.Date(gageR.df$DAILY_DATE, format="%m/%d/%Y"),"%Y")
  gageR.df$MONTH = format(as.Date(gageR.df$DAILY_DATE, format="%m/%d/%Y"),"%m")
  gageR.df$DAILY_DATE = as.Date(format(as.Date(gageR.df$DAILY_DATE, 
                                               format="%m/%d/%Y"),"%Y-%m-%d"))
  rainGages<- 
    unique(gageR.df[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
  rownames(rainGages) <- NULL
  rainGages<-na.omit(rainGages)
}

#---------------------------------------------------
# Retrieve points from futures promise
#---------------------------------------------------
PixelCoords <-future::value(f)

#---------------------------------------------------------
# Define raster mapping extents 
#   and
# Create dataSet: "oneYr" as framework for ECSM pixels
# based upon WMM cell spacing plus rain data South of WWM
#---------------------------------------------------------
cat(paste('Creating Rain pixel structure','\n'))

oneYr<-NRD
i=0

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
if (nrow(oneYr)<3060){
  ymin = floor(min(oneYr[c('Ycoord')])-halfPixel - (15*(pixSz)))
} else {
  ymin = floor(min(oneYr[c('Ycoord')])-halfPixel )
}
ymax = ceiling(max(oneYr[c('Ycoord')])+halfPixel)

rasRows <- floor((ymax - ymin) /(pixSz))
rasCols <- floor((xmax - xmin) /(pixSz))

ras <- raster(nrow=rasRows,ncol=rasCols,xmn=xmin,xmx=xmax,
              ymn=ymin,ymx=ymax,crs=HARNSP17ft)
rasExt <- extent(ras)
clpBnds2 <- gClip(WMDbnd, rasExt)

oneYr<-rename(oneYr,wmmRoCo=variable)
oneYr$Annual<-NULL
oneYr$X<-NULL
dateList<-names(oneYr)[-c(1,2,3)]

#---------------------------------------------
#  If oneYr is smaller than 3060 the source
#  data needs additional rows to cover the Keys
#---------------------------------------------
if (nrow(oneYr)<3060){
  SouthMost<-min(oneYr$Ycoord)-(pixSz) 
  
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
}
oneYr$Xcoord <- as.numeric(oneYr$Xcoord)
oneYr$Ycoord <- as.numeric(oneYr$Ycoord)
rainGages<- na.omit(rainGages)
closest<- nn2(oneYr[,2:3],rainGages[,3:4],1)
index= closest[[1]]
dist=closest[[2]]

nnRowCo<-as.data.frame(cbind(rainGages$DBKEY,
                             as.character(oneYr[index,]$wmmRoCo),dist))
names(nnRowCo)<-c('DBKEY', 'wmmRoCo', 'dist')
nnRowCo$DBKEY <- as.character(nnRowCo$DBKEY)
nnRowCo$dist <- as.numeric(as.character(nnRowCo$dist))

coordinates(oneYr) =  ~ Xcoord + Ycoord
proj4string(oneYr) = HARNSP17ft
oneYr$XHARN <- coordinates(oneYr)[, 1]
oneYr$YHARN <- coordinates(oneYr)[, 2]
oneYr <- spTransform(oneYr,HARNSP17ft)

#---------------------------------------------
# Process calcRainStats:
#---------------------------------------------
cat(paste('Calculating Rain statistics by year for',yr,'\n'))
codeFilter=c('X','M','N','PT','?')
RG<-gageR.df[gageR.df$YEAR==yr &
                 gageR.df$CODE %!in% codeFilter &
                 gageR.df$VALUE >= 0.0,]
iy=1
calcRainMP[[iy]] <- future({calcRainStats(basePath,yr,wideXtra,RG)})
nyr = iy
rainStatData<-list()
# Retrieve values from futures promise as a list of 2 dataframes
for (i in seq(1:nyr)){
  rainStatData[[i]] <-future::value(calcRainMP[[i]])
}

#-------------------------------------------------------
# unlist results returned from FUNCTION "calcRainStats"
#-------------------------------------------------------
Rain <-rainStatData[[1]]$Rain
RainStats <- rainStatData[[1]]$RainStats
if (nyr > 1){
  for (iy in seq(2,nyr)){
    Rain <-rbind(Rain,rainStatData[[iy]]$Rain)
    RainStats <- rbind(RainStats,rainStatData[[iy]]$RainStats)
  }
}

#-------------------------------------------------------
# Create MonthlyRain.pnts from RainStats
#-------------------------------------------------------
allStats<-RainStats 
names(allStats)<- c('RainGage','year','month','XCOORD','YCOORD',
                    'count','sum_Rain','sum_WMM','adjF')

MonthlyRain.pnts <-SpatialPointsDataFrame(coords = 
                                            allStats[as.numeric(allStats$month)==monthNum, 
                                                     c("XCOORD", "YCOORD")],
                                          data = allStats[as.numeric(allStats$month)==monthNum,],
                                          proj4string = HARNSP17ft)
ras<-setValues(ras,0)
mergeRas<-calc(ras,fun=function(x) ifelse(is.na(x[1]),ifelse(is.na(x[2]),0,x[2]),x[1]))
DailyNRD.grid <- as(mergeRas, "SpatialGridDataFrame")
#-------------------------------------------------
#   Krige
#-------------------------------------------------
surf <- autoKrige(formula=sum_Rain ~ 1,input_data=MonthlyRain.pnts,new_data = DailyNRD.grid)
rainRas <- raster(surf$krige_output)
#-------------------------------------------------
#  IDW raster code:
#-------------------------------------------------
gs <- gstat(formula=sum_Rain~1, locations=MonthlyRain.pnts)
idwRas <- interpolate(mergeRas, gs)
#-------------------------------------------------
# Natural Neighbor or Voronoi/Theissen Polygons
#-------------------------------------------------
theisPoly <-   voronoi(MonthlyRain.pnts)
TheisRas <- rasterize(theisPoly, ras, theisPoly$sum_Rain, fun = mean)

#---------------------------------------------
# define list of dates for year being processed
#---------------------------------------------
date_range <-seq(as.Date(paste(yr,'01','01',sep='-')),
                 as.Date(paste(yr,'12','31',sep='-')) ,by=1)
all_dates<-as.data.frame(date_range)
names(all_dates)<-c('dates')
days<-names(oneYr)
oneyrIndex<-cbind(all_dates,format(all_dates$dates,format='%b'),days[-c(1,length(days)-1,length(days))])
names(oneyrIndex)<-c('dates','month','days')
monthOfDays <- oneyrIndex[oneyrIndex$month==month,]$days

#---------------------------------------------
# Create points from monthly Rain Gage summaries
#---------------------------------------------
coordinates(allStats) =  ~ XCOORD + YCOORD
proj4string(allStats) = HARNSP17ft
i=1
dailyPlt=list()

#---------------------------------------------
# Define color ramp and breaks
#---------------------------------------------
top=max(ceiling(MonthlyRain.pnts$sum_WMM))
mid = mean(ceiling(MonthlyRain.pnts$sum_WMM))
atVals <-c(seq(0,top,length=25))

#---------------------------------------------
# Create plot of Monthly raster
#---------------------------------------------
#divTheme = rasterTheme(region = brewer.pal('PiYG', n = 11))
rainTheme = rasterTheme(region = brewer.pal('GnBu', n = 9))
#dailyPlt[[1]]<- plt1ras(idwRas, 'Monthly IDW', MonthlyRain.pnts, rainTheme, atVals,clpBnds2,circles=TRUE)
dailyPlt[[1]]<- plt1ras(rainRas, 'Monthly Krige', MonthlyRain.pnts, rainTheme, atVals,clpBnds2,circles=TRUE)
#plt1ras <- function(rainRas, rainRasName, pnts, rainTheme, atVals,clpBnds2,circles=FALSE){

#dailyPlt[[1]] <- plt1ras(TheisRas, 'Monthly Theis', MonthlyRain.pnts, rainTheme, atVals,clpBnds2,circles=TRUE)

#---------------------------------------------
# Daily color ramp is defined differently 
# with range of values providing better 
# visibility to daily rain
#---------------------------------------------
atVals <-c(seq(0,mid,length=25),seq(mid+1,top,length=2))

#---------------------------------------------
# Produce daily plots
#---------------------------------------------
i=0
for (m in monthOfDays) {
  cat(m,sep='\n')
  i=i+1
  oneYr@data[,c(m)]<-as.numeric(oneYr@data[,c(m)])
  WMMras <-rasterize(oneYr, ras, m,fun=mean) 
  rainRasName =  oneyrIndex[oneyrIndex$days==m,]$dates
  dailyPlt[[i]]<-plt1ras (WMMras, paste(rainRasName, 'WMM'), MonthlyRain.pnts, rainTheme, atVals,clpBnds2)
}


par.main.text=trellis.par.get('par.main.text')
par.main.text$just = "left"
par.main.text$cex = .5
par.main.text$font = 1
par.main.text$x = grid::unit(.4, "in")
trellis.par.set('par.main.text',par.main.text)

# ------------------------------------------  
# combine the plot panels to print on 1 page
# ------------------------------------------   
fileName = paste(outPath,"/dailyWMMrain_",yr,month,".png",sep="")

trellis.device(device="png", filename=fileName, 
               width=7500,height=6000,units="px",res=300)

ncols = 8

for (x in seq(1, i-1)){
  r= floor((x-1)/ncols)
  c= x-(ncols*r)
  #cat(paste(r+1,' ',c,'\n'))
  print(dailyPlt[[x]], split= c(c,r+1,8,4),more=TRUE)
}
print(dailyPlt[[x+1]], split= c(c+1,r+1,8,4),more=FALSE)

dev.off()

cat("PNG file is located:\n")
cat(gsub("/",'\\',fileName,fixed=TRUE),sep='\n')

cat("Red circles indicate under-estimate of Rain at Pixel compared to reported Rain Gage\n")

suppressWarnings(dir.create(paste0(outPath,"/DayByDay"), recursive=TRUE))

#---------------------------------------------
# Export individual maps
#---------------------------------------------
for (x in seq(2, i)){
  fileName = paste(outPath,"/DayByDay/dailyWMMrain_",yr,month,"-",x-1,".png",sep="")
  trellis.device(device="png", filename=fileName,width=1600,height=2000,units="px",res=300)
  print(dailyPlt[[x]])
  dev.off()
}

fileName = paste(outPath,"/DayByDay/WMMrain_",yr,month,".png",sep="")
trellis.device(device="png", filename=fileName,width=1600,height=2000,units="px",res=300)
print(dailyPlt[[1]])
dev.off()
gageRaindf<-dcast(Rain[as.numeric(Rain$MONTH)==monthNum,],
                  wmmRoCo+DBKEY+STATION+AGENCY+XCOORD+YCOORD ~ DAILY_DATE,
                  value.var="Gage",mean)
names(gageRaindf)<-c(names(gageRaindf[1:6]),c(paste0(month,seq(2,i)-1)))
write.csv(gageRaindf,na='',row.names=F,
          file=paste(outPath,"/DayByDay/GageRain_",yr,month,".csv",sep=""))

