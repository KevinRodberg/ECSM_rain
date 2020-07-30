library(utils)
library(dplyr)
library(anytime)
library(reshape2)
NAVD_TideData = data.frame(stn=character(), 
                              Date.Time=as.Date(character()),
                              Prediction=double(), 
                              Type=character(),
                              stringsAsFactors = FALSE)
stn=8722357
# for (stn in c(8722670, 8721604, 8722125, 8722208, 8722219, 8722213)) {
EOM=c('31','28','31','30','31','30','31','31','30','31','30','31')
for (yr in seq(1985, 2019)) {
  for (mon in c('01','02','03','04','05','06','07','08','09','10','11','12')) {
    cat(paste("Station:", stn, " Year:", yr, " Month", mon,'\n'))

    if (mon == '02' & (yr %% 4)== 0){
      eom = '29'
    } else {
      eom =EOM[as.integer(mon)]
    }
    NAVD_TideData <-
      rbind(NAVD_TideData, cbind(stn, 
            read.csv(paste0("https://tidesandcurrents.noaa.gov/api/datagetter",
                            "?begin_date=",yr,mon,"01&end_date=",yr,mon,eom,"&station=",stn,
                            "&product=predictions&datum=MLLW&units=english&time_zone=gmt",
                            "&interval=hilo&application=NOS.COOPS.TAC.WL&format=csv")
                     )
            )
      )
  }
}


# 
NAVD_TideData$Date.Time=anytime(NAVD_TideData$Date.Time)
dailyTides<- NAVD_TideData %>%
  mutate(daily=as.Date(NAVD_TideData$Date.Time, format = "%Y-%m-%d")) %>%
  group_by(stn,daily) %>%
  summarise(meanTide=mean(Prediction)+1.36)

wideMeanTides<-dcast(dailyTides,daily~stn)
write.csv(wideMeanTides,'m:/ECSM/TidalData/MeanTideAt_8722357.csv')

