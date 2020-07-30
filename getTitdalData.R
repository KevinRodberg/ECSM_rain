library(utils)
library(dplyr)
NAVD_TideData = data.frame(stn=character(), 
                              Date.Time=as.Date(character()),
                              Prediction=double(), 
                              stringsAsFactors = FALSE)
for (stn in c(8722670, 8721604, 8722125, 8722208, 8722219, 8722213)) {
  for (yr in seq(1985, 2019)) {
    cat(paste("Station:",stn," Year:",yr,'\n'))
    NAVD_TideData <-rbind(NAVD_TideData,cbind(stn,read.csv(paste0("https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=",
          yr,"0101&end_date=",yr,"1231&station=",stn,
          "&product=predictions&datum=NAVD&units=english&time_zone=gmt&application=NOS.COOPS.TAC.WL&format=csv"
          )))
      )
  }
}

wideTides<-dcast(NAVD_TideData,Date.Time~stn)
write.csv(wideTides,'m:/ECSM/TidalData/ExtraTideStations.csv')

# Oops I goofed up NAVD_TideData  this put it back without requerying
# Tides=melt(wideTides)
# NAVD_TideData=Tides
# names(NAVD_TideData)<-c('Date.Time','stn','Prediction')

NAVD_TideData$Date.Time=anytime(Tides$Date.Time)
dailyTides<- NAVD_TideData %>%
  mutate(daily=as.Date(NAVD_TideData$Date.Time, format = "%Y-%m-%d")) %>%
  group_by(stn,daily) %>%
  summarise(meanTide=mean(Prediction))

wideMeanTides<-dcast(dailyTides,daily~stn)
write.csv(wideMeanTides,'m:/ECSM/TidalData/ExtraMeanTideStations.csv')

