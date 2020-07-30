library(readr)
library(reshape2)
library(tidyr)

library(dplyr)

for (yr in seq(2000,2018)){
  cat(yr,sep='/n')
  ECSM_NRD<- read_csv(paste0("G:/ECSM/Data/ECSM_NRD_",yr,".csv"))
  LongNRD<-melt(ECSM_NRD,id=c("ROWnum","Pixel_id","X","Y"))
  LongNRD$variable = as.character(LongNRD$variable)
  LongWdates<-cbind(LongNRD[LongNRD$variable < "A",c("Pixel_id","X","Y")], LongNRD[LongNRD$variable < "A",] %>% separate(variable, c("year", "month", "day")))
  LongWdates<- LongNRD[LongNRD$variable < "A",] %>% separate(variable, c("year", "month", "day"))
  
  oneYr=dcast(LongWdates,Pixel_id+X+Y~year+month,fun.aggregate = sum)
  
  for (months in seq(1,12)){
    names(oneYr)[months+3]<-paste0("month",months)
  }
  oneYr$Annual <- rowSums(oneYr[,-c(1:3)])
  write.csv(oneYr,paste0('G:/ECSM/Data/NRDrainMonthly',yr,'.csv'))

}