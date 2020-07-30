library(readr)
library(reshape2)

# In Linux, get needed environment to run cell_cat:
# source ~cadavid/.cshrc
# Here is the path to the rainfall binary file (not accessible from PC):
#   /nw/hesm_san/data/rain/rain_v4.7_1914_2018_sfwmd.bin

#  split ecsm_roco...txt into 3 files because cell_cat is ilmited to a 1000 cells.
#  Data gets saved to pc accissible H:/ecsmRain...dat
#  Backup copy made to \\ad.sfwmd.gov\dfsroot\data\wsd\sup\devel\source\R\ECSM_rain

# cell_cat -i ~/ecsm_roco1.txt -o ~/ecsmRain1.dat -y 19850101,20181231 rain_v4.7_1914_2018_sfwmd.bin
# cell_cat -i ~/ecsm_roco2.txt -o ~/ecsmRain2.dat -y 19850101,20181231 rain_v4.7_1914_2018_sfwmd.bin
# cell_cat -i ~/ecsm_roco3.txt -o ~/ecsmRain3.dat -y 19850101,20181231 rain_v4.7_1914_2018_sfwmd.bin

ecsmRain1 <- read_table2("H:/ecsmRain1.dat", 
                         col_types = cols(Da = col_integer(), 
                                          Mo = col_integer(), 
                                          Year = col_integer()))
ecsmRain2 <- read_table2("H:/ecsmRain2.dat", 
                         col_types = cols(Da = col_integer(), 
                                          Mo = col_integer(), 
                                          Year = col_integer()))
ecsmRain3 <- read_table2("H:/ecsmRain3.dat", 
                         col_types = cols(Da = col_integer(), 
                                          Mo = col_integer(), 
                                          Year = col_integer()))

ECSM_WMM_rocoCoords <- read_csv("H:/ECSM_WMM_rocoCoords.csv", 
                                col_types = cols(Col = col_integer(), 
                                                 District_1 = col_skip(), FID = col_skip(), 
                                                 Id = col_skip(), LOK_Active = col_skip(), 
                                                 NSM_Active = col_skip(), NSM_Col = col_skip(), 
                                                 NSM_ROCO = col_skip(), OBJECTID = col_skip(), 
                                                 OBJECTID_1 = col_skip(), Row = col_integer(), 
                                                 SFWMM_Acti = col_skip(), SFWMM_Col = col_skip(), 
                                                 SFWMM_ROCO = col_skip(), Shape_Area = col_skip(), 
                                                 Shape_Leng = col_skip(), cell_area = col_skip(), 
                                                 rf_active = col_skip(), rowco = col_integer()))

Rain1<-(melt(ecsmRain1, id=c("Year","Mo","Da")))
Rain2<-(melt(ecsmRain2, id=c("Year","Mo","Da")))
Rain3<-(melt(ecsmRain3, id=c("Year","Mo","Da")))
LongRain<-do.call(rbind, list(Rain1, Rain2, Rain3))
LongRain$variable <- gsub('roco', '', LongRain$variable)
LongRain$variable <- gsub('_', '0', LongRain$variable)
LongRain$variable <- as.integer(LongRain$variable)
RainwithCoords=merge(x = LongRain, y = ECSM_WMM_rocoCoords, 
      by.x = "variable", by.y="rowco", all.x = TRUE)
RainwithCoords=RainwithCoords[c("variable", "Xcoord", "Ycoord",
                                "Year","Mo","Da","value")]
for (yr in unique(RainwithCoords$Year)){
  #oneYr=dcast(RainwithCoords[RainwithCoords$Year == yr,],variable+Xcoord+Ycoord~Year+Mo+Da)
  # for (days in seq(1,length(names(oneYr))-3)){
  #   names(oneYr)[days+3]<-paste0("day",days)
  # }
  # oneYr$Annual <- rowSums(oneYr[,-c(1:3)])
  # write.csv(oneYr,paste0('G:/ECSM/Data/wmmRain',yr,'.csv'))
  oneYr=dcast(RainwithCoords[RainwithCoords$Year == yr,],variable+Xcoord+Ycoord~Year+Mo,fun.aggregate = sum)
  for (months in seq(1,12)){
     names(oneYr)[months+3]<-paste0("month",months)
  }
  oneYr$Annual <- rowSums(oneYr[,-c(1:3)])
  write.csv(oneYr,paste0('G:/ECSM/Data/wmmRainMonthly',yr,'.csv'))
}
