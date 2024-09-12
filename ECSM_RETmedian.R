library(readr)
library(reshape2)

# In Linux, get needed environment to run cell_cat:
# source ~cadavid/.cshrc
# Here is the path to the RETfall binary file (not accessible from PC):
#   /nw/hesm_san/data/RET/RET_v4.7_1914_2018_sfwmd.bin

#  split ecsm_roco...txt into 3 files because cell_cat is limited to a 1000 cells.
#  Data gets saved to pc accessible H:/ecsmRET...dat
#  Backup copy made to \\ad.sfwmd.gov\dfsroot\data\wsd\sup\devel\source\R\ECSM_RET

# cell_cat -i ~/ecsm_roco1.txt -o ~/ecsmRET1.dat -y 19850101,20181223 /nw/hesm_nas/data/et/ret_48_16.bin
# cell_cat -i ~/ecsm_roco2.txt -o ~/ecsmRET2.dat -y 19850101,20181231 /nw/hesm_nas/data/et/ret_48_16.bin
# cell_cat -i ~/ecsm_roco3.txt -o ~/ecsmRET3.dat -y 19850101,20181231 /nw/hesm_nas/data/et/ret_48_16.bin

ecsmRET1 <- read_table("H:/ecsmRET1.dat", 
                         col_types = cols(Da = col_integer(), 
                                          Mo = col_integer(), 
                                          Year = col_integer()))
ecsmRET2 <- read_table("H:/ecsmRET2.dat", 
                         col_types = cols(Da = col_integer(), 
                                          Mo = col_integer(), 
                                          Year = col_integer()))
ecsmRET3 <- read_table("H:/ecsmRET3.dat", 
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

RET1<-(melt(ecsmRET1, id=c("Year","Mo","Da")))
RET2<-(melt(ecsmRET2, id=c("Year","Mo","Da")))
RET3<-(melt(ecsmRET3, id=c("Year","Mo","Da")))
LongRET<-do.call(rbind, list(RET1, RET2, RET3))
# reformat rocos that look like (01,42) as 01042
LongRET$variable <- gsub('\\(', '', LongRET$variable)
LongRET$variable <- gsub(',', '0', LongRET$variable)
LongRET$variable <- gsub('\\)', '', LongRET$variable)

LongRET$variable <- as.integer(LongRET$variable)
RETwithCoords=merge(x = LongRET, y = ECSM_WMM_rocoCoords, 
      by.x = "variable", by.y="rowco", all.x = TRUE)
RETwithCoords=RETwithCoords[c("variable", "Xcoord", "Ycoord",
                                "Year","Mo","Da","value")]
allSums=data.frame(matrix(ncol=2,nrow=0))
colnames(allSums) <- c('Year', 'RET_meanAnnual')
x=0
for (yr in unique(RETwithCoords$Year)){
  x=x+1
  oneYr=dcast(RETwithCoords[RETwithCoords$Year == yr,],variable+Xcoord+Ycoord~Year+Mo+Da)
  for (days in seq(1,length(names(oneYr))-3)){
    names(oneYr)[days+3]<-paste0("day",days)
  }
  oneYr$Annual <- rowSums(oneYr[,-c(1:3)])
  allSums <- rbind(allSums,cbind(Year=yr,RET_meanAnnual=mean(oneYr[oneYr$Annual>0,]$Annual) ))
#  allSums[x,]$average=  mean(oneYr$Annual)
#  write.csv(oneYr,paste0('H:/wmmDailyRET',yr,'.csv'))
  # oneYr=dcast(RETwithCoords[RETwithCoords$Year == yr,],variable+Xcoord+Ycoord~Year+Mo,fun.aggregate = sum)
  # for (months in seq(1,12)){
  #    names(oneYr)[months+3]<-paste0("month",months)
  # }
  # oneYr$Annual <- rowSums(oneYr[,-c(1:3)])
  # write.csv(oneYr,paste0('H:/wmmRETMonthly',yr,'.csv'),row.names = F)
}
#allSums$Year[which(allSums$RET_meanAnnual==round(median(allSums$RET_meanAnnual)))]
allSums[which(round(allSums$RET_meanAnnual,1)==round(median(allSums$RET_meanAnnual),1)),]

