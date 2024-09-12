library(readr)
filePath = "\\\\ad.sfwmd.gov\\dfsroot\\data\\wsd\\GIS\\GISP_2012\\DistrictAreaProj\\ECSM\\Data\\BiasPlotsWMMb.99-1.01\\"
prefix= "DailybiasXWMM"
suffix = ".csv"

allSums=data.frame(matrix(ncol=2,nrow=0))
colnames(allSums) <- c('Year', 'Rain_meanAnnual')
x=0
for (yr in seq (1985,2018)){
  x=x+1
  filename<- paste0(filePath,prefix,yr,suffix)
  cat(filename,sep='\n')
  tempDF <-read_csv(filename)
  allSums <- rbind(allSums,cbind(Year=yr,Rain_meanAnnual=mean(tempDF$Annual[tempDF$Annual> 0]) ))
}
median(allSums$Rain_meanAnnual)  
allSums[which(round(allSums$Rain_meanAnnual)==round(median(allSums$Rain_meanAnnual))),]
