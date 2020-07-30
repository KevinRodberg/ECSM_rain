rainGageData<-read.csv("//ad.sfwmd.gov/dfsroot/data/wsd/GIS/GISP_2012/DistrictAreaProj/ECSM/Data/ECSM_RainGageV2.csv")
rainGages<- unique(rainGageData[,c("STATION","AGENCY","XCOORD","YCOORD","DBKEY" )])
