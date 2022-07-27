### Calculate GBIF and OBIS-based STIs ######

args = commandArgs(trailingOnly=TRUE)

startspno<-1000*as.numeric(args[1])

# Batch file version

#library(data.table)
library(dismo)
library(rgbif)
library(robis)
library(raster)

setwd("W:/sa01mb/BIOTIME/")

#Hello - test  commit in r
#

# Only frequent species and Including realm
speciesinci2<-read.csv("speciesinci2.csv")

# Maxent predictors

oceanpreds<-brick("oceanpreds.tif") # 1 is SST
landpreds<-brick("landpreds.tif")   # 1 is LST

## Create bias set - Terrestrial ###

biasras<-brick("GBIFeffort3.tif")
xbias<-biasras[[4]]
xbias[is.na(landpreds[[1]])]<-0
bgt <- xyFromCell(xbias, sample(ncell(xbias), 10000, prob=values(xbias)/255))

## Create bias set - Marine ###

xbiasocean<-resample(biasras[[4]],oceanpreds)
xbiasocean[is.na(oceanpreds[[1]])]<-0
bgo <- xyFromCell(xbiasocean, sample(ncell(xbiasocean), 10000, prob=values(xbiasocean)/255))

### Fit Maxent models, predict distribution and get temperature across predicted distbn ######

df<-data.frame(speciesname=character(),X0.=numeric(),X10.=numeric(),X25.=numeric(),X50.=numeric(),
               X75.=numeric(),X90.=numeric(),X100.=numeric(),nrecs=numeric())
#names(df)<-c("speciesname","X0.","X10.","X25.","X50.","X75.","X90.","X100.")

for (i in startspno+1:1000) {
  speciesname<-with(speciesinci2[i,],paste(GENUS,SPECIES,sep=" "))
  realm<-as.character(speciesinci2$REALM[i])
  print(paste(c(i,speciesname,realm),sep=" "))
  
  gotdata=F
  if (realm=="Marine") {
    gbifdata<-robis::occurrence(scientificname=speciesname)# ,limit=1000) # rgbif
    if (dim(gbifdata)[1]>0) {
      gotdata=T
      gbifxy<-with(gbifdata,cbind(decimalLongitude,decimalLatitude))
      gbifxy<-gbifxy[gbifxy[,1]!=0 & gbifxy[,2]!=0,]}
    } else {  # Terrestrial
      gbifdata<-rgbif::occ_data(scientificName=speciesname,limit=5000) # rgbif
      if (length(gbifdata)>0) {
        gotdata=T
        gbifxy<-with(gbifdata$data,cbind(decimalLongitude,decimalLatitude))
        gbifxy<-gbifxy[gbifxy[,1]!=0 & gbifxy[,2]!=0,]}
    }

  if(gotdata==T) {
    if (realm=="Marine") {
      try(mxmodel<-maxent(oceanpreds,gbifxy,a=bgo,removeDuplicates=T))
      specmap<-predict(mxmodel,oceanpreds)
      specmap[specmap<0.4]<-NA
      #    plot(specmap)
      habtemps<-raster::mask(oceanpreds[[1]],specmap) }
    else {    # Terrestrial
      try(mxmodel<-maxent(landpreds,gbifxy,a=bgt,removeDuplicates=T))
      specmap<-predict(mxmodel,landpreds)
      specmap[specmap<0.4]<-NA
      #    plot(specmap)
      habtemps<-raster::mask(landpreds[[1]],specmap) 
    }
    qtemps<-quantile(habtemps[],probs=c(0,0.1,0.25,0.5,0.75,0.9,1),na.rm=T)
    print(qtemps)
    nrecs<-length(gbifxy)
    df<-rbind(df,data.frame(speciesname,t(qtemps),nrecs)) 
    write.csv(df,paste0("STIvals",startspno,".csv"))
  }
}


