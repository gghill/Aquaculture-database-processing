# Terrestrial SDM predictors

bio1<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_01.tif")
bio12<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_12.tif")
bio2<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_02.tif")
bio5<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_05.tif")
bio6<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_06.tif")
bio13<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_13.tif")
bio14<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_14.tif")

#bio4<-raster("D:/Documents/Proposals/2019/NERCHighlightCTIGlobal/Analysis/wc2.0_10m_bio/wc2.0_bio_10m_01.tif")

landpreds<-brick(bio1,bio12,bio2)#,bio5,bio6,bio13,bio14)
writeRaster(landpreds,"landpreds.tif")

# Ocean SDM predictors

sst<-raster("D:/Documents/Proposals/2011/VelocityProposalNERCDecember2011/OISSTV2HR/DailySST/OISSTHRMean8218.tif")
minsst<-raster("D:/Documents/Proposals/2011/VelocityProposalNERCDecember2011/OISSTV2HR/DailySST/OISSTHR05pctMean8218.tif")
maxsst<-raster("D:/Documents/Proposals/2011/VelocityProposalNERCDecember2011/OISSTV2HR/DailySST/OISSTHR95pctMean8218.tif")
depth<-raster("D:/Documents/ArcView Themes/ETOPO/etop0pt25dg.tif")

seasrange<-maxsst-minsst

plot(seasrange)

sst1<-focal(sst,w=matrix(1,3,3),fun=mean,na.rm=T)
seasrange1<-focal(seasrange,w=matrix(1,3,3),fun=mean,na.rm=T)
depth1<-depth
depth1[is.na(depth)]<-0

plot(sst1)
plot(seasrange1)
plot(depth1)

oceanpreds<-brick(sst1,seasrange1,depth1)

writeRaster(oceanpreds,"oceanpreds.tif")

i=1

## Create bias set - Terrestrial ###

biasras<-brick("D:\\Documents\\Proposals\\2019\\NERCHighlightCTIGlobal\\Analysis\\GBIFeffort3.tif")
xbias<-biasras[[4]]

xbias[is.na(landpreds[[1]])]<-0
bgt <- xyFromCell(xbias, sample(ncell(xbias), 10000, prob=values(xbias)/255))
write.csv(bgt,"bgt.csv")

## Create bias set - Marine ###

xbiasocean<-resample(biasras[[4]],sst)

xbiasocean[is.na(oceanpreds[[1]])]<-0

plot(xbiasocean)
bgo <- xyFromCell(xbiasocean, sample(ncell(xbiasocean), 10000, prob=values(xbiasocean)/255))

plot(sst1)
points(bgo,pch=3,cex=0.1)
points(gbifxy,pch=21,cex=0.3,col="red")
