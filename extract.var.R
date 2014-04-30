
env<-read.csv("environment.csv")

head(environment)

library(raster)
library(maptools)
library(tree)

BRmap<-shapefile("~/Documents1/R/3.Atlantic_Forest2014/Shapefiles/BR_Contorno.shp")

BRmap<-shapefile("~/R/3.Atlantic_Forest2014/Shapefiles/BR_Contorno.shp")

plot(BRmap,xlim=c(-61.3,-61),ylim=c(1.5,1.4))

plot(BRmap)

points(environment$LONG,environment$LAT,col=environment$Location)

points(environment$LONG,environment$LAT,col=environment$Location)

plot(environment[,"LONG"],environment[,"LAT"],col=environment$Location)

head(environment)

plot(tapply(environment$LONG,environment$Location,mean,na.rm=TRUE),tapply(environment$LAT,environment$Location,mean,na.rm=TRUE))

text(tapply(environment$LONG,environment$Location,mean,na.rm=TRUE),tapply(environment$LAT,environment$Location,mean,na.rm=TRUE),levels(environment$Location),cex=.8)


isoptera<-read.csv("Isoptera.csv")

attach(isoptera)

loc<-factor(paste(Location,Grid,Trail,Plot),levels=paste(env$Location,env$Grid,env$Trail,env$Plot))

ver<-tapply(Frequency,list(loc,Taxon),sum,na.rm=T)

detach(isoptera)


rich<-rowSums(ver>0,na.rm=T)

environment$LAT

data.frame(environment$Location,environment$Plot,rownames(ver))

summary(lm(rich~environment$LAT))
summary(lm(rich~environment$LONG))


summary(lm(rich~environment$Clay))

plot(rich~environment$Clay,col=env$Location)

termite.grid<-aggregate(termite.plot,list(loc.split$env.Location),sum)

pca.grid<-prcomp(vegdist(termite.grid[,-1],"jaccard"))

plot(pca.grid$x,col=termite.grid[,1])
text(pca.grid$x[,1],pca.grid$x[,2],termite.grid[,1])

library(vegan)

max(termite.plot)

ver<-termite.plot[rich>4,]
environment<-environment[rich>4,]

dist<-dist(t(t(environment[,c("LAT","LONG","Altitude")])*c(1,1,.000)))

jac<-vegdist(ver,"jaccard")
bray<-vegdist(ifelse(is.na(ver),0,ver),"bray")
morisita<-vegdist(ifelse(is.na(ver),0,ver),"morisita")

plot(1-jac~dist,pch=21,col=0,bg=adjustcolor(2,alpha=.1))
plot(1-bray~dist,pch=21,col=0,bg=adjustcolor(2,alpha=.1),cex.lab=1.5,xlim=c(.4,14))
plot(1/(1-morisita)~dist,pch=21,col=0,bg=adjustcolor(2,alpha=.1))

lines(smooth.spline(as.numeric(1-bray[!is.na(dist)])~as.numeric(dist[!is.na(dist)])))

library(mgcv)
model1<-gam(as.numeric(1-bray[!is.na(dist)])~s(ver),data=list(ver=as.numeric(dist[!is.na(dist)])))

lines(seq(0,14,length=100),predict(model1,newdata=list(ver=seq(0,14,length=100))))


dist.clay<-dist(environment$Clay)
dist.P<-dist(environment$P)
dist.PH<-dist(environment$PH_H2O)

mantel(dist,jac,na.rm=TRUE)
mantel(dist.clay,jac,na.rm=TRUE)
mantel(dist.P,jac,na.rm=TRUE)
mantel(dist.PH,jac,na.rm=TRUE)

mantel.partial(jac,dist,dist.clay,na.rm=TRUE)

mantel(dist,bray,na.rm=TRUE)
mantel(dist.clay,bray,na.rm=TRUE)
mantel(dist.P,bray,na.rm=TRUE)


lines(smooth.spline(jac~dist))

plot(scores(cmdscale(jac))[,1]~environment$LAT,col=environment$Location)
plot(scores(cmdscale(jac))[,2]~environment$LAT,col=environment$Location)

plot(scores(cmdscale(bray))[,1]~environment$LAT,col=environment$Location)
plot(scores(cmdscale(bray))[,2]~environment$LAT,col=environment$Location)


plot(scores(cmdscale(jac))[,],col=environment$Location)
ordihull(scores(cmdscale(jac))[,],environment$Location)

plot(scores(cmdscale(jac))[,],col=environment$Location)
ordihull(scores(cmdscale(jac))[,],environment$Location)


NMDSjac<-metaMDS(jac,k=4)

plot(scores(NMDSjac)[,],col=environment$Location)
ordihull(scores(NMDSjac)[,],environment$Location)

NMDSd<-metaMDS(dist2.PDB,k=4)

plot(scores(NMDSd)[,],col=environment$Location)
ordihull(scores(NMDSd)[,],environment$Location)


plot(1-bray~dist(environment[,c("LAT","LONG")]))

lines(smooth.spline(1-bray~dist(environment[,c("LAT","LONG")])))

plot(scores(cmdscale(bray))[,2]~environment$LAT,col=environment$Location)

plot(scores(cmdscale(bray))[,],col=environment$Location)
ordihull(cmdscale(bray),environment$Location)



NMDSbray<-metaMDS(bray,k=4)

plot(scores(NMDSbray)[,],col=environment$Location)
ordihull(scores(NMDSbray)[,],environment$Location)


regtree<-tree(scores(NMDSbray)[,1]~LONG+LAT+P+Clay+Altitude,data=environment)

summary(regtree)


plot(regtree)
text(regtree)

plot(prune.tree(regtree))


plot(scores(NMDSbray)[,1]~environment$LAT)

tobepred<-data.frame(LAT=seq(-10,10,length=100),P=seq(0,0,length=100),Clay=seq(0,0,length=100),Altitude=seq(0,0,length=100))

lines(tobepred$LAT,predict(regtree,newdata=tobepred))





library(rgl)

plot(scores(NMDSbray)[,],col=environment$Location)
ordihull(scores(NMDSbray)[,],environment$Location)

plot3d(scores(NMDSjac)[,1:3],col=as.integer(factor(environment$Location)),size=3,type="s")

plot3d(scores(NMDSbray)[,1:3],col=as.integer(factor(environment$Location)),size=3,type="s")



plot(scores(metaMDS(bray,k=4))[,1]~environment$LAT,col=environment$Location)











library(rgl)

plot(scores(metaMDS(jac,k=4))[,1]~environment$LAT,col=environment$Location)












