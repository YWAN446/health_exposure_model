library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
library(ggmap)
library(geosphere)

setwd("~/stat/Vellore_study/Old_Town_sewage_ODfield_waterpipe_streets")
#sewage
sewage1 <- readOGR(".","BhNagar_sewagechannels_44N")
sewage2 <- readOGR(".","FBD_sewage_44N")
sewage3 <- readOGR(".","MGRNagar_sewageline_44N")
sewage4 <- readOGR(".","SSKMan_sewage_44N")
sewage5 <- readOGR(".","UthK_sewage_44N")

sewage <- gUnion(gUnion(gUnion(gUnion(sewage1,sewage2),sewage3),sewage4),sewage5)
sewage <- spTransform(sewage, CRS("+proj=longlat +datum=WGS84"))

#streets;
street1 <- readOGR(".","BhNagar_streets_44N")
street2 <- readOGR(".","FBD_streets_44N")
street3 <- readOGR(".","MGRNagar_streets_44N")
street4 <- readOGR(".","SSKMan_streets_44N")
street5 <- readOGR(".","UthK_streets_44N")

street <- gUnion(gUnion(gUnion(gUnion(street1,street2),street3),street4),street5)
street <- spTransform(street, CRS("+proj=longlat +datum=WGS84"))

#fecal area;
fecal1 <- readOGR(".","BhNagar_fecalpolygons_44N")
fecal2 <- readOGR(".","FBD_fecalpolygons_44N")
fecal3 <- readOGR(".","MGRNagar_fecalpolygons_44N")
fecal4 <- readOGR(".","SSKMan_fecalpolygons_44N")
fecal5 <- readOGR(".","UthK_fecalpolygons_44N")

fecal <- gUnion(gUnion(gUnion(gUnion(fecal1,fecal2),fecal3),fecal4),fecal5)
fecal <- spTransform(fecal, CRS("+proj=longlat +datum=WGS84"))

#the.sewage.projected <- spTransform(sewage, CRS( "+init=epsg:32635" ))
#sew <- raster(extent(sewage), crs=projection(sewage))
#the.points.sp<-SpatialPointsDataFrame(dat0[, c("LON", "LAT")], data.frame(ID=seq(1:nrow(dat0))),proj4string=CRS("+proj=longlat +datum=WGS84"))
#the.points.projected <- spTransform(the.points.sp[1, ], CRS( "+init=epsg:32635" ))
#the.circles.projected <- gBuffer(the.points.projected, width=100, byid=TRUE)
#the.circles.sp <- spTransform(the.circles.projected, CRS("+proj=longlat +datum=WGS84"))
#the.point.sp <- spTransform(the.circles.projected, CRS("+proj=longlat +datum=WGS84"))

#plot(the.sewage.projected)
#points(the.points.projected,col="red",cex=1)
#lines(the.circles.projected)

#function to calculate the length of line shape;
line.length <- function(the.lines.projected, the.circles.projected) {
  if (gIntersects(the.lines.projected, the.circles.projected)) {
    lines_crp <- crop(the.lines.projected, the.circles.projected)
    lines_crp_length <- gLength(lines_crp)
    return(lines_crp_length)
  } else {
    return(0)
  }
}

#sewage length
the.points.sp<-SpatialPointsDataFrame(dat0[, c("LON", "LAT")], data.frame(ID=seq(1:nrow(dat0))),proj4string=CRS("+proj=longlat +datum=WGS84"))
the.sewage.projected <- spTransform(sewage, CRS( "+init=epsg:32635" ))
sewage.length <- c()
for (i in 1:length(the.points.sp[,1])) {
  the.points.projected <- spTransform(the.points.sp[i, ], CRS( "+init=epsg:32635" ))
  the.circles.projected <- gBuffer(the.points.projected, width=100, byid=TRUE)
  sewage.length[i]<-line.length(the.sewage.projected, the.circles.projected)
}

#street length
the.points.sp<-SpatialPointsDataFrame(dat0[, c("LON", "LAT")], data.frame(ID=seq(1:nrow(dat0))),proj4string=CRS("+proj=longlat +datum=WGS84"))
the.street.projected <- spTransform(street, CRS( "+init=epsg:32635" ))
street.length <- c()
for (i in 1:length(the.points.sp[,1])) {
  the.points.projected <- spTransform(the.points.sp[i, ], CRS( "+init=epsg:32635" ))
  the.circles.projected <- gBuffer(the.points.projected, width=100, byid=TRUE)
  street.length[i]<-line.length(the.street.projected, the.circles.projected)
}

#distance to sewage
dist2sewage <- c()
for (i in 1:length(dat0[,1])) {
  dist2sewage[i]<-dist2Line(dat0[i,c(3,2)], sewage)[1]
}

#distance to street
dist2street <- c()
for (i in 1:length(dat0[,1])) {
  dist2street[i]<-dist2Line(dat0[i,c(3,2)], street)[1]
}


#distance to fecal area
dist2fecal <- c()
for (i in 1:length(dat0[,1])) {
  dist2fecal[i]<-dist2Line(dat0[i,c(3,2)], fecal)[1]
}





#############################################################
gnomic.buffer <- function(p, r) {
  stopifnot(length(p) == 1)
  gnom <- sprintf("+proj=gnom +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                  p@coords[[2]], p@coords[[1]])
  projected <- spTransform(p, CRS(gnom))
  buffered <- gBuffer(projected, width=r, byid=TRUE)
  spTransform(buffered, p@proj4string)
}

custom.buffer <- function(p, r) {
  stopifnot(length(p) == 1)
  cust <- sprintf("+proj=tmerc +lat_0=%s +lon_0=%s +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", 
                  p@coords[[2]], p@coords[[1]])
  projected <- spTransform(p, CRS(cust))
  buffered <- gBuffer(projected, width=r, byid=TRUE)
  spTransform(buffered, p@proj4string)
}

test.1 <- gnomic.buffer(the.points.sp[1,], 100)
test.2 <- custom.buffer(the.points.sp[1,], 100)

#plot(sewage)
#lines(test.1,col="red")
#lines(test.2,col="blue")
#points(the.points.sp)
##############################################################

cnt.child <- c()
for (i in 1:length(the.points.sp[,1])) {
  cnt.child[i] <- sum(over(the.points.sp,gnomic.buffer(the.points.sp[i,], 100))==i,na.rm = TRUE)
}

f.dat<-dat0
f.dat$n_neg_bact<-f.dat$n_test-f.dat$n_pos_bact
f.dat$n_neg_virus<-f.dat$n_test-f.dat$n_pos_virus
f.dat$sewage.length<-sewage.length
f.dat$street.length<-street.length
f.dat$dist2fecal<-dist2fecal
f.dat$dist2sewage<-dist2sewage
f.dat$dist2street<-dist2street
f.dat$cnt.child<-cnt.child

#remove those far away from the SaniPath samples
#plot(dat0$LON,dat0$LAT)
#points(dat_all$ev_long,dat_all$ev_lat,col="red")
f.dat <- f.dat[which(f.dat$LAT>12.910 & f.dat$LAT<12.9151),]



