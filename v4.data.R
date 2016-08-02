#Vellore environmental sample data;
setwd("~/stat/CMC/data/exposure data/")
re_dr<-read.csv("Drain_Sample_Database_092214_with_ec_conc.csv",header=TRUE)
hh_dat<-read.csv("Ind_data_102715.csv",header=TRUE)
pd_dat<-read.csv("phase2 public domain samples.csv",header=TRUE)

library(dplyr)
dat1<-re_dr %>% 
  select(SampleID,HouseID,ev_lat,ev_long,neighbor,ec_conc) %>% 
  mutate(study="re",SampleType="Drain")
dat2<-hh_dat %>% 
  select(SampleID,HouseID,ev_lat,ev_long,neighbor,SampleType,ec_conc) %>% 
  mutate(study="hh")
dat3<-pd_dat %>% 
  select(sample_name,GPS_latitude,GPS_longitude,neighborhood,sample_type_name,ec_conc) %>% 
  mutate(study="hh",HouseID="NA") %>% 
  rename(SampleID=sample_name,ev_lat=GPS_latitude,ev_long=GPS_longitude,neighbor=neighborhood,SampleType=sample_type_name)
#remove duplicated samples
dat3<-dat3[which(!duplicated(dat3[,1:6])),]

dat_all0<-rbind(dat1,dat2,dat3)
dat_all<-dat_all0[which(!is.na(dat_all0$ec_conc) & dat_all0$ev_long>=78 & dat_all0$ev_long<=79.16 & dat_all0$ev_lat>=12.8),]
dat_all$ec_lnconc<-log10(dat_all$ec_conc)
#replace log10(0) to a small value (-3);
dat_all$ec_lnconc[which(dat_all$ec_conc==0)]<--3

#Old Town = 1, Cinna Allapuram = 2;
dat_all$neighbor[which(dat_all$neighbor=="Old Town")]<-"1"
dat_all$neighbor[which(dat_all$neighbor=="Cinna Allapuram")]<-"2"

dat_all<-dat_all[which(dat_all$neighbor=="1"),]

########################################################################################
#Vellore Health outcome data
setwd("~/stat/CMC/data")
cmc<-read.csv("data_cmc.csv",header=TRUE)
cmc$sample_date<-as.Date(as.character(cmc$sdate),"%m/%d/%y")

#####################################################
cmc<-cmc[c(1:3752),]
cmc$last=0
count=0
for (i in 1:(length(cmc$IDNO)-1)){
  if (cmc$IDNO[i]==cmc$IDNO[i+1] & cmc$stooltype[i]=="D1"){
    count=count+1;
    cmc$n_illness[i]=count;
  }
  if (cmc$IDNO[i]!=cmc$IDNO[i+1] & cmc$stooltype[i]=="D1"){
    cmc$last[i]=1;
    count=count+1;
    cmc$n_illness[i]=count;
    count=0;
  }
  if (cmc$IDNO[i]==cmc$IDNO[i+1] & cmc$stooltype[i]=="M1"){
    cmc$n_illness[i]=count;
  }
  if (cmc$IDNO[i]!=cmc$IDNO[i+1] & cmc$stooltype[i]=="M1"){
    cmc$last[i]=1;
    cmc$n_illness[i]=count;
    count=0;
  }  
}
#hard code the last line;
cmc$last[length(cmc$last)]=1;
cmc$n_illness[length(cmc$last)]<-1;

cmc$bact<-cmc$bacteria
cmc$bact[which(cmc$stooltype=="D1")]=0
cmc$vir<-cmc$virus
cmc$vir[which(cmc$stooltype=="D1")]=0
cmc$test<-1
library(plyr)
n_pos_bact <- ddply ( cmc, .(cmc$IDNO), function(x) cumsum(x[55]))[,2]
n_pos_virus <- ddply ( cmc, .(cmc$IDNO), function(x) cumsum(x[56]))[,2]
n_test <- ddply ( cmc, .(cmc$IDNO), function(x) cumsum(x[57]))[,2]
cmc$n_pos_bact<-n_pos_bact
cmc$n_pos_virus<-n_pos_virus
cmc$n_test<-n_test

cmc_clu0<-cmc[which(cmc$last==1),]

start_date<-aggregate(sample_date ~ IDNO, cmc, function(x) min(x))
end_date<-aggregate(sample_date ~ IDNO, cmc, function(x) max(x))
study_date<-merge(start_date,end_date,by="IDNO")
colnames(study_date)[2:3]<-c("start_date","end_date")
study_date$duration<-study_date$end_date-study_date$start_date+1
cmc_clu<-merge(cmc_clu0,study_date)

#link each child with the closest environmental samples
# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

samtype<-c("Bathing Water", "Child Hand Rinse", "Drain", "Particulate", "Piped Water", "Produce", "Swabs", "Toy Feeding Spoon Rinse")

distance<-list()
sampleID<-list()
ec_lnconc<-list()
for (s in 1:length(samtype)){
  sp_dat<-dat_all[which(dat_all$SampleType==samtype[s]),]
  can<-array(NA,dim=c(length(cmc_clu$IDNO),length(sp_dat$SampleID)))
  id<-array(NA,dim=c(length(cmc_clu$IDNO),length(sp_dat$SampleID)))
  conc<-array(NA,dim=c(length(cmc_clu$IDNO),length(sp_dat$SampleID)))
  for (i in 1:length(cmc_clu$IDNO)){
    for (j in 1:length(sp_dat$SampleID)){
      can[i,j]<-gcd.hf(cmc_clu$LON[i],cmc_clu$LAT[i],sp_dat$ev_long[j],sp_dat$ev_lat[j])
      id[i,j]<-levels(sp_dat$SampleID)[sp_dat$SampleID[j]]
      conc[i,j]<-sp_dat$ec_lnconc[j]
    }
  }
  distance[[s]]<-can
  sampleID[[s]]<-id
  ec_lnconc[[s]]<-conc
}

#sampleID[[1]][1,which(distance[[1]][1,]==min(distance[[1]][1,]))]
#distance[[1]][1,which(distance[[1]][1,]==min(distance[[1]][1,]))]
#ec_lnconc[[1]][1,which(distance[[1]][1,]==min(distance[[1]][1,]))]

dat0<-cmc_clu[,c(2,5,6,7,12,16,54,58,59,60,61,62,63)]
for(x in (length(dat0[1,])+1):(length(dat0[1,])+24)){
  dat0[,x] <- NA
}
colnames(dat0)<-c("IDNO","LAT","LON","ALTITUDE","sanipath_pid","neighborhood","n_illness",
                  "n_pos_bact","n_pos_virus","n_test","start_date","end_date","duration",
                  "Bath_ID","Bath_Dis","Bath_lnconc","CHR_ID","CHR_Dis","CHR_lnconc",
                  "Drain_ID","Drain_Dis","Drain_lnconc","Part_ID","Part_Dis","Part_lnconc",
                  "Piped_ID","Piped_Dis","Piped_lnconc","Pro_ID","Pro_Dis","Pro_lnconc",
                  "Swab_ID","Swab_Dis","Swab_lnconc","TFSR_ID","TFSR_Dis","TFSR_lnconc")

for (s in 1:length(samtype)){
  for (i in 1:length(cmc_clu$IDNO)){
    dat0[i,3*(s-1)+length(dat0[1,])-23]<-paste0(sampleID[[s]][i,which(distance[[s]][i,]==min(distance[[s]][i,]))],collapse = ',')
    dat0[i,3*(s-1)+length(dat0[1,])-22]<-min(distance[[s]][i,which(distance[[s]][i,]==min(distance[[s]][i,]))])
    dat0[i,3*(s-1)+length(dat0[1,])-21]<-mean(ec_lnconc[[s]][i,which(distance[[s]][i,]==min(distance[[s]][i,]))],na.rm=TRUE)
  }
}
