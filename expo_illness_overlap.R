#Vellore environmental sample analysis;
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

dat_all0<-rbind(dat1,dat2,dat3)
dat_all<-dat_all0[which(!is.na(dat_all0$ec_conc) & dat_all0$ev_long>=78 & dat_all0$ev_long<=79.16 & dat_all0$ev_lat>=12.8),]
dat_all$ec_lnconc<-log10(dat_all$ec_conc)
#replace log10(0) to a small value (-3);
dat_all$ec_lnconc[which(dat_all$ec_conc==0)]<--3

#Old Town = 1, Cinna Allapuram = 2;
dat_all$neighbor[which(dat_all$neighbor=="Old Town")]<-"1"
dat_all$neighbor[which(dat_all$neighbor=="Cinna Allapuram")]<-"2"

dat_all<-dat_all[which(dat_all$neighbor=="1"),]

range01 <- function(x){((x-min(x))/(max(x)-min(x)))*2+0.1}

dat_all1<-dat_all[which(dat_all$SampleType=="Bathing Water"),]
plot(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
     col="blue",xlim=c(min(dat_all1$ev_long),max(dat_all1$ev_long)),
     ylim=c(min(dat_all1$ev_lat),max(dat_all1$ev_lat)),xlab="Longtitude",ylab="Latitude",main="Standardized log10 E.coli concentration")
dat_all1<-dat_all[which(dat_all$SampleType=="Drain"),]
points(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
       col="black")
dat_all1<-dat_all[which(dat_all$SampleType=="Child Hand Rinse"),]
points(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
       col="green")
dat_all1<-dat_all[which(dat_all$SampleType=="Piped Water"),]
points(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
       col="purple")
dat_all1<-dat_all[which(dat_all$SampleType=="Particulate"),]
points(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
       col="yellow")
dat_all1<-dat_all[which(dat_all$SampleType=="Produce"),]
points(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
       col="red")
dat_all1<-dat_all[which(dat_all$SampleType=="Swabs"),]
points(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
       col="grey")
dat_all1<-dat_all[which(dat_all$SampleType=="Toy Feeding Spoon Rinse"),]
points(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
       col="orange")

legend("topright",legend=c("Bathing Water","Drain","Child Hand Rinse","Piped Water",
                              "Particulate","Produce","Swabs","Toy Feeding Spoon Rinse"),
       col=c("blue","black","green","purple","yellow","red","grey","orange"),pch=1,cex=0.5)





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
library(plyr)
n_pos_bact <- ddply ( cmc, .(cmc$IDNO), function(x) cumsum(x[55]))[,2]
n_pos_virus <- ddply ( cmc, .(cmc$IDNO), function(x) cumsum(x[56]))[,2]
cmc$n_pos_bact<-n_pos_bact
cmc$n_pos_virus<-n_pos_virus

cmc_clu<-cmc[which(cmc$last==1),]

#illness
for (i in 1:length(unique(cmc_clu$n_illness))){
  points(cmc_clu$LON[which(cmc_clu$n_illness==i)],cmc_clu$LAT[which(cmc_clu$n_illness==i)],pch=16,
         col=heat.colors(length(unique(cmc_clu$n_illness)))[length(unique(cmc_clu$n_illness))-i+1])
}



#virus infection
for (i in 1:length(unique(cmc_clu$n_pos_virus))){
  points(cmc_clu$LON[which(cmc_clu$n_pos_virus==i)],cmc_clu$LAT[which(cmc_clu$n_pos_virus==i)],pch=16,
         col=heat.colors(length(unique(cmc_clu$n_pos_virus)))[length(unique(cmc_clu$n_pos_virus))-i+1])
}