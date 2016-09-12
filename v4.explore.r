#plot
plot(fecal,col="gold4",xlim=c(79.13718, 79.15012),ylim=c(12.9105,12.9151))
lines(sewage)
lines(street,col="brown")
points(the.points.sp)
lines(test.1,col="red")
legend(79.146,12.917,c("street","drain","open defecation area","child","100m circle"),bty='n',
       col=c("brown","black","gold4","black","red"),lty=c(1,1,-1,-1,-1),pch=c(NA,NA,15,1,1),cex=0.8)

#pdf(file="~/stat/CMC/v1/output/std_conc_plot_v1.pdf",height=8.5,width=14)
range01 <- function(x){((x-min(x))/(max(x)-min(x)))*2+0.1}
dat_all1<-dat_all[which(dat_all$SampleType=="Bathing Water"),]
plot(dat_all1$ev_long,dat_all1$ev_lat,cex=(range01(dat_all1$ec_lnconc)),
     xlim=c(79.13718, 79.15012),ylim=c(12.9105,12.9151),
     col="blue",
     xlab="Longtitude",ylab="Latitude",main="Standardized log10 E.coli concentration")
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
points(the.points.sp,pch=4)

legend(79.147,12.915,legend=c("Bathing Water","Drain","Child Hand Rinse","Piped Water",
                              "Particulate","Produce","Swabs","Toy Feeding Spoon Rinse","Child"),
       col=c("blue","black","green","purple","yellow","red","grey","orange","black"),pch=c(1,1,1,1,1,1,1,1,4),
       cex=0.8,bty='n')
#dev.off()

#Univariate logistic regression
res.bact<-cbind(f.dat$n_pos_bact,f.dat$n_neg_bact)
res.virus<-cbind(f.dat$n_pos_virus,f.dat$n_neg_virus)
res.ill<-round(f.dat$n_illness/as.numeric(f.dat$duration)*365)

#Bacteria
f.dat$w.Drain_lnconc<-f.dat$Drain_lnconc-log10(f.dat$Drain_Dis)
glm.res.bact <- glm(res.bact ~ f.dat$w.Drain_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)

glm.res.bact <- glm(res.bact ~ f.dat$Drain_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)

plot(f.dat$Bath_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$CHR_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$Drain_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$Part_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$Piped_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$Pro_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$Swab_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$TFSR_lnconc,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$sewage.length,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$street.length,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$dist2fecal,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$dist2sewage,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$dist2street,res.bact[,1]/(res.bact[,1]+res.bact[,2]))
plot(f.dat$cnt.child,res.bact[,1]/(res.bact[,1]+res.bact[,2]))

plot(f.dat$Bath_lnconc,jitter(res.ill))
plot(f.dat$CHR_lnconc,jitter(res.ill))
plot(f.dat$Drain_lnconc,jitter(res.ill))
plot(f.dat$Part_lnconc,jitter(res.ill))
plot(f.dat$Piped_lnconc,jitter(res.ill))
plot(f.dat$Pro_lnconc,jitter(res.ill))
plot(f.dat$Swab_lnconc,jitter(res.ill))
plot(f.dat$TFSR_lnconc,jitter(res.ill))
plot(f.dat$sewage.length,jitter(res.ill))
plot(f.dat$street.length,jitter(res.ill))
plot(f.dat$dist2fecal,jitter(res.ill))
plot(f.dat$dist2sewage,jitter(res.ill))
plot(f.dat$dist2street,jitter(res.ill))
plot(f.dat$cnt.child,jitter(res.ill))
