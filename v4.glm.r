#logistic regression model
res.bact<-cbind(f.dat$n_pos_bact,f.dat$n_neg_bact)
glm.res.bact <- glm(res.bact ~ f.dat$Bath_lnconc + f.dat$CHR_lnconc + f.dat$Drain_lnconc + f.dat$Part_lnconc + f.dat$Piped_lnconc +
                        f.dat$Pro_lnconc + f.dat$Swab_lnconc + f.dat$TFSR_lnconc + f.dat$sewage.length + f.dat$street.length +
                        f.dat$dist2fecal + f.dat$dist2sewage + f.dat$dist2street + f.dat$cnt.child, family=binomial(link='logit'))

summary(glm.res.bact)
#hist(glm.res.bact$residuals)

res.virus<-cbind(f.dat$n_pos_virus,f.dat$n_neg_virus)
glm.res.virus <- glm(res.virus ~ f.dat$Bath_lnconc + f.dat$CHR_lnconc + f.dat$Drain_lnconc + f.dat$Part_lnconc + f.dat$Piped_lnconc +
                      f.dat$Pro_lnconc + f.dat$Swab_lnconc + f.dat$TFSR_lnconc + f.dat$sewage.length + f.dat$street.length +
                      f.dat$dist2fecal + f.dat$dist2sewage + f.dat$dist2street + f.dat$cnt.child, family=binomial(link='logit'))

summary(glm.res.virus)

res.ill<-round(f.dat$n_illness/as.numeric(f.dat$duration)*365)
glm.res.ill <- glm(res.ill ~ f.dat$Bath_lnconc + f.dat$CHR_lnconc + f.dat$Drain_lnconc + f.dat$Part_lnconc + f.dat$Piped_lnconc +
                     f.dat$Pro_lnconc + f.dat$Swab_lnconc + f.dat$TFSR_lnconc + f.dat$sewage.length + f.dat$street.length +
                     f.dat$dist2fecal + f.dat$dist2sewage + f.dat$dist2street + f.dat$cnt.child, family=poisson(link='log'))

summary(glm.res.ill)

#Univariate logistic regression
res.bact<-cbind(f.dat$n_pos_bact,f.dat$n_neg_bact)
res.virus<-cbind(f.dat$n_pos_virus,f.dat$n_neg_virus)
res.ill<-round(f.dat$n_illness/as.numeric(f.dat$duration)*365)

#Bacteria
glm.res.bact <- glm(res.bact ~ f.dat$Bath_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$CHR_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$Drain_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$Part_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$Piped_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$Pro_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$Swab_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$TFSR_lnconc, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$sewage.length, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$street.length, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$dist2fecal, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$dist2sewage, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$dist2street, family=binomial(link='logit'))
summary(glm.res.bact)
glm.res.bact <- glm(res.bact ~ f.dat$cnt.child, family=binomial(link='logit'))
summary(glm.res.bact)

#Virus
glm.res.virus <- glm(res.virus ~ f.dat$Bath_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$CHR_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$Drain_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$Part_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$Piped_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$Pro_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$Swab_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$TFSR_lnconc, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$sewage.length, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$street.length, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$dist2fecal, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$dist2sewage, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$dist2street, family=binomial(link='logit'))
summary(glm.res.virus)
glm.res.virus <- glm(res.virus ~ f.dat$cnt.child, family=binomial(link='logit'))
summary(glm.res.virus)
#Virus seems not association with log10 E. coli concentration for all pathways

#Illness
glm.res.ill <- glm(res.ill ~ f.dat$Bath_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$CHR_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$Drain_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$Part_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$Piped_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$Pro_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$Swab_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$TFSR_lnconc, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$sewage.length, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$street.length, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$dist2fecal, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$dist2sewage, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$dist2street, family=poisson(link='log'))
summary(glm.res.ill)
glm.res.ill <- glm(res.ill ~ f.dat$cnt.child, family=poisson(link='log'))
summary(glm.res.ill)

#model selection
#Bacteria
#AIC
glm.res.bact <- glm(res.bact ~ f.dat$Bath_lnconc + f.dat$CHR_lnconc + f.dat$Drain_lnconc + f.dat$Part_lnconc + f.dat$Piped_lnconc +
                      f.dat$Pro_lnconc + f.dat$Swab_lnconc + f.dat$TFSR_lnconc + f.dat$sewage.length + f.dat$street.length +
                      f.dat$dist2fecal + f.dat$dist2sewage + f.dat$dist2street + f.dat$cnt.child, family=binomial(link='logit'))
step(glm.res.bact)

glm.res.bact.aic <- glm(res.bact ~ f.dat$CHR_lnconc + f.dat$Drain_lnconc + 
                        f.dat$Part_lnconc + f.dat$TFSR_lnconc + f.dat$sewage.length + 
                        f.dat$cnt.child, family = binomial(link = "logit"))
summary(glm.res.bact.aic)

#BIC
step(glm.res.bact,k=log(165))

#Virus
#AIC
glm.res.virus <- glm(res.virus ~ f.dat$Bath_lnconc + f.dat$CHR_lnconc + f.dat$Drain_lnconc + f.dat$Part_lnconc + f.dat$Piped_lnconc +
                      f.dat$Pro_lnconc + f.dat$Swab_lnconc + f.dat$TFSR_lnconc + f.dat$sewage.length + f.dat$street.length +
                      f.dat$dist2fecal + f.dat$dist2sewage + f.dat$dist2street + f.dat$cnt.child, family=binomial(link='logit'))
step(glm.res.virus)

glm.res.virus.aic <- glm(res.virus ~ f.dat$Drain_lnconc + f.dat$Pro_lnconc + 
                           f.dat$Swab_lnconc + f.dat$street.length, family = binomial(link = "logit"))
summary(glm.res.virus.aic)

#BIC
step(glm.res.virus,k=log(165))

#illness
#AIC
glm.res.ill <- glm(res.ill ~ f.dat$Bath_lnconc + f.dat$CHR_lnconc + f.dat$Drain_lnconc + f.dat$Part_lnconc + f.dat$Piped_lnconc +
                       f.dat$Pro_lnconc + f.dat$Swab_lnconc + f.dat$TFSR_lnconc + f.dat$sewage.length + f.dat$street.length +
                       f.dat$dist2fecal + f.dat$dist2sewage + f.dat$dist2street + f.dat$cnt.child, family=poisson(link='log'))
step(glm.res.ill)

#BIC
step(glm.res.ill,k=log(165))


