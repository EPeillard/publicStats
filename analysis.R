setwd("~/ExoAR/dataAnalysis/Rstuff")

library(car)
library(ggplot2)

library("afex")     # needed for ANOVA functions.
library("emmeans")  # emmeans must now be loaded explicitly for follow-up tests.
library("multcomp") # for advanced control for multiple testing/Type 1 errors.

afex_options(emmeans_model = "multivariate") # use multivariate model for all follow-up tests.

d <- read.table('trimmedData.csv',header = TRUE, sep = ',')

d$ans <- d$ans*100

#trim data to delete errors
d$relativeAns <- d$ans/d$e
d<-d[d$relativeAns>0.35,]

rawD <- d

#full factorial model
d$realIsRight <- factor(d$realIsRight)
d$rep <- factor(d$rep)
d$e <- factor(d$e)

# mean between the two rep
d<-aggregate(ans~realIsRight+e+type+id,data=d,FUN=mean)

#ggplot(data=d, aes(x=e, y=relativeAns, colour=factor(type)))+ stat_smooth(method=lm, fullrange=FALSE) + geom_point()

# test if the VR/RV factor is relevent
dVR <- d[which(d$type=="RV"),]
boxplot(dVR$ans ~ dVR$realIsRight)
fitRV <- lm(ans ~ e + realIsRight + (1|id), data=dVR)
summary(fitRV)
anova(fitRV)

boxVR0 <- dVR[which(dVR$realIsRight==0),]
boxVR1 <- dVR[which(dVR$realIsRight==1),]
par(mfrow=c(2,1))
boxplot(boxVR0$ans ~ boxVR0$id)
boxplot(boxVR1$ans ~ boxVR1$id)

par(mfrow=c(1,1))
boxplot(boxVR$ans ~ boxVR$realIsRight + boxVR$id)

# complete model
fitAll <- lm(ans ~ e * type * realIsRight + (1|id), data=d)
summary(fitAll)
anova(fitAll)


### previous
# since the anova is not significant, we delete this model and average the results
d<-aggregate(ans~e+type+id,data=d,FUN=mean)


#some plots
d$relativeAns <- d$ans/d$e

boxplot(d$ans ~ d$id)
boxplot(d$relativeAns ~ d$id)

boxplot(d$ans ~ d$e)
boxplot(d$relativeAns ~ d$e)

boxplot(d$ans ~ d$type)
boxplot(d$relativeAns ~ d$type)

### Old model, not relevent
fit1 <- lm(ans ~ e + type*realIsRight, data=d)
summary(fit1)
ano1<-aov(fit1)
TukeyHSD(ano1)

fit2 <- lm(ans ~ e + type, data=d)
summary(fit2)
ano2<-aov(fit2)
summary(ano2)

# layout(matrix(1:4,2,2))
# plot(fit1)

### plot the linear fit
par(mfrow=c(1,1))
plot <- ggplot(data=d, aes(x=e, y=ans, colour=factor(type)))
plot + stat_smooth(method=lm, fullrange=FALSE) + geom_point()
plot + stat_smooth(method=lm, fullrange=FALSE) + geom_point() + facet_wrap(~as.factor(id), nrow=5)

### compute anova
a1 <- aov_ez("id", "ans", d, between = NULL, within = c("e","type"),anova_table = list(es = "pes"))
nice(a1$anova_table)

### do post-hoc
m1 <- emmeans(a1$aov, ~ type)
summary(as.glht(pairs(m1)), test=adjusted("bonferroni"))
m2 <- emmeans(a1$aov, ~ e)
summary(as.glht(pairs(m2)), test=adjusted("bonferroni"))
m3 <- emmeans(a1$aov, ~ type+e)
summary(as.glht(pairs(m3)), test=adjusted("bonferroni"))

### agregate data for plot
aggD<-aggregate(ans~e+type, rawD ,function(x) c(mean=mean(x),sd = sd(x)))
aggD$sd<-aggD$ans[,"sd"]
aggD$ans<-aggD$ans[,"mean"]
aggD<-aggD[order(e,type),]

n<-length(unique(d$id))
x0<-aggD[aggD$type=="RR",]$e
x1<-aggD[aggD$type=="RV",]$e-1
x2<-aggD[aggD$type=="VV",]$e+1
y0<-aggD[aggD$type=="RR",]$ans
y1<-aggD[aggD$type=="RV",]$ans
y2<-aggD[aggD$type=="VV",]$ans
se0<-aggD[aggD$type=="RR",]$sd/sqrt(n)
se1<-aggD[aggD$type=="RV",]$sd/sqrt(n)
se2<-aggD[aggD$type=="VV",]$sd/sqrt(n)
plot(y0 ~ x0,
     ylim=c(20,100),
     xlim=c(20,65),
     main=" ",
     ylab="Perceived distance (cm)",
     xlab="Real distance (cm)",
     col="red")
points(y1 ~ x1, col="green")
points(y2~x2, col="blue")
arrows(x0,y0-se0/2,x0,y0+se0/2, length=0.05, angle=90, code=3, col="red")
arrows(x1,y1-se1/2,x1,y1+se1/2, length=0.05, angle=90, code=3, col="green")
arrows(x2,y2-se2/2,x2,y2+se2/2, length=0.05, angle=90, code=3, col="blue")
lines(c(0,120),c(0,120),lty=2,col="black")
legend(57,35,
       c("RR","RV","VV"),
       pch = 1,
       col=c("red","green","blue"))



### NEW

d <- rawD

dRR <- d[which(d$type=="RR"),]
dRR$type <- "RR"
dRV <- d[which(d$type=="RV" & d$realIsRight==0),]
dRV$type <- "RV"
dVR <- d[which(d$type=="RV" & d$realIsRight==1),]
dVR$type <- "VR"
dVV <- d[which(d$type=="VV"),]
dVV$type <- "VV"

d<- rbind(dRR,dRV,dVR,dVV)
d$type <- factor(d$type)

fit3 <- lm(ans ~ e + type, data=d)
summary(fit3)
ano3<-aov(fit3)
summary(ano3)


### agregate data for plot
aggD<-aggregate(ans~e+type, d ,function(x) c(mean=mean(x),sd = sd(x)))
aggD$sd<-aggD$ans[,"sd"]
aggD$ans<-aggD$ans[,"mean"]
aggD<-aggD[order(e,type),]

n<-length(unique(d$id))
x0<-aggD[aggD$type=="RR",]$e-2
x1<-aggD[aggD$type=="RV",]$e
x3<-aggD[aggD$type=="VR",]$e-1
x2<-aggD[aggD$type=="VV",]$e+1
y0<-aggD[aggD$type=="RR",]$ans
y1<-aggD[aggD$type=="RV",]$ans
y3<-aggD[aggD$type=="VR",]$ans
y2<-aggD[aggD$type=="VV",]$ans
se0<-aggD[aggD$type=="RR",]$sd/sqrt(n)
se1<-aggD[aggD$type=="RV",]$sd/sqrt(n)
se3<-aggD[aggD$type=="VR",]$sd/sqrt(n)
se2<-aggD[aggD$type=="VV",]$sd/sqrt(n)
plot(y0 ~ x0,
     ylim=c(20,100),
     xlim=c(20,65),
     main=" ",
     ylab="Perceived distance (cm)",
     xlab="Real distance (cm)",
     col="red")
points(y1 ~ x1, col="green")
points(y2~x2, col="blue")
points(y3~x3, col="orange")
arrows(x0,y0-se0/2,x0,y0+se0/2, length=0.05, angle=90, code=3, col="red")
arrows(x1,y1-se1/2,x1,y1+se1/2, length=0.05, angle=90, code=3, col="green")
arrows(x2,y2-se2/2,x2,y2+se2/2, length=0.05, angle=90, code=3, col="blue")
arrows(x3,y3-se3/2,x3,y3+se3/2, length=0.05, angle=90, code=3, col="orange")
lines(c(0,120),c(0,120),lty=2,col="black")
legend(57,35,
       c("RR","RV","VR","VV"),
       pch = 1,
       col=c("red","green","orange","blue"))



# realtive 

d$relativeAns <- d$ans / d$e

### agregate data for plot
aggD<-aggregate(relativeAns~e+type, d ,function(x) c(mean=mean(x),sd = sd(x)))
aggD$sd<-aggD$relativeAns[,"sd"]
aggD$relativeAns<-aggD$relativeAns[,"mean"]
aggD<-aggD[order(e,type),]

n<-length(unique(d$id))
x0<-aggD[aggD$type=="RR",]$e-2
x1<-aggD[aggD$type=="RV",]$e-1
x3<-aggD[aggD$type=="VR",]$e
x2<-aggD[aggD$type=="VV",]$e+1
y0<-aggD[aggD$type=="RR",]$relativeAns
y1<-aggD[aggD$type=="RV",]$relativeAns
y3<-aggD[aggD$type=="VR",]$relativeAns
y2<-aggD[aggD$type=="VV",]$relativeAns
se0<-aggD[aggD$type=="RR",]$sd/sqrt(n)
se1<-aggD[aggD$type=="RV",]$sd/sqrt(n)
se3<-aggD[aggD$type=="VR",]$sd/sqrt(n)
se2<-aggD[aggD$type=="VV",]$sd/sqrt(n)
plot(y0 ~ x0,
     ylim=c(1,1.5),
     xlim=c(20,65),
     main=" ",
     ylab="Perceived distance (cm)",
     xlab="Real distance (cm)",
     col="red")
points(y1 ~ x1, col="green")
points(y2~x2, col="blue")
points(y3~x3, col="orange")
arrows(x0,y0-se0/2,x0,y0+se0/2, length=0.05, angle=90, code=3, col="red")
arrows(x1,y1-se1/2,x1,y1+se1/2, length=0.05, angle=90, code=3, col="green")
arrows(x2,y2-se2/2,x2,y2+se2/2, length=0.05, angle=90, code=3, col="blue")
arrows(x3,y3-se3/2,x3,y3+se3/2, length=0.05, angle=90, code=3, col="orange")
#lines(c(0,120),c(0,120),lty=2,col="black")
legend(57,35,
       c("RR","RV","VR","VV"),
       pch = 1,
       col=c("red","green","orange","blue"))


# gaucher

d <- rawD

dRR <- d[which(d$type=="RR"),]
dRR$type <- "RR"
dRV <- d[which(d$type=="RV" & d$realIsRight==0),]
dRV$type <- "RV"
dVR <- d[which(d$type=="RV" & d$realIsRight==1),]
dVR$type <- "VR"
dVV <- d[which(d$type=="VV"),]
dVV$type <- "VV"

d<- rbind(dRR,dRV,dVR,dVV)
d$type <- factor(d$type)

d<- d[which(d$id == 11),]

boxplot(d$ans~d$type+d$e)

## participants

d <- rawD

dRR <- d[which(d$type=="RR"),]
dRR$type <- "RR"
dRV <- d[which(d$type=="RV" & d$realIsRight==0),]
dRV$type <- "RV"
dVR <- d[which(d$type=="RV" & d$realIsRight==1),]
dVR$type <- "VR"
dVV <- d[which(d$type=="VV"),]
dVV$type <- "VV"

d<- rbind(dRR,dRV,dVR,dVV)

p <- read.table('participants.csv',header=TRUE, sep = ',')

d<-merge(d,p,by.x = "id",by.y = "id")

dVdom1 <- d[which(d$type=="RV" & d$eye=="Droit"),]
dVdom1$type <- "Vdom"
dVdom2 <- d[which(d$type=="VR" & d$eye=="Gauche"),]
dVdom2$type <- "Vdom"
dRdom1 <- d[which(d$type=="RV" & d$eye=="Gauche"),]
dRdom1$type <- "Rdom"
dRdom2 <- d[which(d$type=="VR" & d$eye=="Droit"),]
dRdom2$type <- "Rdom"

dRR$eye <- NA
dVV$eye <- NA

d<-rbind(dVdom1,dVdom2,dRdom1,dRdom2,dRR,dVV)

d$type<-factor(d$type)

fit4 <- lm(ans ~ e + type, data=d)
summary(fit4)
ano4<-aov(fit4)
summary(ano4)


### agregate data for plot
aggD<-aggregate(ans~e+type, d ,function(x) c(mean=mean(x),sd = sd(x)))
aggD$sd<-aggD$ans[,"sd"]
aggD$ans<-aggD$ans[,"mean"]
aggD<-aggD[order(e,type),]

n<-length(unique(d$id))
x0<-aggD[aggD$type=="RR",]$e-2
x1<-aggD[aggD$type=="Vdom",]$e
x3<-aggD[aggD$type=="Rdom",]$e-1
x2<-aggD[aggD$type=="VV",]$e+1

y0<-aggD[aggD$type=="RR",]$ans
y2<-aggD[aggD$type=="VV",]$ans
se0<-aggD[aggD$type=="RR",]$sd/sqrt(n)
se2<-aggD[aggD$type=="VV",]$sd/sqrt(n)
y1<-aggD[aggD$type=="Vdom",]$ans
y3<-aggD[aggD$type=="Rdom",]$ans
se1<-aggD[aggD$type=="Vdom",]$sd/sqrt(n)
se3<-aggD[aggD$type=="Rdom",]$sd/sqrt(n)
plot(y0 ~ x0,
     ylim=c(0,100),
     xlim=c(20,65),
     main=" ",
     ylab="Perceived distance (cm)",
     xlab="Real distance (cm)",
     col="red")
points(y1 ~ x1, col="green")
points(y2~x2, col="blue")
points(y3~x3, col="orange")
arrows(x0,y0-se0/2,x0,y0+se0/2, length=0.05, angle=90, code=3, col="red")
arrows(x1,y1-se1/2,x1,y1+se1/2, length=0.05, angle=90, code=3, col="green")
arrows(x3,y3-se3/2,x3,y3+se3/2, length=0.05, angle=90, code=3, col="orange")
arrows(x2,y2-se2/2,x2,y2+se2/2, length=0.05, angle=90, code=3, col="blue")
lines(c(0,120),c(0,120),lty=2,col="black")
legend(57,35,
       c("RR","Rdom","Vdom","VV"),
       pch = 1,
       col=c("red","orange","green","blue"))

