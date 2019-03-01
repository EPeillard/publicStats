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
fitRV <- lm(ans ~ e + realIsRight, data=d)
summary(fitRV)
anova(fitRV)

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
#fit1 <- lm(ans ~ e + type*realIsRight, data=d)
#summary(fit1)
#ano1<-aov(fit1)
#TukeyHSD(ano1)

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