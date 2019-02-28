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

# fit2 <- lm(ans ~ e + type, data=d)
# summary(fit2)
# ano2<-aov(fit2)
# summary(ano2)

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
