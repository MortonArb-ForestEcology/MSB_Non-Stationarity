
library(mgcv)
library(ggplot2)
library(car)

# GAM script for non-stationarity paper
# Reading in Harvard Forest data
### Only red oak from two different plots at Harvard Forest
dat.raw <- read.csv("../input_raw/tree_ring_input_data.csv", header=T)
summary(dat.raw)
################################################### 
# HERE'S THE GAMM!!!
################################################### 
# Illustrating non-stationarity as an issue
# Trying an LM using year as a factor to see if effect of tmean varies by year
lm.temp <- lm(log(BA.inc) ~ tmean*as.factor(Year) + 
                Site.Code + PlotID  + TreeID, 
              data=dat.raw)
anova(lm.temp)
summary(lm.temp)$r.squared

# Running several different GAMs to look at the sensitivities of the different forms

# gam.time
# Size only GAM -- no cliamte effects, just equivalent to standard tree ring a priori detrending with a spline
gam.time <- gam(log(BA.inc) ~ s(Year, k=3, by=TreeID) +
                  Site.Code + PlotID  + TreeID,
                # random=list(Site=~1, PlotID=~1, TreeID=~1),
                data=dat.raw)

summary(gam.time)$r.sq # R-squared
summary(gam.time)$dev.expl # explained deviance
# anova(gam.time)
AIC(gam.time)


# gam.temp
# Temperature only -- temperature effect if we don't detrend for time
gam.temp <- gam(log(BA.inc)~ s(tmean, k=3) +
                  Site.Code + PlotID  + TreeID,
                # random=list(Site=~1, PlotID=~1, TreeID=~1),
                data=dat.raw)
summary(gam.temp)$r.sq # R-squared
summary(gam.temp)$dev.expl # explained deviance
anova(gam.temp) 
AIC(gam.temp)


# gam.time.temp
# Temperature and Size GAM
gam.time.temp <- gam(log(BA.inc)~ s(tmean, k=3) +
                       # s(dbh.recon, k=3) +
                       s(Year, k=3, by=TreeID)+
                       Site.Code + PlotID  + TreeID,
                     # random=list(Site=~1, PlotID=~1, TreeID=~1),
                     data=dat.raw)
# Look at the R-squared and explained deviance
summary(gam.time.temp)$r.sq # R-squared
summary(gam.time.temp)$dev.expl # explained deviance
# anova(gam.time.temp) 
AIC(gam.time.temp)
plot(resid(gam.time.temp) ~ dat.raw$Year, ylim=c(-3,2))

# Trying a model with temperature on the residuals of the time-only model just for kicks
dat.raw$time.pred <- predict(gam.time, newdata=dat.raw)
dat.raw$time.resid <- resid(gam.time)
plot(time.resid ~ Year, data=dat.raw, ylim=c(-3,2))
summary(dat.raw)
gam.temp2 <- gam(time.resid~ s(tmean, k=3) +
                  Site.Code + PlotID  + TreeID,
                # random=list(Site=~1, PlotID=~1, TreeID=~1),
                data=dat.raw)
summary(gam.temp2)$r.sq # R-squared
summary(gam.temp2)$dev.expl # explained deviance
anova(gam.temp2) 
AIC(gam.temp2)
par(mfrow=c(2,1))
plot(gam.temp)
plot(gam.temp2)
par(mfrow=c(1,1))



# gam.full
# 'full' model looking at effects of size, temperature, and precipitation
gam.full <- gam(log(BA.inc)~ s(tmean, k=3) +
                  s(precip, k=3) +
                  s(Year, k=3, by=TreeID)+
                  Site.Code + PlotID  + TreeID,
                # random=list(Site=~1, PlotID=~1, TreeID=~1),
                data=dat.raw)
summary(gam.full)$r.sq # R-squared
summary(gam.full)$dev.expl # explained deviance
anova(gam.full) 
AIC(gam.full)

png("../figures/Resids_TempModels.png", height=8, width=11, units="in", res=180)
par(mfrow=c(3,1))
plot(resid(gam.temp) ~ dat.raw$Year, ylim=c(-3,2), main="Resids: Temperature Only"); abline(h=0, col="red")
plot(resid(gam.time.temp) ~ dat.raw$Year, ylim=c(-3,2), main="Resids: Temperature + Time"); abline(h=0, col="red")
plot(resid(gam.full) ~ dat.raw$Year, ylim=c(-3,2), main="Resids: Temperature + Time + Precip"); abline(h=0, col="red")
dev.off(); par(mfrow=c(1,1))


save(gam.temp, file="../output_derived/gam.temp_temp_only.Rdata") 
save(gam.time, file="../output_derived/gam.time_time_only.Rdata")
save(gam.time.temp, file="../output_derived/gam.time.temp_time_temp.Rdata") 
save(gam.full, file="../output_derived/gam.full_full_model.Rdata") 

#----------------------------------------------
# Creating a dataframe that will be used to calculate the weights and the sensitivities
n <- 100
data <- dat.raw


n.out = n

new.dat <- data.frame(Model="nonstationarity",
                      Extent=as.factor(paste(min(data$Year), max(data$Year), sep="-")))

# Figure out which vars are numeric vs. factor
vars.num <- vector()

for(v in predictors.all){
  if(class(data[,v]) %in% c("numeric", "integer")) vars.num <- c(vars.num, v)
}
vars.num

# Getting the unique values of our factor variables and adding them to the data frame
predictors.all <- c("tmean", "precip", "Site.Code", "Year", "PlotID", "TreeID") 
spline.by=c("TreeID") # the "by" terms in the models you're running

for(v in predictors.all[!predictors.all %in% vars.num & !(predictors.all %in% c("Site"))]){
  # if v is a factor, merge all unique values into the dataframe
  if(!(v %in% spline.by)){ # Only pull the full range of values for whatever the "by" term was by, otherwise everythign should have the same shape, just different intercepts
    var.temp <- data.frame(x=unique(data[,v])[1]) 
  }else {
    var.temp <- data.frame(x=unique(data[,v])) 
  }
  names(var.temp) <- v
  new.dat <- merge(new.dat, var.temp, all.x=T, all.y=T)
}
# Separate out Plot & Site
# new.dat$PlotID <- as.factor(ifelse(substr(new.dat$TreeID, 1, 3)=="HOW", substr(new.dat$TreeID, 1, 4), substr(new.dat$TreeID, 1, 3)))

# Matching the site for the plot
for(p in unique(new.dat$PlotID)){
  new.dat[new.dat$PlotID==p, "Site"] <- unique(dat.raw[dat.raw$PlotID==p, "Site.Code"])
}
new.dat$Site <- as.factor(new.dat$Site)
summary(new.dat)


# Putting the numerical variables into an array and adding it in 
var.temp <- data.frame(array(dim=c(n.out, length(vars.num))))
names(var.temp) <- vars.num
for(v in vars.num){
  var.temp[,v] <- seq(min(data[,v], na.rm=T), max(data[,v], na.rm=T), length.out=n.out)
}								
new.dat <- merge(new.dat, var.temp, all.x=T, all.y=T)
summary(new.dat)


write.csv(new.dat, file="../output_derived/sensitivity_extaction_dataframe.csv", row.names=F)								


#############################################
# Plotting sensitivity curves
#############################################
#----------------------------------------------
# Temperature effect only
load("../output_derived/gam.temp_temp_only.Rdata")		

new.dat.temp <- new.dat
vars.fac <- c("Site.Code", "PlotID", "TreeID")
var.smooth <- "TreeID"
for(v in vars.fac){
  if(v %in% var.smooth) next # keep all levels for our "by" variable
  # Get rid of unimportant levels for everything else
  l1 <- unique(new.dat.temp[,v])[1]
  new.dat.temp <- new.dat.temp[new.dat.temp[,v]==l1,]
}
summary(new.dat.temp)


source("0_Calculate_GAMM_Posteriors.R")
temp.ci.terms.pred <- post.distns(model.gam=gam.temp, model.name="Site_level", n=n, newdata=new.dat.temp, vars=c("tmean"), terms=T)

temp.ci.out <- temp.ci.terms.pred$ci # separting out the confidence interval 
temp.ci.out[,predictors.all[!predictors.all %in% vars.num]] <- new.dat.temp[,predictors.all[!predictors.all %in% vars.num]] # copying over our factor labels
temp.ci.out$x <- as.numeric(temp.ci.out$x) # making x numeric; will make factors NA; NA's are ok here
summary(temp.ci.out)

# Convert mean, lwr, upr to BAI units
temp.ci.out[,c("mean.bai", "lwr.bai", "upr.bai")] <- exp(temp.ci.out[,c("mean", "lwr", "upr")])
summary(temp.ci.out)


save(temp.ci.out, file="../output_derived/gam.temp_response_graph.Rdata")

pdf("../figures/prelim_figures/gam.temp_sensitivities_tmean_only.pdf", width= 13, height = 8.5)
ggplot(data=temp.ci.out[temp.ci.out$Effect %in% "tmean", ]) + 
  geom_ribbon(aes(x=x, ymin=exp(lwr), ymax=exp(upr)), alpha=0.5) +
  geom_line(aes(x=x, y=exp(mean)))+
  geom_hline(yintercept=1, linetype="dashed") +
  theme_bw()+
  labs(x = "Temperature", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")"))))
dev.off()


# ---------------
# Year effects only
load("../output_derived/gam.time_time_only.Rdata")

n <- 100						
# SOurce & run the function
source("0_Calculate_GAMM_Posteriors.R")
# Make things run faster by reducing dimensions
new.dat.time <- new.dat
vars.fac <- c("Site.Code", "PlotID", "TreeID")
var.smooth <- "TreeID"
for(v in vars.fac){
  if(v == var.smooth) next # keep all levels for our "by" variable
  # Get rid of unimportant levels for everything else
  l1 <- unique(new.dat.time[,v])[1]
  new.dat.time <- new.dat.time[new.dat.time[,v]==l1,]
}
summary(new.dat.time)

# ------
# Fixing species and canopy class names so that they make sense
# ------
summary(new.dat.time)
# ------


time.ci.terms.pred <- post.distns(model.gam=gam.time, model.name="time only", n=n, newdata=new.dat.time, vars="Year", terms=T)

time.ci.out <- time.ci.terms.pred$ci # separting out the confidence interval 
time.ci.out[,predictors.all[!predictors.all %in% vars.num]] <- new.dat.time[,predictors.all[!predictors.all %in% vars.num]] # copying over our factor labels
time.ci.out$x <- as.numeric(time.ci.out$x) # making x numeric; will make factors NA
summary(time.ci.out)

# The model was on log(BAI), so we need to re-transform into real units
time.ci.out[,c("mean.bai", "lwr.bai", "upr.bai")] <- exp(time.ci.out[,c("mean", "lwr", "upr")])

# Truncating to observed range --> not necessary
# Year
# for(TREE in unique(dat.raw$TreeID)){
#   yr.min <- min(dat.raw[dat.raw$TreeID==TREE & !is.na(dat.raw$BA.inc), "Year"])
#   yr.max <- max(dat.raw[dat.raw$TreeID==TREE & !is.na(dat.raw$BA.inc), "Year"])
#   
#   time.ci.out$x <- ifelse(time.ci.out$TreeID!=TREE | time.ci.out$Effect!="year" | (time.ci.out$Effect==time.ci.out$x>=dbh.min & time.ci.out$x<=dbh.max), time.ci.out$x, NA)
#   
# }

save(time.ci.out,file="../output_derived/gam.time_response_graph.Rdata")

# For the year-only model, show the spline for some random trees
set.seed(1300)
trees.random <- sample(unique(time.ci.out$TreeID),25)

pdf("../figures/prelim_figures/gam.time_sensitivities_year.pdf", width= 13, height = 8.5)		
ggplot(data=time.ci.out[time.ci.out$TreeID %in% trees.random, ]) + 
  facet_wrap(~TreeID, scales="free_y") +
  geom_ribbon(aes(x=x, ymin=exp(lwr), ymax=exp(upr)), alpha=0.5) +
  geom_line(aes(x=x, y=exp(mean))) +
  theme_bw()+
  labs(x = "Year", y = expression(paste("Effect on BAI (mm"^"2", "y"^"-1",")")))
dev.off()


#----------------------------------------------
# Year + temperature effects
load("../output_derived/gam.time.temp_time_temp.Rdata")

new.dat.time.temp <- new.dat
vars.fac <- c("Site.Code", "PlotID", "TreeID")
var.smooth <- "TreeID"
for(v in vars.fac){
  if(v == var.smooth) next # keep all levels for our "by" variable
  # Get rid of unimportant levels for everything else
  l1 <- unique(new.dat.time.temp[,v])[1]
  new.dat.time.temp <- new.dat.time.temp[new.dat.time.temp[,v]==l1,]
}
summary(new.dat.time.temp)


# SOurce & run the function
source("0_Calculate_GAMM_Posteriors.R")
time.temp.ci.terms.pred2 <- post.distns(model.gam=gam.time.temp, model.name="spp_only", n=n, newdata=new.dat.time.temp, vars=c("Year", "tmean"), terms=T)

time.temp.ci.out2 <- time.temp.ci.terms.pred2$ci # separting out the confidence interval 
time.temp.ci.out2[,predictors.all[!predictors.all %in% vars.num]] <- new.dat.time.temp[,predictors.all[!predictors.all %in% vars.num]] # copying over our factor labels
time.temp.ci.out2$x <- as.numeric(time.temp.ci.out2$x) # making x numeric; will make factors NA
summary(time.temp.ci.out2)

time.temp.ci.out2[,c("mean.bai", "lwr.bai", "upr.bai")] <- exp(time.temp.ci.out2[,c("mean", "lwr", "upr")])




save(time.temp.ci.out2, file="../output_derived/gam.time.temp_response_time_temp.R")

png("../figures/prelim_figures/gam.time.temp_sensitivities_tmean.png", width= 8, height = 8, units="in", res=180)		
ggplot(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ]) + 
  geom_hline(yintercept=1, linetype="dashed")+
  geom_ribbon(aes(x=x, ymin=lwr.bai, ymax=upr.bai), alpha=0.5) +
  geom_line(aes(x=x, y=mean.bai)) +
  theme_bw()+
  labs(x = "Climate Variable", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")")))) 
dev.off()	


png("../figures/prelim_figures/Sensitivity_Comparison_tmean.png", width= 8, height = 8, units="in", res=180)		
ggplot() + 
  geom_hline(yintercept=1, linetype="dashed")+
  geom_ribbon(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp Only"), alpha=0.5) +
  geom_line(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp Only")) +
  geom_ribbon(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp + Time"), alpha=0.5) +
  geom_line(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp + Time")) +
  guides(color=F) +
  theme_bw() +
  theme(legend.position=c(0.75,0.25)) +
  labs(x = "Climate Variable", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")")))) 
dev.off()

png("../figures/prelim_figures/gam.time.temp_sensitivities_year.png", width= 11, height = 8, units="in", res=180)		
ggplot(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("Year") & time.temp.ci.out2$TreeID %in% trees.random, ]) + 
  facet_wrap(~TreeID, scales="free_y") +
  # geom_line(aes(x=x, y=0), linetype="dashed")+
  geom_ribbon(aes(x=x, ymin=lwr.bai, ymax=upr.bai), alpha=0.5) +
  geom_line(aes(x=x, y=mean.bai))+
  theme_bw()+
  labs(x = "Climate Variable", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")"))))+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+
  scale_y_continuous(expand=c(0,0))
dev.off()	





#----------------------------------------------
# Full Model: Temp + Precip + Year

load("../output_derived/gam.full_full_model.Rdata")		

new.dat.full <- new.dat
vars.fac <- c("Site.Code", "PlotID", "TreeID")
var.smooth <- c("TreeID")
for(v in vars.fac){
  if(v %in% var.smooth) next # keep all levels for our "by" variable
  # Get rid of unimportant levels for everything else
  l1 <- unique(new.dat.full[,v])[1]
  new.dat.full <- new.dat.full[new.dat.full[,v]==l1,]
}
# new.dat.full$spp.cc <- as.factor(paste(new.dat.full$Species, new.dat.full$Canopy.Class, sep="."))
summary(new.dat.full)

source("0_Calculate_GAMM_Posteriors.R")
full.ci.terms.pred <- post.distns(model.gam=gam.full, model.name="Additive", n=n, newdata=new.dat.full, vars=c("Year", "tmean", "precip"), terms=T)

full.ci.out <- full.ci.terms.pred$ci # separting out the confidence interval 
full.ci.out[,predictors.all[!predictors.all %in% vars.num]] <- new.dat.full[,predictors.all[!predictors.all %in% vars.num]] # copying over our factor labels
full.ci.out$x <- as.numeric(full.ci.out$x) # making x numeric; will make factors NA; NA's are ok here
summary(full.ci.out)

# Convert mean, lwr, upr to BAI units
full.ci.out[,c("mean.bai", "lwr.bai", "upr.bai")] <- exp(full.ci.out[,c("mean", "lwr", "upr")])
summary(full.ci.out)


summary(full.ci.out)
save(full.ci.out, file="../output_derived/gam.full_response_graph.Rdata")

png("../figures/prelim_figures/gam.full_sensitivities_observed_tmean_precip.pdf", width= 8, height = 8, units="in", res=180)
ggplot(data=full.ci.out[full.ci.out$Effect %in% c("tmean", "precip"), ]) + 
  facet_wrap(~Effect, scales="free_x") +
  geom_ribbon(aes(x=x, ymin=exp(lwr), ymax=exp(upr)), alpha=0.5) +
  geom_line(aes(x=x, y=exp(mean)))+
  geom_hline(yintercept=1, linetype="dashed") +
  coord_cartesian(ylim=c(0.5, 1.5)) +
  theme_bw()+
  labs(x = "Climate Variable", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")"))))
dev.off()

png("../figures/prelim_figures/Sensitivity_Comparison_tmean.png", width= 8, height = 8, units="in", res=180)		
ggplot() + 
  geom_hline(yintercept=1, linetype="dashed")+
  geom_ribbon(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp Only"), alpha=0.5) +
  geom_line(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp Only")) +
  geom_ribbon(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp + Time"), alpha=0.5) +
  geom_line(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp + Time")) +
  geom_ribbon(data=full.ci.out[full.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp + Time + Precip"), alpha=0.5) +
  geom_line(data=full.ci.out[full.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp + Time + Precip")) +
  guides(color=F) +
  theme_bw() +
  theme(legend.position=c(0.75,0.25)) +
  labs(x = "Temperature", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")")))) 
dev.off()


png("../figures/prelim_figures/gam.full_sensitivities_year.png", width= 11, height = 8, units="in", res=180)		
ggplot(data=full.ci.out[full.ci.out$Effect %in% c("Year") & full.ci.out$TreeID %in% trees.random, ]) + 
  facet_wrap(~TreeID, scales="free_y") +
  # geom_line(aes(x=x, y=0), linetype="dashed")+
  geom_ribbon(aes(x=x, ymin=lwr.bai, ymax=upr.bai), alpha=0.5) +
  geom_line(aes(x=x, y=mean.bai))+
  theme_bw()+
  labs(x = "Climate Variable", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")"))))+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+
  scale_y_continuous(expand=c(0,0))
dev.off()	






################################################################################
# Some stats on the models
################################################################################
# These probably need to be cleaned up
dat.raw$gam.time.pred <- predict(gam.time, newdata=dat.raw)
dat.raw$gam.time.temp.pred <- predict(gam.time.temp, newdata=dat.raw)
dat.raw$gam.temp.pred <- predict(gam.temp, newdata=dat.raw)
dat.raw$gam.full.pred <- predict(gam.full, newdata=dat.raw)

# Looking at the obs vs. full predicted fixed and mixed effects
meow1 <- lm(log(BA.inc)~gam.time.pred, data=dat.raw)
plot(log(BA.inc)~gam.time.pred, data=dat.raw)
abline(meow1, col="red", lwd=3)

# Residuals through time
dat.raw$gam.time.res <- resid(gam.time)
plot(dat.raw$gam.time.res~dat.raw$Year)


meow2 <- lm(log(BA.inc)~gam.time.temp.pred, data=dat.raw)
plot(log(BA.inc)~gam.time.temp.pred, data=dat.raw)
abline(meow2, col="red", lwd=3)

# Residuals through time
dat.raw$gam.time.temp.res <- resid(gam.time.temp)
plot(dat.raw$gam.time.temp.res~dat.raw$Year)


meow4 <- lm(log(BA.inc)~gam.temp.pred, data=dat.raw)
plot(log(BA.inc)~gam.temp.pred, data=dat.raw)
abline(meow4, col="red", lwd=3)

# Residuals through time
dat.raw$gam.temp.res <- resid(gam.temp)
plot(dat.raw$gam.temp.res~dat.raw$Year)


meow5 <- lm(log(BA.inc)~gam.full.pred, data=dat.raw)
plot(log(BA.inc)~gam.full.pred, data=dat.raw)
abline(meow5, col="red", lwd=3)


# Residuals through time
dat.raw$gam.full.res <- resid(gam.full)
plot(dat.raw$gam.full.res~dat.raw$Year)

