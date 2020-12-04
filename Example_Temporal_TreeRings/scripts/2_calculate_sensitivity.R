################################################### 
# Copied over from 0_process_gamm.R
# This will give us the sensitivities of RW in a pretty format.

# Set up a dummy dataset for the script to run correctly
# number of simulations to run

load("../output_derived/gam.full_full_model.Rdata")


source("helper_functions/0_Calculate_GAMM_Posteriors.R")
# Fitting our model to the data to see if we're doing a pasable job of capturing the variance
# If things don't match, we shoudl take our sensitiivty curves with a grain of salt

n <- 100
model.pred <- post.distns(model.gam=gam.full, model.name="full_model", newdata=dat.raw, vars=predictors.all, n=n, terms=F)
summary(model.pred$ci)

# Need help dealing with the list that is set up here.  Need help with the aggregation
model.pred2 <- model.pred$ci
summary(model.pred2)
dim(model.pred2)

# Change predicted CI to BAI units
model.pred2[,c("mean.bai", "lwr.bai", "upr.bai")] <- exp(model.pred2[,c("mean", "lwr", "upr")])
summary(model.pred2)

# Aggregating to the group level in teh same way we do below with the ring widths
# We can then compare our modeled RW to our measured RW and see how things look
# Sanity Check #1

mean.model <- aggregate(model.pred2$mean, by = model.pred2[, c("Site.Code", "Year")], FUN=mean, na.rm=T)
names(mean.model)[names(mean.model)=="x"] <- c("BAI.mean") 
mean.model[,"BAI.lwr"] <- aggregate(model.pred2$mean, by=model.pred2[,c("Site.Code", "Year")], FUN=quantile, probs=0.025, na.rm=T)[,"x"]
mean.model[,"BAI.upr"] <- aggregate(model.pred2$mean, by=model.pred2[,c("Site.Code", "Year")], FUN=quantile, probs=0.975, na.rm=T)[,"x"]
head(mean.model)

mean.modelB <- aggregate(model.pred2$mean.bai, by = model.pred2[, c("Site.Code", "Year")], FUN=mean, na.rm=T)
names(mean.modelB)[names(mean.modelB)=="x"] <- c("BAI.mean") 
mean.modelB[,"BAI.lwr"] <- aggregate(model.pred2$mean.bai, by=model.pred2[,c("Site.Code", "Year")], FUN=quantile, probs=0.025, na.rm=T)[,"x"]
mean.modelB[,"BAI.upr"] <- aggregate(model.pred2$mean.bai, by=model.pred2[,c("Site.Code", "Year")], FUN=quantile, probs=0.975, na.rm=T)[,"x"]
head(mean.modelB)

# aggregating the raw data for graphing
mean.rw <- aggregate(dat.raw$BA.inc, by=dat.raw[, c("Site.Code", "Year")], FUN=mean, na.rm=T)
names(mean.rw)[names(mean.rw)=="x"]<- c("BAI.mean")                  


mean.rw[,"BAI.lwr"] <- aggregate(dat.raw$BA.inc, by=dat.raw[,c("Site.Code", "Year")], FUN=quantile, probs=0.025, na.rm=T)[,"x"]
mean.rw[,"BAI.upr"] <- aggregate(dat.raw$BA.inc, by=dat.raw[,c("Site.Code", "Year")], FUN=quantile, probs=0.975, na.rm=T)[,"x"]


head(mean.rw)


# Setting up lm between the modeled values and the observed	
dim(mean.rw)
dim(mean.model)


summary(mean.rw); 
summary(mean.model)

mean.model2 <- mean.model
names(mean.model2)[3:5] <- c("mod.mean", "mod.lwr", "mod.up")

mean.model2$exp.mod.mean <- exp(mean.model2$mod.mean)

gam.full.resid <- resid(gam.full)
hist(gam.full.resid)
dat.raw$resid.gam.full <- resid(gam.full)
plot(log(BA.inc) ~ resid.gam.full, data=dat.raw)
abline(a=0, b=1, col="red")

summary(mean.rw)

# Sanity Check #1 graph
pdf("../figures/prelim_figures/gam.full_sanitycheck1_mod_vs_real_rw.pdf", width= 13, height = 8.5)
 ggplot(data=mean.rw) + facet_grid(Site.Code~., scales="fixed") + theme_bw() +
          # plot the data
          geom_ribbon(aes(x=Year, ymin=BAI.lwr, ymax=BAI.upr), alpha=0.5) +
          geom_line(aes(x=Year, y=BAI.mean), size=1) +
          # Plot our model
          geom_ribbon(data=mean.model, aes(x=Year, ymin=exp(BAI.lwr), ymax=exp(BAI.upr)), fill="red3", alpha=0.3) +
          geom_line(data=mean.model, aes(x=Year, y=exp(BAI.mean)), color="red3", alpha=0.8, size=1) +
          labs(title="Gamm Model vs. Data", x="Year", y="BAI")

dev.off()


# Sanity Check #2
# Pulling random trees from both the data.use and the model.pred2 to see how they compare

n <- 10
data.use2 <- dat.raw


sanity2.trees <- sample(dat.raw$TreeID, size=n, replace=F) 
summary(sanity2.trees)

summary(mean.rw)
summary(model.pred2)

# Sanity Check #2 graph
pdf("../figures/prelim_figures/gam.full_sanitycheck2_tree_level.pdf", width= 13, height = 8.5)
ggplot(data=dat.raw[dat.raw$TreeID %in% sanity2.trees,]) + facet_wrap(~TreeID, scales="fixed") + theme_bw() +
  # plot the data
  #geom_ribbon(aes(x=Year, ymin=rw.lwr, ymax=rw.upr), alpha=0.5) +
  geom_line(aes(x=Year, y=BA.inc), size=1) +
  # Plot our model
  #geom_ribbon(data=model.pred2[model.pred2$TreeID %in% sanity2.trees,], aes(x=Year, ymin=rw.lwr, ymax=rw.upr), fill="red3", alpha=0.3) +
  geom_line(data=model.pred2[model.pred2$TreeID %in% sanity2.trees,], aes(x=Year, y=exp(mean)), color="red3", alpha=0.8, size=1) +
  labs(title="Gamm Model vs. Data Indiv. Trees", x="Year", y="RW")
dev.off()


# lm2 for sanity check 2
summary(model.pred2)
summary(dat.raw)

model.pred2$RW <- dat.raw$RW
model.pred2$BAI <- dat.raw$BA.inc

# LM for indiv. trees
sanity.lm2 <- lm(log(BAI) ~ mean, data=model.pred2)
summary(sanity.lm2)



# running scripts to get the weights
source("helper_functions/0_Calculate_GAMM_Weights.R")



predictors.all


vars <- c("tmean",  "precip", "Site.Code", "Year", "PlotID", "TreeID") # This should be your splines for those splines
gam.full.weights <- factor.weights(model.gam = gam.full, model.name = "size_temp", newdata = dat.raw, extent = "", vars = vars, limiting=T)

summary(gam.full.weights)
summary(dat.raw)
gam.full.weights[,c("BA.inc", "Site.Code")] <- dat.raw[,c("BA.inc", "Site.Code")] # Adding in factors we forgot

gam.full.weights[,c("fit.full.bai", "tmean.bai","precip.bai", "Year.bai")] <- exp(gam.full.weights[,c("fit.full", "fit.tmean","fit.precip", "fit.Year")])
summary(gam.full.weights)

# Just the weights of tmean and Precip, ignoring size
vars2 <- c("tmean.bai", "precip.bai", "Year.bai")
fit.spline2 <- rowSums(abs(gam.full.weights[,vars2]), na.rm=T)


for(v in vars2){
  gam.full.weights[,paste("weight", v, "2", sep=".")] <- gam.full.weights[,v]/fit.spline2
}
summary(gam.full.weights)


cols.weights <- c("weight.tmean.bai.2","weight.precip.bai.2","weight.Year.bai.2")
gam.full.weights$factor.max2 <- NA
for(i in 1:nrow(gam.full.weights)){
  fweight <- abs(gam.full.weights[i,cols.weights])
  gam.full.weights[i,"max2"] <- max(fweight, na.rm=T)
  gam.full.weights[i,"factor.max2"] <- c("tmean", "precip","Year")[which(fweight==max(fweight))]
}
gam.full.weights$factor.max2 <- as.factor(gam.full.weights$factor.max2)
summary(gam.full.weights)

save(gam.full.weights, file="../output_derived/gam.full_weights.Rdata")

