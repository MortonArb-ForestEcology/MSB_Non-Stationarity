dat.raw <- read.csv("~/Desktop/Research/MSB_Non-Stationarity/Example_Temporal_TreeRings/input_raw/tree_ring_input_data.csv")

summary(dat.raw)

dat.agg <- aggregate(dat.raw[,c("tmean", "precip", "vpd.min", "vpd.max")], by=dat.raw[,c("Site.Code", "Year")], FUN=mean, na.rm=T)
summary(dat.agg)

lm.warm <- lm(tmean ~ Year, data=dat.agg)
sum.warm <- summary(lm.warm)
sum.warm$coefficients[2,1]*10

plot(tmean ~ Year, data=dat.agg, pch=19)
abline(lm.warm, col="red")

lm.precip <- lm(precip ~ Year, data=dat.agg)
sum.precip <- summary(lm.precip)
sum.precip
sum.precip$coefficients[2,1]*10

plot(precip ~ Year, data=dat.agg, pch=19)
abline(lm.precip, col="red")