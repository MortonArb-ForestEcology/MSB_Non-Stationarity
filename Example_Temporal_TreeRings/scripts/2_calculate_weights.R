library(ggplot2)
library(car)
load(file="gam.full_weights.Rdata")

se <- function(x){
  sd(x, na.rm=TRUE) / sqrt((length(!is.na(x))))}

summary(gam.full.weights)
factors.fits <- c("fit.tmean", "fit.precip", "fit.Year", "fit.full", "BA.inc")

factors.weights <- c("weight.tmean", "weight.Year", "weight.precip")
gam.full.weights[,c("weight.tmean", "weight.precip", "weight.Year")] <- gam.full.weights[,c("weight.tmean", "weight.precip", "weight.Year")]/rowSums(gam.full.weights[,c("weight.tmean", "weight.precip", "weight.Year")],na.rm=T)

# Transforming things back to BA.inc rather than log
gam.full.weights[,which(substr(names(gam.full.weights),1,3)=="fit")] <- exp(gam.full.weights[,which(substr(names(gam.full.weights),1,3)=="fit")] )

othervars <- c("Year", "Site.Code", "Model")

data.graph1 <- aggregate(gam.full.weights[,factors.fits], by = gam.full.weights[,othervars], FUN= mean, na.rm=T)

data.graph1[,paste(factors.fits, "upr", sep=".")] <- aggregate(gam.full.weights[,factors.fits], by = gam.full.weights[,othervars], FUN= quantile, prob= 0.975, na.rm=T)[,factors.fits]

data.graph1[,paste(factors.fits, "lwr", sep=".")] <- aggregate(gam.full.weights[,factors.fits], by = gam.full.weights[,othervars], FUN= quantile, prob= 0.05, na.rm=T)[,factors.fits]



summary(data.graph1)

data.graph <- aggregate(abs(gam.full.weights[,factors.weights]), by = gam.full.weights[,othervars], FUN= mean, na.rm=T)

data.graph[,paste(factors.weights, "upr", sep=".")] <- aggregate(abs(gam.full.weights[,factors.weights]), by = gam.full.weights[,othervars], FUN= quantile, prob= 0.975, na.rm=T)[,factors.weights]

data.graph[,paste(factors.weights, "lwr", sep=".")] <- aggregate(abs(gam.full.weights[,factors.weights]), by = gam.full.weights[,othervars], FUN= quantile, prob= 0.05, na.rm=T)[,factors.weights]

data.graph[,paste(factors.weights, "SD", sep=".")] <- aggregate(abs(gam.full.weights[,factors.weights]), by = gam.full.weights[,othervars], FUN= sd, na.rm=T)[,factors.weights]

data.graph[,paste(factors.weights, "SE", sep=".")] <- aggregate(abs(gam.full.weights[,factors.weights]), by = gam.full.weights[,othervars], FUN= se)[,factors.weights]


summary(data.graph)

data.graph <- merge(data.graph1, data.graph, all.x=T, all.y=T)

# data.graph <- gam.full.weights[gam.full.weights$TreeID== "MMA014",]
summary(data.graph)
gam.full.weights$wts.check <- rowSums(abs(gam.full.weights[,factors.weights]))
data.graph$wts.check <- rowSums(abs(data.graph[,factors.weights]))

summary(gam.full.weights)
save(gam.full.weights, file="gam.full_weights_processed.Rdata")
save(data.graph, file="gam.full_data_graph.Rdata")
summary(data.graph)

#  PLotting only one site because they should be identical lines  due to the model structure
gam.effects <- ggplot(data.graph[data.graph$Site.Code=="LF",]) + 
                  	geom_hline(aes(yintercept=1), linetype="dashed") +
                    geom_ribbon(aes(x=Year, ymin=fit.tmean.lwr, ymax=fit.tmean.upr), fill="red", alpha=0.5) +
                  	geom_ribbon(aes(x=Year, ymin=fit.precip.lwr, ymax=fit.precip.upr), fill="blue", alpha=0.5) +
                  
                  	
                  	geom_line(aes(x=Year, y=fit.tmean), size=1, color="red") +
                  	geom_line(aes(x=Year, y=fit.precip), size=1, color="blue") +
                  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank())+
                  	theme(axis.line.x = element_line(color="black", size = 0.5),
                          axis.line.y = element_line(color="black", size = 0.5))+
                  	scale_x_continuous(breaks = c(1890, 1900, 1910, 1920,  1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020)) +
                    labs(x=expression(bold(paste("Year"))), y = expression(bold(paste("Relative Effect Size"))))
                  
save(gam.effects, file = "gam.effects.graph.Rdata")

pdf("HF_QURU_effect_size.pdf", width= 13, height = 8.5)
  gam.effects
dev.off()
