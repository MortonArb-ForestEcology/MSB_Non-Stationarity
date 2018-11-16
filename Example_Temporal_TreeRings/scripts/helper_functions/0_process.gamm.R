process.gamm <- function(gamm.model, data, model.name, vars, resolution="t.001", extent=c(850,2010), 
						  	fweights=T, sites=T, ci.model=T, ci.terms=T, n=250,
						  	write.out=T, outdir, control=list()){
	# data     = data frame with data in it
	# model.name    = which model.name to subset
	# k        = number of knots in the spline
	# n        = number of simulations for generating confidence intervals
	# outdir   = where to save the .Rdata file
	library(mgcv)
	source("R/0_Calculate_GAMM_Weights.R")
	source("R/0_Calculate_GAMM_Posteriors.R")

	out <- list(data=data, gamm=gamm.model)
    # -----------
	# Calculating the Factor Weights through time
	# -----------
	if(fweights==T){	
		f.weights <- factor.weights(model.gam=gamm.model, model.name=model.name, newdata=data, extent=extent, vars=vars); 
		out[["weights"]] <- f.weights 
	}	
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.model==T){
		ci.response <- post.distns(model.gam=gamm.model, model.name=model.name, n=n, newdata=data, vars=vars, terms=F)
		out[["ci.response"]]  <- ci.response$ci 
		out[["sim.response"]] <- ci.response$sims 
	}
	# -----------
	
    # -----------
	# Calculating the CI around our response prediction
	# -----------
	if(ci.terms==T){
		n.out = n

		new.dat <- data.frame(Model=model.name,
							  Extent=as.factor(paste(extent[1], extent[2], sep="-")), 
							  Resolution=resolution)
		
		# Common extra variables used in fitting the models
		vars.extra <- c("Site", "PlotID", "TreeID", "PFT")
		for(v in vars.extra){
			if(v %in% names(data)){
			  if(v %in% vars){
			 	# if v is a variable used to fit curves, merge it into the data frame
			 	var.temp <- data.frame(x=unique(data[,v])) 
			 	names(var.temp) <- v
			 	new.dat <- merge(new.dat, var.temp)
			  } else {
			 	# if v is not a sensitivity variable, just pick the first value 
			 	# because it doesn't matter
			 	new.dat[,v] <- unique(data[,v])[1]
			  }				
			} 
		}

		# Figure out which vars are numeric vs. factor
		vars.num <- vector()
		for(v in vars){
			if(class(data[,v]) %in% c("numeric", "integer")) vars.num <- c(vars.num, v)
		}

		var.temp <- data.frame(array(dim=c(n.out, length(vars.num))))
		names(var.temp) <- vars.num
		for(v in vars.num){
			var.temp[,v] <- seq(min(data[,v], na.rm=T), max(data[,v], na.rm=T), length.out=n.out)
		}								
		new.dat <- merge(new.dat, var.temp, all.x=T, all.y=T)
								
		ci.terms.pred <- post.distns(model.gam=gamm.model, model.name=model.name, n=n, newdata=new.dat, vars=vars, terms=T)

		out[["ci.terms"]]  <- ci.terms.pred$ci 
		out[["sim.terms"]] <- ci.terms.pred$sims 
	}	
	# -----------
	if(write.out==T) save(out, file=file.path(outdir, paste("gamm", model.name, resolution, "AllSites", "Rdata", sep=".")))
	return(out)
	# -----------

}