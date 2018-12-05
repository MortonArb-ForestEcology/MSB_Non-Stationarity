# composite figure script

# --------------------------------------
# Loading in data to make figures from script 8

load("../output_derived/gam.temp_temp_only.Rdata") # gam.temp
load("../output_derived/gam.time_time_only.Rdata") # gam.time
load("../output_derived/gam.time.temp_time_temp.Rdata") # gam.time.temp 
load("../output_derived/gam.full_full_model.Rdata") # gam.full

#################################################################
# Important!!!
# Run Lines 110+ in 8_GAM_script_clean.R to ensure that data are loaded to make this figure
#################################################################

# needs reorganized to plot in gglpot.
resid.temp <- data.frame(resids = resid(gam.temp),
                        Year = dat.raw$Year,
                        type = as.factor("Temp Only"))

resid.time.temp <- data.frame(resids = resid(gam.time.temp),
                              Year = dat.raw$Year,
                              type = as.factor("Temp + Time"))

resid.full <- data.frame(resids = resid(gam.full),
                              Year = dat.raw$Year,
                              type = as.factor("Temp + Time + Precip."))

resid.graph <- rbind(resid.temp, resid.time.temp, resid.full)


#  Residuals
residual.plot <- ggplot(data=resid.graph[resid.graph$Year<2013,]) + facet_grid(type~.) +
                    geom_point(aes(x=Year, y=resids), alpha=0.2, stroke=0) +
                    geom_hline(aes(yintercept=0), col="red") +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank())+
                    theme(axis.line.x = element_line(color="black", size = 0.5),
                          axis.line.y = element_line(color="black", size = 0.5),
                          strip.text.y=element_text(face="bold")) +
                    scale_x_continuous(breaks = c(1900, 1920, 1940, 1960, 1980, 2000, 2020)) +
                    labs(x=expression(bold(paste("Year"))), y = expression(bold(paste("Residual Value"))))
  


# Sensitivity curves
load("../output_derived/gam.full_response_graph.Rdata")
load("../output_derived/gam.time.temp_response_time_temp.Rdata")
load("../output_derived/gam.temp_response_graph.Rdata")

sens.curves <- ggplot() + 
                  geom_hline(yintercept=1, linetype="dashed")+
                  geom_ribbon(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp Only"), alpha=0.5) +
                  geom_line(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp Only")) +
                  geom_ribbon(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp + Time"), alpha=0.5) +
                  geom_line(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp + Time")) +
                  geom_ribbon(data=full.ci.out[full.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai, ymax=upr.bai, fill="Temp + Time + Precip"), alpha=0.5) +
                  geom_line(data=full.ci.out[full.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai, color="Temp + Time + Precip")) +
                  guides(color=F, fill=guide_legend(title=NULL)) +
                  scale_fill_manual(values=c("#0072B2", "#009E73", "#E69F00")) +
                  scale_color_manual(values=c("#0072B2", "#009E73", "#E69F00")) +
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank())+
                  theme(axis.line.x = element_line(color="black", size = 0.5),
                        axis.line.y = element_line(color="black", size = 0.5)) +
                  theme(legend.position=c(0.4,0.25)) +
                  labs(x = "Temperature", y = expression(bold(paste("Effect on BAI (mm"^"2","y"^"-1",")")))) 

#---------------------------------------
# Loading in figure from script 9

load("gam.full_data_graph.Rdata")

# Effects Curve

gam.effects <- ggplot(data.graph[data.graph$Site.Code=="LF" & data.graph$Year<2013,]) + 
                    geom_hline(aes(yintercept=1), linetype="dashed") +
                    geom_ribbon(aes(x=Year, ymin=fit.tmean.lwr, ymax=fit.tmean.upr, fill="Temp"), alpha=0.5) +
                    geom_ribbon(aes(x=Year, ymin=fit.precip.lwr, ymax=fit.precip.upr, fill="Precip"), alpha=0.5) + 
                    
                    
                    
                    geom_line(aes(x=Year, y=fit.tmean, color="Temp"), size=1) +
                    geom_line(aes(x=Year, y=fit.precip, color="Precip"), size=1) +
                    
                    scale_color_manual(values=c("blue", "red"), labels=c("Precip", "Temp")) +
                    scale_fill_manual(values=c("blue", "red"), labels=c("Precip", "Temp"), name="") +
                    guides(color=F) +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank())+
                    theme(axis.line.x = element_line(color="black", size = 0.5),
                          axis.line.y = element_line(color="black", size = 0.5),
                          legend.position = c(0.75, 0.25))+
                    
                    labs(x=expression(bold(paste("Year"))), y = expression(bold(paste("Relative Effect Size"))))


library(cowplot)

combo.curves.effects <- plot_grid(sens.curves, gam.effects, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels = c("B)", "C)"))

comp.plot <- plot_grid(residual.plot, combo.curves.effects, ncol = 2, rel_heights = c(1/2, 1/2), labels = c("A)", ""))

png(filename="../figures/composite_nonstationarity_TR_graph.png", height=8, width=11, unit="in", res=300)
comp.plot
dev.off()


