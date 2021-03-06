library(ggplot2)
# composite figure script
path.ms <- "/Volumes/GoogleDrive/My Drive/Non-Stationarity_MSB/FEE_Copyedits/figures/"
# --------------------------------------
# Loading in data to make figures from script 8
dat.raw <- read.csv("../input_raw/tree_ring_input_data.csv", stringsAsFactors = T)

load("../output_derived/gam.temp_temp_only.Rdata") # gam.temp
load("../output_derived/gam.time_time_only.Rdata") # gam.time
load("../output_derived/gam.time.temp_time_temp.Rdata") # gam.time.temp 
load("../output_derived/gam.full_full_model.Rdata") # gam.full

#################################################################
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
resid.simple <- ggplot(data=resid.graph[resid.graph$Year<2013 & resid.graph$type %in% c("Temp Only", "Temp + Time"),], aes(x=Year, y=resids)) + facet_grid(type~.) +
  geom_point(alpha=0.2, stroke=0, size=0.5) +
  geom_hline(aes(yintercept=0), col="red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        strip.text.y=element_text(face="bold", size=10),
        axis.text=element_text(color="black", size=9),
        axis.title=element_text(face="bold", size=12)) +
  scale_x_continuous(breaks = c(1900, 1925, 1950, 1975, 2000)) +
  coord_cartesian(ylim=quantile(resid.graph$resids, c(0.001, 0.999), na.rm=T)) +
  labs(x=expression(bold(paste("Year"))), y = expression(bold(paste("Residual Value"))))



# Sensitivity curves
load("../output_derived/gam.full_response_graph.Rdata")
load("../output_derived/gam.time.temp_response_time_temp.Rdata")
load("../output_derived/gam.temp_response_graph.Rdata")


sens.curves.simple <- ggplot() + 
  geom_hline(yintercept=100, linetype="dashed")+
  geom_ribbon(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai*100, ymax=upr.bai*100, fill="Temp Only"), alpha=0.5) +
  geom_line(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai*100, color="Temp Only")) +
  geom_ribbon(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai*100, ymax=upr.bai*100, fill="Temp + Time"), alpha=0.5) +
  geom_line(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai*100, color="Temp + Time")) +
  guides(color=F, fill=guide_legend(title=NULL)) +
  scale_fill_manual(values=c("#0072B2", "#E69F00")) +
  scale_color_manual(values=c("#0072B2", "#E69F00")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text=element_text(color="black", size=12),
        axis.title=element_text(face="bold", size=12)) +
  theme(legend.position=c(0.65,0.25),
        legend.background = element_blank(),
        legend.key = element_rect(fill=NA),
        legend.text = element_text(size=11)) +
  labs(x = "Temperature", y = "Relativized BAI (%)")


# Saving the manuscript figure
png(filename=file.path(path.ms, "Figure4_composite_nonstationarity_TR_graph.png"), height=7.5, width=11.5, unit="cm", res=300)
cowplot::plot_grid(resid.simple, sens.curves.simple, ncol = 2, labels = c("A)", "B)"))
dev.off()

jpeg(filename=file.path(path.ms, "Figure4_composite_nonstationarity_TR_graph.jpeg"), height=7.5, width=11.5, unit="cm", res=300)
cowplot::plot_grid(resid.simple, sens.curves.simple, ncol = 2, labels = c("A)", "B)"))
dev.off()

tiff(filename=file.path(path.ms, "Figure4_composite_nonstationarity_TR_graph.tiff"), height=7.5, width=11.5, unit="cm", res=300)
cowplot::plot_grid(resid.simple, sens.curves.simple, ncol = 2, labels = c("A)", "B)"))
dev.off()



resid.simple2 <- resid.simple + geom_point(size=1, stroke=0, alpha=0.2) +  geom_hline(aes(yintercept=0), col="red") + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), strip.text.y=element_text(size=16))

tiff(filename=file.path(path.ms, "Figure4a_composite_nonstationarity_TR_graph.tiff"), height=7.5*1.5, width=11.5, unit="cm", res=300)
cowplot::plot_grid(resid.simple2, ncol = 1, labels = c("A)"))
dev.off()


sens.curves.simple2 <- sens.curves.simple + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text =element_text(size=16))
tiff(filename=file.path(path.ms, "Figure4b_composite_nonstationarity_TR_graph.tiff"), height=7.5*1.5, width=11.5, unit="cm", res=300)
cowplot::plot_grid(sens.curves.simple2, ncol = 1, labels = c("B)"))
dev.off()

#---------------------------------------
# Loading in figure from script 9

# load("gam.full_data_graph.Rdata")
load("../output_derived/gam.full_data_graph.Rdata") # gam.full


residual.plot <- ggplot(data=resid.graph[resid.graph$Year<2013,]) + facet_grid(type~.) +
  geom_point(aes(x=Year, y=resids), alpha=0.2, stroke=0, size=1.25) +
  geom_hline(aes(yintercept=0), col="red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        strip.text.y=element_text(face="bold", size=15),
        axis.text=element_text(color="black", size=18),
        axis.title=element_text(face="bold", size=18)) +
  scale_x_continuous(breaks = c(1900, 1920, 1940, 1960, 1980, 2000, 2020)) +
  labs(x=expression(bold(paste("Year"))), y = expression(bold(paste("Residual Value"))))

sens.curves <- ggplot() + 
  geom_hline(yintercept=100, linetype="dashed")+
  geom_ribbon(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai*100, ymax=upr.bai*100, fill="Temp Only"), alpha=0.5) +
  geom_line(data=temp.ci.out[temp.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai*100, color="Temp Only")) +
  geom_ribbon(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai*100, ymax=upr.bai*100, fill="Temp + Time"), alpha=0.5) +
  geom_line(data=time.temp.ci.out2[time.temp.ci.out2$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai*100, color="Temp + Time")) +
  geom_ribbon(data=full.ci.out[full.ci.out$Effect %in% c("tmean"), ], aes(x=x, ymin=lwr.bai*100, ymax=upr.bai*100, fill="Temp + Time + Precip"), alpha=0.5) +
  geom_line(data=full.ci.out[full.ci.out$Effect %in% c("tmean"), ], aes(x=x, y=mean.bai*100, color="Temp + Time + Precip")) +
  guides(color=F, fill=guide_legend(title=NULL)) +
  scale_fill_manual(values=c("#0072B2", "#009E73", "#E69F00")) +
  scale_color_manual(values=c("#0072B2", "#009E73", "#E69F00")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text=element_text(color="black", size=16),
        axis.title=element_text(face="bold", size=16),
        legend.text=element_text(size=14)) +
  theme(legend.position=c(0.6,0.25)) +
  labs(x = "Temperature", y = "Relativized BAI (%)")



# Effects Curve
gam.effects <- ggplot(data.graph[data.graph$Site.Code=="LF" & data.graph$Year<2013,]) + 
                    geom_hline(aes(yintercept=100), linetype="dashed") +
                    geom_ribbon(aes(x=Year, ymin=fit.tmean.lwr*100, ymax=fit.tmean.upr*100, fill="Temp"), alpha=0.5) +
                    geom_ribbon(aes(x=Year, ymin=fit.precip.lwr*100, ymax=fit.precip.upr*100, fill="Precip"), alpha=0.5) + 
                    
                    
                    
                    geom_line(aes(x=Year, y=fit.tmean*100, color="Temp"), size=1) +
                    geom_line(aes(x=Year, y=fit.precip*100, color="Precip"), size=1) +
                    
                    scale_color_manual(values=c("blue", "red"), labels=c("Precip", "Temp")) +
                    scale_fill_manual(values=c("blue", "red"), labels=c("Precip", "Temp"), name="") +
                    guides(color=F) +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank())+
                    theme(axis.line.x = element_line(color="black", size = 0.5),
                          axis.line.y = element_line(color="black", size = 0.5),
                          axis.text=element_text(color="black", size=16),
                          axis.title=element_text(color="black", face="bold", size=16),
                          legend.text=element_text(color="black", size=14),
                          legend.position = c(0.75, 0.25))+
                    labs(x=expression(bold(paste("Year"))), y = "Relativized BAI (%)")


library(cowplot)

combo.curves.effects <- plot_grid(sens.curves, gam.effects, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels = c("B)", "C)"))

comp.plot <- plot_grid(residual.plot, combo.curves.effects, ncol = 2, rel_heights = c(1/2, 1/2), labels = c("A)", ""))

tiff(filename=file.path(path.ms, "WebFigure2_composite_nonstationarity_TR_graph_full.tiff"), height=8, width=11, unit="in", res=300)
comp.plot
dev.off()

png(filename=file.path(path.ms, "WebFigure2_composite_nonstationarity_TR_graph_full.png"), height=8, width=11, unit="in", res=300)
comp.plot
dev.off()

