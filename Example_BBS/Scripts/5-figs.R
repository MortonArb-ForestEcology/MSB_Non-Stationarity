####################
# 5 - create figs
#
# -NOCA range with boxes on map
# -Slope estimates for BIO15 covariate on map
# -Ensemble predictions on map
# -Global - ensemble predictions on map
# 
# Creates:
# Figures
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(rgeos)
library(ggplot2)
library(sp)
library(RColorBrewer)


# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Example_BBS/Data/'))


# Medium box RDS ----------------------------------------------------------

#global and regional fit summaries
global_fit_df <- readRDS('global_fit_df_M.rds')
regional_fit_df <- readRDS('regional_fit_df_M.rds')

#global and regional fit objects
predict_df <-readRDS('predict_df_M.rds')

#box spdf
box_spdf <- readRDS('BOX_M_spdf.rds')

#BBS sites
BBS_bio <- readRDS('bbs_bio_BOX_M_5930_2018-10-30.rds')


# range map ---------------------------------------------------------------

#shp file from BirdLife International
setwd(paste0(dir, 'Example_BBS/Data/shp_files'))
sp_rng <- rgdal::readOGR('Cardinalis_cardinalis_22723819.shp', verbose = FALSE)

sp_rng@data$id <- rownames(sp_rng@data)
sp_rng_pts <- ggplot2::fortify(sp_rng, region = 'id')
sp_rng_df <- merge(sp_rng_pts, sp_rng@data, by = 'id')


usa_m <- maps::map('usa', plot = FALSE, fill = TRUE)
IDs <- sapply(strsplit(usa_m$names, ":"), function(x) x[1])
usa <- maptools::map2SpatialPolygons(usa_m, IDs = IDs, proj4string = sp::CRS('+init=epsg:4326'))


spydf_states <- rgeos::gBuffer(usa, byid=TRUE, width=0)

# process fit data ------------------------------------------------------------

#load US map
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])

#get box centroids and add to regional_fit_df
b_cen <- rgeos::gCentroid(box_spdf, byid=TRUE)

regional_fit_df[c('cen_x', 'cen_y')] <- cbind(b_cen@coords[,1], b_cen@coords[,2])

#add proportion of surveys detected
regional_fit_df <- dplyr::mutate(regional_fit_df, prop = n_pres/n_survey)

#detections must be X% of total surveys for box
plt_rf <- regional_fit_df[which(!is.na(regional_fit_df$gamma_beta_bio1)),]

#box number used in slope plotting
bn <- plt_rf$box_num



# Process ensemble predict data -------------------------------------------

#col indices with 'box_num' used for slope plotting
cind <- which(colnames(predict_df) %in% paste0('predict_box_num_', bn))

#function to see how many cells created each ensemble
sna_fun <- function(x)
{
  tt <- sum(!is.na(x[cind]))
  return(tt)
}

nc_en <- apply(predict_df, 1, sna_fun)

mdfun <- function(x)
{
  tt <- median(as.numeric(x[cind]), na.rm = TRUE)
  return(tt)
}

ensemble <- apply(predict_df, 1, mdfun)
to.rm <- which(is.na(ensemble))

predict_ensemble_df_p <- data.frame(predict_df[,1:8], ensemble, 
                                    diff = predict_df[,8] - ensemble)

predict_ensemble_df <- predict_ensemble_df_p[-to.rm,]




# Panel a -------------------------------------------------

#convert boxes to df
box_spdf@data$id <- rownames(box_spdf@data)
box_pts <- ggplot2::fortify(box_spdf, region = 'id')
box_df <- merge(box_pts, box_spdf@data, by = 'id')

#sample 75 boxes
set.seed(5)
box_id_s <- sample(box_df$id, size = 75)
box_df_s <- dplyr::filter(box_df, id %in% box_id_s)



box_map <- ggplot() +
  #geom_path(data = usamap, aes(x = x, y = y), color = 'black', alpha = 0.7) +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  geom_polygon(data = sp_rng_df,
               aes(x = long, y = lat, group=group), fill = 'orange', alpha = 0.8) +
  geom_path(data = box_df_s, aes(x = long, y = lat, group = group), color = 'black',
            alpha = 0.5) +
  geom_point(data = BBS_bio, aes(x = lng, y = lat), color = 'black',
            alpha = 0.4, size = 0.5) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12))


box_map



# Panel b ----------------------------------------

hist(plt_rf$gamma_beta_bio15)
BREAKS <- c(-1, -0.2, -0.1, -0.025, 0.025, 0.1, 0.2, 1)
cols <- colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(BREAKS)-1)


bio15_map <- ggplot() +
  #geom_path(data = usamap, aes(x = x, y = y), color = 'black', alpha = 0.7) +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  geom_point(data = plt_rf, inherit.aes = FALSE,
             aes(x = cen_x, y = cen_y, color = cut(gamma_beta_bio15, breaks = BREAKS)), 
             alpha = 0.5, size = 4) +
  scale_color_manual('Slope - BIO15', values = rev(cols), 
                     labels = c('< -0.2', '-0.2 to -0.1', '-0.1 to -0.025',
                                '-0.025 to 0.025', 
                                '0.025 to 0.1', '0.1 to 0.2', '> 0.2')) +
  guides(color = guide_legend(override.aes = list(size = 10), reverse = TRUE)) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))

bio15_map


# Panel c ---------------------------------------------

hist(predict_ensemble_df$ensemble)
BREAKS <- c(-500, 10, 20, 30, 40, 50, 60, 500)
cols <- colorRampPalette(brewer.pal(9, 'Reds'))(length(BREAKS) - 1)

ens_pred_map <- ggplot() +
  #geom_path(data = usamap, aes(x = x, y = y), color = 'black', alpha = 0.7) +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  geom_point(data = predict_ensemble_df, inherit.aes = FALSE,
             aes(x = x, y = y, color = cut(ensemble, breaks = BREAKS)), 
             alpha = 0.8, size = 0.2) +
  scale_color_manual('Predicted abundance', values = cols, 
                     labels = c('< 10', '10 to 20', '20 to 30', 
                                '30 to 40', '40 to 50', '50 to 60', '> 60')) +
  guides(color = guide_legend(override.aes = list(size = 10), reverse = TRUE)) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))

ens_pred_map



# Panel d -----------------


BREAKS <- c(100, 18, 10, 3, -3, -10, -18, -100)
cols <- colorRampPalette(brewer.pal(11, 'PiYG'))(length(BREAKS) - 1)


ge_pred_map <- ggplot() +
  #geom_path(data = usamap, aes(x = x, y = y), color = 'black', alpha = 0.7) +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  geom_point(data = predict_ensemble_df, inherit.aes = FALSE,
             aes(x = x, y = y, color = cut(diff, breaks = BREAKS)), 
             alpha = 0.8, size = 0.2) +
  scale_color_manual(bquote(Delta ~ bold('predicted abundance')), values = cols, 
                     labels = c('< -18', '-18 to -10', '-10 to -3', 
                                '-3 to 3', '3 to 10', '10 to 18', '> 18')) +
  guides(color = guide_legend(override.aes = list(size = 10), reverse = TRUE)) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))

ge_pred_map
