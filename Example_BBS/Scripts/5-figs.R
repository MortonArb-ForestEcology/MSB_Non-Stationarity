####################
# 5 - create figs
#
# -Calculate residuals
# -NOCA range with boxes on map
# -Slope estimates for BIO15 covariate on map
# -Ensemble predictions on map
# -Global - ensemble predictions on map
# 
# Creates:
# Figures
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/Project_archive/'


# Load packages -----------------------------------------------------------

library(rgeos)
library(ggplot2)
library(sp)
library(RColorBrewer)
library(raster)
library(stringr)

# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'BBS_ns/Data/'))


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
setwd(paste0(dir, 'BBS_ns/Data/shp_files'))
sp_rng <- rgdal::readOGR('Cardinalis_cardinalis_22723819.shp', verbose = FALSE)
usa_4326 <- rgdal::readOGR('USA_4326.shp', verbose = FALSE)

sp_rng_usa <- raster::intersect(usa_4326, sp_rng)

sp_rng_usa@data$id <- rownames(sp_rng_usa@data)
sp_rng_usa_pts <- ggplot2::fortify(sp_rng_usa, region = 'id')
sp_rng_usa_df <- merge(sp_rng_usa_pts, sp_rng_usa@data, by = 'id')


usa_m <- maps::map('usa', plot = FALSE, fill = TRUE)
IDs <- sapply(strsplit(usa_m$names, ":"), function(x) x[1])
usa <- maptools::map2SpatialPolygons(usa_m, IDs = IDs, proj4string = sp::CRS('+init=epsg:4326'))




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





# Calculate residuals -----------------------------------------------------

cell_cnt <- data.frame(event_id = BBS_bio$event_id, 
                       bio1 = BBS_bio$bio1, bio4 = BBS_bio$bio4,
                       bio12 = BBS_bio$bio12, bio15 = BBS_bio$bio15,
                       lng = BBS_bio$lng, lat = BBS_bio$lat, 
                       target_count = BBS_bio$target_count, cell = NA,
                       predict_global = NA, ensemble = NA)


#BBS spatial points
BBS_sp_pts <- sp::SpatialPoints(cbind(cell_cnt$lng, cell_cnt$lat))

#Bioclim spatial points
BC_sp_pts <- sp::SpatialPoints(cbind(predict_ensemble_df_p$x, 
                                     predict_ensemble_df_p$y))

#specify CRS
sp::proj4string(BBS_sp_pts) <- sp::CRS('+init=epsg:4326')
sp::proj4string(BC_sp_pts) <- sp::CRS('+init=epsg:4326')

#transform to projected
tr_BBS <- sp::spTransform(BBS_sp_pts, sp::CRS("+init=epsg:3857"))
tr_BC <- sp::spTransform(BC_sp_pts, sp::CRS("+init=epsg:3857"))

#index of closest bioclim cell to BBS site
cl_bc <- apply(rgeos::gDistance(tr_BBS, tr_BC, byid = TRUE), 2, which.min)

#fill df
cell_cnt$cell <- predict_ensemble_df_p[cl_bc, ]$cell
cell_cnt$predict_global <- predict_ensemble_df_p[cl_bc, ]$predict_global
cell_cnt$ensemble <- predict_ensemble_df_p[cl_bc, ]$ensemble

#residuals df
resid_df <- data.frame(cell_cnt, 
                       resid_global = cell_cnt$target_count - cell_cnt$predict_global, 
                       resid_ensemble = cell_cnt$target_count - cell_cnt$ensemble)

#calculate RSS for global and ensemble
RSS_global <- sum((resid_df$resid_global^2))
RSS_ensemble <- sum((resid_df$resid_ensemble^2))

#amount of reduction in RSS by using ensemble model compared to global model
(RSS_global - RSS_ensemble) / RSS_global



# difference in residuals -------------------------------------------------

#positive # is ensemble better, neg # is global model better
err <- abs(resid_df$resid_global) - abs(resid_df$resid_ensemble)

#krig difference in resid
fit1 <- fields::spatialProcess(resid_df[,c('lng', 'lat')], err)
predict_ensemble_df$pred_fit1 <- predict(fit1, predict_ensemble_df[,c('x','y')])




# Panel a -------------------------------------------------

#convert boxes to df
box_spdf@data$id <- rownames(box_spdf@data)
box_pts <- ggplot2::fortify(box_spdf, region = 'id')
box_df <- merge(box_pts, box_spdf@data, by = 'id')

#sample 75 boxes
set.seed(7)
box_id_s <- sample(box_df$id, size = 50)
box_df_s <- dplyr::filter(box_df, id %in% box_id_s)



box_map <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  geom_polygon(data = sp_rng_usa_df,
               aes(x = long, y = lat, group=group), fill = 'orange', alpha = 0.8) +
  geom_path(data = usa_m, aes(x = long, y = lat, group = group), 
            color = 'black', alpha = 0.7) +
  geom_path(data = box_df_s, aes(x = long, y = lat, group = group), color = 'black',
            alpha = 0.5, lwd = 0.8) +
  geom_point(data = BBS_bio, aes(x = lng, y = lat), color = 'black',
            alpha = 0.4, size = 1) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12))


#box_map

setwd(paste0(dir, 'BBS_ns/Figs/'))
ggsave('panel_1.png', box_map)


# Panel b ----------------------------------------

hist(plt_rf$gamma_beta_bio15)
BREAKS <- c(-1, -0.2, -0.1, -0.03, 0.03, 0.1, 0.2, 1)
cols <- colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(BREAKS)-1)


bio15_map <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  geom_point(data = plt_rf, inherit.aes = FALSE,
             aes(x = cen_x, y = cen_y, color = cut(gamma_beta_bio15, breaks = BREAKS)), 
             alpha = 0.5, size = 4) +
  # scale_color_manual('', values = rev(cols),
  #                    labels = c('< -0.2', '-0.2 to -0.1', '-0.1 to -0.03',
  #                               '-0.03 to 0.03',
  #                               '0.03 to 0.1', '0.1 to 0.2', '> 0.2')) +
  scale_color_manual('Effect of precipitation seasonality', values = rev(cols),
                     labels = c('< -0.2', '-0.2 to -0.1', '-0.1 to -0.03',
                                '-0.03 to 0.03',
                                '0.03 to 0.1', '0.1 to 0.2', '> 0.2')) +
  guides(color = guide_legend(override.aes = list(size = 10), reverse = TRUE)) +
  geom_path(data = usa_m, aes(x = long, y = lat, group = group), 
            color = 'black', alpha = 0.7) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))

#bio15_map

ggsave('panel_2.png', bio15_map)


# Panel c ---------------------------------------------

hist(predict_ensemble_df$ensemble)
BREAKS <- c(-500, 10, 20, 30, 40, 50, 60, 500)
cols <- colorRampPalette(brewer.pal(9, 'Reds'))(length(BREAKS) - 1)

ens_pred_map <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  geom_point(data = predict_ensemble_df, inherit.aes = FALSE,
             aes(x = x, y = y, color = cut(ensemble, breaks = BREAKS)), 
             alpha = 0.8, size = 0.2) +
  scale_color_manual(str_wrap('Predicted abundance', 10), values = cols,
                     labels = c('< 10', '10 to 20', '20 to 30',
                                '30 to 40', '40 to 50', '50 to 60', '> 60')) +
  guides(color = guide_legend(override.aes = list(size = 10), reverse = TRUE)) +
  geom_path(data = usa_m, aes(x = long, y = lat, group = group), 
            color = 'black', alpha = 0.7) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))

#ens_pred_map

ggsave('panel_3.png', ens_pred_map)



# Panel d -----------------


BREAKS <- c(-1000, -3, -2, -1, 1, 2, 3, 1000)
cols <- colorRampPalette(brewer.pal(11, 'PiYG'))(length(BREAKS) - 1)

#positive # is ensemble better, neg # is global model better

d_err_map <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = rgb(0,0,0, 0.2)) +
  coord_map("ortho", orientation = c(35, -100, 0), 
            xlim = c(-125, -70), ylim = c(23, 50)) +
  geom_point(data = predict_ensemble_df, inherit.aes = FALSE,
             aes(x = x, y = y, color = cut(pred_fit1, breaks = BREAKS)), 
             alpha = 0.8, size = 0.2) +
  scale_color_manual(expression(atop(bold('|global err| - |ensemble err|'), ), 10),
                     # scale_color_manual(expression(atop(bold(''), ), 10), 
                     values = cols,
                     labels = c('< -3', '-3 to -2',
                                '-2 to -1', '-1 to 1', 
                                '1 to 2', '2 to 3', '> 3')) +
  guides(color = guide_legend(override.aes = list(size = 10), reverse = TRUE)) +
  geom_path(data = usa_m, aes(x = long, y = lat, group = group), 
            color = 'black', alpha = 0.7) +
  xlab('Longitude') +
  ylab('Latitude') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))

ggsave('panel_4.png', d_err_map)
