####################
# 2 - Process data
#
# -Merge BBS and bioclim data (just range over NOCA)
# -Use particle filtering and dggridR to create hex grid
# -Select points in each hex proportial to area that contains US
# -Create boxes 'spatial support sets' for each point (Make sure each box has 25 BBS locations in it)
# -Get bioclim variables for each pixel within each box (AKA 'spatial support set')
#
# Creates: 
# * bbs_bio_BOX_5930_DATE.rds -> bbs data merged with bioclim data at survey points and box #
# * BOX_spdf.rds -> Spatial Polygon Dataframe for boxes
# * cell_bioclim.rds -> bioclim measures at each cell for every box
####################



# specify box size and N BBS surveys in box ----------------------------------

# Medium support set box from Fink et al. 2010 - 6 degrees lat x 9 degrees lon
lat_dim <- 6
lon_dim <- 9

# Mininum number of BBS sites for each box
NBBS <- 25


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dggridR)
library(sp)
library(maptools)
library(maps)
library(dplyr)
library(rgdal)
library(raster)


# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Example_BBS/Data/'))



# read in BBS data --------------------------------------------------------


BBS_query <- readRDS('bbs_query_5930_2018-10-30.rds')



# create spatial polygons object of US ------------------------------------

usa_m <- maps::map('usa', plot = FALSE, fill = TRUE)
IDs <- sapply(strsplit(usa_m$names, ":"), function(x) x[1])
usa <- maptools::map2SpatialPolygons(usa_m, IDs = IDs, proj4string = sp::CRS('+init=epsg:4326'))



# merge with bioclim ------------------------------------------------------

# BIO1 = Annual Mean Temperature
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO12 = Annual Precipitation
# BIO15 = Precipitation Seasonality (Coefficient of Variation)

#BIO 1, BIO 4, BIO 12, BIO 15
bioclim_full <- raster::getData('worldclim', var = 'bio', res = 2.5)
bioclim_sub <- bioclim_full[[c(1,4,12,15)]]

#read in shp file of cardinal range
NOCA_rng <- rgdal::readOGR('shp_files/Cardinalis_cardinalis_22723819.shp')

#extract cell values at bbs surveys from rasters
bbs_sp_pts_p <- sp::SpatialPoints(cbind(BBS_query$lng, BBS_query$lat),
                              proj4string = sp::CRS('+init=epsg:4326'))

NOCA_rng_sp <- sp::SpatialPolygons(NOCA_rng@polygons)
sp::proj4string(NOCA_rng_sp) <- sp::CRS('+init=epsg:4326')

#data just over range of NOCA
mrg_ind <- which(!is.na(sp::over(bbs_sp_pts_p, NOCA_rng_sp, fn = NULL)))
bbs_sp_pts_p2 <- bbs_sp_pts_p[mrg_ind]
bbs_ev_id_p <- BBS_query$event_id[mrg_ind]

#data just over continental US
mrg_ind2 <- which(!is.na(sp::over(bbs_sp_pts_p2, usa, fn = NULL)))
bbs_sp_pts <- bbs_sp_pts_p2[mrg_ind2]
bbs_ev_id_p2 <- bbs_ev_id_p[mrg_ind2]
bbs_d <- filter(BBS_query, event_id %in% bbs_ev_id_p2)

bio1_bbs <- raster::extract(bioclim_sub[[1]], bbs_sp_pts)
bio4_bbs <- raster::extract(bioclim_sub[[2]], bbs_sp_pts)
bio12_bbs <- raster::extract(bioclim_sub[[3]], bbs_sp_pts)
bio15_bbs <- raster::extract(bioclim_sub[[4]], bbs_sp_pts)

BBS_bio <- data.frame(bbs_d, 
                      bio1 = bio1_bbs, 
                      bio4 = bio4_bbs, 
                      bio12 = bio12_bbs, 
                      bio15 = bio15_bbs)


# plot(bioclim_sub[[3]], xlim = c(-120, -50), ylim = c(20, 60))
# points(BBS_bio$lng, BBS_bio$lat, col = 'red', pch ='.')
# plot(usa, add = TRUE)



# distribute points uniformly on sphere -----------------------------------

#from: http://mathworld.wolfram.com/SpherePointPicking.html
set.seed(1)

N <- 10000000    #How many cells to sample
u     <- runif(N)
v     <- runif(N)
theta <- 2*pi*u      * 180/pi
phi   <- acos(2*v-1) * 180/pi
lon   <- theta-180
lat   <- phi-90
df    <- data.frame(lat = lat, lon = lon)



# convert to spatial points - WGS84 ---------------------------------------

pts <- sp::SpatialPoints(cbind(df$lon, df$lat), 
                         proj4string = sp::CRS('+init=epsg:4326'))




# particle filter to select points inside US ------------------------------

nn <- which(!is.na(sp::over(pts, usa)))
npts <- pts[nn]
ndf <- df[nn,]



# construct hex cells using points in US ----------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)
cells <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                  in_lon_deg = ndf$lon, in_lat_deg = ndf$lat)[[1]]
cell_grid <- dggridR::dgcellstogrid(hexgrid6, cells)



# sample points within each hex cell --------------------------------------

#number of points to sample from each box, scaled by area (i.e., # of boxes to create from each hex)
N <- 15

#write grid to shp file and read back in
dggridR::dgcellstogrid(hexgrid6, cells, savegrid = getwd())
dgg <- rgdal::readOGR('dggrid.shp')

raster::crs(dgg) <- raster::crs(npts)

# plot(dgg, add = TRUE)

#randomly sample N points from each hex cell - store in list
dgg2 <- as(dgg, 'SpatialPolygons')

#determine which fraction of points are in each cell (psuedo-area in each cell)
sizes <- c()
for (i in 1:length(dgg2))
{
  temp <- sp::over(npts, dgg2[i,])
  t_npts <- npts[!is.na(temp),]
  
  #proportion total points
  sizes <- c(sizes, length(t_npts)/length(npts))
}

#scale between 0 and 1
range01 <- function(x){(x - min(x))/(max(x) - min(x))}
#psuedo-area
ps_area <- range01(sizes)


#sample points from cells uniformly
set.seed(1)
box_centers <- c()
for (i in 1:length(dgg2))
{
  #i <- 10
  temp <- sp::over(npts, dgg2[i,])
  t_npts <- npts[!is.na(temp),]
  
  #scale sampling by area
  N_area <- ceiling(ps_area[i]*10)
  
  if (NROW(t_npts) >= N_area)
  {
    N_SAMPLE <- N_area
  } else {
    N_SAMPLE <- NROW(t_npts)
  }
  
  s_pts <- dplyr::sample_n(data.frame(t_npts@coords), N_SAMPLE, replace = FALSE)
  box_centers <- rbind(box_centers, s_pts)
  # points(s_pts, col = 'red', pch = '.')
}



# create boxes for each sampled point -------------------------------------

#get box coordinates

#lon
b_left <- round(box_centers[,1] - lon_dim/2, digits = 4)
b_right <- round(box_centers[,1] + lon_dim/2, digits = 4)

#lat
b_top <- round(box_centers[,2] - lat_dim/2, digits = 4)
b_bottom <- round(box_centers[,2] + lat_dim/2, digits = 4)

p_df <- data.frame(b_left, b_right, b_top, b_bottom)


box_col_names <- c()

poly_list <- list()
counter <- 1
for (i in 1:NROW(p_df))
{
  #i <- 1
  x <- p_df[i,]
  x_coords <- c(rep(x$b_top, 2), rep(x$b_bottom, 2))
  y_coords <- c(x$b_left, x$b_right, x$b_right, x$b_left)
  
  xym <- cbind(y_coords, x_coords)
  names(xym) <- c('lat', 'lon')
  
  p <- sp::Polygon(xym)
  ps <- sp::Polygons(list(p), 1)
  sps <- sp::SpatialPolygons(list(ps))
  proj4string(sps) <- sp::CRS('+init=epsg:4326')
  
  
  BBS_box_ind <- which(!is.na(over(bbs_sp_pts, sps)))
  if (length(BBS_box_ind) >= NBBS)
  {
    bcn <- paste0('box_num_', i)
    BBS_bio[BBS_box_ind, bcn] <- TRUE
    poly_list[[counter]] <- sps
    counter <- counter + 1
    #plot(sps, add = TRUE, col = rgb(1,0,0,0.2))
    
    box_col_names <- c(box_col_names, bcn)
  }
}

#merge polygons together into one Spatial Polygon
polygon_sp <- do.call(bind, poly_list)

#support set spdf
box_spdf <- sp::SpatialPolygonsDataFrame(Sr = polygon_sp, 
                                          data = data.frame(box = box_col_names))



# save as RDS objects -----------------------------------------------------

saveRDS(BBS_bio, 'bbs_bio_BOX_M_5930_2018-10-30.rds')
saveRDS(box_spdf, 'BOX_M_spdf.rds')



# bioclim variables inside each box ---------------------------------------

#mask for only US (study area)
bioclim_sub_mask <- raster::mask(bioclim_sub, usa)
bioclim_NOCA_sub_mask <- raster::mask(bioclim_sub_mask, NOCA_rng_sp)

#extract bioclim values for all cells in study area
extr_bio_p <- raster::extract(bioclim_NOCA_sub_mask, usa, cellnumbers = TRUE)
#bind list together (created bc usa has multiple features)
extr_bio <- do.call('rbind', extr_bio_p)
#get cell coordinates
extr_coord <- data.frame(extr_bio, xyFromCell(bioclim_NOCA_sub_mask[[1]], extr_bio[,1]))


#add box name columns
cell_bioclim_p <- extr_coord
cell_bioclim_p[box_col_names] <- FALSE


for (i in 1:length(box_spdf))
{
  # i <- 1
  #get cell values with associated cell_numbers
  cell_nums <- raster::extract(bioclim_NOCA_sub_mask, box_spdf[i,], cellnumbers = TRUE)[[1]][,1]
  
  cell_ind <- which(cell_bioclim_p$cell %in% cell_nums)
  
  #add box columns and fill TRUE if on that box
  cell_bioclim_p[cell_ind, box_col_names[i]] <- TRUE
}

#remove rows where bioclim values are NA (out of NOCA range)
cell_bioclim <- cell_bioclim_p[which(!is.na(cell_bioclim_p$bio1)),]

# raster::plot(box_spdf[1,], add = TRUE)

saveRDS(cell_bioclim, 'cell_bioclim_M.rds')

