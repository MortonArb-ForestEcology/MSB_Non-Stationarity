####################
# 4 - predict counts at each cell based on global and regional model fits
#
# -Predict # for each BIOclim cell across US and within each box
#
# Creates:
# * predict_df.rds -> predicted counts for each cell based on global and regional models
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Example_BBS/Data/'))


# read in RDS -------------------------------------------------------------

#bioclim data for each cell
cell_bioclim <- readRDS('cell_bioclim_M.rds')

#global and regional fit summaries
global_fit_df <- readRDS('global_fit_df_M.rds')
regional_fit_df <- readRDS('regional_fit_df_M.rds')

#global and regional fit objects
global_bin_fit_obj <- readRDS('global_bin_fit_obj_M.rds')
global_gamma_fit_obj <- readRDS('global_gamma_fit_obj_M.rds')
regional_bin_fit_obj <- readRDS('regional_bin_fit_obj_M.rds')
regional_gamma_fit_obj <- readRDS('regional_gamma_fit_obj_M.rds')




# create prediction df ----------------------------------------------------

#same format as bioclim df
predict_df <- cell_bioclim[,1:7]

#get col names for boxes
cn <- paste0('predict_', colnames(cell_bioclim)[grep('box_num', colnames(cell_bioclim))])

#add slots for prediction for each box and global
predict_df['predict_global'] <- NA
predict_df[cn] <- NA



# predict count in each cell (GLOBAL) -------------------------------------

#use bioclim data in each cell to predict count and add predictions to predict df

#predictions for binomial model
pred_bin <- predict(global_bin_fit_obj, 
                    newdata = predict_df,
                    type = "response")

#predictinos for gamma model
pred_gamma <- predict(global_gamma_fit_obj, 
                      newdata = predict_df,
                      type = "response")

#multiply to get mean
predict_df$predict_global <- pred_bin * pred_gamma




# predict count in each cell (REGIONAL) -----------------------------------

col_box_num <- colnames(cell_bioclim)[grep('box_num', colnames(cell_bioclim))]

for (i in 1:length(cn))
{
  #i <- 39
  #which indices correspond to that box (i.e., which cells are in which box)
  t_ind <- which(cell_bioclim[col_box_num[i]] == TRUE)
  
  cb_filt <- cell_bioclim[t_ind, 1:7]
  
  if (length(regional_gamma_fit_obj[[i]]) > 1)
  {
    #predictions for binomial model
    t_pred_bin <- predict(regional_bin_fit_obj[[i]], 
                          newdata = cb_filt,
                          type = "response")
    
    #predictinos for gamma model
    t_pred_gamma <- predict(regional_gamma_fit_obj[[i]], 
                          newdata = cb_filt,
                          type = "response")
    
    #multiply to get mean
    predict_df[t_ind, cn[i]] <- t_pred_bin * t_pred_gamma
  } else {
    predict_df[t_ind, cn[i]] <- NA
  }
}



# save to RDS -------------------------------------------------------------

saveRDS(predict_df, 'predict_df_M.rds')

