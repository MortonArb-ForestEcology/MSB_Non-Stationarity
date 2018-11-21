####################
# 3 - Fit global and regional models
#
# -Fit global glm
# -Fit regional glms
#
# Creates:
# * global_fit_df.rds -> fit summary for global model
# * regional_fit_df.rds -> fit summary for regional models
# * global_fit_obj.rds -> model obj for global model
# * regional_fit_obj -> model obj for regional model
####################



# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(MASS)


# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Example_BBS/Data/'))


# read in RDS -------------------------------------------------------------

BBS_bio <- readRDS('bbs_bio_BOX_M_5930_2018-10-30.rds')


# Compare models ------------------------------------------------------------

#distribution is zero inflated and overdispersed

#Actual data
hist(BBS_bio$target_count, xlim = c(0, 200), breaks = 20, ylim = c(0, 500))

#poisson
pois_fit_obj <- glm(target_count ~ bio1 + bio4 + bio12 + bio15, family = poisson, data = BBS_bio)

#negative binomial
nb_fit_obj <- MASS::glm.nb(target_count ~ bio1 + bio4 + bio12 + bio15, data = BBS_bio)

#hurdle gamma
non_zero <- ifelse(BBS_bio$target_count > 0, 1, 0)
d <- data.frame(BBS_bio, non_zero)
#binomial model
fit_bin <- glm(non_zero ~ bio1 + bio4 + bio12 + bio15, 
               data = d, family = binomial(link = logit))
#gamma model
fit_gamma <- glm(target_count ~ bio1 + bio4 + bio12 + bio15, 
                 data = subset(d, non_zero == 1), family = Gamma(link = log))

hg_LL <- logLik(fit_bin)[1] + logLik(fit_gamma)[1]

#AIC pois
-2 * logLik(pois_fit_obj)[1] + (2 * 5)
#AIC negbinom
-2 * logLik(nb_fit_obj)[1] + (2 * 5)
#AIC hurdle gamma - LOWEST
-2 * hg_LL + (2 * 11)



# Global model - hurdle gamma ------------------------------------------------------------

non_zero <- ifelse(BBS_bio$target_count > 0, 1, 0)
d <- data.frame(BBS_bio, non_zero)
#binomial model
fit_bin <- glm(non_zero ~ bio1 + bio4 + bio12 + bio15, 
          data = d, family = binomial(link = logit))
#gamma model
fit_gamma <- glm(target_count ~ bio1 + bio4 + bio12 + bio15, 
          data = subset(d, non_zero == 1), family = Gamma(link = log))

sum_bin <- summary(fit_bin)
sum_gamma <- summary(fit_gamma)



#slopes from count model
global_fit_df <- data.frame(n_survey = NROW(BBS_bio), 
                            n_pres = sum(BBS_bio$target_count > 0),
                            bin_int = sum_bin$coefficients[1,1], 
                            bin_beta_bio1 = sum_bin$coefficients[2,1], 
                            bin_beta_bio4 = sum_bin$coefficients[3,1],
                            bin_beta_bio12 = sum_bin$coefficients[4,1],
                            bin_beta_bio15 = sum_bin$coefficients[5,1],
                            gamma_int = sum_gamma$coefficients[1,1], 
                            gamma_beta_bio1 = sum_gamma$coefficients[2,1], 
                            gamma_beta_bio4 = sum_gamma$coefficients[3,1],
                            gamma_beta_bio12 = sum_gamma$coefficients[4,1],
                            gamma_beta_bio15 = sum_gamma$coefficients[5,1],
                            gamma_int_p = sum_gamma$coefficients[1,4], 
                            gamma_beta_bio1_p = sum_gamma$coefficients[2,4], 
                            gamma_beta_bio4_p = sum_gamma$coefficients[3,4],
                            gamma_beta_bio12_p = sum_gamma$coefficients[4,4],
                            gamma_beta_bio15_p = sum_gamma$coefficients[5,4])



# Regional models ---------------------------------------------------------

full_box_names <- grep('box_num', colnames(BBS_bio))

#make sure boxes have BBS data in them
relevant_boxes <- names(which(apply(BBS_bio[,full_box_names], 2, 
                                                 function(x) sum(x, na.rm = TRUE)) > 0))

N <- length(relevant_boxes)

regional_fit_df <- data.frame(box_num = rep(NA, N), 
                              n_survey = rep(NA, N), 
                              n_pres = rep(NA, N),
                              bin_int = rep(NA, N), 
                              bin_beta_bio1 = rep(NA, N), 
                              bin_beta_bio4 = rep(NA, N),
                              bin_beta_bio12 = rep(NA, N),
                              bin_beta_bio15 = rep(NA, N),
                              gamma_int = rep(NA, N), 
                              gamma_beta_bio1 = rep(NA, N), 
                              gamma_beta_bio4 = rep(NA, N),
                              gamma_beta_bio12 = rep(NA, N),
                              gamma_beta_bio15 = rep(NA, N),
                              gamma_int_p = rep(NA, N), 
                              gamma_beta_bio1_p = rep(NA, N), 
                              gamma_beta_bio4_p = rep(NA, N),
                              gamma_beta_bio12_p = rep(NA, N),
                              gamma_beta_bio15_p = rep(NA, N))

regional_bin_fit_obj <- list()
regional_gamma_fit_obj <- list()
nd_gamma <- c()
for (i in 1:length(relevant_boxes))
{
  #i <- 2
  
  #filter data for box
  t_col <- relevant_boxes[i]
  #t_col <- 'box_num_651'
  t_ind <- which(BBS_bio[t_col] == TRUE)
  t_BBS_bio <- BBS_bio[t_ind, 1:17]
  non_zero <- ifelse(t_BBS_bio$target_count > 0, 1, 0)
  temp_d <- data.frame(t_BBS_bio, non_zero)
  
  #proportion of detections
  prop <- round(sum(t_BBS_bio$target_count > 0)/length(t_BBS_bio$target_count), digits = 2)
  
  #binomial model
  temp_fit_bin <- glm(non_zero ~ bio1 + bio4 + bio12 + bio15,
                      data = temp_d, family = binomial(link = logit),
                      control = list(maxit = 500))

  #gamma model
  # need more than 5 data points, bc there are 5 params - 6 to be safe
  if (sum(t_BBS_bio$target_count > 0) <= 6)
  {
    temp_fit_gam <- NA
    nd_gamma <- c(nd_gamma, i)
  } else {
    temp_fit_gam <- glm(target_count ~ bio1 + bio4 + bio12 + bio15,
                   data = subset(temp_d, non_zero == 1), family = Gamma(link = log),
                   control = list(maxit = 500))
  }

  #insert model fit into list
  regional_bin_fit_obj[[i]] <- temp_fit_bin
  regional_gamma_fit_obj[[i]] <- temp_fit_gam

  #extract and save parameter estimates
  regional_fit_df$box_num[i] <- substr(t_col, start = 9, stop = nchar(t_col))
  regional_fit_df$n_survey[i] <- NROW(t_BBS_bio)
  regional_fit_df$n_pres[i] <- sum(t_BBS_bio$target_count > 0)

  if (length(temp_fit_bin) > 1)
  {
    temp_sm_bin <- summary(temp_fit_bin)

    regional_fit_df$bin_int[i] <- temp_sm_bin$coefficients[1,1]
    regional_fit_df$bin_beta_bio1[i] <- temp_sm_bin$coefficients[2,1]
    regional_fit_df$bin_beta_bio4[i] <- temp_sm_bin$coefficients[3,1]
    regional_fit_df$bin_beta_bio12[i] <- temp_sm_bin$coefficients[4,1]
    regional_fit_df$bin_beta_bio15[i] <- temp_sm_bin$coefficients[5,1]
  } else {
    regional_fit_df$bin_int[i] <- NA
    regional_fit_df$bin_beta_bio1[i] <- NA
    regional_fit_df$bin_beta_bio4[i] <- NA
    regional_fit_df$bin_beta_bio12[i] <- NA
    regional_fit_df$bin_beta_bio15[i] <- NA
  }

  if (length(temp_fit_gam) > 1)
  {
    temp_sm_gam <- summary(temp_fit_gam)

    regional_fit_df$gamma_int[i] <- temp_sm_gam$coefficients[1,1]
    regional_fit_df$gamma_beta_bio1[i] <- temp_sm_gam$coefficients[2,1]
    regional_fit_df$gamma_beta_bio4[i] <- temp_sm_gam$coefficients[3,1]
    regional_fit_df$gamma_beta_bio12[i] <- temp_sm_gam$coefficients[4,1]
    regional_fit_df$gamma_beta_bio15[i] <- temp_sm_gam$coefficients[5,1]
    regional_fit_df$gamma_int_p[i] <- temp_sm_gam$coefficients[1,4]
    regional_fit_df$gamma_beta_bio1_p[i] <- temp_sm_gam$coefficients[2,4]
    regional_fit_df$gamma_beta_bio4_p[i] <- temp_sm_gam$coefficients[3,4]
    regional_fit_df$gamma_beta_bio12_p[i] <- temp_sm_gam$coefficients[4,4]
    regional_fit_df$gamma_beta_bio15_p[i] <- temp_sm_gam$coefficients[5,4]
  } else {
    regional_fit_df$gamma_int[i] <- NA
    regional_fit_df$gamma_beta_bio1[i] <- NA
    regional_fit_df$gamma_beta_bio4[i] <- NA
    regional_fit_df$gamma_beta_bio12[i] <- NA
    regional_fit_df$gamma_beta_bio15[i] <- NA
    regional_fit_df$gamma_int_p[i] <- NA
    regional_fit_df$gamma_beta_bio1_p[i] <- NA
    regional_fit_df$gamma_beta_bio4_p[i] <- NA
    regional_fit_df$gamma_beta_bio12_p[i] <- NA
    regional_fit_df$gamma_beta_bio15_p[i] <- NA
  }
}



# save RDS objects --------------------------------------------------------

#write fit summary to RDS
saveRDS(global_fit_df, 'global_fit_df_M.rds')
saveRDS(regional_fit_df, 'regional_fit_df_M.rds')

#write fit obj to RDS
saveRDS(fit_bin, 'global_bin_fit_obj_M.rds')
saveRDS(fit_gamma, 'global_gamma_fit_obj_M.rds')
saveRDS(regional_bin_fit_obj, 'regional_bin_fit_obj_M.rds')
saveRDS(regional_gamma_fit_obj, 'regional_gamma_fit_obj_M.rds')

