library(dplyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(pROC)
library(ggplot2)
library(corrplot)
run <- c("big_creek","salmon_river", "elk_river")[3] # CHANGE ME!

cross_validation <- FALSE
ocean <- readRDS("ocean_data/fall_chinook_data.rds")
ocean <- ocean[complete.cases(ocean),] # we want no missing values

if(run == "big_creek") {
  results <- readRDS("results/fall_chinook_results_big_creek.rds")
  cfs <- readRDS("river_data/columbia_dalles_CFS.rds")
  cfs <- dplyr::filter(cfs, month == 4) %>%
    dplyr::rename(cfs_mean = CFS.mean)
}
if(run == "salmon_river") {
  results <- readRDS("results/fall_chinook_results_salmon_river.rds")
  cfs <- readRDS("river_data/siletz_august_CFS.rds") # no temp gauge in salmon use siletz (nestucca not a solid time series)
  cfs <-  dplyr::rename(cfs, cfs_mean = CFS.mean) # would like to have temperature too but not available until 2023
}

if(run == "elk_river") {
  results <- readRDS("results/fall_chinook_results_elk_river.rds")
  #No CFS for elk river fish
}

# results - filter out frowny faces
results <- dplyr::filter(results,
                         direction == "down")
ocean <- ocean[,c(1, which(names(ocean) %in% results$var))]

if(run == "big_creek" | run == "salmon_river") cfs$year <- as.numeric(cfs$year)

# bring in hatchery data

fish_dat <- readRDS("fish_data/chinook_by_survival_toage2.rds")

#filter by hatchery
if(run == "big_creek") fish <- dplyr::filter(fish_dat, stock_location_name == "BIG CR HATCHERY")
if(run == "salmon_river") fish <- dplyr::filter(fish_dat, stock_location_name == "SALMON R")
if(run == "elk_river") fish <- dplyr::filter(fish_dat, stock_location_name == "ELK R (ELK R HT)")

# define run year
if(run == "elk_river") fish$year <- fish$brood_year + 2
if(run == "big_creek" | run == "salmon_river") fish$year <- fish$brood_year + 1

fish <- dplyr::filter(fish, year %in% c(ocean$year))
fish <- dplyr::select(fish, brood_year, year, total_release, total_return,
                      avg_weight, avg_yday)

dat <- dplyr::left_join(fish, ocean)
dat<-dat[complete.cases(dat), ]

# make sure total return is an integer
dat$total_return <- round(dat$total_return)
dat$not_return <- dat$total_release - dat$total_return

# join cfs only for big creek and salmon river (not elk river)
if(run == "big_creek" | run == "salmon_river") dat <- dplyr::left_join(dat, cfs[,c("year","cfs_mean")])

# Simulated data and model fitting.
# One decision: what's the maximum number of covariates that any model can have?
if(run == "big_creek" | run == "salmon_river") covariates <- names(dat)[-which(names(dat) %in% c("brood_year","year","year_L1","number_cwt_estimated","total_return","cfs_mean","not_return"))]#setdiff(names(df), "y")
if(run == "elk_river") covariates <- names(dat)[-which(names(dat) %in% c("brood_year","year","year_L1","number_cwt_estimated","total_return","not_return"))]#setdiff(names(df), "y")

combinations <- lapply(1:3, function(i) {
  combn(covariates, i, simplify = FALSE)
})

combinations <- unlist(combinations, recursive = FALSE)
# Sort each combination and remove duplicates
combinations <- unique(lapply(combinations, function(x) sort(x)))

models <- list()
results <- data.frame()

for (i in seq_along(combinations)) {
  # ps here represents a P-spline / penalized regression spline
  # k represent the number of parameters / knots estimating function at, should be small
  smooth_terms <- paste("s(", combinations[[i]], ", k = 3)", collapse = " + ")
  if(run == "big_creek" | run == "salmon_river") smooth_terms <- paste(smooth_terms, " + s(cfs_mean,k=3)") # add cfs to model for big creek and salmon river
  formula_str <- paste("cbind(total_return, not_return) ~ ", smooth_terms)

  if(cross_validation) {
    predictions <- numeric(nrow(dat))
    n_year <- length(unique(dat$year))
    # Loop over each observation
    for (j in 1:n_year) {
      train_index <- setdiff(1:n_year, j)  # All indices except the j-th
      test_index <- j                 # The j-th index

      # Fit model on n-1 observations
      gam_model <- gam(as.formula(formula_str),
                       family = binomial(),
                       #weights = number_cwt_estimated,
                       data = dat[which(dat$year != unique(dat$year)[j]), ])

      # Predict the excluded observation
      predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
    }
  } else {
    # just fit the model once
    gam_model <- gam(as.formula(formula_str),
                     family = binomial(),
                     #weights = number_cwt_estimated,
                     data = dat)
    complete_dat <- which(complete.cases(dat)==TRUE)
    dat$predictions <- NA
    dat$predictions[complete_dat] <- as.numeric(predict(gam_model, type="response"))
  }
  dat$obs_survival <- dat$total_return / (dat$total_release)
  logit <- function(x) {return(log(x/(1-x)))}
  rmse <- sqrt(mean((logit(dat$obs_survival) - logit(dat$predictions))^2, na.rm=T))
  #pred_ROCR <- prediction(predictions, dat$jack)
  #auc_ROCR <- performance(pred_ROCR, measure = "auc")
  #auc <- auc_ROCR@y.values[[1]]

  # Extract variable names
  var_names <- gsub("s\\(([^,]+),.*", "\\1", combinations[[i]])
  # Store results with variable names padded to ensure there are always 3 columns
  padded_vars <- c(var_names, rep(NA, 3 - length(var_names)))

  # re-fit the model -- primarily for cross - validation case
  if(cross_validation==TRUE) {
    gam_model <- gam(as.formula(formula_str),
                     family = binomial(),
                     #weights = number_cwt_estimated,
                     data = dat)
  }
  # Store results
  models[[i]] <- gam_model
  results <- rbind(results, data.frame(
    ModelID = i,
    AIC = AIC(gam_model),
    RMSE = rmse,
    #AUC = auc,
    var1 = padded_vars[1],
    var2 = padded_vars[2],
    var3 = padded_vars[3]
  ))
  print(i)
}

# View the results dataframe
print(results)


# Create a baseline model with only the intercept
if(cross_validation) {
  predictions <- numeric(nrow(dat))
  n_year <- length(unique(dat$year))

  # Loop over each observation
  for (j in 1:n_year) {
    train_index <- setdiff(1:n_year, j)  # All indices except the j-th
    test_index <- j                 # The j-th index

    # Fit model on n-1 observations
    if(run == "big_creek" | run == "salmon_river") gam_model <- gam(cbind(total_return, not_return) ~ s(cfs_mean,k=3),
                     family = binomial(),
                     #weights = number_cwt_estimated,
                     data = dat[which(dat$year != unique(dat$year)[j]), ])
    if(run == "elk_river") gam_model <- gam(cbind(total_return, not_return) ~ 1,
                                            family = binomial(),
                                            #weights = number_cwt_estimated,
                                            data = dat[which(dat$year != unique(dat$year)[j]), ])
    # Predict the excluded observation
    predictions[which(dat$year == unique(dat$year)[j])] <- predict(gam_model, newdata = dat[which(dat$year == unique(dat$year)[j]), ])
  }
} else {
  if(run == "big_creek" | run == "salmon_river") gam_model <- gam(cbind(total_return, not_return) ~ s(cfs_mean,k=3),
                   family = binomial(),
                   #weights = number_cwt_estimated,
                   data = dat)
  if(run == "elk_river") gam_model <- gam(cbind(total_return, not_return) ~ 1,
                   family = binomial(),
                   #weights = number_cwt_estimated,
                   data = dat)

  complete_dat <- which(complete.cases(dat)==TRUE)
  dat$predictions <- NA
  dat$predictions[complete_dat] <- as.numeric(predict(gam_model, type="response"))
}
dat$obs_survival <- dat$total_return / (dat$total_release)
logit <- function(x) {return(log(x/(1-x)))}
baseline_rmse  <- sqrt(mean((logit(dat$obs_survival) - logit(dat$predictions))^2, na.rm=T))

# we can calculate the marginal improvement for each covariate
results$n_cov <- ifelse(!is.na(results$var1), 1, 0) + ifelse(!is.na(results$var2), 1, 0) +
  ifelse(!is.na(results$var3), 1, 0)
marginals <- data.frame(cov = covariates, "rmse_01" = NA, "rmse_12" = NA, "rmse_23" = NA,
                        "aic_01" = NA, "aic_12" = NA, "aic_23" = NA)
for(i in 1:length(covariates)) {
  sub <- dplyr::filter(results, n_cov == 1,
                       var1 == covariates[i])
  marginals$rmse_01[i] <- sub$RMSE / baseline_rmse
  marginals$aic_01[i] <- AIC(gam_model) - sub$AIC

  # next look at all values of models that have 2 covariates and include this model
  sub1 <- dplyr::filter(results, n_cov == 1)
  sub2 <- dplyr::filter(results, n_cov == 2) %>%
    dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0)) %>%
    dplyr::filter(keep == 1) %>% dplyr::select(-keep)
  # loop over every variable in sub2, and find the simpler model in sub1 that just represents the single covariate
  sub2$rmse_diff <- 0
  sub2$AIC_diff <- 0
  for(j in 1:nrow(sub2)) {
    vars <- sub2[j,c("var1", "var2")]
    vars <- vars[which(vars != covariates[i])]
    indx <- which(sub1$var1 == as.character(vars))
    sub2$rmse_diff[j] <- sub2$RMSE[j] / sub1$RMSE[indx]
    sub2$AIC_diff[j] <- sub1$AIC[indx] - sub2$AIC[j]
  }

  # Finally compare models with 3 covariates to models with 2 covariates
  sub2_all <- dplyr::filter(results, n_cov == 2)
  # Apply a function across the rows to sort the values in var1 and var2
  sorted_names <- t(apply(sub2_all[, c("var1", "var2")], 1, function(x) sort(x)))
  # Replace the original columns with the sorted data
  sub2_all$var1 <- sorted_names[, 1]
  sub2_all$var2 <- sorted_names[, 2]

  sub3 <- dplyr::filter(results, n_cov == 3) %>%
    dplyr::mutate(keep = ifelse(var1 == covariates[i],1,0) + ifelse(var2 == covariates[i],1,0) + ifelse(var3 == covariates[i],1,0)) %>%
    dplyr::filter(keep == 1) %>% dplyr::select(-keep)
  sub3$rmse_diff <- 0
  sub3$AIC_diff <- 0
  for(j in 1:nrow(sub3)) {
    vars <- sub3[j,c("var1", "var2","var3")]
    vars <- sort(as.character(vars[which(vars != covariates[i])]))
    # find the same in sub2
    indx <- which(paste(sub2_all$var1, sub2_all$var2) == paste(vars, collapse=" "))
    sub3$rmse_diff[j] <- sub3$RMSE[j] / sub2_all$RMSE[indx]
    sub3$AIC_diff[j] <- sub2_all$AIC[indx] - sub3$AIC[j]
  }

  # Fill in summary stats
  marginals$rmse_12[i] <- mean(sub2$rmse_diff)
  marginals$aic_12[i] <- mean(sub2$AIC_diff)
  marginals$rmse_23[i] <- mean(sub3$rmse_diff)
  marginals$aic_23[i] <- mean(sub3$AIC_diff)
}

# Calculate avergages of averages
marginals$total_rmse <- apply(marginals[,c("rmse_01","rmse_12", "rmse_23")], 1, mean)
marginals$total_aic <- apply(marginals[,c("aic_01","aic_12", "aic_23")], 1, mean)

marginals$null_rmse <- marginals[,c("rmse_01")]
table <- dplyr::arrange(marginals, null_rmse) %>%
  dplyr::select(cov, total_rmse)

# CREATE TABLE OF RMSE/AIC------------------------------------

table <- dplyr::arrange(marginals, total_rmse) %>%
  dplyr::select(cov, total_rmse, total_aic)

table <- dplyr::arrange(marginals[,c("rmse_01")]) %>%
  dplyr::select(cov, total_rmse, total_aic)

# Fit a model that includes multiple variables from the top model
if(run == "big_creek" | run == "salmon_river") combos = c("cfs_mean", table$cov[1], table$cov[2], table$cov[3], table$cov[4], table$cov[5])
if(run == "elk_river") combos = c(table$cov[1], table$cov[2], table$cov[3], table$cov[4], table$cov[5], table$cov[6])

smooth_terms <- paste("s(", combos, ", k = 3)", collapse = " + ")
formula_str <- paste("cbind(total_return, not_return) ~ ", smooth_terms)
gam_model <- gam(as.formula(formula_str),
                 family = binomial(),
                 #weights = number_cwt_estimated,
                 data = dat)

# Make a plot of the marginal effects from this best model
partial_effects <- draw(gam_model, transform = "response")

ggsave(partial_effects, filename = paste0("plots/fall_chinook/fall_chinook_",run,"_gam_partial_effects.png"), height = 7, width = 10)

# USE Z-SCORE TO CALCULATE COVARIATE EFFECT SIZE --------------------------

table$avg_absolute_effect <- 0
if(run == "big_creek" | run == "salmon_river") n_cov <- 5
if(run == "elk_river") n_cov <- 6
for(i in 1:n_cov) {
  # prediction data frame
  pred_df <- 0 * dat[1:3,which(names(dat) %in% table$cov[1:n_cov])]
  if(run == "big_creek" | run == "salmon_river") pred_df$cfs_mean <- mean(dat$cfs_mean,na.rm=T) # hold at mean, not standardized
  if("avg_weight" %in% names(pred_df)) {
    pred_df[,which(names(pred_df) == "avg_weight")] <- mean(dat$avg_weight,na.rm=T)
  }
  if("avg_yday" %in% names(pred_df)) {
    pred_df[,which(names(pred_df) == "avg_yday")] <- mean(dat$avg_yday,na.rm=T)
  }
  # throught the loop change which variable we're calculating effect sizes for
  pred_df[,which(names(pred_df) == table$cov[i])] = c(-1,0,1)
  # predictions
  pred_df$est_survival <- predict(gam_model, newdata = pred_df, type = "response")

  # average absolute change in survival given a 1 change sd in the covariate
  table$avg_absolute_effect[i] <- mean(abs(diff(pred_df$est_survival)))

}
if(run == "big_creek") saveRDS(table, "results/fall_chinook_results_table_big creek.rds")
if(run == "salmon_river") saveRDS(table, "results/fall_chinook_results_table_salmon_river.rds")
if(run == "elk_river") saveRDS(table, "results/fall_chinook_results_table_elk_river.rds")

# MAKE CORRELATION PLOT OF TOP AIC COVARIATES -----------------------------

# Check out correlation for variables with expected relationship
v <- table$cov[1:6]
ocean_simp <- ocean[which(names(ocean) %in% v)]
cc<-cor(ocean_simp)


if(run == "big_creek") file_path= "plots/fall_chinook/fall_chinook_big_creek_corrplot.png"
if(run == "salmon_river") file_path= "plots/fall_chinook/fall_chinook_salmon_river_corrplot.png"
if(run == "elk_river") file_path= "plots/fall_chinook/fall_chinook_elk_river_corrplot.png"

png(height=500, width=500, file=file_path)
corrplot(cc, method = 'number', type = 'lower', diag = FALSE, tl.col="black", number.cex=1)
dev.off()


