## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE, warning=FALSE, out.width="800px"-----------------------
knitr::include_graphics("https://github.com/mrecos/klrfome/blob/master/README_images/KLRfome_dataflow.png?raw=true")

## ----packages, message=FALSE, warning=FALSE, paged.print=FALSE-----------
library("ggplot2")   # for plotting results
library("NLMR")      # for creating simulated landscapes
library("rasterVis") # for plotting simulated lan
library("pROC")      # for evaluation of model AUC metric
library("dplyr")     # for data manipulation
library("knitr")     # for printing tables in this document
library("klrfome")   # for modeling

## ----params--------------------------------------------------------------
#Parameters
set.seed(232)
sigma = 0.5
lambda = 0.1
dist_metric = "euclidean"

## ----sim_data------------------------------------------------------------
### Simulate Training Data
sim_data <- get_sim_data(site_samples = 800, N_site_bags = 75)
formatted_data <- format_site_data(sim_data, N_sites=10, train_test_split=0.8,
                                   sample_fraction = 0.9, background_site_balance=1)
train_data <- formatted_data[["train_data"]]
train_presence <- formatted_data[["train_presence"]]
test_data <- formatted_data[["test_data"]]
test_presence <- formatted_data[["test_presence"]]

## ----fit_model-----------------------------------------------------------
##### Logistic Mean Embedding KRR Model
#### Build Kernel Matrix
K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric, progress = FALSE)
#### Train KLR model
train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 2)
#### Predict KLR model on test data
test_log_pred <- KLR_predict(test_data, train_data, dist_metric = dist_metric,
                             train_log_pred[["alphas"]], sigma, progress = FALSE)

### Plot K Matrix
K_corrplot(K,train_data,clusters=4)

### Plot Test Set Prediction
predicted_log <- data.frame(pred = test_log_pred, obs = test_presence)
ggplot(predicted_log, aes(x = as.factor(obs), y = pred, color = as.factor(obs))) +
  geom_jitter(width = 0.1) +
  theme_bw() +
  ylim(c(0,1)) +
  labs(y = "Predicted Probability", x = "Site Presence",
       title = "Kernel Logistic Regression",
       subtitle = "test set predictions; simulated data") +
  theme(
    legend.position = "none"
  )

### Save parameters for later prediction
params <- list(train_data = train_data,
               alphas_pred = train_log_pred[["alphas"]],
               sigma = sigma,
               lambda = lambda,
               means = formatted_data$means,
               sds = formatted_data$sds)

## ----echo=FALSE, warning=FALSE, out.width="800px"------------------------
knitr::include_graphics("https://github.com/mrecos/klrfome/blob/master/README_images/KLRfome_prediction.png?raw=true")

## ----predict_rasters-----------------------------------------------------
### width and hieght of roving focal window (required)
ngb = 5
### Number of rows and columns in prediction rasters
## needed for making simulated rasters, as well as for predicting real-world rasters
cols = 100
rows = 100

### Create simulated environmental rasters  (sim data only) ####
s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
s_var1 <- rescale_sim_raster(s_var1r, 50, 10) 
s_var2 <- rescale_sim_raster(s_var1r, 3, 2) 
b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
b_var1 <- rescale_sim_raster(b_var1r, 100, 20) 
b_var2 <- rescale_sim_raster(b_var1r, 6, 3) 
### Create a site-present trend surface  (sim data only)
trend_coords <- sim_trend(cols, rows, n = 3)
coords <- trend_coords$coords
trend <- trend_coords$trend
inv_trend <- abs(1-trend)
var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#### end simulated data creation ####

### Create raster stack of predictor variables
pred_var_stack <- raster::stack(var1, var2)
names(pred_var_stack) <- c("var1","var2")
### scale rasters to training data
pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
### Predict raster (single chunk, not in parallel) 
pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
                                progress = FALSE, parallel = FALSE)
### plot with simulated sites
rasterVis::levelplot(pred_rast, margin = FALSE, par.settings=viridisTheme()) +
  layer(sp.points(sp.points(SpatialPoints(coords), pch=15, cex = 2.25, col = "red")), columns=1)


## ----multi-proc----------------------------------------------------------
library("doParallel")
### create and register parallel backend
cl <- makeCluster(detectCores())
doParallel::registerDoParallel(cl)

### Use same KLR_raster_predict function with parallel = TRUE
pred_rast_list <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = TRUE, ppside = 5,
                                   progress = FALSE, parallel = TRUE, output = "list",
                                   save_loc = NULL, overwrite = TRUE, cols = cols, rows = rows)
### Merge list back to a single raster
pred_rast <-  do.call(merge, pred_rast_list)
### plot with simulated sites
rasterVis::levelplot(pred_rast, margin = FALSE, par.settings=viridisTheme()) +
  layer(sp.points(sp.points(SpatialPoints(coords), pch=15, cex = 2.25, col = "red")), columns=1)

### Or set output = "save" to save each prediction block out to a folder as a GeoTiff # not run
# pred_rast_list <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = TRUE, ppside = 5,
#                                    progress = FALSE, parallel = TRUE, output = "save",
#                                    save_loc = "c:/Temp/tif", overwrite = TRUE)

stopCluster(cl)

## ----model_evla----------------------------------------------------------
### Make some polygons around the simulated site points.
### If all you have is points for sites, site radius can be an assumption
site_pnts  <- SpatialPoints(coords)
site_polys <- rgeos::gBuffer(site_pnts, width = 6, byid = FALSE)

### extract sensitivity raster values to site areas
site_sample <- raster::extract(pred_rast, site_polys, weights = FALSE, 
                               small = TRUE, df = TRUE) %>%
  rename(pred = layer) %>%
  mutate(presence = 1)
### sample for an environmental background of sensitivity values. (e.g. n = 500)
bkg_sample <- data.frame(ID = 0, pred = sampleRandom(pred_rast, 500),
                         presence = 0)
model_pred <- rbind(site_sample, bkg_sample)

### A vector of the sensitivity thresholds that you want to evaluate the model at
threshold <- seq(0,1,0.1)
### Compute True Positive, True Negative, False Positive, and False Negative values at each threshold
kstats <- CM_quads(model_pred, threshold)

### use the pROC::auc and klrfome::metrics functions to compute the metrics of choice at each threshold
Test_area_metrics <- kstats %>%
  group_by(Threshold) %>%
  dplyr::mutate(AUC = round(pROC::auc(model_pred$presence, model_pred$pred, type = "linear"),3),
                YoudensJ = round(metrics(TP,TN,FP,FN)$Informedness,3),
                KG       = round(metrics(TP,TN,FP,FN)$KG,3),
                Sensitivity = round(metrics(TP,TN,FP,FN)$Sensitivity,3),
                FPR = round(metrics(TP,TN,FP,FN)$FPR,3),
                FNR = round(metrics(TP,TN,FP,FN)$FNR,3)) %>%
  data.frame()

knitr::kable(Test_area_metrics)

