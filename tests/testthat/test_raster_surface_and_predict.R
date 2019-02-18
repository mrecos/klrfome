context("Test raster data prediction")
library(klrfome)


test_that("sim_trend and pred_var_stack", {
  set.seed(717)
  sigma = 0.5
  lambda = 0.1
  dist_metric = "euclidean"
  sim_data <- get_sim_data(site_samples = 800, N_site_bags = 75)
  formatted_data <- format_site_data(sim_data, N_sites=10, train_test_split=0.8,
                                     sample_fraction = 0.9, background_site_balance=1)
  train_data <- formatted_data[["train_data"]]
  train_presence <- formatted_data[["train_presence"]]
  test_data <- formatted_data[["test_data"]]
  test_presence <- formatted_data[["test_presence"]]
  
  K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric, progress = FALSE)
  train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 0)
  test_log_pred <- KLR_predict(test_data, train_data, dist_metric = dist_metric,
                               train_log_pred[["alphas"]], sigma, progress = FALSE)
  params <- list(train_data = train_data,
                 alphas_pred = train_log_pred[["alphas"]],
                 sigma = sigma,
                 lambda = lambda,
                 means = formatted_data$means,
                 sds = formatted_data$sds)
  ### width and hieght of roving focal window (required)
  ngb = 5
  ### Number of rows and columns in prediction rasters
  ## needed for making simulated rasters, as well as for predicting real-world rasters
  cols = 50
  rows = 50
  ### Create simulated environmental rasters  (sim data only) ####
  s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
  s_var1 <- rescale_sim_raster(s_var1r, 50, 10) 
  s_var2 <- rescale_sim_raster(s_var1r, 3, 2) 
  b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
  b_var1 <- rescale_sim_raster(b_var1r, 100, 20) 
  b_var2 <- rescale_sim_raster(b_var1r, 6, 3) 
  ### Create a site-present trend surface  (sim data only)
  sim_sites_n = 3
  trend_coords <- sim_trend(cols, rows, n = sim_sites_n)
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
  
  # check to make sure all sim site coordinates are reported
  expect_true(nrow(trend_coords$coords) == sim_sites_n)
  # check to make sure trend_coords has a raster
  expect_is(trend_coords$trend, "RasterLayer")
  # check to see that pred_var_stack is a raster stack
  expect_is(pred_var_stack, "RasterStack")
  # check to see that pred_var_stack_scaled is a raster stack
  expect_is(pred_var_stack_scaled, "RasterStack")
  # check to see that pred_rast is a raster
  expect_is(pred_rast, "RasterLayer")
  # check to make sure raster predictions are greater than or equal to zero
  expect_true(all(range(pred_rast@data@values) >= 0))
})
