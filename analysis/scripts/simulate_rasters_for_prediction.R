rescale_raster <- function(rast,mean,sd){
  r_vals <- rast@data@values
  rast[] <- mean + (r_vals - mean(r_vals)) * (sd/sd(r_vals))
  return(rast)
}
sim_trend <- function(cols, rows, n = 2, size = 2){
  trend <- raster(ncols = cols,nrows = rows,res = 1,
                  xmn=0, xmx = cols, ymn=0, ymx=rows, vals = 0)
  coords <- matrix(ncol=2,nrow=n)
  colnames(coords) <- c("x","y")
  for(i in seq_len(n)){
    crd <- c(sample(1:(rows-size),1),sample(1:(cols-size),1))
    coords[i,1] <- crd[1]
    coords[i,2] <- crd[2]
    # had to adjust origin to match expected from crds.
    site_i <- NLMR::nlm_distancegradient(cols,rows,
                                   origin = c(rows-crd[2],rows-crd[2]+size,
                                              crd[1],crd[1]-size))^0.5
    trend <- trend + site_i
  }
  mn <- cellStats(trend, min)
  mx <- cellStats(trend, max)
  trend <- abs(1-((trend-mn)/(mx-mn))) # standardize and invert
  return(list(trend = trend, coords = coords))
}
scale_prediction_rasters <- function(pred_var_stack, params, verbose = 1){
  pred_var_stack_scaled <- pred_var_stack
  for(i in seq_len(dim(pred_var_stack)[3])){
    var_name <- names(pred_var_stack[[i]])
    if(verbose == 1){
      cat("Normalizing:", var_name, "\n")
    }
    pred_var_stack_scaled[[var_name]] <- scale(pred_var_stack[[var_name]],
                                               center = params$means[var_name],
                                               scale  = params$sds[var_name])
  }
  return(pred_var_stack_scaled)
}
KLR_raster_predict <- function(rast_stack, ngb, params, progress = TRUE){
  preds_focal <- getValuesFocal(rast_stack, ngb = ngb, names = TRUE)
  # turn neighborhoods into list for KLR_predict
  test_data_raw <- list()
  for(z in 1:nrow(preds_focal[[1]])){
    preds_focal_row_z <- lapply(preds_focal, `[`,z,) # this way of subsetting only works with lapply()
    preds_focal_element_z <- sapply(preds_focal_row_z, cbind)
    if(ngb == 1){
      test_data_raw[[z]] <- t(as.matrix(preds_focal_element_z))
    } else {
      test_data_raw[[z]] <- preds_focal_element_z
    }
  }
  # get only the neighborhoods that have data, ignore the NA hoods.
  is_na <- which(sapply(test_data_raw, sum, na.rm=TRUE) == 0)
  has_data_index <- which(sapply(test_data_raw, sum, na.rm=TRUE) != 0)
  has_data <- test_data_raw[has_data_index]
  # send t0 KLR_predict
  tpred <- KLR_predict(has_data, params$train_data, params$alphas_pred, params$sigma, progress = progress)
  # assign predicted values to raster, and NA to the NA spots
  pred_rast <- raster(rast_stack)
  pred_rast[is_na] <- NA
  pred_rast[has_data_index] <- tpred
  return(pred_rast)
}

library("NLMR")
library("rasterVis")
cols = 100
rows = 100
ngb = 3

params <- list(train_data = train_data,
               alphas_pred = train_log_pred[["alphas"]],
               sigma = sigma,
               lambda = lambda,
               vars = vars,
               means = formatted_data$means,
               sds = formatted_data$sds)

# could be a function...
s_var1r <- nlm_gaussianfield(cols,rows, autocorr_range = 20)
s_var1 <- rescale_raster(s_var1r, 50, 10) 
s_var2 <- rescale_raster(s_var1r, 3, 2) 
b_var1r <- nlm_gaussianfield(cols,rows,autocorr_range = 20)
b_var1 <- rescale_raster(b_var1r, 100, 20) 
b_var2 <- rescale_raster(b_var1r, 6, 3) 
# get site-present trend and its inverse, multiply site by trend and add to anti-trend
trend_coords <- sim_trend(cols, rows, n = 3)
coords <- trend_coords$coords
trend <- trend_coords$trend
inv_trend <- abs(1-trend)
var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
pred_var_stack <- raster::stack(var1, var2)
names(pred_var_stack) <- c("var1","var2")
plot(pred_var_stack)
# scale rasters to training data
pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
# Predict raster (single chunk) 
pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb =15, params, progress = TRUE)
# plot with simulated sites
rasterVis::levelplot(pred_rast, margin = FALSE, ,par.settings=viridisTheme()) +
  layer(sp.points(sp.points(SpatialPoints(coords), pch=15, cex = 3, col = "red")), columns=1)




