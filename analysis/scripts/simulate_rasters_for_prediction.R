# rescale_raster <- function(rast,mean,sd){
#   r_vals <- rast@data@values
#   rast[] <- mean + (r_vals - mean(r_vals)) * (sd/sd(r_vals))
#   return(rast)
# }
# sim_trend <- function(cols, rows, n = 2, size = 2){
#   trend <- raster(ncols = cols,nrows = rows,res = 1,
#                   xmn=0, xmx = cols, ymn=0, ymx=rows, vals = 0)
#   coords <- matrix(ncol=2,nrow=n)
#   colnames(coords) <- c("x","y")
#   for(i in seq_len(n)){
#     crd <- c(sample(1:(rows-size),1),sample(1:(cols-size),1))
#     coords[i,1] <- crd[1]
#     coords[i,2] <- crd[2]
#     # had to adjust origin to match expected from crds.
#     site_i <- NLMR::nlm_distancegradient(cols,rows,
#                                    origin = c(rows-crd[2],rows-crd[2]+size,
#                                               crd[1],crd[1]-size))^0.5
#     trend <- trend + site_i
#   }
#   mn <- cellStats(trend, min)
#   mx <- cellStats(trend, max)
#   trend <- abs(1-((trend-mn)/(mx-mn))) # standardize and invert
#   return(list(trend = trend, coords = coords))
# }
# scale_prediction_rasters <- function(pred_var_stack, params, verbose = 1){
#   pred_var_stack_scaled <- pred_var_stack
#   for(i in seq_len(dim(pred_var_stack)[3])){
#     var_name <- names(pred_var_stack[[i]])
#     if(verbose == 1){
#       cat("Normalizing:", var_name, "\n")
#     }
#     pred_var_stack_scaled[[var_name]] <- scale(pred_var_stack[[var_name]],
#                                                center = params$means[var_name],
#                                                scale  = params$sds[var_name])
#   }
#   return(pred_var_stack_scaled)
# }
# 
# split_raster_stack <- function(rast_stack, ppside, split=TRUE){
#   # modified from - https://stackoverflow.com/questions/29784829/r-raster-package-split-image-into-multiples
#   h        <- ceiling(ncol(rast_stack)/ppside)
#   v        <- ceiling(nrow(rast_stack)/ppside)
#   agg      <- aggregate(rast_stack,fact=c(h,v))
#   agg_out  <- agg
#   r_list <- vector(mode="list", length = length(agg)/dim(agg)[3])
#   if(isTRUE(split)){
#     agg[]    <- 1:ncell(agg)
#     agg_poly <- rasterToPolygons(agg)
#     names(agg_poly) <- rep("polis", length(names(agg_poly)))
#     r_list <- list()
#     for(i in 1:ncell(agg)){
#       e1          <- extent(agg_poly[agg_poly$polis==i,])
#       crop_rasters_i <- crop(rast_stack,e1)
#       r_list[[i]] <- crop_rasters_i
#     }
#   }
#   return(r_list)
# }
# 
# KLR_raster_predict <- function(rast_stack, ngb, params, split = FALSE, ppside = NULL, 
#                                progress = TRUE, parallel = FALSE){
#   if(!isTRUE(split) & is.null(ppside)){
#     if(isTRUE(parallel)){
#       message("Parallel only works with split == TRUE and ppside > 1. Continuing without parallel","\n")
#     }
#     pred_rast <- KLR_predict_each(rast_stack, ngb, params, progress)
#     return(pred_rast)
#   } else if(isTRUE(split) & !is.null(ppside)){
#       cat("Splitting rasters into blocks","\n")
#       split_stack <- split_raster_stack(rast_stack, ppside)
#       if(isTRUE(parallel)){
#         if(getDoParRegistered() == FALSE){
#           stop("You must make and register a parallel cluster with' doParallel' package","\n")
#         }
#         cat("Predicting splits in parallel on",getDoParWorkers(),"cores","\n")
#         pred_rast_list = foreach(i=seq_along(split_stack), .inorder = TRUE,
#                                   .packages=c('klrfome','raster')) %dopar% {
#                                     klrfome::KLR_predict_each(split_stack[[i]], ngb, params, progress)
#                                   }
#       } else if(!isTRUE(parallel)){
#         cat("You should really consider using doParallel package to use your multiple core!","\n")
#         pred_rast_list <- vector(mode = "list", length = length(split_stack))
#       for(i in seq_along(split_stack)){
#         cat("Predicting for block", i, "of", length(split_stack),"\n")
#         pred_rast_i <- KLR_predict_each(split_stack[[i]], ngb, params, progress)
#         pred_rast_list[[i]] <- pred_rast_i
#         }
#       }
#       return(pred_rast_list)
#   } else if(isTRUE(split) & is.null(ppside)){
#     message("If you want to split the raster predictions, please provide a ppside > 1","\n")
#   }
# }
# 
# KLR_predict_each <- function(rast_stack, ngb, params, progress = progress){
#   preds_focal <- getValuesFocal(rast_stack, ngb = ngb, names = TRUE)
#   # turn neighborhoods into list for KLR_predict
#   test_data_raw <- list()
#   for(z in 1:nrow(preds_focal[[1]])){
#     preds_focal_row_z <- lapply(preds_focal, `[`,z,) # this way of subsetting only works with lapply()
#     preds_focal_element_z <- sapply(preds_focal_row_z, cbind)
#     if(ngb == 1){
#       test_data_raw[[z]] <- t(as.matrix(preds_focal_element_z))
#     } else {
#       test_data_raw[[z]] <- preds_focal_element_z
#     }
#   }
#   # get only the neighborhoods that have data, ignore the NA hoods.
#   is_na <- which(sapply(test_data_raw, sum, na.rm=TRUE) == 0)
#   has_data_index <- which(sapply(test_data_raw, sum, na.rm=TRUE) != 0)
#   has_data <- test_data_raw[has_data_index]
#   # send t0 KLR_predict
#   tpred <- KLR_predict(has_data, params$train_data, params$alphas_pred, params$sigma, progress = progress)
#   # assign predicted values to raster, and NA to the NA spots
#   pred_rast <- raster(rast_stack)
#   pred_rast[is_na] <- NA
#   pred_rast[has_data_index] <- tpred
#   return(pred_rast)
# }


library("NLMR")
library("rasterVis")
library("klrfome")

cols = 50
rows = 50
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
pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, progress = TRUE)
# plot with simulated sites
rasterVis::levelplot(pred_rast, margin = FALSE, ,par.settings=viridisTheme()) +
  layer(sp.points(sp.points(SpatialPoints(coords), pch=15, cex = 3, col = "red")), columns=1)


pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, progress = TRUE)

pred_rast2 <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params,
                                 split = TRUE, ppside = 2, progress = TRUE)

xx <-  do.call(merge, pred_rast1)
plot(xx)

