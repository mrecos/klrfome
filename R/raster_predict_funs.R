#' rescale_sim_raster
#' 
#' Helper function to rescale a simulated landscape to match the mean and standard deviation of the simulated data used to fit a model.
#' `rescale_sim_raster()` - Function that rescales simulated rasters from `NLMR::nlm_gaussianfield` (Sciaini, Fritsch, and Simpkins 2017) or whatever you want to use, to the mean and standard deviation of the simulated data used to fit the klr model. You will have to add the mean and sd arguments manually based on what you put into the get_sim_data function. The example in the code above inputs the default mean and sd values from the defualts of the `get_sim_data()` function. Returned is a raster scaled too your simualted training data.
#'
#' @param rast [raster] input raster object
#' @param mean [numeric] value of mean used in `get_sim_data()` function. 
#' @param sd   [numeric] value of standard deviation used in `get_sim_data()` function. 
#'
#' @return [raster] a raster object scaled to the input data
#' @export
#' 
#' @examples 
#' \dontrun{
#' ### width and hieght of roving focal window (required)
#' ngb = 5
#' ### Number of rows and columns in prediction rasters
#' ## needed for making simulated rasters, as well as for predicting real-world rasters
#' cols = 100
#' rows = 100
#' 
#' ### Create simulated environmental rasters  (sim data only) ####
#' s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
#' s_var1 <- rescale_sim_raster(s_var1r, 50, 10) 
#' s_var2 <- rescale_sim_raster(s_var1r, 3, 2) 
#' b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
#' b_var1 <- rescale_sim_raster(b_var1r, 100, 20) 
#' b_var2 <- rescale_sim_raster(b_var1r, 6, 3) 
#' ### Create a site-present trend surface  (sim data only)
#' trend_coords <- sim_trend(cols, rows, n = 3)
#' coords <- trend_coords$coords
#' trend <- trend_coords$trend
#' inv_trend <- abs(1-trend)
#' var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
#' var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#' #### end simulated data creation ####
#' 
#' ### Create raster stack of predictor variables
#' pred_var_stack <- raster::stack(var1, var2)
#' names(pred_var_stack) <- c("var1","var2")
#' ### scale rasters to training data
#' pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
#' ### Predict raster (single chunk, not in parallel) 
#' pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
#'                                 progress = FALSE, parallel = FALSE)
#' }
#'
rescale_sim_raster <- function(rast,mean,sd){
  r_vals <- rast@data@values
  rast[] <- mean + (r_vals - mean(r_vals)) * (sd/sd(r_vals))
  return(rast)
}


#' sim_trend
#' 
#' A helper function to create a trend surface that when applied to a simulated landscape raster correlates it with simulated site locations
#' 
#' `sim_trend()` - Function is used to take create `n` number of simulated site locations of `size` cell dimensions on a rows by cols raster. The latter two arguments should match the size of your simulated rasters. The function randomly locates the sites and then creates a distance gradient (trend) from the site locations outward. The trend is a value 1 at the sites and reduces to 0 at the maximum combined distance from all sites. The output of this function is a list of a matrix of simulated site x/y coordinates (centers) and a raster of the trend surface. The point of the trend is to then combine it with the simulated rasters (as down in the code above) such that the raster depicting site-likely conditions is multiplied by the trend to output a raster retaining site-likely conditions near simulated site locations. Conversely, the site-unlikely simulated raster is multiplied by the inverse of the trend to result in a raster retaining site-unlikely characteristics away from the site locations. When those two rasters are added you get a simulated environment that is more preferable to site locations near site locations. It is a bit involved for something that had nothing to do with the actual KLRfome model, but it is needed to produce actual correlated environments for model testing.
#'
#' @param cols [integer] the number of columns in the simulated landscape raster
#' @param rows [integer] the number of rows in the simulated landscape raster
#' @param n    [integer] the number of columns in the simulated landscape raster
#' @param size [integer] the pixel dimensions for simualted sites
#'
#' @importFrom NLMR nlm_distancegradient
#' @return [list] a list with two elements, `coords` are the center coordiantes of the simulated sites, and `trend` is the raster trend surface to be applied to the simulated landscape.
#' @export
#'
#' @examples 
#' \dontrun{
#' ### width and hieght of roving focal window (required)
#' ngb = 5
#' ### Number of rows and columns in prediction rasters
#' ## needed for making simulated rasters, as well as for predicting real-world rasters
#' cols = 100
#' rows = 100
#' 
#' ### Create simulated environmental rasters  (sim data only) ####
#' s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
#' s_var1 <- rescale_sim_raster(s_var1r, 50, 10) 
#' s_var2 <- rescale_sim_raster(s_var1r, 3, 2) 
#' b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
#' b_var1 <- rescale_sim_raster(b_var1r, 100, 20) 
#' b_var2 <- rescale_sim_raster(b_var1r, 6, 3) 
#' ### Create a site-present trend surface  (sim data only)
#' trend_coords <- sim_trend(cols, rows, n = 3)
#' coords <- trend_coords$coords
#' trend <- trend_coords$trend
#' inv_trend <- abs(1-trend)
#' var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
#' var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#' #### end simulated data creation ####
#' 
#' ### Create raster stack of predictor variables
#' pred_var_stack <- raster::stack(var1, var2)
#' names(pred_var_stack) <- c("var1","var2")
#' ### scale rasters to training data
#' pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
#' ### Predict raster (single chunk, not in parallel) 
#' pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
#'                                 progress = FALSE, parallel = FALSE)
#' }
#'
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
  mn <- raster::cellStats(trend, min)
  mx <- raster::cellStats(trend, max)
  trend <- abs(1-((trend-mn)/(mx-mn))) # standardize and invert
  return(list(trend = trend, coords = coords))
}

#' scale_prediction_rasters
#' 
#' Center and scale the prediction raster stack to the parameters of the training data used to fit a model
#' 
#'  `scale_prediction_rasters()` - Function scales your predictor rater stack based on the params list created in the model fitting process. This script simply loops over the rasters in the stack and centers and scales based on mean and sd of the training data used to fit the klr model. The function outputs a raster stack. 
#'
#' @param pred_var_stack [raster stack] stack of input landscape rasters for prediction
#' @param params [list] list of mean and standard deviation parameters returned from formatting data. (See example below) 
#' @param verbose [integer] a value of `1` to print out the landscape variable name as it is being normalized.
#'
#' @importFrom raster scale
#' @return [raster stack] a stack of the prediction rasters scaled to the the input training parameters.
#' @export
#' 
#' @examples 
#' \dontrun{
#' ### Create param list
#' params <- list(train_data = train_data,
#' alphas_pred = train_log_pred[["alphas"]],
#' sigma = sigma,
#' lambda = lambda,
#' means = formatted_data$means,
#' sds = formatted_data$sds)
#' 
#' ### width and hieght of roving focal window (required)
#' ngb = 5
#' ### Number of rows and columns in prediction rasters
#' ## needed for making simulated rasters, as well as for predicting real-world rasters
#' cols = 100
#' rows = 100
#' 
#' ### Create simulated environmental rasters  (sim data only) ####
#' s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
#' s_var1 <- rescale_sim_raster(s_var1r, 50, 10) 
#' s_var2 <- rescale_sim_raster(s_var1r, 3, 2) 
#' b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
#' b_var1 <- rescale_sim_raster(b_var1r, 100, 20) 
#' b_var2 <- rescale_sim_raster(b_var1r, 6, 3) 
#' ### Create a site-present trend surface  (sim data only)
#' trend_coords <- sim_trend(cols, rows, n = 3)
#' coords <- trend_coords$coords
#' trend <- trend_coords$trend
#' inv_trend <- abs(1-trend)
#' var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
#' var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#' #### end simulated data creation ####
#' 
#' ### Create raster stack of predictor variables
#' pred_var_stack <- raster::stack(var1, var2)
#' names(pred_var_stack) <- c("var1","var2")
#' ### scale rasters to training data
#' pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
#' ### Predict raster (single chunk, not in parallel) 
#' pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
#'                                 progress = FALSE, parallel = FALSE)
#' }
#'
scale_prediction_rasters <- function(pred_var_stack, params, verbose = 1){
  pred_var_stack_scaled <- pred_var_stack
  for(i in seq_len(dim(pred_var_stack)[3])){
    var_name <- names(pred_var_stack[[i]])
    if(verbose == 1){
      cat("Normalizing:", var_name, "\n")
    }
    pred_var_stack_scaled[[var_name]] <- raster::scale(pred_var_stack[[var_name]],
                                               center = params$means[var_name],
                                               scale  = params$sds[var_name])
  }
  return(pred_var_stack_scaled)
}


#' split_raster_stack
#' 
#' Internal function to used by `KLR_raster_predict()` split input rasters into overlapping blocks for prediction.
#' 
#' `split_raster_stack()` takes a raster stack, creates a grid of `ppside` pixels, adds a collar, clips out each raster, and returns a list  
#'
#' @param rast_stack 
#' @param ppside 
#' @param split 
#'
#' @importFrom raster rasterToPolygons extent crop ncell aggregate
#' @importFrom rgeos gBuffer
#' @return [list] a list of raster chips
#'
#'
split_raster_stack <- function(rast_stack, ppside, ngb, split=TRUE){
  # modified from - https://stackoverflow.com/questions/29784829/r-raster-package-split-image-into-multiples
  h        <- ceiling(ncol(rast_stack)/ppside)
  v        <- ceiling(nrow(rast_stack)/ppside)
  agg      <- raster::aggregate(rast_stack,fact=c(h,v))
  agg_out  <- agg
  r_list <- vector(mode="list", length = length(agg)/dim(agg)[3])
  if(isTRUE(split)){
    agg[]    <- 1:ncell(agg)
    agg_poly <- raster::rasterToPolygons(agg)
    names(agg_poly) <- rep("polis", length(names(agg_poly)))
    ## adds collar to allow focal window function to be accuracte at the margins
    exp_agg <- rgeos::gBuffer(agg_poly, byid = TRUE, 
                              width = ngb, joinStyle = "mitre", mitreLimit = ngb)
    r_list <- list()
    for(i in 1:raster::ncell(agg)){
      e1             <- raster::extent(exp_agg[exp_agg$polis==i,])
      crop_rasters_i <- raster::crop(rast_stack,e1)
      r_list[[i]]    <- crop_rasters_i
    }
  }
  return(r_list)
}


#' crop_raster_collar
#' 
#' Utility function to crop a raster by a neighborhood. Not exported
#'
#' @param pred_rast [raster] predicted raster
#' @param ngb [integer] square neighborhodd size in cells
#' @param cols [integer] cols of crop
#' @param rows [integer] rows of crop
#' 
#' @importFrom raster crop extent
#' @return [raster] a cropped raster
#'
crop_raster_collar <- function(pred_rast, ngb, cols, rows){
  ext_i    <- extent(pred_rast)
  xmin_i   <- ifelse(ext_i[1] == 0, 0, ext_i[1]+ngb)
  xmax_i   <- ifelse(ext_i[2] == cols, cols, ext_i[2]-ngb)
  ymin_i   <- ifelse(ext_i[3] == 0, 0, ext_i[3]+ngb)
  ymax_i   <- ifelse(ext_i[4] == rows, rows, ext_i[4]-ngb)
  crop_ext <- raster::extent(xmin_i,xmax_i,ymin_i,ymax_i)
  pred_rast_crop <- raster::crop(pred_rast, crop_ext)
}

#' KLR_raster_predict
#' 
#' `KLR_raster_predict()` predicts a fit KLR model to new data.
#' 
#' `KLR_raster_predict()` - Function predicts the probability of site-presence based on a `rast_stack` raster stack of center/scaled predictor rasters, a focal neighborhood size in cells as `ngb`, and the params list of model fit parameters. Finally, the function also needs to the the number of columns as cols and rows as rows of the study areas raster stack. The rest of the arguments default to predicting the entire raster stack in one pass and only on a single core (not in parallel). The rest of the argument control whether the study area is split (`split` = TRUE) into a grid of blocks. The `ppside` positive integer controls the number of blocks along each axis of the study area. If you wish to compute the prediction in `parallel`, you will need to split it into blocks so that each block can be sent to a different processor core. The final set of optional arguments control how the prediction is returned, either as a `output` = "list" or `output` = "save" for returning a list of rasters of saving each out as a GeoTiff, and then arguments for the GeoTiff location with `save_loc` and whether to overwrite existing GeoTiffs with `overwrite.` The function contains a bit of logic to try and assist the user in which arguments go with what. Perhaps future versions will streamline this a bit.
#'
#' @param rast_stack [raster stack] of prediction rasters (usually scaled)
#' @param ngb [integer] pixel dimension of square neighborhood used to predict model on
#' @param params [list] list of important model parameters. (See example)
#' @param split [logical] TRUE/FALSE for whether to split raster into chunks. If `parallal` == TRUE, then `split` must == TRUE
#' @param ppside [integer] the number of blocks to split the study area into if `split` == TRUE
#' @param progress [logical] if TRUE, displays progress of function
#' @param parallel [logical] if TRUE, prediction will be excuted as parallel computations. Requires a parallel backend, `split` == TRUE, and number for `ppside`
#' @param output [string] either 'list' or 'save'. If 'list', returns list of predicted blocks. If 'save', blocks are saved to save_loc as GeoTiff
#' @param save_loc [string] Location to save raster blocks (GeoTiff). Uses getwd() is NULL
#' @param overwrite [logical] TRUE will overwrite saved raster GeoTiffs in save_loc
#'
#' @import foreach
#' @import doParallel
#' @return [raster], [list], or [nothing] - If `parallel` == FALSE, a raster object is returned. If `parallel` == TRUE then either a list of raster objects if `output` == "list", or files are saved as images if `output` == "save"
#' @export
#'
#' @examples 
#' \dontrun{
#' ### Create param list
#' params <- list(train_data = train_data,
#' alphas_pred = train_log_pred[["alphas"]],
#' sigma = sigma,
#' lambda = lambda,
#' means = formatted_data$means,
#' sds = formatted_data$sds)
#' 
#' ### width and hieght of roving focal window (required)
#' ngb = 5
#' ### Number of rows and columns in prediction rasters
#' ## needed for making simulated rasters, as well as for predicting real-world rasters
#' cols = 100
#' rows = 100
#' 
#' ### Create simulated environmental rasters  (sim data only) ####
#' s_var1r <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 20)
#' s_var1 <- rescale_sim_raster(s_var1r, 50, 10) 
#' s_var2 <- rescale_sim_raster(s_var1r, 3, 2) 
#' b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
#' b_var1 <- rescale_sim_raster(b_var1r, 100, 20) 
#' b_var2 <- rescale_sim_raster(b_var1r, 6, 3) 
#' ### Create a site-present trend surface  (sim data only)
#' trend_coords <- sim_trend(cols, rows, n = 3)
#' coords <- trend_coords$coords
#' trend <- trend_coords$trend
#' inv_trend <- abs(1-trend)
#' var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
#' var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#' #### end simulated data creation ####
#' 
#' ### Create raster stack of predictor variables
#' pred_var_stack <- raster::stack(var1, var2)
#' names(pred_var_stack) <- c("var1","var2")
#' ### scale rasters to training data
#' pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
#' ### Predict raster (single chunk, not in parallel) 
#' pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
#'                                 progress = FALSE, parallel = FALSE)
#' ### plot with simulated sites
#' rasterVis::levelplot(pred_rast, margin = FALSE, par.settings=viridisTheme()) +
#'layer(sp.points(sp.points(SpatialPoints(coords), pch=15, cex = 2.25, col = "red")), columns=1)
#'
#'    ## Or with parallel backend ##
#' library("doParallel")
#' ### create and register parallel backend
#' cl <- makeCluster(detectCores())
#' doParallel::registerDoParallel(cl)
#' 
#' ### Use same KLR_raster_predict function with parallel = TRUE
#' pred_rast_list <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = TRUE, ppside = 5,
#'                                      progress = FALSE, parallel = TRUE, output = "list",
#'                                      save_loc = NULL, overwrite = TRUE, cols = cols, rows = rows)
#' #> Splitting rasters into blocks 
#' #> Predicting splits in parallel on 8 cores
#' ### Merge list back to a single raster
#' pred_rast <-  do.call(merge, pred_rast_list)
#' ### plot with simulated sites
#' rasterVis::levelplot(pred_rast, margin = FALSE, par.settings=viridisTheme()) +
#'   layer(sp.points(sp.points(SpatialPoints(coords), pch=15, cex = 2.25, col = "red")), columns=1)    
#'                          
#' }
#'
KLR_raster_predict <- function(rast_stack, ngb, params, cols, rows, split = FALSE, ppside = NULL, 
                               progress = TRUE, parallel = FALSE, output = "list",
                               save_loc = NULL, overwrite = FALSE){
  if(isTRUE(parallel) & foreach::getDoParRegistered() == FALSE){
    stop("You must make and register a parallel cluster with' doParallel' package","\n")
  }
  if(!isTRUE(split) & is.null(ppside)){ # default settings
    if(isTRUE(parallel)){
      message("Parallel only works when split == TRUE and ppside > 1. Continuing without parallel","\n")
    }
    pred_rast <- klrfome:::KLR_predict_each(rast_stack, ngb, params, progress)
    return(pred_rast)
  } else if(isTRUE(split) & !is.null(ppside)){ # split & ppside == TRUE
    if(!(output %in% c("list","save"))){ # check for output to be set
      stop("In order to use split the raster, you need to set output = to 'list' or 'save' & 'save_loc","\n")
    }
    cat("Splitting rasters into blocks","\n")
    split_stack <- klrfome:::split_raster_stack(rast_stack, ppside, ngb, split=TRUE)
    if(isTRUE(parallel)){ # split, ppside, and parallel = TRUE
      cat("Predicting splits in parallel on",getDoParWorkers(),"cores","\n")
      if(output == "list"){ # split, ppside, parallel, and list output
        pred_rast_list = foreach(i=seq_along(split_stack), .inorder = TRUE, .verbose = progress,
                                 .packages=c('klrfome','raster')) %dopar% {
                                   pred_rast_i <- klrfome:::KLR_predict_each(split_stack[[i]], ngb, params, progress)
                                   crop_raster_collar(pred_rast_i, ngb, cols, rows)
                                 }
        return(pred_rast_list)
      } else if(output == "save"){ # split, ppside, parallel, and save output
        message("Predicting blocks in parallel and saving to 'save_loc'")
        if(is.null(save_loc)){
          message("save_loc not set, saving to working directory")
          save_loc <- getwd()
        }
        foreach(i=seq_along(split_stack), .inorder = TRUE, .verbose = progress,
                .packages=c('klrfome','raster')) %dopar% {
                  pred_rast_i <- klrfome:::KLR_predict_each(split_stack[[i]], ngb, params, progress)
                  pred_rast_i <- crop_raster_collar(pred_rast_i, ngb, cols, rows)
                  writeRaster(pred_rast_i, filename=file.path(save_loc, paste("prediction_block_",i,".tif",sep="")),
                              format="GTiff",datatype="FLT4S",overwrite = overwrite)
                }
      }
    } else if(!isTRUE(parallel)){ # split, ppside == TRUE, no parallel
      cat("You should really consider using doParallel package to use your multiple cores!","\n")
      if(output == "save" & is.null(save_loc)){
        message("save_loc not set, saving to working directory")
        save_loc <- getwd()
      }
      pred_rast_list <- vector(mode = "list", length = length(split_stack))
      for(i in seq_along(split_stack)){
        cat("Predicting for block", i, "of", length(split_stack),"\n")
        pred_rast_i <- klrfome:::KLR_predict_each(split_stack[[i]], ngb, params, progress)
        # deal with edge effects of blocking
        pred_rast_i <- crop_raster_collar(pred_rast_i, ngb, cols, rows)
        if(output == "list"){ # split, ppside, no parallel, list output
          pred_rast_list[[i]] <- pred_rast_i
        } else if(output == "save"){ # split, ppside, no parallel, save output
          writeRaster(pred_rast_i, filename=file.path(save_loc, paste("prediction_block_",i,".tif",sep="")),
                      format="GTiff",datatype="FLT4S",overwrite = overwrite)
        }
      }
      if(output == "list"){
        return(pred_rast_list)
      }
    }
    # return(pred_rast_list)
  } else if(isTRUE(split) & is.null(ppside)){ # split, ppside is NULL
    message("If you want to split the raster predictions, please provide a ppside > 1","\n")
  }
}

#' KLR_predict_each
#' 
#' `KLR_predict_each()` predict probability of presence on a raster or neighborhood. Not exported
#' 
#' This is an internal function to bridge between `KLR_raster_predict()` and `KLR_predict()`. This takes a raster stack of the study area of blocks, formats that data into the tabular format for prediction, and send the data to `KLR_predict()`.
#'
#' @param rast_stack [raster stack] either a full study area stack or blocks from within a split stack
#' @param ngb [integer] pixel dimension of focal neighborhood
#' @param params [list] params from fitting process, see `KLR_raster_predict()` examples
#' @param progress [logical] TRUE to show function progress
#'
#' @importFrom raster getValuesFocal raster
#' @return [raster] a raster of predicted presence
#'
KLR_predict_each <- function(rast_stack, ngb, params, progress = progress){
  preds_focal <- raster::getValuesFocal(rast_stack, ngb = ngb, names = TRUE)
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
  # send to KLR_predict
  tpred <- klrfome::KLR_predict(has_data, params$train_data, params$alphas_pred, params$sigma, progress = progress)
  # assign predicted values to raster, and NA to the NA spots
  pred_rast <- raster::raster(rast_stack)
  pred_rast[is_na] <- NA
  pred_rast[has_data_index] <- tpred
  return(pred_rast)
}
