#' rescale_raster
#'
#' @param rast 
#' @param mean 
#' @param sd 
#'
#' @return [raster]
#' @export
#'

rescale_raster <- function(rast,mean,sd){
  r_vals <- rast@data@values
  rast[] <- mean + (r_vals - mean(r_vals)) * (sd/sd(r_vals))
  return(rast)
}

#' sim_trend
#'
#' @param cols 
#' @param rows 
#' @param n 
#' @param size 
#'
#' @importFrom NLMR nlm_distancegradient
#' @return [raster]
#' @export
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
  mn <- cellStats(trend, min)
  mx <- cellStats(trend, max)
  trend <- abs(1-((trend-mn)/(mx-mn))) # standardize and invert
  return(list(trend = trend, coords = coords))
}

#' scale_prediction_rasters
#'
#' @param pred_var_stack 
#' @param params 
#' @param verbose 
#'
#' @importFrom raster scale
#' @return [raster]
#' @export
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
#' @param rast_stack 
#' @param ppside 
#' @param split 
#'
#' @importFrom raster rasterToPolygons extent crop ncell aggregate
#' @return [raster]
#' @export
#'

split_raster_stack <- function(rast_stack, ppside, split=TRUE){
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
    r_list <- list()
    for(i in 1:raster::ncell(agg)){
      e1          <- raster::extent(agg_poly[agg_poly$polis==i,])
      crop_rasters_i <- raster::crop(rast_stack,e1)
      r_list[[i]] <- crop_rasters_i
    }
  }
  return(r_list)
}

#' KLR_raster_predict
#'
#' @param rast_stack 
#' @param ngb 
#' @param params 
#' @param split 
#' @param ppside 
#' @param progress 
#' @param parallel 
#' @param output [string] either 'list' or 'save'. If 'list', returns list of predicted blocks. If 'save', blocks are saved to save_loc as GeoTiff
#' @param save_loc [string] Location to save raster blocks (GeoTiff). Uses getwd() is NULL
#' @param overwrite [logical] TRUE will overwrite saved raster GeoTiffs in save_loc
#'
#' @import foreach
#' @import doParallel
#' @return [raster]
#' @export
#'

KLR_raster_predict <- function(rast_stack, ngb, params, split = FALSE, ppside = NULL, 
                               progress = TRUE, parallel = FALSE, output = "list",
                               save_loc = NULL, overwrite = FALSE){
  if(isTRUE(parallel) & foreach::getDoParRegistered() == FALSE){
    stop("You must make and register a parallel cluster with' doParallel' package","\n")
  }
  if(!isTRUE(split) & is.null(ppside)){ # default settings
    if(isTRUE(parallel)){
      message("Parallel only works when split == TRUE and ppside > 1. Continuing without parallel","\n")
    }
    pred_rast <- klrfome::KLR_predict_each(rast_stack, ngb, params, progress)
    return(pred_rast)
  } else if(isTRUE(split) & !is.null(ppside)){ # split & ppside == TRUE
    if(!(output %in% c("list","save"))){ # check for output to be set
      stop("In order to use split the raster, you need to set output = to 'list' or 'save' & 'save_loc","\n")
    }
    cat("Splitting rasters into blocks","\n")
    split_stack <- klrfome::split_raster_stack(rast_stack, ppside)
    if(isTRUE(parallel)){ # split, ppside, and parallel = TRUE
      cat("Predicting splits in parallel on",getDoParWorkers(),"cores","\n")
      if(output == "list"){ # split, ppside, parallel, and list output
        pred_rast_list = foreach(i=seq_along(split_stack), .inorder = TRUE, .verbose = progress,
                                 .packages=c('klrfome','raster')) %dopar% {
                                   klrfome::KLR_predict_each(split_stack[[i]], ngb, params, progress)
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
                  pred_rast_i <- KLR_predict_each(split_stack[[i]], ngb, params, progress)
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
        pred_rast_i <- klrfome::KLR_predict_each(split_stack[[i]], ngb, params, progress)
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
#' @param rast_stack 
#' @param ngb 
#' @param params 
#' @param progress 
#'
#' @importFrom raster getValuesFocal raster
#' @return [raster]
#' @export
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
  # send t0 KLR_predict
  tpred <- klrfome::KLR_predict(has_data, params$train_data, params$alphas_pred, params$sigma, progress = progress)
  # assign predicted values to raster, and NA to the NA spots
  pred_rast <- raster::raster(rast_stack)
  pred_rast[is_na] <- NA
  pred_rast[has_data_index] <- tpred
  return(pred_rast)
}
