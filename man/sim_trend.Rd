% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raster_predict_funs.R
\name{sim_trend}
\alias{sim_trend}
\title{sim_trend}
\usage{
sim_trend(cols, rows, n = 2, size = 2)
}
\arguments{
\item{cols}{[integer] the number of columns in the simulated landscape raster}

\item{rows}{[integer] the number of rows in the simulated landscape raster}

\item{n}{[integer] the number of columns in the simulated landscape raster}

\item{size}{[integer] the pixel dimensions for simualted sites}
}
\value{
[list] a list with two elements, `coords` are the center coordiantes of the simulated sites, and `trend` is the raster trend surface to be applied to the simulated landscape.
}
\description{
A helper function to create a trend surface that when applied to a simulated landscape raster correlates it with simulated site locations
}
\details{
`sim_trend()` - Function is used to take create `n` number of simulated site locations of `size` cell dimensions on a rows by cols raster. The latter two arguments should match the size of your simulated rasters. The function randomly locates the sites and then creates a distance gradient (trend) from the site locations outward. The trend is a value 1 at the sites and reduces to 0 at the maximum combined distance from all sites. The output of this function is a list of a matrix of simulated site x/y coordinates (centers) and a raster of the trend surface. The point of the trend is to then combine it with the simulated rasters (as down in the code above) such that the raster depicting site-likely conditions is multiplied by the trend to output a raster retaining site-likely conditions near simulated site locations. Conversely, the site-unlikely simulated raster is multiplied by the inverse of the trend to result in a raster retaining site-unlikely characteristics away from the site locations. When those two rasters are added you get a simulated environment that is more preferable to site locations near site locations. It is a bit involved for something that had nothing to do with the actual KLRfome model, but it is needed to produce actual correlated environments for model testing.
}
\examples{
\dontrun{
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
}

}
