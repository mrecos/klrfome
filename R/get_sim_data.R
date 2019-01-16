#' get_sim_data
#' 
#' `get_sim_data()` is a function to simulate archaeological site data for use in this package. Beacuse archaeological site data is often a protected class of data or because this package requires data in a particular format; this function helps to test run models without neededing real data.
#' 
#' This function requires only the setting of parameters and does not rely on input data. The function generated as data.frame of with four columns including `presence` as a 1/0 variable for site presence/absence, `SITENO` that is a group level label for each site (presence == 1) and for each collection of background data points (presence == 0), and `var1` and `var2` that are the simulated environmental variables used as regression features. The number of sites and background samples to generate is controlled by the `site_samples` and `N_site_bags` arguments respectively. The `background_site_balance` controls the ratio of site to background groups in the output data.frame. The remaining arguments are used to control how similar the two simulated environmental features are within and between site presence/absence groups. The features are simulated from a normal distribution and can be controlled by adjusting the `mean` and `sd` (standard deviation) arguments for each. Similar distribution parameters between site and background approximate a more homogeneous environment where sites are found on less distinct landforms. Alternatively, greater difference in the distribution parameters between sites and background approximate an environment where sites are found on distinct landforms. Adjusting this can give insight into how the model works for a range of proxy environments.
#'
#' @param site_samples - [scaler] Number of sites to include in simulated data
#' @param N_site_bags - [scaler] Number of background groups to use in simulated data
#' @param sites_var1_mean - [scaler] mean of variable 1 for sites
#' @param sites_var1_sd - [scaler] sd of variable 1 for sites
#' @param sites_var2_mean - [scaler] mean of variable 2 for sites
#' @param sites_var2_sd - [scaler] sd of variable 2 for sites
#' @param backg_var1_mean - [scaler] mean of variable 1 for background
#' @param backg_var1_sd - [scaler] sd of variable 1 for background
#' @param backg_var2_mean - [scaler] mean of variable 2 for background
#' @param backg_var2_sd - [scaler] sd of variable 2 for background
#' @param background_site_balance - [scaler] Ratio of site to bacground groups
#'
#' @return - data.frame of simulated training and test set data
#' @import dplyr
#' @export
#' 
#' @examples
#'\dontrun{
#' sim_data <- get_sim_data(site_samples = 800, N_site_bags = 75,
#' sites_var1_mean = 80, sites_var1_sd   = 10,
#' sites_var2_mean = 5,  sites_var2_sd   = 2,
#' backg_var1_mean = 100,backg_var1_sd   = 20,
#' backg_var2_mean = 6,  backg_var2_sd   = 3)
#' formatted_data <- format_site_data(sim_data, N_sites=10, train_test_split=0.8,
#'                                    sample_fraction = 0.9, background_site_balance=1)
#' train_data <- formatted_data[["train_data"]]
#' train_presence <- formatted_data[["train_presence"]]
#' test_presence <- formatted_data[["test_presence"]]
#'
#' ##### Logistic Mean Embedding KLR Model
#' #### Build Kernel Matrix
#' K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric)
#' #### Train
#' train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 2)
#' #### Predict
#' test_log_pred <- KLR_predict(test_data, train_data, dist_metric = dist_metric,
#'                             train_log_pred[["alphas"]], sigma)
#'}
#'
#'
get_sim_data <- function(site_samples, N_site_bags,
                          sites_var1_mean = 50, sites_var1_sd   = 10,
                          sites_var2_mean = 3,  sites_var2_sd   = 2,
                          backg_var1_mean = 100,backg_var1_sd   = 20,
                          backg_var2_mean = 6,  backg_var2_sd   = 3,
                          background_site_balance = 1){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  back_samples <- site_samples*background_site_balance
  sites <- data.frame(presence = 1,
                      SITENO = sample(paste0("Site",1:N_site_bags),site_samples, replace = TRUE),
                      var1 = rnorm(site_samples,sites_var1_mean,sites_var1_sd),
                      var2 = rnorm(site_samples,sites_var2_mean,sites_var2_sd),
                      stringsAsFactors = FALSE) %>%
    arrange(SITENO)
  backs <- data.frame(presence = 0,
                      SITENO = "background",
                      var1 = rnorm(site_samples,backg_var1_mean,backg_var1_sd),
                      var2 = rnorm(site_samples,backg_var2_mean,backg_var2_sd))
  sim_data   <- rbind(sites,backs)
  return(sim_data)
}
