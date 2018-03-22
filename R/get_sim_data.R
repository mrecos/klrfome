#' get_sim_data
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
#' @return - data.frame and list of simulated training and test set data
#' @import dplyr
#' @export
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
