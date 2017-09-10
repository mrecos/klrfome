get_sim_data <- function(sites_var1_mean = 50,
                         sites_var1_sd   = 10,
                         sites_var2_mean = 3,
                         sites_var2_sd   = 2,
                         backg_var1_mean = 100,
                         backg_var1_sd   = 20,
                         backg_var2_mean = 6,
                         backg_var2_sd   = 3,
                         site_samples    = 400,
                         N_site_bags     = 0,
                         background_site_balance = 1,
                         test_train_split = 0.75){
  back_samples <- site_samples*background_site_balance
  N_site_bags <- site_samples/10
  sites <- data.frame(var1 = rnorm(site_samples,sites_var1_mean,sites_var1_sd),
                      var2 = rnorm(site_samples,sites_var2_mean,sites_var2_sd),
                      SITENO = "Site")
  backs <- data.frame(var1 = rnorm(site_samples,backg_var1_mean,backg_var1_sd),
                      var2 = rnorm(site_samples,backg_var2_mean,backg_var2_sd),
                      SITENO = "Back")
  x_data   <- rbind(sites,backs)
  x_data_norm   <- data.frame(apply(x_data[,-3],2,scale))
  x_data_norm$SITENO = x_data$SITENO
  ### Split out sites
  site_list <- dplyr::filter(x_data_norm, SITENO == "Site") %>%
    split(sample(N_site_bags, nrow(.), replace=T))
  names(site_list) <- sample(paste0("Site",1:N_site_bags),length(site_list))
  ### Split out background
  N_back_bags <- N_site_bags * background_site_balance
  back_list <- dplyr::filter(x_data_norm, SITENO == "Back") %>%
    split(sample(N_back_bags, nrow(.), replace=T))
  names(back_list) <- sample(paste0("Back",1:N_back_bags),length(back_list))
  # Merge site and background bags together
  x_data_norm <- c(site_list, back_list)
  # Shuffle list
  x_data_norm <- sample(x_data_norm,length(x_data_norm))
  train_index <- sample(1:length(x_data_norm), length(x_data_norm) * test_train_split)
  train_dat <- x_data_norm[train_index]
  train_presence <- ifelse(grepl("Site",names(train_dat)),1,0)
  test_dat <- x_data_norm[-train_index]
  test_presence <- ifelse(grepl("Site",names(test_dat)),1,0)
  return(list(train_data = train_dat, train_presence = train_presence,
              test_data = test_dat, test_presence = test_presence))
}
