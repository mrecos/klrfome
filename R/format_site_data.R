#' format_site_data
#'
#'`format_site_data()` formats a data.frame of a specific format into a list suitable for use with the klrfome model.
#'
#'The function takes a data.frame that at a mimimum must include a column for site presence/absence, a column for the site identifier, and one or more columns of covariates to be used in the regression model. There are three additional required parameters to this function and they are; `N_sites` which can be used to limit the number of sites (groups of the presence == 1 class) returned in the results; `train_test_split` that is used to split the site present data into training and testing data sets; and `background_site_balance` that is the ratio of background observation groups to include for each site group. The choice of `N_sites` and `background_site_balance` argument values are influenced by a number of factors including the amount of data, length of computation, and site prevelance. 
#'
#'This function returns a list object which is a reformatting of the input data where each list element if a grouping of present or absence observations. This is needed because the klrfome model works on groups of observations.
#'
#'
#' @param dat - [data.frame] A data.frame of presence and absence records. Column "presence" must contain presence/absence as 1/0, and column "SITENO" contains the grouping variable.
#' @param N_sites - [scalar] The number of sites to randomly select for analysis
#' @param train_test_split - [scalar] a float from 0 to 1 indicating the percent of N_sites to be used as training dataset vs testing dataset.
#' @param background_site_balance - [scalar] Integer > 0 indicating how many background groups per site group
#'
#' @return - list of various ways to arrange site data and mean/sd of data
#' @import dplyr
#' @export
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
format_site_data <- function(dat, N_sites, train_test_split, background_site_balance, sample_fraction){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!("SITENO" %in% names(dat) & "presence" %in% names(dat))) {
    stop("Data must contain colums named 'presence' and 'SITENO'. See ?format_site_data",
         call. = FALSE)
  }
  variables <- setdiff(colnames(dat), c("presence", "SITENO"))
  means  <- sapply(dat[variables], mean, na.rm=T)
  sds    <- sapply(dat[variables], sd, na.rm=T)
  dat_c   <- data.frame(apply(dat[, variables],2,scale))  # should be scaled on only training data mean/sd
  dat   <- cbind(dat_c, dat[,c("presence","SITENO")])
  N_back_bags <- N_sites * background_site_balance
  ## Reduce number of sites to N_sites
  sites <- filter(dat, presence == 1)
  site_names <- unique(sites$SITENO)
  N_sites_index <- sample(site_names, N_sites)
  sites <- filter(sites, SITENO %in% N_sites_index)
  ### Split Sites Data
  sites_train_index <- sample(N_sites_index, length(N_sites_index) * train_test_split)
  train_sites <- filter(sites, SITENO %in% sites_train_index)
  test_sites  <- filter(sites, !SITENO %in% sites_train_index)
  ### Split Background Data
  train_background <- filter(dat, presence == 0) %>%
    sample_n(nrow(train_sites) * background_site_balance, replace = TRUE) %>%
    mutate(presence = 0)
  test_background <- filter(dat, presence == 0) %>%
    sample_n(nrow(test_sites) * background_site_balance, replace = TRUE) %>%
    mutate(presence = 0)
  ### Tablular data - REDUCE BY [sample_fraction]
  tbl_train_data  <- rbind(train_sites, train_background) %>%
    sample_frac(size = sample_fraction)
  tbl_train_presence <- dplyr::select(tbl_train_data,presence)
  tbl_train_presence <- as.numeric(tbl_train_presence$presence)
  tbl_test_data   <- rbind(test_sites, test_background) %>%
    sample_frac(size = sample_fraction)
  tbl_test_presence <- dplyr::select(tbl_test_data,presence)
  tbl_test_presence <- as.numeric(tbl_test_presence$presence)
  ### Split out background - Still wonky, but works for now # NEEDS ATTENTION
  train_back_list <- dplyr::filter(tbl_train_data, SITENO == "background") %>%
    dplyr::select(-presence,-SITENO) %>%
    split(sample(N_back_bags, nrow(.), replace=T))
  names(train_back_list) <- sample(paste0("background",1:N_back_bags),length(train_back_list))
  train_site_list <- dplyr::filter(tbl_train_data, SITENO != "background") %>%
    split(f = .$SITENO ) %>%
    lapply(., function(x) x[!(names(x) %in% c("presence", "SITENO"))])
  test_site_list <-  group_by(tbl_test_data, SITENO) %>%
    mutate(id = paste0(SITENO, "_", seq_len(n()))) %>%
    split(f = .$id) %>%
    lapply(., function(x) x[!(names(x) %in% c("presence", "SITENO", "id"))])
  # Merge site and background bags together
  train_data <- c(train_site_list, train_back_list)
  # don't need to split background and site lists, so no need to c()
  test_data <- test_site_list
  # Shuffle list
  train_data <- sample(train_data,length(train_data))
  test_data  <- sample(test_data,length(test_data))
  train_presence <- ifelse(grepl("background", names(train_data)),0,1)
  test_presence  <- ifelse(grepl("background", names(test_data)),0,1)
  return(list(train_data = train_data,
              test_data = test_data,
              train_presence = train_presence,
              test_presence = test_presence,
              tbl_train_data = tbl_train_data,
              tbl_train_presence = tbl_train_presence,
              tbl_test_data = tbl_test_data,
              tbl_test_presence = tbl_test_presence,
              means = means,
              sds = sds))
}
