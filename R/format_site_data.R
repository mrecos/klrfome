#' format_site_data
#'
#' @param dat - [data.frame] A data.frame of presence and absence records. Column "presence" must contain presence/absence as 1/0, and column "SITENO" contains the grouping variable.
#' @param N_sites - [scaler] The number of sites to randomly select
#' @param train_test_split - [scaler] a float from 0 to 1 indicating the percent of N_sites to be used as training data.
#' @param background_site_balance - [scaler] Integer > 0 indicating how many background groups per site group
#'
#' @return - list of various ways to arrange site data and mean/sd of data
#' @import dplyr
#' @export
#'
format_site_data <- function(dat, N_sites, train_test_split, background_site_balance){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  variables <- setdiff(colnames(dat1), c("presence", "SITENO"))
  means  <- sapply(dat1[variables], mean, na.rm=T)
  sds    <- sapply(dat1[variables], sd, na.rm=T)
  dat1   <- data.frame(apply(dat1[, variables],2,scale))  # should be scaled on only training data mean/sd
  dat1   <- cbind(dat1, dat[,c("presence","SITENO")])
  ## Reduce number of sites to N_sites
  sites <- filter(dat1, presence == 1)
  site_names <- unique(sites$SITENO)
  N_sites_index <- sample(site_names, N_sites)
  sites <- filter(sites, SITENO %in% N_sites_index)
  ### Split Sites Data
  sites_train_index <- sample(N_sites_index, length(N_sites_index) * train_test_split)
  train_sites <- filter(sites, SITENO %in% sites_train_index)
  test_sites  <- filter(sites, !SITENO %in% sites_train_index)
  ### Split Background Data
  train_background <- filter(dat1, presence == 0) %>%
    sample_n(nrow(train_sites) * background_site_balance, replace = TRUE) %>%
    mutate(presence = 0)
  test_background <- filter(dat1, presence == 0) %>%
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
