context("Data creation and formatting")
library(klrfome)

test_that("get_sim_data returns proper data", {
  stes_sample_n = 100
  balance = 1
  sim_data_test <- get_sim_data(site_samples = stes_sample_n, N_site_bags = 75, background_site_balance = balance)
  # check to make sure the number of rows is the number of sites when balance == 1
  expect_equal(nrow(sim_data_test), stes_sample_n*2)
  # check to make sure presence column has both and only c(1,0)
  expect_equal(unique(sim_data_test$presence), c(1,0))
  # check to see if both "background" and "Site*" are present in SITENO
  expect_true(TRUE %in% grepl("background", sim_data_test$SITENO) & TRUE %in% grepl("Site", sim_data_test$SITENO))
})


test_that("format_site_data errors and returns proper data", {
  set.seed(717)
  N_sites = 5
  dat_test_good <- data.frame(presence = rep(c(0,1),times=N_sites), 
                              var1     = rnorm(N_sites*2,0,1)) %>% 
    rowwise() %>% 
    mutate(SITENO   = if_else(presence == 1,
                              sample(LETTERS,1),
                              sample(paste0("background",sample(1:100)),1)))
  # column name "presence" missing
  dat_test_bad  <- data.frame(BAD_NAME = rep(c(0,1),times=N_sites), 
                              var1     = rnorm(N_sites*2,0,1)) %>% 
    rowwise() %>% 
    mutate(SITENO   = if_else(BAD_NAME == 1,
                              sample(LETTERS,1),
                              sample(paste0("background",sample(1:100)),1)))
  train_test_split = 0.5
  background_site_balance = 1
  sample_fraction = 1
  real_data_test <- format_site_data(dat_test_good, N_sites, 
                                     train_test_split, background_site_balance, 
                                     sample_fraction)
  # check to make sure its a list
  expect_is(real_data_test, "list")
  # check to make sure list length is 10
  expect_length(real_data_test,10)
  # expect error on data with bad field names
  expect_error(format_site_data(dat_test_bad, N_sites, 
                                train_test_split, background_site_balance, 
                                sample_fraction))
})
