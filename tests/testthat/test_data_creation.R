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
