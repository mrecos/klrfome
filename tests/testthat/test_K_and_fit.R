context("Modeling functions")
library(klrfome)

test_that("build_k returns correct matrix", {
  set.seed(717)
  sigma = 0.5
  lambda = 0.1
  dist_metric = "euclidean"
  sim_data <- get_sim_data(site_samples = 800, N_site_bags = 75)
  formatted_data <- format_site_data(sim_data, N_sites=10, train_test_split=0.8,
                                     sample_fraction = 0.9, background_site_balance=1)
  train_data <- formatted_data[["train_data"]]
  train_presence <- formatted_data[["train_presence"]]
  test_data <- formatted_data[["test_data"]]
  test_presence <- formatted_data[["test_presence"]]
  
  K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric, progress = FALSE)
  
  # check to make sure K is a matrix
  expect_is(K, "matrix")
  # check to make sure K is square
  expect_true(nrow(K) == ncol(K))
  # check to see if all values of K are zero or greater
  expect_true(all(K>=0))
})

test_that("KLR outputs predictions and weights", {
  set.seed(717)
  sigma = 0.5
  lambda = 0.1
  dist_metric = "euclidean"
  sim_data <- get_sim_data(site_samples = 800, N_site_bags = 75)
  formatted_data <- format_site_data(sim_data, N_sites=10, train_test_split=0.8,
                                     sample_fraction = 0.9, background_site_balance=1)
  train_data <- formatted_data[["train_data"]]
  train_presence <- formatted_data[["train_presence"]]
  test_data <- formatted_data[["test_data"]]
  test_presence <- formatted_data[["test_presence"]]
  
  K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric, progress = FALSE)
  train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 0)
  
  # check to make sure train_log_pred is a list
  expect_is(train_log_pred, "list")
  # check to make sure train_log_pred is length two
  expect_length(train_log_pred, 2)
  # check to see that predictions are greater or equal to zero
  expect_true(all(train_log_pred[["pred"]] >= 0))
  # check to see elements of train_log_pred are same length
  expect_true(length(train_log_pred[[1]]) == length(train_log_pred[[2]]))
})

test_that("KLR_predict outputs", {
  set.seed(717)
  sigma = 0.5
  lambda = 0.1
  dist_metric = "euclidean"
  sim_data <- get_sim_data(site_samples = 800, N_site_bags = 75)
  formatted_data <- format_site_data(sim_data, N_sites=10, train_test_split=0.8,
                                     sample_fraction = 0.9, background_site_balance=1)
  train_data <- formatted_data[["train_data"]]
  train_presence <- formatted_data[["train_presence"]]
  test_data <- formatted_data[["test_data"]]
  test_presence <- formatted_data[["test_presence"]]
  
  K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric, progress = FALSE)
  train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 0)
  test_log_pred <- KLR_predict(test_data, train_data, dist_metric = dist_metric,
                               train_log_pred[["alphas"]], sigma, progress = FALSE)
  
  # check to make sure test_log_pred is numeric
  expect_is(test_log_pred, "numeric")
  # check to make sure test_log_pred length = test_data length
  expect_true(length(test_log_pred) == length(test_data))
  # check to see that predictions are greater or equal to zero
  expect_true(all(test_log_pred >= 0))
})



