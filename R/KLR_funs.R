#' KLR
#' 
#' `KLR()` Is the primary model fitting function in the `klrfome` package. This function fits a Kernel Logistic Regression (KLR) model to a similarity or distance matrix.
#' 
#' the `KLR` function takes the similarity kernel matrix `K`, a vector `presnece` of presence/absence coded as 1 or 0, and a scalar value for the `lambda` regularizing hyperparameter; optionally values for maximum iterations and threshold. This function performs Kernel Logistic Regression (KLR) via Iterative Re-weighted Least Squares (IRLS). The objective is to approximate a set of parameters that minimize the negative likelihood of the parameters given the data and response. This function returns a list of `pred`, the estimated response (probability of site-presence) for the training data, and `alphas`, the approximated parameters from the IRLS algorithm.
#'
#' @param K - [NxN] Mean embedding kernel matrix
#' @param presence - [vector] or presence = 1 or absence = 0
#' @param lambda - [scaler] Ridge regularization parameter
#' @param maxiter - [integer] Maximum iterations for IRLS algorithm
#' @param tol - [double] The convergence tolerance
#' @param verbose - [scaler] 0 = No notice; 1 = Reports the number of steps until convergence; 2 = Reports each iteration
#'
#' @return list: `pred` - predicted probabiity of positive class, `alphas` - estimated coefficients
#' @importFrom Matrix Diagonal
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
KLR <- function(K, presence, lambda, maxiter = 100, tol = 0.01, verbose=1){
  if(is.vector(K)){
    N = length(K)
  } else if(is.matrix(K)){
    N = nrow(K)
  }
  alpha = rep(1/N, N) # initial value
  iter = 1
  while(TRUE) {
    Kalpha = as.vector(K %*% alpha)
    spec = 1 + exp(-Kalpha)
    pi = 1 / spec
    diagW = pi * (1 - pi)
    # same as: Zhu & Hastie (2002:19)
    # W^-1(y-pi); W = diag[pi(1 − pi)]N×N
    # e = solve(as.matrix(Matrix::Diagonal(x=diagW))) %*% (presence - pi)
    z = Kalpha + ((presence - pi) / diagW)
    # same as: Zhu & Hastie (2002:19)
    # z = (K_a alpha(k−1) + W^−1(y − pi)) or (y-pi)/W
    # z = as.vector(K %*% alpha) + ((presence - pi) / diagW)
    alpha_new = try(solve(K + lambda * Matrix::Diagonal(x=1/diagW), z))
    if (class(alpha_new) == "try-error") {
      cat("Error in calculating solution.","\n")
      break
    }
    alphan = as.vector(alpha_new)
    if(verbose == 2){
      cat("Step ", iter, ". Absolute Relative Approximate Error = ",
          round(abs(sum(alphan - alpha)/sum(alphan))*100, 4), "\n", sep = "")
      
    }
    if (any(is.nan(alphan)) || all(abs(alphan - alpha) <= tol)) {
      if(verbose %in% c(1,2)){
        cat("Found solution in", iter, "steps.","\n")
      }
      break
    }
    else if (iter > maxiter) {
      cat("Maximum iterations for KRR Logit!", "\n",
          "May be non-optimum solution.","\n")
      break
    }
    else {
      alpha = alphan
      iter = iter + 1
    }
  }
  log_pred <-  1 / (1 + exp(-as.vector(K %*% alpha_new)))
  return(list(pred = log_pred, alphas = alpha_new))
}


#' get_k
#' 
#' `get_k()` is an internal function used by `build_k()` to calcualte similarity kernel for each pair of observations
#' 
#' This function takes two data.frames or matricies, calculates the cross distance between them, and then applies a radial basis function to the resulting distance matrix.
#'
#' @param y1 - [NxP] Matrix of data from bag i
#' @param y2 - [NxP] Matrix of data from bag j
#' @param sigma - [scaler] smoothing hyperparameters for RBF kernel
#' @param dist_metric [character] One of the distance methods from rdist::cdist. Default = "euclidean". see ?rdist::cdist
#'
#' @return Matrix G
#' @importFrom rdist cdist
#'
get_k <- function(y1,y2,sigma, dist_metric = "euclidean"){
  g = rdist::cdist(as.matrix(y1),as.matrix(y2), metric=dist_metric)
  g = exp(-g^2/(2*sigma^2)) # equal to exp(-sigma2*g^2) ; sigma2 = 1/(2*sigma^2)
  return(g)
}


#' build_k
#' 
#' `build_k()` is a primary package function that takes in the formatted list of site/background data and builds a similarity matrix suitable for computation with the `KLR()` function or direct study.
#' 
#' This function takes list of training data, scalar value for `sigma` hyperparameter, and a distance method to compute a mean embedding similarity kernel. This kernel is a pair-wise (N x N) matrix of the mean similarity between the attributes describing each site location and background group. Optional inouts are `progress` for a progress bar and `dist_metric` for the distance computation. By default, the distance metric is euclidean and should likely stay as such unless you have explored other distances and know why/how you want to use them.
#'
#' @param y1 - [list] List of site/background data formatted by `format_site_data()`
#' @param y2 - [list] Typically left blank as y2 == y1.
#' @param sigma - [scaler] smoothing hyperparameters for RBF kernel
#' @param progress - [logical] False = no progress bar; 1 = show progress bar
#' @param dist_metric [character] One of the distance methods from rdist::cdist. Default = "euclidean". see ?rdist::cdist
#'
#' @return - matrix K
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
build_K <- function(y1, y2=y1, sigma, progress = TRUE, dist_metric = "euclidean"){
  K <- matrix(nrow = length(y1), ncol = length(y2))
  if(isTRUE(progress)){
    total_iter <- sum(seq(length(y1)-1))
    pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
  }
  iter <- 0
  for(i in 1:length(y1)){
    for(j in i:length(y2)){ 
      g <- get_k(y1[[i]], y2[[j]], sigma, dist_metric = dist_metric)
      k <- round(mean(g, na.rm = TRUE),3)
      K[i,j] <- k
      if(isTRUE(progress)){setTxtProgressBar(pb, iter)}
      iter <- iter + 1
    }
  }
  if(isTRUE(progress)){close(pb)}
  K <- tri_swap(K)
  return(K)
}


#' tri_swap
#'
#' `tri_swap()` is an internal helper function used to flip the sides of a matrix
#'
#' @param m - Matrix
#'
#' @return - Matrix
#'
tri_swap <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}


#' KLR_predict
#' 
#' `KLR_predict()` is a function to predict the probability of site presence to a new list of data based on the fitted alpha parameters returned from the `KLR()` funtion.
#' 
#' This function takes a list of the `test_data`, a list of the `train_data`, a vector of the approximated alpha parameters as `alpha_pred`, a scalar value for the `sigma` kernel hyperparameter, and a distance method (deafult = "Euclidean"). This function predicts the probability of site presence for new observations based on the training data and `alphas` parameters. This is accomplished by building the `k*k` kernel matrix as the similarity between the training test data then computing the inverse logit of `k*k %*% alphas`. The output is the predicted probability of site presence for each training data example.
#'
#' @param test_data - [list] Training data used to create similarity kernel matrix
#' @param train_data - [list] Testing data to predict class
#' @param alphas_pred - [vector] Numeric vector of alpha parameters from KLR function
#' @param sigma - [scaler] Smoothing parameter for RBF kernel
#' @param progress - [logical] False = no progress bar; 1 = show progress bar
#' @param dist_metric [character] One of the distance methods from rdist::cdist. Default = "euclidean". see ?rdist::cdist
#'
#' @return - [vector] - predicted probabiity of positive class
#' @export
#'
#' @examples
#' \dontrun{
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
KLR_predict <- function(test_data, train_data, alphas_pred, sigma, progress = TRUE, dist_metric = "euclidean"){
  kstark <- matrix(nrow = length(test_data), ncol = length(train_data))
  if(isTRUE(progress)){
    total_iter <- length(test_data) * length(train_data)
    pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
  }
  iter <- 0
  for(j in 1:length(test_data)){
    for(i in 1:length(train_data)){
      g_i <- get_k(train_data[[i]],
                   test_data[[j]], sigma, dist_metric = dist_metric)
      k_i <- round(mean(g_i, na.rm = TRUE),3)
      kstark[j,i] <- k_i
      if(isTRUE(progress)){setTxtProgressBar(pb, iter)}
      iter <- iter + 1
    }
  }
  if(isTRUE(progress)){close(pb)}
  pred <- 1 / (1 + exp(-as.vector(kstark %*% alphas_pred)))
  return(pred)
}
