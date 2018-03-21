#' KLR
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
#'
KLR <- function(K, presence, lambda, maxiter = 100, tol = 0.01, verbose=1){
  # LOGISTIC - optimize alpha
  if(is.vector(K)){
    N = length(K)# NOT SURE IF THIS WORKS WITH REST OF FUNCTION
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
    # e = (presence - pi) / diagW
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
      cat("Step ", iter, ". Change in alpha parameters = ",
          round(sum(abs(alphan - alpha)),4), "\n", sep = "")
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
#' @param y1 - [NxP] Matrix of data from bag i
#' @param y2 - [NxP] Matrix of data from bag j
#' @param sigma - [scaler] smoothing hyperparameters for RBF kernel
#' @param dist_method - [object] distance method from proxy package, e.g. proxy::pr_DB$get_entry("Euclidean")
#'
#' @return Matrix G
#' @importFrom proxy dist
#' @export
#'
get_k <- function(y1,y2,sigma, dist_method = dist_method){
  g = proxy::dist(as.matrix(y1),as.matrix(y2), method = dist_method,
                  by_rows = TRUE, auto_convert_data_frames = FALSE) # speed bottle neck according to profvis
  ## my version
  g = exp(-g^2/(2*sigma^2)) # equal to exp(-sigma2*g^2) ; sigma2 = 1/(2*sigma^2)
  # g = exp(-(1/(2*sigma^2))*g^2) # equal to exp(-g^2/(2*sigma^2))
  return(g)
}


#' build_k
#'
#' @param y1 - [NxP] Matrix of data from bag i
#' @param y2 - [NxP] Matrix of data from bag j
#' @param sigma - [scaler] smoothing hyperparameters for RBF kernel
#' @param progress - [logical] False = no progress bar; 1 = show progress bar
#' @param dist_method - [object] distance method from proxy package, e.g. proxy::pr_DB$get_entry("Euclidean")
#'
#' @return - matrix K
#' @export
#'
build_K <- function(y1,y2=y1,sigma, progress = TRUE, dist_method){
  # example: K <- build_K(x_data_norm, x_data_norm, sigma)
  K <- matrix(nrow = length(y1), ncol = length(y2))
  if(isTRUE(progress)){
    total_iter <- sum(seq(length(y1)-1)) # only works for symetrical matricies
    pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
  }
  iter <- 0
  for(i in 1:length(y1)){
    for(j in i:length(y2)){  ## index to i so that only does upper triangle (only for square matrix!)
      # print(paste0(i, " : ", j))
      g <- get_k(y1[[i]], y2[[j]], sigma, dist_method)
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
#' @param test_data - [list] Training data used to create similarity kernel matrix
#' @param train_data - [list] Testing data to predict class
#' @param alphas_pred - [vector] Numeric vector of alpha parameters from KLR function
#' @param sigma - [scaler] Smoothing parameter for RBF kernel
#' @param progress - [logical] False = no progress bar; 1 = show progress bar
#' @param dist_method - [object] distance method from proxy package, e.g. proxy::pr_DB$get_entry("Euclidean")
#'
#' @return - [vector] - predicted probabiity of positive class
#' @export
#'
KLR_predict <- function(test_data, train_data, alphas_pred, sigma, dist_method = dist_method, progress = TRUE){
  kstark <- matrix(nrow = length(test_data), ncol = length(train_data))
  if(isTRUE(progress)){
    total_iter <- length(test_data) * length(train_data)
    pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
  }
  iter <- 0
  for(j in 1:length(test_data)){
    for(i in 1:length(train_data)){
      g_i <- get_k(train_data[[i]],
                   test_data[[j]], sigma, dist_method)
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
