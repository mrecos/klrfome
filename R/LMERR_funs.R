# active Developments
KRR_logit <- function(K,y, lambda){
  #### Logistic KRR
  # KRR_logit(K,presence)
  if(is.vector(K)){
    N = length(K)
  } else if(is.matrix(K)){
    N = nrow(K)
  }
  alpha = rep(1/N, N)
  Kalpha = as.vector(K %*% alpha)
  spec = 1 + exp(-Kalpha)
  pi = 1 / spec
  diagW = pi * (1 - pi)
  e = (y - pi) / diagW
  q = Kalpha + e
  ident.N <- diag(rep(1,N)) # added by me
  theSol = solve(K + lambda * ident.N, q)
  log_pred <- 1 / (1 + exp(-as.vector(K %*% theSol)))
  return(list(pred = log_pred, alphas = theSol))
}
KRR_logit_optim <- function(K, presence, lambda, maxiter = 100, tol = 0.01){
  # LOGISTIC - optimize alpha
  # example: KRR_logit_optim(K,presence, 100, 0.01)
  if(is.vector(K)){
    N = length(K)# NOT SURE IF THIS WORKS WITH REST OF FUNCTION
  } else if(is.matrix(K)){
    N = nrow(K)
  }
  alpha = rep(1/N, N) # initial value
  iter = 1
  while (TRUE) {
    Kalpha = as.vector(K %*% alpha)
    spec = 1 + exp(-Kalpha)
    pi = 1 / spec
    diagW = pi * (1 - pi)
    e = (presence - pi) / diagW
    q = Kalpha + e
    theSol = try(solve(K + lambda * Diagonal(x=1/diagW), q))
    if (class(theSol) == "try-error") {
      break
    }
    alphan = as.vector(theSol)
    if (any(is.nan(alphan)) || all(abs(alphan - alpha) <= tol)) {
      break
    }
    else if (iter > maxiter) {
      cat("klogreg:maxiter!")
      break
    }
    else {
      alpha = alphan
      iter = iter + 1
    }
  }
  log_pred <-  1 / (1 + exp(-as.vector(K %*% theSol)))
  return(list(pred = log_pred, alphas = theSol))
}
get_k <- function(y1,y2,sigma, dist_method = "Euclidean"){
  g = proxy::dist(as.matrix(y1),as.matrix(y2), method = dist_method,
                  by_rows = TRUE, auto_convert_data_frames = FALSE) # speed bottle neck accordingn to profvis
  g = exp(-g^2 / (2 * sigma^2))
  return(g)
}
build_K <- function(y1,y2,sigma, progress = TRUE, ...){
  # example: K <- build_K(x_data_norm, x_data_norm, sigma)
  K <- matrix(nrow = length(y1), ncol = length(y2))
  if(isTRUE(progress)){
    total_iter <- sum(seq(length(y1)-1)) # only works for symetrical matri
    pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
  }
  iter <- 0
  for(i in 1:length(y1)){
    for(j in i:length(y2)){  ## index to i so that only does upper triangle
      # print(paste0(i, " : ", j))
      g <- get_k(y1[[i]],
                 y2[[j]], sigma, ...)
      k <- round(mean(g),3)
      K[i,j] <- k
      if(isTRUE(progress)){setTxtProgressBar(pb, iter)}
      iter <- iter + 1
    }
  }
  if(isTRUE(progress)){close(pb)}
  K <- tri_swap(K)
  return(K)
}
tri_swap <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

KRR_logit_predict <- function(test_data, train_data, alphas_pred, sigma, dist_method = "Euclidean", progress = TRUE){
  # example: KRR_logit_predict(test_dat, train_dat, theSol, sigma)
  pred_yhat <- matrix(nrow = length(test_data), ncol = length(train_data))
  if(isTRUE(progress)){
    total_iter <- length(test_data) * length(train_data)
    pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
  }
  iter <- 0
  for(j in 1:length(test_data)){
    for(i in 1:length(train_data)){
      g_i <- get_k(train_data[[i]],
                   test_data[[j]], sigma, dist_method = dist_method)
      k_i <- round(mean(g_i),3)
      pred_yhat[j,i] <- k_i
      if(isTRUE(progress)){setTxtProgressBar(pb, iter)}
      iter <- iter + 1
    }
  }
  if(isTRUE(progress)){close(pb)}
  pred <- 1 / (1 + exp(-as.vector(pred_yhat %*% alphas_pred)))
  return(pred)
}
