# K <- build_K(train_data, train_data, sigma, dist_method = method_object)
y1 = train_data
y2 = train_data

progress <- TRUE
K <- matrix(nrow = length(y1), ncol = length(y2))
total_iter <- sum(seq(length(y1)-1)) # only works for symetrical matri
pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
iter <- 0
for(i in 1:length(y1)){
  for(j in i:length(y2)){  ## index to i so that only does upper triangle
    # print(paste0(i, " : ", j))
    g <- get_k(y1[[i]],
               y2[[j]], sigma,  dist_method = dist_method)
    k <- round(mean(g),3)
    K[i,j] <- k
    if(isTRUE(progress)){setTxtProgressBar(pb, iter)}
    iter <- iter + 1
  }
}

sigma = 1
sigmal = 1/(2*sigma^2)
g1 = proxy::dist(as.matrix(y1[[i]]),
                 as.matrix(y2[[j]]), dist_method = dist_method,
                by_rows = TRUE, auto_convert_data_frames = FALSE) # speed bottle neck accordingn to profvis
g1b = mydist(as.matrix(y1[[i]]),as.matrix(y2[[j]]))
g1_mat <- matrix(as.numeric(g1), nrow = nrow(g1))
g2 = exp(-g1^2 / (2 * sigma^2))
round(mean(g2),3)

library("kernlab")
library("mmpp")
rbf <- rbfdot(sigma = sigmal)
g2_k <- kernelMatrix(x = g1_mat, kernel = rbf)
g2_r <- rdetools::rbfkernel(g1_mat, sigma = sigma)
g2_f <- exp(-sigmal*g1_mat^2)
g2_m <- k2d(g1_mat, direction = "d2k", method = "exp", scale = sigma, pos = TRUE)

mydist <- function(a,b){
  la <- nrow(a)
  lb <- nrow(b)
  csscrrm <- matrix(nrow=la,ncol=lb)
  for(i in 1:la)
    for(j in 1:lb)
      csscrrm[i,j] <- sqrt(sum((a[i,] - b[j,])^2))
  return(csscrrm)
}
euclidean_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}
library("compiler")
cmp_dist = cmpfun(mydist)

# Testing equality between results
library("rbenchmark")
reps <- 10000
benchmark(mydist(as.matrix(y1[[i]]),as.matrix(y2[[j]])),
          proxy::dist(as.matrix(y1[[i]]),as.matrix(y2[[j]]),
                      method = dist_method, by_rows = TRUE,
                      auto_convert_data_frames = FALSE),
          cmp_dist(as.matrix(y1[[i]]),as.matrix(y2[[j]])),
          outer(as.matrix(y1[[i]]),as.matrix(y2[[j]]),Vectorize(euclidean_distance)),
          pdist::pdist(as.matrix(y1[[i]]),as.matrix(y2[[j]])),
          wordspace::dist.matrix(as.matrix(y1[[i]]),as.matrix(y2[[j]]),
                                method="euclidean") ,
          columns=c("test", "elapsed", "relative"),
          order="relative", replications=reps)

## Lookign into purr
get_mean_sim <- function(y1,y2, sigma, dist_method){
  g1 = proxy::dist(as.matrix(y1),as.matrix(y2), dist_method = dist_method,
                   by_rows = TRUE, auto_convert_data_frames = FALSE)
  g2 = exp(-(1/(2*sigma^2))*g1^2)
  k_ij <- round(mean(g2),3)
  return(k_ij)
}
purrr_K <- function(y1,y2,sigma,dist_method){
  cx_list <- purrr::cross2(y1,y2)
  l1 <- lapply(cx_list, '[[', 1)
  l2 <- lapply(cx_list, '[[', 2)
  K2 <- purrr::map2(l1, l2, ~get_mean_sim(.x,.y,sigma,dist_method)) %>%
    unlist() %>%
    matrix(., nrow =length(y1))
  return(K2)
}

K_1 <- build_K(y1,y2,sigma,TRUE,dist_method)
K_2 <- purrr_K(y1,y2,sigma,dist_method)

reps <- 1000
benchmark(build_K(y1,y2,sigma,FALSE,dist_method),
          purrr_K(y1,y2,sigma,dist_method),
          columns=c("test", "elapsed", "relative"),
          order="relative", replications=reps)
