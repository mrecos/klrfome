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

###### FROM my METRICS FUNCTIONS
#### NEEDS REFRESH WITH UDATED KG CALC!!
cohens_kappa <- function(TP,TN,FP,FN){
  A <- TP
  B <- FP
  C <- FN
  D <- TN
  Po <- (A+D)/(A+B+C+D)
  Pe_a <- ((A+B)*(A+C))/(A+B+C+D)
  Pe_b <- ((C+D)*(B+D))/(A+B+C+D)
  Pe <- (Pe_a+Pe_b)/(A+B+C+D)
  k <- (Po-Pe)/(1-Pe)
  return(k)
}
get_metric <- function(TP,TN,FP,FN,metric_type){
  m <- metrics(TP,TN,FP,FN)
  m <- m[[metric_type]]
  return(m)
}
##### NEED TO UPDATE WITH NEW KG CALCs
metrics <- function(TP,TN,FP,FN){
  # metrics derived from TP,TN,FP, and FN
  ### "Summary" scores
  # Sens, Spec, Precision, Recall calculated here and reused below
  Sensitivity <- TP/(TP+FN) # TPR, Recall
  Specificity <- TN/(FP+TN) # TNR
  Precision   <- TP/(TP+FP) # PPV
  Recall      <- TP/(TP+FN) # Sensitivity, TPR
  pred        <- c(rep(1,TP),rep(1,FP),rep(0,FN),rep(0,TN))
  obs         <- c(rep(1,TP),rep(0,FP),rep(1,FN),rep(0,TN))
  metrics <- list(
    Sensitivity = Sensitivity, # TPR, Recall
    Specificity = Specificity, # TNR
    Prevalence  = (TP+FN)/(TP+TN+FP+FN),
    Accuracy    = (TP+TN)/(TP+TN+FP+FN),
    Err_Rate    = (FP+FN)/(TP+TN+FP+FN),
    other       = (TN+FP)/(TP+TN+FP+FN), # % all area with no sites
    Precision   = Precision, # PPV
    Recall      = Recall, # Sensitivity, TPR
    F_Measure   = (2*Precision*Recall)/(Precision+Recall),
    Geo_Mean    = sqrt(TP*TN),
    FPR         = 1 - Specificity, #Fall-Out
    FNR         = 1 - Sensitivity,
    TPR         = Sensitivity, # Sensitivity, Recall
    TNR         = Specificity, # Specificity
    FOR         = FN/(FN+TN),
    FDR         = FP/(TP+FP),
    Power       = 1-(1-Sensitivity),
    LRP         = Sensitivity/(1-Specificity), # TPR/FPR
    log_LRP     = log10(Sensitivity/(1-Specificity)),
    LRN         = (1-Sensitivity)/Specificity, # FNR/TNR
    PPV         = TP/(TP+FP), # Precision
    NPV         = TN/(FN+TN),
    KG          = 1-((1-Specificity)/Sensitivity), # 1-(FPR/TPR)
    DOR         = (Sensitivity/(1-Specificity)/((1-Sensitivity)/Specificity)), # LRP/LRN
    log_DOR     = log10((Sensitivity/(1-Specificity)/((1-Sensitivity)/Specificity))),
    # D & S - https://en.wikipedia.org/wiki/Diagnostic_odds_ratio
    D           = boot::logit(Sensitivity) - boot::logit(1-Specificity),
    S           = boot::logit(Sensitivity) + boot::logit(1-Specificity),
    Kappa       = cohens_kappa(TP,TN,FP,FN),
    # Kappa agrees with http://terpconnect.umd.edu/~dchoy/thesis/Kappa/#
    # http://aircconline.com/ijdkp/V5N2/5215ijdkp01.pdf
    Opp_Precision = ((TP+TN)/(TP+TN+FP+FN))-(abs(Specificity-Sensitivity)/(Specificity+Sensitivity)),
    # https://en.wikipedia.org/wiki/Precision_and_recall
    # MCC         = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),
    Informedness  = Sensitivity+Specificity-1, #TSS, Younden's J
    Markedness    = (TP/(TP+FP))+(TN/(FN+TN))-1,
    # http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2006.01214.x/full
    # TSS         = Sensitivity+Specificity-1 #Informedness
    # TEST = 1-(FNR/Spec) or 1-(FNR/TNR)
    # one minus % sites incorrect divided by % site not-likely background
    # also, 1 minus (% misclassifications / % non-site area)
    # Reach = TEST
    Reach         = 1-((1-Sensitivity)/Specificity),
    # Reach2        = 1-((1-(TP/(TP+FN)))/(TN/(FP+TN)))
    ## From Verhagen(2007; 121)
    AFK           = suppressWarnings(sqrt(Sensitivity*((Sensitivity-(1-Specificity))/((TN+FP)/(TP+TN+FP+FN))))),
    Indicative      = Sensitivity-(1-Specificity),
    Indicative_norm = (Sensitivity-(1-Specificity))/((TN+FP)/(TP+TN+FP+FN)),
    Brier         = mean((obs-pred)^2), # MSE for binary class problems
    # adding more
    MAE           = mean(abs(pred-obs)),
    RMSE          = sqrt(mean((pred-obs)^2))

  )
  return(metrics)
}

format_site_data <- function(dat, N_sites, train_test_split, background_site_balance){
  library("dplyr")
  variables <- setdiff(colnames(dat1), c("presence", "SITENO"))
  dat1   <- data.frame(apply(dat1[, variables],2,scale))
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
              tbl_test_presence = tbl_test_presence))
}

make_xstats <- function(results){
  library("pROC")
  xstats <- group_by(results, rep, model) %>%
    summarise(TP = sum(pred_cat == 1 & obs == 1, na.rm = TRUE),
              FP = sum(pred_cat == 1 & obs == 0, na.rm = TRUE),
              TN = sum(pred_cat == 0 & obs == 0, na.rm = TRUE),
              FN = sum(pred_cat == 0 & obs == 1, na.rm = TRUE),
              auc = pROC::auc(obs,pred, type = "linear")) %>%
    group_by(rep) %>%
    dplyr::mutate(Reach = get_metric(TP,TN,FP,FN,"Reach"),
                  KG = get_metric(TP,TN,FP,FN,"KG"),
                  Sensitivity = get_metric(TP,TN,FP,FN,"Sensitivity"),
                  `1-Specificity` = 1-get_metric(TP,TN,FP,FN,"Specificity"),
                  avg_metric = (KG + Reach)/2) %>%
    data.frame()
  return(xstats)
}
