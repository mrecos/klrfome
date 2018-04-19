get_metrics <- function(dat){
  TP <- sum(dat$pred_cat == 1 & dat$obs == 1, na.rm = TRUE)
  FP <- sum(dat$pred_cat == 1 & dat$obs == 0, na.rm = TRUE)
  TN <- sum(dat$pred_cat == 0 & dat$obs == 0, na.rm = TRUE)
  FN <- sum(dat$pred_cat == 0 & dat$obs == 1, na.rm = TRUE)
  metrics   <- suppressWarnings(metrics(TP,TN,FP,FN))
  return(metrics)
}
data_helper <- function(params,dat){
  cat("Getting Data...","\n")
  dat_split <- format_site_data(dat,
                                params$N_sites,
                                params$train_test_split,
                                params$background_site_balance)
}
logreg_helper <- function(dat_split){
  lr_data <- dat_split$tbl_train_data %>%
    dplyr::select(.,-SITENO)
  lr <- glm(presence ~ ., data = lr_data, family = "binomial")
}
SVM_helper <- function(dat_split, cost = 0.001){
  svm_data <- dat_split$tbl_train_data %>%
    dplyr::select(.,-SITENO)
  svm <- svm(presence ~ .,
             data = svm_data, cost = cost,
             family = "binomial")
}
KRR_helper <- function(dat_split, sigma, lambda, dist_metric){
  krr_data <- dat_split$train_data
  K <- build_K(krr_data, krr_data, sigma, dist_method = dist_metric, progres = TRUE)
  #### Train
  krr <- KLR(K, dat_split$train_presence, lambda, 100, 0.001, verbose = 0)
  return(krr)
}
predict_helper <- function(model,splits,model_type,
                           confusion_matrix_cutoff=0.5,Y_threshold=FALSE){
  if(model_type %in% c("LR","SVM")){
    newdata <- splits$tbl_test_data
    preds <- predict(model, newdata = newdata, type = "response")
    predicted <- data.frame(pred = preds,
                            obs  = splits$tbl_test_presence)
  }
  if(model_type == "KRR"){
    alphas_pred   <- model[["alphas"]]
    newdata    <- splits$test_data
    train_data <- splits$train_data
    preds <- KLR_predict(newdata, train_data, alphas_pred, sigma, dist_method = dist_metric)
    predicted <- data.frame(pred = preds,
                            obs  = splits$test_presence)
  }
  if(isTRUE(Y_threshold)){
    confusion_matrix_cutoff <- get_YJ_max_threshold(predicted)
  }
  predicted$pred_cat <- ifelse(predicted$pred >= confusion_matrix_cutoff,1,0)
  return(predicted)
}
get_AUC <- function(preds){
  roc_obj <- pROC::roc(preds$obs, preds$pred)
  AUC <- pROC::auc(roc_obj)
}
get_YJ_max_threshold <- function(preds){
  roc_obj <- pROC::roc(preds$obs, preds$pred)
  sens <- roc_obj$sensitivities
  spec <- roc_obj$specificities
  YJ   <- sens + spec - 1
  threshold <- roc_obj$thresholds[which.max(YJ)]
  return(threshold)
}


library("corrplot")
library("pROC")
library("data.table")
library("tidyverse")
library("e1071")         # for SVM comparison
library("klrfome")
library("future")

sigma <- 1
lambda <- 0.11

### Data parameters
N_back_bags = 50 # need to figure out better way to measure this
N_sites     = 50
background_site_balance = 1
sample_fraction = 0.50
train_test_split = 0.75
confusion_matrix_cutoff = 0.5
dist_metric = "euclidean"

physioshed_ids <- c(1,2,6,8,12)

for(z in seq_along(physioshed_ids)){
  physioshed_z <- physioshed_ids[z]
  physioshed   <- paste0("r91_all_upland_section_",physioshed_z, "_regression_data_SITENO.csv")
  data_location = file.path("C:/Users/matthew.d.harris/Dropbox/R/PASS_regression",physioshed)
  # data_location = file.path("/home/rstudio/",physioshed)

  ### Load Data
  dat <- fread(data_location)
  dat <- data.frame(dat)
  potential_vars <- c("presence", "SITENO",
                      "ed_h6", "std_32c", "slpvr_32c", "slpvr_16c", "rng_16c",
                      "ed_h2", "cd_conf", "elev_2_drainh")
  dat1 <- dat[,which(names(dat) %in% potential_vars)]
  ### Center and Standardize data
  variables <- setdiff(colnames(dat1), c("presence", "SITENO"))
  dat1   <- data.frame(apply(dat1[, variables],2,scale))
  dat1   <- cbind(dat1, dat[,c("presence","SITENO")])

  ###### parallel fitting and model compare
  message(paste0("starting maps for ",physioshed,"\n"))
  searches <- 2
  model_data <- seq_len(searches) %>%
    data_frame(
      id = .,
      N_sites = N_sites,
      train_test_split = train_test_split,
      background_site_balance = background_site_balance
    ) %>%
    nest(-id, .key = "params") %>%
    mutate(splits  = map(params,  ~future::future(data_helper(.x,dat1)))) %>%
    mutate(splits  = map(splits,  ~future::value(.x)))
  model_fits <- model_data %>%
    mutate(model_LR    = map(splits,     ~future::future(logreg_helper(.x))),
           model_SVM   = map(splits,     ~future::future(SVM_helper(.x, cost = 0.001))),
           model_KRR   = map(splits,     ~future::future(KRR_helper(.x, sigma, lambda, dist_metric)))) %>%
    mutate(model_LR    = map(model_LR,   ~future::value(.x)),
           model_SVM   = map(model_SVM,  ~future::value(.x)),
           model_KRR   = map(model_KRR,  ~future::value(.x)))
  model_preds <- model_fits %>%
    mutate(preds_LR    = map2(model_LR, splits, ~predict_helper(.x,.y,model_type="LR",
                                                                Y_threshold = TRUE)),
           preds_SVM   = map2(model_SVM, splits, ~predict_helper(.x,.y,model_type="SVM",
                                                                 Y_threshold = TRUE)),
           preds_KRR   = map2(model_KRR, splits, ~predict_helper(.x,.y,model_type="KRR",
                                                                 Y_threshold = TRUE))) %>%
    dplyr::select(-splits)
  model_metrics <- model_preds %>%
    mutate(AUC_LR      = map_dbl(preds_LR, get_AUC),
           metric_LR   = map(preds_LR, get_metrics),
           ROC_LR      = map(preds_LR, get_YJ_max_threshold),
           AUC_SVM     = map_dbl(preds_SVM, get_AUC),
           metric_SVM  = map(preds_SVM, get_metrics),
           ROC_SVM     = map(preds_SVM, get_YJ_max_threshold),
           AUC_KRR     = map_dbl(preds_KRR, get_AUC),
           metric_KRR  = map(preds_KRR, get_metrics),
           ROC_KRR     = map(preds_KRR, get_YJ_max_threshold))

  model_preds_strip <- model_preds %>%
    dplyr::select(-starts_with("model"))
  model_metrics_strip <- model_metrics %>%
    dplyr::select(-starts_with("model"), -starts_with("preds"))

  metrics_AUC_long <- model_metrics_strip %>%
    dplyr::select(-starts_with("ROC"),-starts_with("metric"),-params) %>%
    gather(model, value, -id) %>%
    separate(model, into = c("metrics","model"), sep = "_")

  metrics_threshold_long <- model_metrics_strip %>%
    select(-starts_with("ROC"),-starts_with("AUC"),-params)  %>%
    mutate(metrics = map(metric_LR, names)) %>%
    mutate_if(is.list, simplify_all) %>%
    unnest()  %>%
    gather(model,value,-id,-metrics) %>%
    separate(model,into=c("a","model"),sep="_") %>%
    select(-a)

  save_folder <- paste0(searches,"_KRR_LR_SVM_compare_r91U",physioshed_z)
  if(!dir.exists(save_folder)){
    dir.create(save_folder)
  }
  ## save(model_fits,    file = file.path(save_folder,"model_fits.RData"))
  save(model_preds_strip,   file = file.path(save_folder,"model_preds_srtip.RData"))
  save(model_metrics_strip, file = file.path(save_folder,"model_metrics_strip.RData"))
  write.csv(metrics_threshold_long, file = file.path(save_folder,"metrics_threshold_long.csv"))
  write.csv(metrics_AUC_long, file = file.path(save_folder,"metrics_AUC_long.csv"))

  rm(model_fits)
  rm(dat)
  gc()
}
