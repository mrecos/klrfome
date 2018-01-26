get_metrics <- function(dat){
  TP <- sum(dat$pred_cat == 1 & dat$obs == 1, na.rm = TRUE)
  FP <- sum(dat$pred_cat == 1 & dat$obs == 0, na.rm = TRUE)
  TN <- sum(dat$pred_cat == 0 & dat$obs == 0, na.rm = TRUE)
  FN <- sum(dat$pred_cat == 0 & dat$obs == 1, na.rm = TRUE)
  inf   <- suppressWarnings(metrics(TP,TN,FP,FN)$Informedness)
  sens  <- suppressWarnings(metrics(TP,TN,FP,FN)$Sensitivity)
  spec2 <- suppressWarnings(1-metrics(TP,TN,FP,FN)$Specificity)
  list(Informedness = inf, Sensitivity = sens, FPR = spec2)
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
KRR_helper <- function(dat_split, sigma, lambda, dist_method){
  krr_data <- dat_split$train_data
  K <- build_K(krr_data, krr_data, sigma, dist_method = dist_method, progres = TRUE)
  diag(K) <- 1
  #### Train
  krr <- KRR_logit_optim(K, dat_split$train_presence, lambda, 100, 0.001, verbose = 0)
  return(krr)
}
predict_helper <- function(model,splits,model_type,confusion_matrix_cutoff=0.5){
  if(model_type %in% c("LR","SVM")){
    newdata <- splits$tbl_test_data
    preds <- predict(model, newdata = newdata, type = "response")
    predicted <- data.frame(pred = preds,
                            obs  = splits$tbl_test_presence,
                            pred_cat = ifelse(preds >= confusion_matrix_cutoff,1,0))
  }
  if(model_type == "KRR"){
    method_object <- proxy::pr_DB$get_entry("Euclidean")
    alphas_pred   <- model[["alphas"]]
    newdata    <- splits$test_data
    train_data <- splits$train_data
    preds <- KRR_logit_predict(newdata, train_data, alphas_pred, sigma, dist_method = method_object)
    predicted <- data.frame(pred = preds,
                            obs  = splits$test_presence,
                            pred_cat = ifelse(preds >= confusion_matrix_cutoff,1,0))
  }
  return(predicted)
}

library("corrplot")
library("latex2exp")
library("data.table")
library("tidyverse")
library("e1071")         # for SVM comparison
library("DistRegLMERR")
library("future")

method_object <- proxy::pr_DB$get_entry("Euclidean")
sigma <- 1
lambda <- 0.11
searches <- 1

### Data parameters
N_back_bags = 50 # need to figure out better way to measure this
N_sites     = 50
background_site_balance = 1
sample_fraction = 0.50
train_test_split = 0.75
confusion_matrix_cutoff = 0.5
# data_location = "data/r91_all_upland_section_6_regression_data_SITENO.csv"
# data_location = "/Users/mattharris/Dropbox/R/PASS_regression/r91_all_upland_section_12_regression_data_SITENO.csv"
data_location = "C:/Users/matthew.d.harris/Dropbox/R/PASS_regression/r91_all_upland_section_12_regression_data_SITENO.csv"
### Load Data
# dat <- fread(data_location)
dat <- data.frame(dat)
dat1 <- dplyr::select(dat,presence, SITENO,
                      ed_h6, std_32c, ed_h7, slpvr_32c, ed_h2, cd_conf, ed_h2, elev_2_conf)
# elev_2_strm, e_hyd_min, ed_drnh, elev_2_drainh,
# tri_16c, cd_h5)
### Center and Standardize data
variables <- setdiff(colnames(dat1), c("presence", "SITENO"))
dat1   <- data.frame(apply(dat1[, variables],2,scale))
dat1   <- cbind(dat1, dat[,c("presence","SITENO")])

###### parallel fitting and model compare
model_data <- seq_len(searches) %>%
  data_frame(
    id = .,
    N_sites = N_sites,
    train_test_split = train_test_split,
    background_site_balance = background_site_balance
    ) %>%
  nest(-id, .key = "params") %>%
  # mutate(splits = map(params, ~data_helper(.x,dat1)))
  mutate(splits  = map(params,  ~future::future(data_helper(.x,dat1)))) %>%
  mutate(splits  = map(splits,  ~future::value(.x)))
model_fits <- model_data %>%
  mutate(model_LR    = map(splits,     ~future::future(logreg_helper(.x))),
         model_SVM   = map(splits,     ~future::future(SVM_helper(.x, cost = 0.001))),
         model_KRR   = map(splits,     ~future::future(KRR_helper(.x, sigma, lambda, method_object)))) %>%
  mutate(model_LR    = map(model_LR,   ~future::value(.x)),
         model_SVM   = map(model_SVM,  ~future::value(.x)),
         model_KRR   = map(model_KRR,  ~future::value(.x)))
model_preds <- model_fits %>%
  mutate(preds_LR    = map2(model_LR, splits, ~predict_helper(.x,.y,model_type="LR")),
         preds_SVM   = map2(model_SVM, splits, ~predict_helper(.x,.y,model_type="SVM")),
         preds_KRR   = map2(model_KRR, splits, ~predict_helper(.x,.y,model_type="KRR")))
model_metrics <- model_preds %>%
  mutate(metric_LR   = map(preds_LR, get_metrics),
         inform_LR   = map_dbl(metric_LR,"Informedness"),
         metric_SVM  = map(preds_SVM, get_metrics),
         inform_SVM  = map_dbl(metric_SVM,"Informedness"),
         metric_KRR  = map(preds_KRR, get_metrics),
         inform_KRR  = map_dbl(metric_KRR,"Informedness"))

save(model_fits, "model_fits.RData")
save(model_preds, "model_fits.RData")
save(model_metrics, "model_fits.RData")

model_metrics %>%
  dplyr::select(inform_LR, inform_KRR) %>%
  summarise(mean_LR = mean(inform_LR),
            mean_SVM = mean(inform_SVM),
            mean_KRR = mean(inform_KRR),
            sd_LR   = sd(inform_LR),
            sd_SVM  = sd(inform_SVM),
            sd_KRR  = sd(inform_KRR))
