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

logreg2_helper <- function(dat_split, frac = 1){
  lr_data <- dat_split$tbl_train_data
  lr_presence <- lr_data %>%
    filter(presence == 1) %>%
    group_by(SITENO) %>%
    summarise_all(mean)
  lr_absence <- lr_data %>%
    filter(presence == 0) %>%
    sample_n(size = nrow(lr_presence)*frac)
  lr_data <- rbind(lr_presence, lr_absence) %>%
    dplyr::select(-SITENO)
  lr <- glm(presence ~ ., data = lr_data, family = "binomial")
}

SVM_helper <- function(dat_split, cost = 0.001){
  svm_data <- dat_split$tbl_train_data %>%
    dplyr::select(.,-SITENO)
  svm <- svm(presence ~ .,
             data = svm_data, cost = cost,
             family = "binomial")
}
predict_helper <- function(model,splits,model_type,FPR=NULL,confusion_matrix_cutoff=0.5){
  if(model_type %in% c("LR","SVM")){
    newdata <- splits$tbl_test_data
    preds <- predict(model, newdata = newdata, type = "response")
    predicted <- data.frame(pred = preds,
                            obs  = splits$tbl_test_presence)
  }
  if(model_type == "KRR"){
    method_object <- proxy::pr_DB$get_entry("Euclidean")
    alphas_pred   <- model[["alphas"]]
    newdata    <- splits$test_data
    train_data <- splits$train_data
    preds <- KRR_logit_predict(newdata, train_data, alphas_pred, sigma, dist_method = method_object)
    predicted <- data.frame(pred = preds,
                            obs  = splits$test_presence)
  }
  #### CUT OFF CLASSIFICATION IS ITS OWN STEP
  # if(is.numeric(FPR)){
  #   confusion_matrix_cutoff <- get_cm_cutoff(predicted,FPR=FPR)
  # }
  predicted$pred_cat <- ifelse(predicted$pred >= confusion_matrix_cutoff,1,0)
  return(predicted)
}
get_AUC <- function(preds){
  roc_obj <- roc(preds$obs, preds$pred)
  AUC <- auc(roc_obj)
}
get_cm_cutoff <- function(preds,FPR=0.66){ #### <- SHOULD be FPR OF JUST BACKGROUND!!!!
  ############ DOES NOT SEEM TO BE FIXING PREDICTIONS AT FPR - CONTINUE TO WORK!!!!!!!!!!!!
  # FPR cutoff is not implimented at the moment, only AUC and 0.5 threshold measures
  preds_sort <- preds[order(preds$pred, decreasing = T),"pred"]
  cm_cutoff <- preds_sort[floor(length(preds_sort)*FPR)]
}

library("corrplot")
library("data.table")
library("tidyverse")
library("e1071")         # for SVM comparison
library("DistRegLMERR")
library("future")
library("pROC")

# method_object <- proxy::pr_DB$get_entry("Euclidean")
# sigma <- 1
# lambda <- 0.11

### Data parameters
N_back_bags = 50 # need to figure out better way to measure this
N_sites     = 50
background_site_balance = 1
sample_fraction = 0.50
train_test_split = 0.75
confusion_matrix_cutoff = 0.5
#
physioshed_ids <- c(1,2,6,8,12)

for(z in seq_along(physioshed_ids)){
    physioshed_z <- physioshed_ids[z]
    physioshed   <- paste0("r91_all_upland_section_",physioshed_z, "_regression_data_SITENO.csv")
    # data_location = "/home/rstudio/r91_all_upland_section_12_regression_data_SITENO.csv"
    data_location = file.path("C:/Users/matthew.d.harris/Dropbox/R/PASS_regression",physioshed)

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
    searches <- 100
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
             model_LR2   = map(splits,     ~future::future(logreg2_helper(.x, frac = 1)))) %>%
      mutate(model_LR    = map(model_LR,   ~future::value(.x)),
             model_LR2   = map(model_LR2,  ~future::value(.x)))
    model_preds <- model_fits %>%
      mutate(preds_LR    = map2(model_LR,  splits, ~predict_helper(.x,.y,model_type="LR")),
             preds_LR2   = map2(model_LR2, splits, ~predict_helper(.x,.y,model_type="LR"))) %>%
      dplyr::select(-splits)
    model_metrics <- model_preds %>%
      mutate(AUC_LR = map_dbl(preds_LR, get_AUC),
             AUC_LR2 = map_dbl(preds_LR2, get_AUC),
             metric_LR   = map(preds_LR, get_metrics),
             inform_LR   = map_dbl(metric_LR,"Informedness"),
             metric_LR2  = map(preds_LR2, get_metrics),
             inform_LR2  = map_dbl(metric_LR2,"Informedness"))

    model_preds_strip <- model_preds %>%
      dplyr::select(-starts_with("model"))
    model_metrics_strip <- model_metrics %>%
      dplyr::select(-starts_with("model"), -starts_with("preds"))

    metrics_AUC_long <- model_metrics_strip %>%
      dplyr::select(-starts_with("inform"),-starts_with("metric"),-params) %>%
      gather(model, value, -id) %>%
      separate(model, into = c("metrics","model"), sep = "_")

    ggplot(metrics_AUC_long, aes(x=model, y=value, color = model)) +
      geom_boxplot(width=0.2)+
      geom_jitter(width=0.1) +
      theme_bw()

    metrics_threshold_long <- model_metrics_strip %>%
      select(-starts_with("inform"),-starts_with("AUC"),-params)  %>%
      mutate(metrics = map(metric_LR, names)) %>%
      mutate_if(is.list, simplify_all) %>%
      unnest()  %>%
      gather(model,value,-id,-metrics) %>%
      separate(model,into=c("a","model"),sep="_") %>%
      select(-a)

    save_folder <- paste0(searches,"_LR_cells_mean_compare_r91U",physioshed_z)
    if(!dir.exists(save_folder)){
      dir.create(save_folder)
    }
    # save(model_fits,    file = file.path(save_folder,"model_fits.RData"))
    save(model_preds_strip,   file = file.path(save_folder,"model_preds_srtip.RData"))
    save(model_metrics_strip, file = file.path(save_folder,"model_metrics_strip.RData"))
    write.csv(metrics_threshold_long, file = file.path(save_folder,"metrics_threshold_long.csv"))
    write.csv(metrics_AUC_long, file = file.path(save_folder,"metrics_AUC_long.csv"))

    # model_metrics %>%
    #   dplyr::select(inform_LR, inform_LR2) %>%
    #   summarise(median_LR = median(inform_LR),
    #             median_LR2 = median(inform_LR2),
    #             sd_LR   = sd(inform_LR),
    #             sd_LR2  = sd(inform_LR2),
    #             high_LR   = max(inform_LR),
    #             high_LR2  = max(inform_LR2),
    #             low_LR   = min(inform_LR),
    #             low_LR2  = min(inform_LR2)) %>%
    #   gather(metric, value)
    #
    # metrics_threshold_long %>%
    #   rbind(metrics_AUC_long) %>%
    #   filter(metrics %in% c("Informedness","Sensitivity","PPG","FPR","AUC")) %>%
    #   ggplot(.,aes(x=model,y=value,group=model, color = model)) +
    #   geom_boxplot(width = 0.2) +
    #   geom_jitter(width=0.2)+
    #   facet_wrap(~metrics, ncol=1, scales = "free") +
    #   theme_bw()
}
