get_metrics <- function(dat, pred_cat, obs){
  TP <- sum(dat$pred_cat == 1 & dat$obs == 1, na.rm = TRUE)
  FP <- sum(dat$pred_cat == 1 & dat$obs == 0, na.rm = TRUE)
  TN <- sum(dat$pred_cat == 0 & dat$obs == 0, na.rm = TRUE)
  FN <- sum(dat$pred_cat == 0 & dat$obs == 1, na.rm = TRUE)
  inf   <- suppressWarnings(metrics(TP,TN,FP,FN)$Informedness)
  sens  <- suppressWarnings(metrics(TP,TN,FP,FN)$Sensitivity)
  spec2 <- suppressWarnings(1-metrics(TP,TN,FP,FN)$Specificity)
  list(inf, sens, spec2)
}

############### Logistic Mean Embedding KRR (and non-Logistic)
### Works with functions in KRR_Logit_fit_predict_functions.R
### Full process example with REAL data

### IDEA: in predict, limit examples/coefficients to those in neighborhood of test bag; local similaritiy

library("dplyr")
library("corrplot")
library("latex2exp")
library("data.table")
library("ggplot2")
library("e1071")         # for SVM comparison
library("klrfome")

#Parameters
# set.seed(3849)
sigma = 2
lambda = 0.81
dist_metric = "euclidean"

### Data parameters
N_back_bags = 50 # need to figure out better way to measure this
N_sites     = 50
background_site_balance = 1
sample_fraction = 0.5
train_test_split = 0.80
confusion_matrix_cutoff = 0.5
# data_location = "data/r91_all_upland_section_6_regression_data_SITENO.csv"
# data_location = "/Users/mattharris/Dropbox/R/PASS_regression/r91_all_upland_section_12_regression_data_SITENO.csv"
data_location = "C:/Users/matthew.d.harris/Dropbox/R/PASS_regression/r91_all_upland_section_6_regression_data_SITENO.csv"


### Load Data
dat <- fread(data_location)
dat <- data.frame(dat)
dat1 <- dplyr::select(dat,presence, SITENO,
                      ed_h6, std_32c, cd_h7, slpvr_16c, ed_h2, cd_conf, elev_2_conf)
# elev_2_strm, e_hyd_min, ed_drnh, elev_2_drainh,
# tri_16c, cd_h5)
## Reduce number of sites to N_sites
formatted_data <- format_site_data(dat1, N_sites, train_test_split, background_site_balance,
                                   sample_fraction = sample_fraction)
train_data              <- formatted_data[["train_data"]]
test_data               <- formatted_data[["test_data"]]
train_presence          <- formatted_data[["train_presence"]]
test_presence           <- formatted_data[["test_presence"]] # 1/0 from shuffled list
tbl_train_data          <- formatted_data[["tbl_train_data"]]
tbl_train_presence      <- formatted_data[["tbl_train_presence"]]
tbl_test_data           <- formatted_data[["tbl_test_data"]]
tbl_test_presence       <- formatted_data[["tbl_test_presence"]]
# all.equal(test_presence, ifelse(grepl(pattern="background", names(test_data)), 0, 1)) # for testing & reasurnace

## Logistic Mean Embedding KRR Model
#### Build Kernel Matrix
K <- build_K(train_data, train_data, sigma, dist_metric = dist_metric)
# diag(K) <- 1
#### Train
train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 2)
#### Predict
test_log_pred <- KLR_predict(test_data, train_data, train_log_pred[["alphas"]], 
                             sigma, dist_metric = dist_metric)

### Performance data frames
train_log_pred_plot <- data.frame(pred = train_log_pred[["pred"]],
                                  obs = train_presence)
predicted_log <- data.frame(pred = test_log_pred,
                            obs = test_presence,
                            pred_cat = ifelse(test_log_pred >= confusion_matrix_cutoff,1,0))

### Performance Metrics ##
preds <- data.frame(preds = test_log_pred, obs = test_presence)
cm_cutoff <- get_YJ_max_threshold(preds)
cm <- make_quads(ifelse(test_log_pred >= cm_cutoff, 1, 0), preds$obs)
klrfome::metrics(cm["TP"],cm["TN"],cm["FP"],cm["FN"])$Informedness
# Get AUC #
roc_obj <- pROC::roc(test_presence, test_log_pred)
as.numeric(pROC::auc(roc_obj))
group_by(predicted_log, obs) %>%
  summarize(mean_pred = mean(pred))

### compare to logit
lr <- glm(presence ~ ., data = dplyr::select(tbl_train_data, -SITENO), family = "binomial")
lr_pred <- predict(lr, newdata = tbl_test_data, type = "response")
predicted_lr <- data.frame(pred = lr_pred,
                            obs = tbl_test_presence,
                            pred_cat = ifelse(lr_pred >= confusion_matrix_cutoff,1,0))
### LR Performance Metrics
group_by(predicted_lr, obs) %>%
  summarize(mean_pred = mean(pred))
# x <- confusionMatrix(predicted_lr$pred_cat, predicted_lr$obs, positive = "1")
lr_metrics <- get_metrics(predicted_lr)

### compare to SVM
svm_mod <- svm(presence ~ .,
               data = dplyr::select(tbl_train_data, -SITENO), cost = 0.001,
               family = "binomial")
svm_pred <- predict(svm_mod, newdata = tbl_test_data, type = "response")
predicted_svm <- data.frame(pred = svm_pred,
                            obs = tbl_test_presence,
                            pred_cat = ifelse(svm_pred >= confusion_matrix_cutoff,1,0))
### LR Performance Metrics
group_by(predicted_svm, obs) %>%
  summarize(mean_pred = mean(pred))
# x <- confusionMatrix(predicted_lr$pred_cat, predicted_lr$obs, positive = "1")
svm_metrics <- get_metrics(predicted_svm)

data.frame(metric = c("Informedness", "Sensitivity", "1-Specificity"),
           KRR    = c(krr_metrics[[1]], krr_metrics[[2]], krr_metrics[[3]]),
           LR     = c(lr_metrics[[1]], lr_metrics[[2]], lr_metrics[[3]]),
           SVM    = c(svm_metrics[[1]], svm_metrics[[2]], svm_metrics[[3]])) %>%
  mutate_if(is.numeric, round, 2)


##### Plots
### Plot K Matrix
# colnames(K) <- names(train_data)
# rownames(K) <- names(train_data)
# col3 <- colorRampPalette(c("red", "white", "blue"))
# corrplot::corrplot(K,tl.cex = 0.5, tl.col = "black",
#                    order="hclust", col=col3(10), cl.lim=c(0,1),
#                    addrect = 6)

ggplot(predicted_log, aes(x = as.factor(obs), y = pred, color = as.factor(obs))) +
  geom_jitter(width = 0.1, alpha = 0.15) +
  theme_bw() +
  ylim(c(0,1)) +
  labs(y = "Predicted Probability", x = "Site Presence",
       title = "Kernel Logistic Regression",
       subtitle = "test set predictions; simulated data") +
  theme(
    legend.position = "none"
  )



