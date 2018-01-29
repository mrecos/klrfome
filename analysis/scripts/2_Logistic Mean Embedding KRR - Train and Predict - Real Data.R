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
library("DistRegLMERR")

#Parameters
# set.seed(3849)
sigma = 1
lambda = 0.11

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
dat <- fread(data_location)
dat <- data.frame(dat)
dat1 <- dplyr::select(dat,presence, SITENO,
                      ed_h6, std_32c, ed_h7, slpvr_32c, ed_h2, cd_conf, ed_h2, elev_2_conf)
# elev_2_strm, e_hyd_min, ed_drnh, elev_2_drainh,
# tri_16c, cd_h5)
### Center and Standardize data
variables <- setdiff(colnames(dat1), c("presence", "SITENO"))
dat1   <- data.frame(apply(dat1[, variables],2,scale))
dat1   <- cbind(dat1, dat[,c("presence","SITENO")])
## Reduce number of sites to N_sites
formatted_data <- format_site_data(dat1, N_sites, train_test_split, background_site_balance)
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
method_object <- proxy::pr_DB$get_entry("Euclidean")
K <- build_K(train_data, train_data, sigma, dist_method = method_object)
diag(K) <- 1
#### Train
train_log_pred <- KRR_logit_optim(K, train_presence, lambda, 100, 0.001, verbose = 2)
# train_log_pred <- KRR_logit(K, train_presence, lambda)
alphas_pred   <- train_log_pred[["alphas"]]
#### Predict
test_log_pred <- KRR_logit_predict(test_data, train_data, alphas_pred, sigma, dist_method = method_object)

### Performance data frames
train_log_pred_plot <- data.frame(pred = train_log_pred[["pred"]],
                                  obs = train_presence)

predicted_log <- data.frame(pred = test_log_pred,
                            obs = test_presence,
                            pred_cat = ifelse(test_log_pred >= confusion_matrix_cutoff,1,0))

### Performance Metrics
group_by(predicted_log, obs) %>%
  summarize(mean_pred = mean(pred))
# x <- confusionMatrix(predicted_log$pred_cat, predicted_log$obs, positive = "1")
krr_metrics <- get_metrics(predicted_log)

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
group_by(predicted_vsm, obs) %>%
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
# ### Plot Fit
# ggplot(train_log_pred_plot, aes(x = obs, y = pred)) +
#   geom_jitter(width=0.05) +
#   theme_bw() +
#   ylim(c(0,1))
### Plot Prediction
subtitle <- TeX('$1/(1 + exp(-\\frac{1}{n}\\sum{n}^{i=1}(K_g(x,x`)) + \\lambda ||f||^2_H))\n;\n K_g(x,x`) = K(\\mu_x,\\mu_{x`})$')
ggplot(predicted_log, aes(y = as.factor(obs), x = pred,
                          color = as.factor(obs))) +
  geom_point(alpha = 0.5, position = position_jitter(w = 0, h = 0.3)) +
  scale_color_manual(values=c("blue","orange")) +
  theme_bw() +
  xlim(c(0,1)) +
  # labs(y = "Predicted Probability",
  #      x = "Site Presence",
  #      title = "Mean Embedding Logistic Kernel Ridge Regression",
  #      subtitle = subtitle) +
  theme(
    legend.position = "none",
    text=element_text(family="Trebuchet MS", size = 10),
    plot.title = element_text(size = 12, family = "TrebuchetMS-Bold")
  )

inf_plot(predicted_log$pred, predicted_log$obs, threshold = 0.5)


# d(x,y) = sqrt(K(x,x)  + K(y,y) - 2K(x,y))





