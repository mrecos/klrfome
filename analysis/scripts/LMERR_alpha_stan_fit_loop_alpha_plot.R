# options(mc.cores = parallel::detectCores())
library("rstan")
library("bayesplot")
model_dat2 <- list(K = K, y = train_presence, P = nrow(K), N = length(train_presence))
file_loc <- "/Users/mattharris/Dropbox/R/ITE_MERR_Dist_Regression/"
fit_normal2 <- stan(file = paste0(file_loc, "LMERR_lambda_stan.stan"), 
                   data = model_dat2, control = list(adapt_delta = 0.90))

yhat2_summary <- summary(fit_normal2, pars = c("yhat2"), probs = c(0.25,0.5,0.95))$summary
yhat <- as.numeric(yhat2_summary[,"mean"])
mod_results <- data.frame(stan = yhat, analytical = train_log_pred$pred, obs = train_presence)
plot_mod_results <- tidyr::gather(mod_results, model, value, -obs)
ggplot(plot_mod_results, aes(x = obs, y = value, color = model)) +
  geom_jitter(width = 0.15) +
  ylim(c(0,1)) +
  theme_bw()

## look at lambda posterior
lambda_draw <- as.matrix(fit_normal2, par = "lambda")
mcmc_areas(lambda_draw)
summary(fit_normal2, pars = "lambda")$summary

### Post model fit
## Loop predictons across many parts of the alpha_hat posterior distribution to see effect
probs <- seq(0.1,1,0.1)
probs_label <- paste0("P",probs)
params <- summary(fit_normal2, pars = c("alpha_hat"), probs = probs)$summary %>%
  data.frame() %>%
  mutate(SITENO = names(train_data),
         presence = train_presence)
colnames(params) <- c("mean","se_mean","sd",probs_label,"n_eff","Rhat","SITENO","presence")

alpha_results <- NULL
for(i in seq_along(probs_label)){
  print(paste0("iteration: ", i,", alpha: ", probs_label[i]))
  alphas_pred   <- as.matrix(dplyr::select_(params, probs_label[i]) )
  #### Predict
  test_log_pred <- KRR_logit_predict(test_data, train_data, alphas_pred, 
                                     sigma, dist_method = method_object)
  # train_log_pred_plot <- data.frame(pred = train_log_pred[["pred"]],
  #                                   obs = train_presence)
  predicted_log <- data.frame(pred = test_log_pred,
                              obs = test_presence,
                              pred_cat = ifelse(test_log_pred >= confusion_matrix_cutoff,1,0),
                              alpha = probs_label[i])
  alpha_results <- rbind(alpha_results, predicted_log)
}

xstats <- group_by(alpha_results, alpha) %>%
  summarise(TP = sum(pred_cat == 1 & obs == 1, na.rm = TRUE),
            FP = sum(pred_cat == 1 & obs == 0, na.rm = TRUE),
            TN = sum(pred_cat == 0 & obs == 0, na.rm = TRUE),
            FN = sum(pred_cat == 0 & obs == 1, na.rm = TRUE),
            auc = pROC::auc(obs,pred, type = "linear")) %>%
  group_by(alpha) %>%
  dplyr::mutate(Reach = get_metric(TP,TN,FP,FN,"Reach"),
                KG = get_metric(TP,TN,FP,FN,"KG"),
                Sensitvity = get_metric(TP,TN,FP,FN,"Sensitivity"),
                `1-Specificity` = 1-get_metric(TP,TN,FP,FN,"Specificity"),
                mean_met = (Reach+KG)/2) %>%
  data.frame()


ggplot(alpha_results, aes(x = as.factor(obs), y = pred,
                          color = as.factor(obs))) +
  geom_jitter(width = 0.2, alpha = 0.1) +
  scale_color_manual(values=c("blue","orange")) +
  theme_bw() +
  ylim(c(0,1)) +
  labs(y = "Predicted Probability",
       x = "Site Presence",
       title = "Mean Embedding Logistic Kernel Ridge Regression",
       subtitle = "Test Data") +
  facet_wrap(~alpha, ncol = 3) +
  theme(
    legend.position = "none",
    text=element_text(family="Trebuchet MS", size = 10),
    plot.title = element_text(size = 12, family = "TrebuchetMS-Bold")
  )
