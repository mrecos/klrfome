############### Logistic Mean Embedding KRR (and non-Logistic)
### Full process example with simulated data

library("dplyr")
# library("corrplot")
library("latex2exp")
# library("data.table")
library("ggplot2")
library("DistRegLMERR")

#Parameters
set.seed(sample(1:99999,1))
sigma = 0.1
lambda = 0.015

### Simulate Training Data
sim_data <- get_sim_data(sites_var1_mean = 50, sites_var1_sd = 10,
                         sites_var2_mean = 3,   sites_var2_sd   = 2,
                         backg_var1_mean = 60, backg_var1_sd   = 12,
                         backg_var2_mean = 4,   backg_var2_sd   = 2.25,
             site_samples    = 100,
             N_site_bags     = 10,
             background_site_balance = 1,
             test_train_split = 0.50)

train_data <- sim_data[["train_data"]] %>%
  lapply(., function(x) dplyr::select(x,-SITENO)) ## remove SITENO = workaround
train_presence <- sim_data[["train_presence"]]
test_data <- sim_data[["test_data"]] %>%
  lapply(., function(x) dplyr::select(x,-SITENO)) ## remove SITENO = workaround
test_presence <- sim_data[["test_presence"]]

##### Logistic Mean Embedding KRR Model
#### Build Kernel Matrix
method_object <- proxy::pr_DB$get_entry("Euclidean")
K <- build_K(train_data, sigma = sigma, dist_method = method_object)
#### Train
train_log_pred <- KRR_logit_optim(K, train_presence, lambda, 100, 0.01)
alphas_pred   <- train_log_pred[["alphas"]]
#### Predict
test_log_pred <- KRR_logit_predict(test_data, train_data, alphas_pred, sigma, dist_method = method_object)

##### Plots
## response of training data
ggplot(data.frame(x=K[1,],y=train_log_pred$pred),aes(x=x,y=y)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

## response of Test data
ggplot(data.frame(x=K[1,],y=test_log_pred),aes(x=x,y=y)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

### Plot K Matrix
K_corrplot(K,train_data,clusters=4)

### Plot Fit
train_log_pred_plot <- data.frame(pred = train_log_pred[["pred"]],
                                  obs = train_presence)
ggplot(train_log_pred_plot, aes(x = obs, y = pred)) +
  geom_jitter(width=0.05) +
  theme_bw() +
  ylim(c(0,1))
### Plot Prediction
predicted_log <- data.frame(pred = test_log_pred,
                            obs = test_presence)
# group_by(predicted_log, obs) %>%
#   summarize(mean_pred = mean(pred))
subtitle <- TeX('$1/(1 + exp(-\\frac{1}{n}\\sum{n}^{i=1}(K_g(x,x`)) + \\lambda ||f||^2_H))\n;\n K_g(x,x`) = K(\\mu_x,\\mu_{x`})$')
ggplot(predicted_log, aes(x = as.factor(obs), y = pred,
                          color = as.factor(obs))) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values=c("blue","orange")) +
  theme_bw() +
  ylim(c(0,1)) +
  labs(y = "Predicted Probability",
       x = "Site Presence",
       title = "Mean Embedding Logistic Kernel Ridge Regression",
       subtitle = subtitle) +
  theme(
    legend.position = "none",
    text=element_text(family="Trebuchet MS", size = 10),
    plot.title = element_text(size = 12, family = "TrebuchetMS-Bold")
  )
