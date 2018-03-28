library("dplyr")
library("corrplot")
library("ggplot2")
library("klrfome")

#Parameters
set.seed(31337)
sigma = 0.5
lambda = 0.1
dist_metric = "euclidean"

### Simulate Training Data
sim_data <- get_sim_data(site_samples = 800, N_site_bags = 75)
formatted_data <- format_site_data(sim_data, N_sites=10, train_test_split=0.8,
                                   sample_fraction = 0.9, background_site_balance=1)
train_data <- formatted_data[["train_data"]]
train_presence <- formatted_data[["train_presence"]]
test_data <- formatted_data[["test_data"]]
test_presence <- formatted_data[["test_presence"]]

##### Logistic Mean Embedding KRR Model
#### Build Kernel Matrix
K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric)
#### Train
train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 2)
#### Predict
test_log_pred <- KLR_predict(test_data, train_data, dist_metric = dist_metric,
                             train_log_pred[["alphas"]], sigma)

### Plot K Matrix
K_corrplot(K,train_data,clusters=4)

### Plot Prediction
predicted_log <- data.frame(pred = test_log_pred, obs = test_presence)
ggplot(predicted_log, aes(x = as.factor(obs), y = pred, color = as.factor(obs))) +
  geom_jitter(width = 0.1) +
  theme_bw() +
  ylim(c(0,1)) +
  labs(y = "Predicted Probability", x = "Site Presence",
       title = "Kernel Logistic Regression",
       subtitle = "test set predictions; simulated data") +
  theme(
    legend.position = "none"
  )
