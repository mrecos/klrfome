############### Logistic Mean Embedding KRR (and non-Logistic)
### Full process example with simulated data

library("dplyr")
library("corrplot")
library("ggplot2")
library("klrfome")

#Parameters
set.seed(1337)
sigma = 0.5
lambda = 0.1

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
method_object <- proxy::pr_DB$get_entry("Euclidean")
K <- build_K(train_data, sigma = sigma, dist_method = method_object)
#### Train
train_log_pred <- KLR(K, train_presence, lambda, 100, 0.01)
alphas_pred   <- train_log_pred[["alphas"]]
#### Predict
test_log_pred <- KLR_predict(test_data, train_data, alphas_pred, sigma, dist_method = method_object)

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
