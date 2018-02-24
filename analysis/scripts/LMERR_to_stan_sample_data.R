library('rstan')
library("DistRegLMERR")
library("dplyr")

set.seed(717)
iterations = 10000
chains = 3
warmup = 2500
sigma = 1
lambda = 0.01

###### Sample datat for testing linear model to make sure stan and file loc ar working
# data(mtcars)
# dat <-data.frame(x1 = runif(100,-2,2), x2 = runif(100,-2,2))
# X <- model.matrix(~x1+x2,dat)
# betas <- runif(ncol(X),-1,1)
# sigma <- 0.25
# y <- rnorm(100,X%*%betas, sigma)
# model_dat <- list(x1 = dat$x1, x2 = dat$x2, y = y, N = length(y))
# file_loc <- "/Users/mattharris/Dropbox/R/ITE_MERR_Dist_Regression/"
# fit_normal <- stan(file = paste0(file_loc, "LMERR_stan.stan"), data = model_dat)
######## End of example data, hold here for now... ####

### Simulate Training Data
sim_data <- get_sim_data(sites_var1_mean = 50, sites_var1_sd = 10,
                         sites_var2_mean = 3,   sites_var2_sd   = 2,
                         backg_var1_mean = 60, backg_var1_sd   = 12,
                         backg_var2_mean = 4,   backg_var2_sd   = 4,
                         site_samples    = 100,
                         N_site_bags     = 10,
                         background_site_balance = 1,
                         test_train_split = 0.50)

train_data <- sim_data[["train_data"]]
train_data <- lapply(train_data, function(x) x[!(names(x) %in% c("SITENO"))])
train_presence <- sim_data[["train_presence"]]
test_data <- bind_rows(sim_data[["test_data"]]) %>%
              mutate(id = paste0(SITENO, "_", seq_len(n()))) %>%
              split(f = .$id) %>%
              lapply(function(x) x[!(names(x) %in% c("id", "SITENO"))])
test_presence <- ifelse(grepl("Site", names(test_data)), 1, 0)

##### Logistic Mean Embedding KRR Model
method_object <- proxy::pr_DB$get_entry("Euclidean")
K <- build_K(train_data, train_data, sigma, dist_method = method_object)
train_log_pred <- KRR_logit_optim(K, train_presence, lambda, 100, 0.01)
alphas_pred   <- train_log_pred[["alphas"]]
test_log_pred <- KRR_logit_predict(test_data, train_data, alphas_pred,
                                   sigma, dist_method = method_object)
pROC::auc(test_presence, test_log_pred)
p_test <- data.frame(obs = test_presence, pred = test_log_pred)
ggplot(p_test, aes(x = as.factor(obs), y = pred)) +
  geom_jitter()


model_dat <- list(K = K, lambda = lambda, y = train_presence, P = nrow(K),
                  N = length(train_presence))
file_loc <- "C:/R_local/DistRegLMERR/analysis/Stan_models"
fit_normal <- stan(file = file.path(file_loc, "LMERR_stan2.stan"), data = model_dat)

##### plot compare stan to analytical train predictions
yhat1_summary <- summary(fit_normal, pars = c("yhat1"), probs = c(0.05,0.5,0.95))$summary
yhat <- as.numeric(yhat1_summary[,"mean"])
mod_results <- data.frame(stan = yhat, analytical = train_log_pred$pred, obs = train_presence)
plot_mod_results <- tidyr::gather(mod_results, model, value, -obs)
ggplot(plot_mod_results, aes(x = obs, y = value, color = model)) +
  geom_point() +
  ylim(c(0,1)) +
  theme_bw()

