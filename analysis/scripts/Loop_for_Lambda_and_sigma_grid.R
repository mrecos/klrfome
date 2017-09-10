sigma <- seq(0.01,3,length.out = 20)
lambda <- seq(0.002,0.25,length.out = 10)
confusion_matrix_cutoff = 0.5
results <- matrix(ncol=5, nrow = length(lambda)*length(sigma))
colnames(results) <- c("sigma", "lambda", "Reach", "phat Sites", "phat Back")
last_i <- 0
method_object <- pr_DB$get_entry("Euclidean")
for(k in seq_along(sigma)){
  for(i in seq_along(lambda)){
    K <- build_K(train_data, train_data, sigma[k], dist_method = method_object)
    train_log_pred <- KRR_logit_optim(K, train_presence, lambda[i], 1000, 0.001)
    alphas_pred   <- train_log_pred[["alphas"]]
    #### Predict
    test_log_pred <- KRR_logit_predict(test_data, train_data, alphas_pred, sigma[k], dist_method = method_object)

    ### Performance data frames
    train_log_pred_plot <- data.frame(pred = train_log_pred[["pred"]],
                                      obs = train_presence)
    train_log_pred_plot$pred_cat <- ifelse(train_log_pred_plot$pred >= confusion_matrix_cutoff,1,0)

    predicted_log <- data.frame(pred = test_log_pred,
                                obs = test_presence,
                                pred_cat = ifelse(test_log_pred >= confusion_matrix_cutoff,1,0))
    pred_values <- predicted_log
    ### Performance Metrics
    phat <- group_by(pred_values, obs) %>%
      summarize(mean_pred = mean(pred))
    # confusionMatrix(predicted_log$pred_cat, predicted_log$obs, positive = "1")
    TP <- sum(pred_values$pred_cat == 1 & pred_values$obs == 1, na.rm = TRUE)
    FP <- sum(pred_values$pred_cat == 1 & pred_values$obs == 0, na.rm = TRUE)
    TN <- sum(pred_values$pred_cat == 0 & pred_values$obs == 0, na.rm = TRUE)
    FN <- sum(pred_values$pred_cat == 0 & pred_values$obs == 1, na.rm = TRUE)
    results[last_i + i,1] <- round(sigma[k],3)
    results[last_i + i,2] <- round(lambda[i],3)
    results[last_i + i,3] <- metrics(TP,TN,FP,FN)$Reach
    results[last_i + i,4] <- as.numeric(phat[2,2])
    results[last_i + i,5] <- as.numeric(phat[1,2])
    print(paste0("row = ", last_i + i,
                 "; sigma = ", round(sigma[k],3), 
                 "; lambda = ", round(lambda[i],4) , 
                 "; Reach = ", round(metrics(TP,TN,FP,FN)$Reach,3)))
  }
  last_i <- last_i + i
}

print(results)
results[which(results[,"Reach"] == max(results[,"Reach"])),]
ggplot(data.frame(results), aes(x = lambda, y = Reach, 
                                group = sigma, color = as.factor(sigma))) +
  geom_line() + 
  theme_bw() +
  facet_grid(~sigma) +
  theme(
  axis.text.x = element_text(angle = 90, hjust = 1)
  )

library("viridis")
ggplot(data.frame(results), aes(x = lambda, y = sigma, fill = Reach)) +
  geom_raster() + 
  scale_fill_viridis() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

# sigma = 0.575, lambda = 0.002
# sigma = 0.01,  lambda = 0.076
# sigma = 0.507, lambda = 0.051
