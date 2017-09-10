threshold <- seq(0,1,0.1)

model_name <- c("LMERR", "GLM", "XGB", "GLM_mean", "XGB_mean")
model_pred <- list(LMERR = test_log_pred_raw,
                   GLM = glm_pred_raw,
                   XGB = xgb_pred_raw,
                   GLM_mean = glm_pred_mean_raw,
                   XGB_mean = xgb_pred_mean_raw)
model_test_presence <- list(LMERR = test_presence_full,
                            GLM = tbl_test_presence_full,
                            XGB = tbl_test_presence_full,
                            GLM_mean = tbl_test_presence_full,
                            XGB_mean = tbl_test_presence_full)
all_kstats <- NULL
for(j in seq_along(model_name)){
  threshold_class <- NULL
  print(model_name[j])
  pb <- txtProgressBar(min = 0, max = length(threshold), style = 3)
  for(i in seq_along(threshold)){
    threshold_i <- data.frame(pred = model_pred[[j]],
                                obs = model_test_presence[[j]],
                                pred_cat = ifelse(model_pred[[j]] >= threshold[i],1,0),
                                rep = z,
                                model = model_name[j],
                              cutoff = threshold[i])
    threshold_class <- rbind(threshold_class, threshold_i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  kstats <- group_by(threshold_class, cutoff) %>%
    summarise(TP = sum(pred_cat == 1 & obs == 1, na.rm = TRUE),
              FP = sum(pred_cat == 1 & obs == 0, na.rm = TRUE),
              TN = sum(pred_cat == 0 & obs == 0, na.rm = TRUE),
              FN = sum(pred_cat == 0 & obs == 1, na.rm = TRUE),
              auc = pROC::auc(obs,pred, type = "linear"),
              model = unique(model)) %>%
    group_by(cutoff) %>%
    dplyr::mutate(Reach = get_metric(TP,TN,FP,FN,"Reach"),
                  KG = get_metric(TP,TN,FP,FN,"KG"),
                  Sensitivity = get_metric(TP,TN,FP,FN,"Sensitivity"),
                  `1-Specificity` = 1-get_metric(TP,TN,FP,FN,"Specificity"),
                  avg_metric = (KG + Reach)/2,
                  Kappa = get_metric(TP,TN,FP,FN,"Kappa"),
                  FPR = get_metric(TP,TN,FP,FN,"FPR"),
                  TPR = get_metric(TP,TN,FP,FN,"TPR")) %>%
    data.frame()
  
  # kstats_plot <- dplyr::select(kstats, cutoff, Kappa, FPR) %>%
  #   gather(metric, value, -cutoff)
  
  # ggplot(data = kstats, aes(x = FPR, y = TPR)) + # y = TPR or Kappa
  #   geom_line(color = "red") +
  #   geom_line(data = data.frame(x = c(0,1), y = c(0,1)),
  #             aes(x = x, y = y), color = "black", linetype = 3) +
  #   theme_bw() +
  #   scale_x_continuous(breaks=seq(0,1,0.05)) +
  #   scale_y_continuous(breaks=seq(0,1,0.05))
  
  all_kstats <- rbind(all_kstats, kstats)
}    


ggplot(data = all_kstats, aes(x = FPR, y = TPR, 
                          group = model)) + # y = TPR or Kappa
  geom_line(aes(color = model)) +
  geom_line(data = data.frame(x = c(0,1), y = c(0,1), model = "Random"),
            aes(x = x, y = y), color = "black", linetype = 3) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,1,0.05)) +
  scale_y_continuous(breaks=seq(0,1,0.05))


## I made this to check my kappa calcualtion in metrics()
## made it as function to call from mutate()
## My calcualtions are identical and faster to run
## save for some posterity
# library("psych")
# ck <- ck(threshold_class$pred_cat, threshold_class$obs)
# ck <- function(pred_cat, obs){
#   library("psych")
#   ck_data <- data.frame(obs = obs, pred = pred_cat)
#   ckappa <- cohen.kappa(ck_data)
#   kappa <- ckappa$kappa
#   return(kappa)
# }


