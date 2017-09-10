library("tidyr")
library("dplyr")
library("corrplot")
library("latex2exp")
library("data.table")
library("Matrix")
library("ggplot2")
library("proxy")
library("ggplot2")
library("hrbrthemes")

reps <- 1
method_object <- pr_DB$get_entry("Euclidean")
sigma <- 2.5
lambda <- 0.25

N_back_bags = 25 # need to figure out better way to measure this
N_sites     = 50
background_site_balance = 1
sample_fraction = 0.15
train_test_split = 0.75
confusion_matrix_cutoff = 0.5
data_location = "~/Dropbox/R/r91_all_upland_section_6_regression_data_SITENO.csv"

### Load Data
dat <- fread(data_location)
dat <- data.frame(dat)
dat1 <- dplyr::select(dat,presence, SITENO,
                      ed_h6,  std_32c, tpi_250c, elev_2_strm, cd_h4) #,
                      # cd_h7, slpvr_16c, ed_h2, cd_conf, ed_h2, elev_2_conf,
                      # elev_2_strm, e_hyd_min, ed_drnh, elev_2_drainh,
                      # tri_16c, cd_h5)
### Center and Standardize data
variables <- setdiff(colnames(dat1), c("presence", "SITENO"))
dat1   <- data.frame(apply(dat1[, variables],2,scale))
dat1   <- cbind(dat1, dat[,c("presence","SITENO")])

rep_results <- NULL
for(z in seq_len(reps)){
  ####### Data recycle for each rep (this should all be in a function one day...)
  formatted_data <- format_site_data(dat1, N_sites, train_test_split, background_site_balance)
  train_data              <- formatted_data[["train_data"]]
  test_data               <- formatted_data[["test_data"]]
  train_presence          <- formatted_data[["train_presence"]]
  test_presence           <- formatted_data[["test_presence"]] # 1/0 from shuffled list
  tbl_train_data          <- formatted_data[["tbl_train_data"]]
  tbl_train_presence      <- formatted_data[["tbl_train_presence"]]
  tbl_test_data           <- formatted_data[["tbl_test_data"]]
  tbl_test_presence       <- formatted_data[["tbl_test_presence"]]
  
  var_j_results <- NULL
  last_i <- 0
  for(j in 2:length(variables)){
    # set new formula here 
    j_vars <- variables[1:j]
    train_data_j <- lapply(train_data, function(x) x[(names(x) %in% j_vars)])
    test_data_j <- lapply(test_data, function(x) x[(names(x) %in% j_vars)])
    print(paste0("j = ",j, "; ", j_vars))
    animate_results <- NULL # moved to here.  test
    for(k in seq_along(sigma)){
      for(i in seq_along(lambda)){
        K <- build_K(train_data_j, train_data_j, sigma[k], dist_method = method_object)
        train_log_pred <- KRR_logit_optim(K, train_presence, lambda[i], 1000, 0.001)
        alphas_pred   <- train_log_pred[["alphas"]]
        #### Predict
        test_log_pred <- KRR_logit_predict(test_data_j, train_data_j, alphas_pred, 
                                           sigma[k], dist_method = method_object)
        ### Performance data frames
        train_log_pred_plot <- data.frame(pred = train_log_pred[["pred"]],
                                          obs = train_presence)
        train_log_pred_plot$pred_cat <- ifelse(train_log_pred_plot$pred >= confusion_matrix_cutoff,1,0)
        predicted_log <- data.frame(pred = test_log_pred,
                                    obs = test_presence,
                                    pred_cat = ifelse(test_log_pred >= confusion_matrix_cutoff,1,0),
                                    frame = last_i + i,
                                    sigma = sigma[k],
                                    lambda = lambda [i],
                                    rep =   z,
                                    model = "LMERR",
                                    vars = j)
        ### Performance Metrics
        TP <- sum(predicted_log$pred_cat == 1 & predicted_log$obs == 1, na.rm = TRUE)
        FP <- sum(predicted_log$pred_cat == 1 & predicted_log$obs == 0, na.rm = TRUE)
        TN <- sum(predicted_log$pred_cat == 0 & predicted_log$obs == 0, na.rm = TRUE)
        FN <- sum(predicted_log$pred_cat == 0 & predicted_log$obs == 1, na.rm = TRUE)
        print(paste0("row = ", last_i + i, " of ", length(lambda) * length(sigma),
                     "; sigma = ", round(sigma[k],3), 
                     "; lambda = ", round(lambda[i],4) , 
                     "; Reach = ", round(metrics(TP,TN,FP,FN)$Reach,3),
                     "; Rep = ", z, " of ", reps))
        animate_results <- rbind(animate_results, predicted_log)
      }
      last_i <- last_i + i
    }
    ### This is in the J loop, need to move it out one loop and put something here to capture j
    # rep_results <- rbind(rep_results, animate_results)
    var_j_results <- rbind(var_j_results, animate_results)
  }
  rep_results <- rbind(rep_results, var_j_results)
}

# write.csv(rep_results, "stepwiseVars_5rep_by_6params_50sites_analytical_LMERR_results.csv")

######## Analyze results
### this is a differnt output than make_xstats(), should try to make a common functions
xstats <- group_by(rep_results, vars, sigma, lambda, rep) %>%
  summarise(TP = sum(pred_cat == 1 & obs == 1, na.rm = TRUE),
            FP = sum(pred_cat == 1 & obs == 0, na.rm = TRUE),
            TN = sum(pred_cat == 0 & obs == 0, na.rm = TRUE),
            FN = sum(pred_cat == 0 & obs == 1, na.rm = TRUE),
            auc = pROC::auc(obs,pred, type = "linear")) %>%
  group_by(vars, rep, sigma, lambda) %>%
  dplyr::mutate(Reach = get_metric(TP,TN,FP,FN,"Reach"),
                KG = get_metric(TP,TN,FP,FN,"KG"),
                Sensitivity = get_metric(TP,TN,FP,FN,"Sensitivity"),
                `1-Specificity` = 1-get_metric(TP,TN,FP,FN,"Specificity"),
                avg_metric = (KG + Reach)/2,
                Kappa = 1-get_metric(TP,TN,FP,FN,"Kappa")) %>%
  data.frame()

# write.csv(xstats, "1000rep_by_25params_100sites_analytical_LMERR_aggregated_results.csv")

agg_stats <- group_by(xstats, sigma, lambda, vars) %>%
  summarise(auc_mean = mean(auc),
            # auc_sd = sd(auc),
            # n = n()
            # auc_high = auc_mean + 1.645*(auc_sd/n),
            # auc_low = auc_mean - 1.645*(auc_sd/n),
            Reach = mean(Reach),
            KG = mean(KG),
            Sensitivity = mean(Sensitivity),
            X1.Specificity = mean(X1.Specificity),
            avg_metric = mean(avg_metric)) %>%
  data.frame()

plot_agg_stats <-  gather(agg_stats, metric, value, -sigma, -lambda, -vars)
plot_agg_stats$sigma_label <- paste0("sigma ==",  plot_agg_stats$sigma)
plot_agg_stats$lambda_label <- paste0("lambda ==",  plot_agg_stats$lambda)

ggplot(plot_agg_stats, aes(x = factor(vars), y = value, group = metric)) +
  geom_line(aes(color = metric), size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 0.75, linetype = "dashed") +
  labs(x = "Number of Covariates", y="",
       title="Forward Stepwise Variable Selection Over Hyperparameter Grid Search",
       subtitle="Mean Embeddings Logistic Kernel Ridge Regression Model",
       caption="Mean of metric for 25% test sample over 5 repeated split-samples") + 
  theme_ipsum2(plot_title_size = 14,
               panel.spacing.x.unit = 1,
               panel.spacing.y.unit = 1) +
  facet_grid(lambda_label ~ sigma_label, labeller = "label_parsed") +
  scale_color_manual(values = ipsum_palette2) +
  # scale_x_continuous(breaks = lambda, labels = round(lambda,3)) +
  scale_y_continuous(limits=c(0,1),expand = c(0,0,0.1)) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.title=element_blank()
  ) -> gg
# gg_check(gg)
plot(gg)
# ggsave("1000rep_by_25params_100sites_analytical_LMERR_aggregated_results.png",
       # width = 9, height = 5)
