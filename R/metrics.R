###### FROM my METRICS FUNCTIONS
#### NEEDS REFRESH WITH UDATED KG CALC!!
cohens_kappa <- function(TP,TN,FP,FN){
  A <- TP
  B <- FP
  C <- FN
  D <- TN
  Po <- (A+D)/(A+B+C+D)
  Pe_a <- ((A+B)*(A+C))/(A+B+C+D)
  Pe_b <- ((C+D)*(B+D))/(A+B+C+D)
  Pe <- (Pe_a+Pe_b)/(A+B+C+D)
  k <- (Po-Pe)/(1-Pe)
  return(k)
}
get_metric <- function(TP,TN,FP,FN,metric_type){
  m <- metrics(TP,TN,FP,FN)
  m <- m[[metric_type]]
  return(m)
}
##### NEED TO UPDATE WITH NEW KG CALCs
metrics <- function(TP,TN,FP,FN){
  # metrics derived from TP,TN,FP, and FN
  ### "Summary" scores
  # Sens, Spec, Precision, Recall calculated here and reused below
  Sensitivity <- TP/(TP+FN) # TPR, Recall
  Specificity <- TN/(FP+TN) # TNR
  Precision   <- TP/(TP+FP) # PPV
  Recall      <- TP/(TP+FN) # Sensitivity, TPR
  pred        <- c(rep(1,TP),rep(1,FP),rep(0,FN),rep(0,TN))
  obs         <- c(rep(1,TP),rep(0,FP),rep(1,FN),rep(0,TN))
  metrics <- list(
    Sensitivity = Sensitivity, # TPR, Recall
    Specificity = Specificity, # TNR
    Prevalence  = (TP+FN)/(TP+TN+FP+FN),
    Accuracy    = (TP+TN)/(TP+TN+FP+FN),
    Err_Rate    = (FP+FN)/(TP+TN+FP+FN),
    other       = (TN+FP)/(TP+TN+FP+FN), # % all area with no sites
    Precision   = Precision, # PPV
    Recall      = Recall, # Sensitivity, TPR
    F_Measure   = (2*Precision*Recall)/(Precision+Recall),
    Geo_Mean    = sqrt(TP*TN),
    FPR         = 1 - Specificity, #Fall-Out
    FNR         = 1 - Sensitivity,
    TPR         = Sensitivity, # Sensitivity, Recall
    TNR         = Specificity, # Specificity
    FOR         = FN/(FN+TN),
    FDR         = FP/(TP+FP),
    Power       = 1-(1-Sensitivity),
    LRP         = Sensitivity/(1-Specificity), # TPR/FPR
    log_LRP     = log10(Sensitivity/(1-Specificity)),
    LRN         = (1-Sensitivity)/Specificity, # FNR/TNR
    PPV         = TP/(TP+FP), # Precision
    NPV         = TN/(FN+TN),
    KG          = 1-((1-Specificity)/Sensitivity), # 1-(FPR/TPR)
    DOR         = (Sensitivity/(1-Specificity)/((1-Sensitivity)/Specificity)), # LRP/LRN
    log_DOR     = log10((Sensitivity/(1-Specificity)/((1-Sensitivity)/Specificity))),
    # D & S - https://en.wikipedia.org/wiki/Diagnostic_odds_ratio
    D           = boot::logit(Sensitivity) - boot::logit(1-Specificity),
    S           = boot::logit(Sensitivity) + boot::logit(1-Specificity),
    Kappa       = cohens_kappa(TP,TN,FP,FN),
    # Kappa agrees with http://terpconnect.umd.edu/~dchoy/thesis/Kappa/#
    # http://aircconline.com/ijdkp/V5N2/5215ijdkp01.pdf
    Opp_Precision = ((TP+TN)/(TP+TN+FP+FN))-(abs(Specificity-Sensitivity)/(Specificity+Sensitivity)),
    # https://en.wikipedia.org/wiki/Precision_and_recall
    # MCC         = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),
    Informedness  = Sensitivity+Specificity-1, #TSS, Younden's J
    Markedness    = (TP/(TP+FP))+(TN/(FN+TN))-1,
    # http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2006.01214.x/full
    # TSS         = Sensitivity+Specificity-1 #Informedness
    # TEST = 1-(FNR/Spec) or 1-(FNR/TNR)
    # one minus % sites incorrect divided by % site not-likely background
    # also, 1 minus (% misclassifications / % non-site area)
    # Reach = TEST
    Reach         = 1-((1-Sensitivity)/Specificity),
    # Reach2        = 1-((1-(TP/(TP+FN)))/(TN/(FP+TN)))
    ## From Verhagen(2007; 121)
    AFK           = suppressWarnings(sqrt(Sensitivity*((Sensitivity-(1-Specificity))/((TN+FP)/(TP+TN+FP+FN))))),
    Indicative      = Sensitivity-(1-Specificity),
    Indicative_norm = (Sensitivity-(1-Specificity))/((TN+FP)/(TP+TN+FP+FN)),
    Brier         = mean((obs-pred)^2), # MSE for binary class problems
    # adding more
    MAE           = mean(abs(pred-obs)),
    RMSE          = sqrt(mean((pred-obs)^2))

  )
  return(metrics)
}
make_xstats <- function(results){
  library("pROC")
  xstats <- group_by(results, rep, model) %>%
    summarise(TP = sum(pred_cat == 1 & obs == 1, na.rm = TRUE),
              FP = sum(pred_cat == 1 & obs == 0, na.rm = TRUE),
              TN = sum(pred_cat == 0 & obs == 0, na.rm = TRUE),
              FN = sum(pred_cat == 0 & obs == 1, na.rm = TRUE),
              auc = pROC::auc(obs,pred, type = "linear")) %>%
    group_by(rep) %>%
    dplyr::mutate(Reach = get_metric(TP,TN,FP,FN,"Reach"),
                  KG = get_metric(TP,TN,FP,FN,"KG"),
                  Sensitivity = get_metric(TP,TN,FP,FN,"Sensitivity"),
                  `1-Specificity` = 1-get_metric(TP,TN,FP,FN,"Specificity"),
                  avg_metric = (KG + Reach)/2) %>%
    data.frame()
  return(xstats)
}
