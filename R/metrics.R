#' make_quads
#' 
#' `make_quads()` Produces True Positive, False Positive, True Negative, and False Negative quantities from two binary input vectors of predicted and observed classes.
#' 
#' This function takes two binary vectors coded as 1 == presence/positive and 0 == absent/negative. The `pred` vector is the predicted binary class and `obs` vector is the true observed class. The output quantites can be used for additional metric evaluation.
#'
#' @param pred - [vector] Predicted probabilities.
#' @param obs - [vector] Observed presence/absence as 1/0
#'
#' @return [vector] TP, FP, TN, FN names vector
#' @export
#' 
#' @examples 
#' \dontrun{
#' cm <- make_quads(ifelse(test_log_pred >= 0.5, 1, 0), test_presence)
#' metrics(TP = cm[1], TN = cm[3], FP = cm[2], FN = cm[4])$Informedness
#'}
#'
make_quads <- function(pred,obs){
  TP = sum(pred == 1 & obs == 1, na.rm = TRUE)
  FP = sum(pred == 1 & obs == 0, na.rm = TRUE)
  TN = sum(pred == 0 & obs == 0, na.rm = TRUE)
  FN = sum(pred == 0 & obs == 1, na.rm = TRUE)
  return(c("TP"=TP,"FP"=FP,"TN"=TN,"FN"=FN))
}

#' CM_quads
#'
#' `CM_quads()` produces the True Positive, False Positive, True Negative, and False Negative quantities of the Confusion Matrix at one or more threshholds for the binary classification of continuous prediction values.
#' 
#' This function takes two arguments: 1) a data.frame of two columns; `pred` is a column of continuous predicted values and `presence` is a binary observed class value (encoded as `1 == present` or `0 == absent`); and 2) a scalar or numeric vector of threshold values `threshold` at which to classify the continuous values of `pred` into binary 1/0 classes. The results is a data table where each row holds the TP, FP, TN, FN counts for each of the thresholds given in `threshold`.
#'  
#'
#' @param dat - [data.frame] Table with two columns, "pred" and "presence". "pred" is the predicted probability. "presence" is the observed presence/absence as 1/0.
#' @param threshold - [scalar or vector] a scalar or vector of one or more thresholds at which to evaluate the Confusion Matrix quadrants.
#'
#' @return [data.frame] A data.frame of Confusion Matrix quatrants at one or more threshold values.
#'
#' @importFrom dplyr group_by summarise
#' @export
#' 
CM_quads <- function(dat,threshold = 0.5){
  threshold_class <- NULL
  for(i in seq_along(threshold)){
    threshold_i <- data.frame(pred = dat$pred,
                              obs = dat$presence,
                              pred_cat = ifelse(model_pred$pred >= threshold[i],1,0),
                              Threshold = threshold[i])
    threshold_class <- rbind(threshold_class, threshold_i)
  }
  kstats <- threshold_class %>%
    dplyr::group_by(Threshold) %>%
    dplyr::summarise(TP = sum(pred_cat == 1 & obs == 1, na.rm = TRUE),
                     FP = sum(pred_cat == 1 & obs == 0, na.rm = TRUE),
                     TN = sum(pred_cat == 0 & obs == 0, na.rm = TRUE),
                     FN = sum(pred_cat == 0 & obs == 1, na.rm = TRUE))
  return(kstats)
}


#' cohens_kappa
#' 
#' `cohens_kappa()` is a metric scoring function that returns the Cohen's Kappa statistic
#' 
#' This function takes integer counts for True Positive `TP`, True Negative `TN`, False Positive `FP`, and False Negative `FN` values and returns the Cohen's Kappa. This statistic measure the agreement between observations and predictions. See: https://en.wikipedia.org/wiki/Cohen%27s_kappa
#'
#' @param TP - [scalar] True Positives
#' @param TN - [scalar] True Negatives
#' @param FP - [scalar] False Positives
#' @param FN - [scalar] False Negatives
#'
#' @return [scalar] Cohen's Kappa statistic
#' 
#' #'\dontrun{
#' cm <- make_quads(ifelse(test_log_pred >= 0.5, 1, 0), test_presence)
#' cohens_kappa(TP = cm[1], TN = cm[3], FP = cm[2], FN = cm[4])
#'}
#'
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


#' metrics
#' 
#' `metrics()` returns a wide range of binary class evaluation metrics based on the inputs of True Positive, True Negative, Fale Positive, and Fale Negative quantities
#' 
#' This function is a one-stop-shop to compute 50+ metric results based on the input of TP, TN, FP, and FN
#'
#' @param TP - [scalar] True Positives
#' @param TN - [scalar] True Negatives
#' @param FP - [scalar] False Positives
#' @param FN - [scalar] False Negatives
#' @param a - [scalar] alpha level
#'
#' @return [list] - list of all metrics
#' @importFrom boot logit
#' @export
#'
metrics <- function(TP,TN,FP,FN){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("boot needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # metrics derived from TP,TN,FP, and FN
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
    Pm          = (TP+FP)/(TP+TN+FP+FN), # probability of site-likely (Verhagen 2007)
    Pm_prime    = (FN+TN)/(TP+TN+FP+FN), # probability of site-unlikely (Verhagen 2007)
    Precision   = Precision, # PPV
    Recall      = Recall, # Sensitivity, TPR
    F_Measure   = (2*Precision*Recall)/(Precision+Recall),
    Geo_Mean    = sqrt(TP*TN),
    MAE           = mean(abs(pred-obs)),
    RMSE          = sqrt(mean((pred-obs)^2)),
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
    KG2         = 1-(((TP+FP)/(TP+TN+FP+FN))/Sensitivity), # 1-(Pm/sens)
    DOR         = (Sensitivity/(1-Specificity)/((1-Sensitivity)/Specificity)), # LRP/LRN
    log_DOR     = log10((Sensitivity/(1-Specificity)/((1-Sensitivity)/Specificity))),
    # D & S - https://en.wikipedia.org/wiki/Diagnostic_odds_ratio
    D           = boot::logit(Sensitivity) - boot::logit(1-Specificity),
    S           = boot::logit(Sensitivity) + boot::logit(1-Specificity),
    Kappa       = cohens_kappa(TP,TN,FP,FN),
    # http://aircconline.com/ijdkp/V5N2/5215ijdkp01.pdf
    Opp_Precision = ((TP+TN)/(TP+TN+FP+FN))-(abs(Specificity-Sensitivity)/(Specificity+Sensitivity)),
    # https://en.wikipedia.org/wiki/Precision_and_recall
    MCC         = suppressWarnings((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))),
    Informedness  = Sensitivity+Specificity-1, #TSS, Younden's J
    Markedness    = (TP/(TP+FP))+(TN/(FN+TN))-1,
    # http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2006.01214.x/full
    ## From Verhagen(2007; 121)
    AFK           = suppressWarnings(sqrt(Sensitivity*((Sensitivity-(1-Specificity))/((TN+FP)/(TP+TN+FP+FN))))),
    Indicative      = Sensitivity/(1-Specificity),
    Indicative2     = Sensitivity/((TP+FP)/(TP+TN+FP+FN)), # TPR/Pm
    Indicative_norm = (Sensitivity/(1-Specificity))/((TN+FP)/(TP+TN+FP+FN)),
    Indicative_norm2 = (Sensitivity/((TP+FP)/(TP+TN+FP+FN)))/((TN+FP)/(TP+TN+FP+FN)),
    Brier         = mean((obs-pred)^2), # MSE for binary class problems
    X1            = (TP/(TP+FP))/((TP+FN)/(TP+TN+FP+FN)), # PPV/Prev or Prec/Prev
    X2            = (FN/(FN+TN))/((TP+FN)/(TP+TN+FP+FN)), # FOR/Prev
    X3            = (FP/(TP+FP))/((TN+FP)/(TP+TN+FP+FN)), # FDR/other
    X4            = (TN/(FN+TN))/((TN+FP)/(TP+TN+FP+FN)), # NPV/other
    # Oehlert & Shea (2007)
    PPG           = (TP/(TP+FP))/((TP+FN)/(TP+TN+FP+FN)), # PPV/Prev or Prec/Prev or X1
    NPG           = (FN/(FN+TN))/((TP+FN)/(TP+TN+FP+FN)), # FOR/Prev or X2
    # mean of KG+Reach = Balance
    Balance       = ((1-((1-Specificity)/Sensitivity))+
                       (1-((1-Sensitivity)/Specificity)))/2,
    Balance2       = ((1-(((TP+FP)/(TP+TN+FP+FN))/Sensitivity))+
                        (1-((1-Sensitivity)/((FN+TN)/(TP+TN+FP+FN)))))/2
  )
  return(metrics)
}
