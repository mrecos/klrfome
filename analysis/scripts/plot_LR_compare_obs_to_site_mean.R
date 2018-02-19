


library("tidyverse")
library("viridis")
# devtools::install_github("thomasp85/patchwork")
library("patchwork")

physioshed_ids <- c(1,2,6,8,12)
searches = 100

### comparing arbitrary 0.5 thresholds threhsold data
all_dat_LR <- NULL
for(z in seq_along(physioshed_ids)){
  physioshed_z <- physioshed_ids[z]
  # data_loc <- "C:\\Users\\matthew.d.harris\\Dropbox\\R\\ITE_MERR_Dist_Regression\\AWS_RESULTS\\Threshold_05_runs\\100_reps_LR_cell_to_mean_compare_all_physiosheds"
  data_loc <- "C:\\Users\\matthew.d.harris\\Dropbox\\R\\ITE_MERR_Dist_Regression\\AWS_RESULTS\\Threshold_max_Youden_runs\\LR_compare2"
  save_folder <- paste0(searches,"_LR_cells_mean_compare_r91U",physioshed_z)
  metrics_z <- read.csv(file.path(data_loc,save_folder,"metrics_threshold_long.csv"))
  auc_z <- read.csv(file.path(data_loc,save_folder,"metrics_AUC_long.csv"))
  dat_z <- rbind(metrics_z, auc_z) %>%
    mutate(physioshed = physioshed_z,
           threshold  = "max_Youden")
  all_dat_LR <- rbind(all_dat_LR, dat_z)
}
glimpse(all_dat_LR)



plot_dat <- filter(all_dat_LR, metrics %in% c("Informedness", "Sensitivity",
                                              "FPR","PPG","AUC"))
plot_dat %>%
  group_by(model, metrics) %>%
  summarise(median = median(value))
ggplot(plot_dat, aes(x=model, y = value, color = model)) +
  geom_boxplot() +
  scale_color_viridis_d(end = 0.70) +
  geom_jitter(width=0.05, alpha = 0.25)+
  facet_wrap(~metrics, ncol=1, scales = "free") +
  theme_bw()
#
# load(file.path(data_loc,save_folder,"model_metrics_strip.Rdata"))
# load(file.path(data_loc,save_folder,"model_preds_srtip.Rdata"))
