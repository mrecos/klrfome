
library(tidyverse)
physioshed_ids <- c(1,2,6,8,12)
searches = 100
all_dat <- NULL
for(z in seq_along(physioshed_ids)){
  physioshed_z <- physioshed_ids[z]
  data_loc <- "C:\\Users\\matthew.d.harris\\Dropbox\\R\\ITE_MERR_Dist_Regression\\AWS_RESULTS"
  save_folder <- paste0(searches,"_KRR_LR_SVM_compare_r91U",physioshed_z)
  metrics_z <- read.csv(file.path(data_loc,save_folder,"metrics_threshold_long.csv"))
  auc_z <- read.csv(file.path(data_loc,save_folder,"metrics_AUC_long.csv"))
  dat_z <- rbind(metrics_z, auc_z) %>%
    mutate(physioshed = physioshed_z)
  all_dat <- rbind(all_dat, dat_z)
}
glimpse(all_dat)

plot_dat <- filter(all_dat, metrics %in% c("AUC","Informedness"))
ggplot(plot_dat, aes(x=model, y = value, color = model)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha = 0.2)+
  facet_wrap(~metrics, ncol=1, scales = "free") +
  theme_bw()
