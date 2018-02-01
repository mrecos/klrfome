
library("tidyverse")
library("viridis")
# devtools::install_github("thomasp85/patchwork")
library("patchwork")

physioshed_ids <- c(1,2,6,8,12)
searches = 100

### comparing arbitrary 0.5 thresholds threhsold data
all_dat_05 <- NULL
for(z in seq_along(physioshed_ids)){
  physioshed_z <- physioshed_ids[z]
  data_loc <- "C:\\Users\\matthew.d.harris\\Dropbox\\R\\ITE_MERR_Dist_Regression\\AWS_RESULTS\\Threshold_05_runs"
  save_folder <- paste0(searches,"_KRR_LR_SVM_compare_r91U",physioshed_z)
  metrics_z <- read.csv(file.path(data_loc,save_folder,"metrics_threshold_long.csv"))
  auc_z <- read.csv(file.path(data_loc,save_folder,"metrics_AUC_long.csv"))
  dat_z <- rbind(metrics_z, auc_z) %>%
    mutate(physioshed = physioshed_z,
           threshold  = "05")
  all_dat_05 <- rbind(all_dat_05, dat_z)
}
glimpse(all_dat_05)

plot_dat <- filter(all_dat_05, metrics %in% c("AUC","Informedness"))
ggplot(plot_dat, aes(x=model, y = value, color = model)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha = 0.2)+
  facet_wrap(~metrics, ncol=1, scales = "free") +
  theme_bw()


### comparing to Y-stats threhsold data
all_dat_Y <- NULL
for(z in seq_along(physioshed_ids)){
  physioshed_z <- physioshed_ids[z]
  data_loc <- "C:\\Users\\matthew.d.harris\\Dropbox\\R\\ITE_MERR_Dist_Regression\\AWS_RESULTS\\Threshold_max_Youden_runs"
  save_folder <- paste0(searches,"_KRR_LR_SVM_compare_r91U",physioshed_z)
  metrics_z <- read.csv(file.path(data_loc,save_folder,"metrics_threshold_long.csv"))
  auc_z <- read.csv(file.path(data_loc,save_folder,"metrics_AUC_long.csv"))
  dat_z <- rbind(metrics_z, auc_z) %>%
    mutate(physioshed = physioshed_z,
           threshold  = "max_Youden")
  all_dat_Y <- rbind(all_dat_Y, dat_z)
}
glimpse(all_dat_Y)

plot_dat <- filter(all_dat_Y, metrics %in% c("AUC","Informedness"))
ggplot(plot_dat, aes(x=model, y = value, color = model)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha = 0.2)+
  facet_wrap(~metrics, ncol=1, scales = "free") +
  theme_bw()

all_dat_05_sum <- all_dat_05 %>%
  group_by(metrics, model, threshold) %>%
  summarise(median = median(value, na.rm=T))
all_dat_Y_sum <- all_dat_Y %>%
  group_by(metrics, model, threshold) %>%
  summarise(median = median(value, na.rm=T))
all_dat_sum <- cbind(all_dat_05_sum,all_dat_Y_sum)


all_dat_plot <- filter(rbind(all_dat_05, all_dat_Y), metrics %in% c("Informedness", "Sensitivity",
                                                                    "FPR","PPG"))
ggplot(all_dat_plot, aes(x=model, y=value, fill = threshold,
                         group=interaction(model,threshold))) +
  geom_boxplot() +
  facet_wrap(~metrics, scales = "free", nrow=2) +
  theme_bw()

AUC_dat_Y_plot <- filter(all_dat_Y, metrics %in% c("AUC"))
AUC_dat_Y_plot %>%
  group_by(model, physioshed) %>%
  summarise(median = median(value)) %>%
  spread(physioshed,median)
AUC_plot <- ggplot(AUC_dat_Y_plot, aes(x=as.factor(physioshed), y=value,
                           color = as.factor(model),
                           group=interaction(as.factor(physioshed),model))) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(), alpha = 0.2) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "gray50", alpha=0.7) +
  scale_y_continuous(limits=c(0,1),labels=seq(0,1,0.1),breaks=seq(0,1,0.1)) +
  labs(title = "Median AUC value for 100 splits by model type and physioshed",
       y = "Area Under the ROC Curve (AUC)",
       x = "Physioshed") +
  scale_color_viridis_d(name = "model type", option="B", end=0.8) +
  theme_bw()
plot(AUC_plot)

Y_dat_Y_plot <- filter(all_dat_Y, metrics %in% c("Informedness"))
Y_dat_Y_plot %>%
  group_by(model, physioshed) %>%
  summarise(median = median(value)) %>%
  spread(physioshed,median)
YJ_plot <- ggplot(Y_dat_Y_plot, aes(x=as.factor(physioshed), y=value,
                           color = as.factor(model),
                           group=interaction(as.factor(physioshed),model))) +
  geom_boxplot(fill = "gray95") +
  geom_jitter(position=position_jitterdodge(), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50", alpha=0.7) +
  scale_y_continuous(limits=c(-0.5,1),labels=seq(-0.5,1,0.1),breaks=seq(-0.5,1,0.1)) +
  labs(title = "Median Youden's J value for 100 splits by model type and physioshed",
       y = "Younden's J (Informedness)",
       x = "Physioshed") +
  # scale_color_discrete(name = "model type") +
  scale_color_viridis_d(name = "model type", option="B", end=0.8) +
  theme_bw()
plot(YJ_plot)

AUC_plot + YJ_plot + plot_layout(ncol = 1)




