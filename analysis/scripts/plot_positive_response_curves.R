
### could be two plots.... 1) where all train_sites are from 0 to 1 similar to "test" and 2) where all background is the same
### done while holding background or sites (whatever the counter class is) at 0 or 0.5 depending on assumptions
set_p <- function(tr_names,class,p,alt_p){
  if(class == "Site"){
    pattern = ".*Site.*|.*36.*"
    alt_pattern = ".*background.*"
  } else if(class == "background"){
    pattern = ".*background.*"
    alt_pattern = ".*Site.*|.*36.*"
  }
  ks <- gsub(pattern = pattern, replacement = p, tr_names)
  ks <- gsub(pattern = alt_pattern, replacement = alt_p, ks)
  return(as.numeric(ks))
}
get_ks <- function(tr_names, class, alt_p, int = 0.1){
  ks <- NULL
  # if(class == "Site"){
    p_seq = seq(0,1,int)
  # } else if(class == "background"){
  #   p_seq = rev(seq(0,1,int))
  # }
  for(i in seq_along(p_seq)){
    ki <- set_p(tr_names,class,p_seq[i],alt_p)
    ks <- rbind(ks,ki)
  }
  return(ks)
}

tr_names <- names(train_data)
alphas_pred <- params$alphas_pred
# p = 1
alt_p = 0.5
int = 0.1
ks_site <- get_ks(tr_names, "Site", alt_p, int)
ks_background <- get_ks(tr_names, "background", alt_p, int)
Site_pred <- 1 / (1 + exp(-as.vector(ks_site %*% alphas_pred)))
back_pred <- 1 / (1 + exp(-as.vector(ks_background %*% alphas_pred)))
# plot(back_pred, type="l")

plot_dat <- data.frame(pred = c(Site_pred, back_pred), 
                       similarity = rep(seq(0,1,int),2),
                       class = rep(c("Site","Background"),
                                   times=c(length(Site_pred),
                                   length(back_pred))))

ggplot(plot_dat, aes(x=similarity,y=pred,group=class, color=class)) +
  geom_line(size=2) +
  labs(x="Similarity to Sites", y = "Predicted Probability") +
  theme_bw()


### Loop over alt_p
tr_names <- names(train_data)
alphas_pred <- params$alphas_pred
p = 1
int = 0.01
alt_p_seq = seq(0,1,0.1)
Site_pred2 <- NULL
for(i in seq_along(alt_p_seq)){
  ks_site <- get_ks(tr_names, "Site", alt_p_seq[i], int)
  Site_pred <- 1 / (1 + exp(-as.vector(ks_site %*% alphas_pred)))
  pred_i <- data.frame(pred = Site_pred, alt_p = alt_p_seq[i])
  Site_pred2 <- rbind(Site_pred2, pred_i)
}

plot_dat2 <- data.frame(pred = Site_pred2$pred, 
                        alt_p = Site_pred2$alt_p,
                        similarity = rep(seq(0,1,int),length(alt_p_seq)),
                        class = rep(c("Site")))

ggplot(plot_dat2, aes(x=similarity,y=pred,group=alt_p, color = as.factor(alt_p))) +
  geom_line(size=2) +
  scale_colour_viridis_d(option = "A", name = "similarity of\nbackground") +
  labs(x="Similarity to Sites", y = "Predicted Probability",
       title = "Positive Response Curve Over Range of Site Similarities",
       subtitle = "Each curve is the response based on a static similarity for all background areas to all sites") +
  theme_bw()





