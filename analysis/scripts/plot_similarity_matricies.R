library("viridis")
library("ggplot2")

samps <- sample(1:nrow(K),20)

K_high_small <- K_high[samps,samps]
cp_high <- K_corrplot(K_high_small,train_data[samps],clusters=4)
cp_high_tri <- as.data.frame(cp_high) %>%
  mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>%
  gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE) %>%
  mutate(sigma = "sigma == 3")

K_low_small <- K_low[samps,samps]
cp_low <- K_corrplot(K_low_small,train_data[samps],clusters=4)
cp_low_tri <- as.data.frame(cp_low) %>%
  mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>%
  gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE)  %>%
  mutate(sigma = "sigma == 1")

cp_tri <- rbind(cp_high_tri, cp_low_tri)

title <- TeX('Similarity Matrix $K(\\hat{\\mu}_{ij},\\hat{\\mu}_{ij}\\prime)$')
ggplot(data = cp_tri, aes(Var2, Var1, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Similarity", option = "A")+
  guides(fill = guide_colorbar(ticks = FALSE, nbins = 100)) +
  labs(x="",y="",
       title = title) +
  coord_equal() +
  theme_minimal() +
  facet_grid(~sigma, labeller = "label_parsed") +
  theme(
    text = element_text(family="Trebuchet MS"),
    axis.text.x = element_text(angle=90, vjust = 0.5),
    strip.background =element_rect(fill="gray95", color = "white"),
    strip.text = element_text(face="bold", size = 12),
    plot.title = element_text(size = 20, face = "bold", hjust=0.5)
  )
ggsave(file.path("/Users/mattharris/Dropbox/R/SAA_2018_poster/images","similarity_matricies.png"),
       width = 9, height = 5)
