library(tidyverse)
library(arrow)
library(patchwork)

factor_list <- c("Ahr", "Irf4", "Zbtb46", "Maf", "Mafb")

filepath <- "results/all_cells.perturb_scores.parquet"
df <- read_parquet(filepath)

limit <- max(abs(df$score)) * c(-1, 1)
limit <- 2 * c(-1, 1)
p1 <- df %>%
  filter(factor %in% factor_list) %>%
  ggplot(aes(x, y)) +
  geom_point(aes(color = score)) +
  scale_color_distiller(type = "div", limit = limit) +
  facet_wrap(~factor, nrow = 1)

filepath <- "results/all_cells.mc_transitions.500.parquet"
df <- read_parquet(filepath)

p2 <- df %>%
  filter(factor %in% factor_list) %>%
  ggplot(aes(umap_1, umap_2)) +
  geom_hex(aes(fill = after_stat(density)), bins = 20) +
  facet_wrap(~factor, nrow = 1)

p1 / p2 +
  plot_layout(axes = "collect")
