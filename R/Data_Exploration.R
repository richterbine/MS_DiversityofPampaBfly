## data exploration ##

library(ggplot2)
library(tidyverse)
library(reshape)

bfly.pampa <- read.csv(here::here("data/raw/Species_data.csv"), header = T, 
                       sep = ";")
head(bfly.pampa)

length(unique(bfly.pampa$Species))

# percentual of species in each subfamily
round((sum(subset(bfly.pampa, Subfamily == "Satyrinae")$Abundance)*100)/sum(bfly.pampa$Abundance), 1)
round((sum(subset(bfly.pampa, Subfamily == "Charaxinae")$Abundance)*100)/sum(bfly.pampa$Abundance), 1)
round((sum(subset(bfly.pampa, Subfamily == "Biblidinae")$Abundance)*100)/sum(bfly.pampa$Abundance), 1)
round((sum(subset(bfly.pampa, Subfamily == "Nymphalinae")$Abundance)*100)/sum(bfly.pampa$Abundance), 1)

# Evaluation by site/local, note that here, we evaluate the local as a whole, not 
# separating into sampling occasions

df <- bfly.pampa[, c("Local", "SO", "Environment", "Tribe", "Subfamily", "Species")]
head(df)
df$comm <- paste(df$Local, df$SO, sep = " ")
  
grouped.df <- df %>%
  group_by( Environment, Tribe ) %>%
  count( comm ) %>%
  ungroup()

grouped.df$Tribe <- factor(grouped.df$Tribe, levels = c("Nymphalini", "Coeini",
                                                        "Ageroniini", "Biblidini",
                                                        "Callicorini", "Catonephelini",
                                                        "Anaeini", "Preponini", "Brassolini",
                                                        "Morphini", "Satyrini"))

p.tribo <- grouped.df[-which(grouped.df$comm == "EEA 1"),] %>%
  ggplot() +
  geom_col( aes(x = reorder(comm, -n), y = n, fill = Tribe), position = position_dodge2(width = .1)) +
  facet_wrap( ~ Environment ) + scale_y_sqrt() + 
  labs(x = NULL, y = "Abundance", fill = "Tribe", tag = "b)") +
  scale_fill_viridis_d() + theme(legend.position = "bottom")
p.tribo

# by subfamily
grouped.df <- df %>%
  group_by( Environment, Subfamily ) %>%
  count( comm ) %>%
  ungroup()

grouped.df$Subfamily <- factor(grouped.df$Subfamily, levels = c("Nymphalinae", "Biblidinae",
                                                                "Charaxinae", "Satyrinae"))

p.subf <- grouped.df[-which(grouped.df$comm == "EEA 1"),] %>%
  ggplot() +
  geom_col( aes(x = reorder(comm, -n), y = n, fill = Subfamily), position = position_dodge2(width = .1) ) +
  facet_wrap( ~ Environment ) + scale_y_sqrt() +
  labs(x = NULL, y = "Abundance", fill = "Subfamily", tag = "a)") +
  scale_fill_viridis_d() + theme(legend.position = "bottom")
p.subf


grouped.df <- df %>%
  group_by( Environment, Species ) %>%
  count( comm ) %>%
  ungroup()

#grouped.df <- arrange(grouped.df, desc(n), .by_group = T)

p.spp <- grouped.df[-which(grouped.df$comm == "EEA 1"),] %>%
  ggplot() +
  geom_col( aes(x = reorder(comm, -n), y = n, fill = Species), position = position_dodge2(width = .1)) +
  facet_wrap( ~ Environment ) + scale_y_sqrt() + 
  labs(x = "Sampling Locals", y = "Abundance", fill = "Species", tag = "c)") +
  scale_fill_viridis_d() + theme(legend.position = "none")
p.spp

p.fim <- cowplot::plot_grid(p.subf, p.tribo, p.spp, ncol = 1)
cowplot::save_plot(filename = here::here("output/figures/Final_figures/FigS1_summary.png"), plot = p.fim, 
                   base_height = 10, base_width = 8)
