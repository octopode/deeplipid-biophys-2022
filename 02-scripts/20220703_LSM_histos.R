library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)

file_lsm_gp = here("01-rawdata", "20220530_LSM_GP_uncorrd.tsv")
data_lsm_gp = read_tsv(file_lsm_gp)

data_lsm_norm = data_lsm_gp %>% 
  group_by(sp, temp_incu) %>% 
  # normalize
  mutate(count = count/sum(count))

data_lsm_norm %>% 
  filter(gp >= 0) %>% 
  arrange(temp_incu) %>% 
  mutate(temp_incu = temp_incu %>% factor(., levels = unique(.))) %>% 
  ggplot(
    aes(
      x = gp,
      y = count,
      color = temp_incu,
      group = temp_incu
    )
  ) +
  facet_grid(rows = vars(sp)) +
  geom_point(alpha = 0.4) +
  geom_smooth(span = 0.1) +
  theme_pubr() +
  scale_color_manual(
    values = c(
      "2"  = "#1A2A6C",
      "10" = "#FDBB2C",
      "20" = "#B2201F"
    )
  ) +
  labs(
    x = "C-laurdan GP",
    y = "Normalized pixel count",
    color = "Incubation temperature (°C)"
  )
ggsave(here("03-mainfigs", "depthdists", "20220703_clau_histos.pdf"), width=6, height=3)

