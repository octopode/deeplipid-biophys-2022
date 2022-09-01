# Plotting for MD of low-complexity synthetic lipid systems

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)

source(here("02-scripts", "20220502_lipidomics_helpers.R"))

simdata_summ = read_tsv(
  here("01-rawdata", "20220523_simpars.tsv"),
  col_types = cols(comp = col_character(), ph = col_character(), .default = col_double())
) %>% 
  # come up with the delta APL parameter (1 bar-normalized)
  group_by(comp) %>% 
  arrange(press_sim) %>% 
  mutate(
    dapl_mean = apl_mean / apl_mean[[1]],
    dapl_serr = apl_serr / apl_mean[[1]]
  ) %>% 
  pivot_longer(matches("mean|serr"), names_to = "var", values_to = "val") %>% 
  separate(var, into = c("par", "stat"), sep = '_') %>% 
  group_by(across(-c("stat", "val"))) %>% 
  pivot_wider(names_from = stat, values_from = val) %>% 
  separate(comp, into = c("chains", "class"), sep = ' ', remove = FALSE)

# some really simple facets

barwid = 50
simdata_summ %>% 
  filter(
    (press_sim <= 1000) |
      (chains == "18:0/20:4")
  ) %>% 
  #mutate(comp = str_glue(comp, ' ', temp_sim, "°C")) %>% 
  ggplot(
    aes(
      x = press_sim,
      y = mean,
      ymin = mean - serr,
      ymax = mean + serr,
      color = class,
      group = class
    )
  ) +
  facet_grid(rows = vars(par), cols = vars(chains), scales = "free_y", switch = 'y') +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  geom_errorbar(position = position_dodge(width = 2*barwid)) +
  geom_point(position = position_dodge(width = 2*barwid)) +
  scale_color_manual(values = chroma_cl) +
  theme_pubr() +
  labs(
    title = "MD results summaries",
    x = "simulation pressure",
    y = "parameter value"
  ) +
  guides(color = "none")
ggsave(here("04-suppfigs", "20220523_sims_1.pdf"), width=8, height=6)

# turn it on end!
# no dioleoyl, for main text
simdata_summ %>% 
  filter(
    #(ph == "la") &
      (chains != "18:1/18:1")
  ) %>% 
  arrange(par,desc(chains)) %>% 
  rowwise() %>% 
  mutate(
    chainpar = paste(chains, par),
    chainpar = chainpar %>% factor(., levels = unique(.))
  ) %>%
  ggplot(
    aes(
      y = press_sim,
      x = mean,
      xmin = mean - serr,
      xmax = mean + serr,
      color = class,
      group = class
    )
  ) +
  facet_grid(cols = vars(chainpar), scales = "free_x") +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  geom_errorbarh(height = barwid) +
  geom_point() +
  scale_color_manual(values = chroma_cl) +
  scale_y_reverse() +
  theme_pubr() +
  labs(
    title = "MD results summaries",
    y = "simulation pressure",
    x = "parameter value"
  ) +
  guides(color = "none")
ggsave(here("03-mainfigs", "synsim", "20220524_sims_1.pdf"), width=10, height=3)
