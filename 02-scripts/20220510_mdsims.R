# Plotting for MD of low-complexity synthetic lipid systems

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)

source(here("02-scripts", "20220502_lipidomics_helpers.R"))

simdata_summ = read_tsv(here("01-rawdata", "20220510_simpars.tsv")) %>% 
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
ggsave(here("04-suppfigs", "20220510_sims_2.pdf"), width=8, height=6)

# how about testing the differences in slope?

# SCD analysis

scds = read_tsv(here("01-rawdata", "20220609_Scd.tsv")) %>% 
  pivot_longer(cols = c(sn1, sn2), names_to = "tail", values_to = "scd")# %>% 
  #drop_na(scd)

# FIG. 3 PANELS

# bind in SCD data

simdata_all = simdata_summ %>% 
  bind_rows(
    scds %>% 
      separate(comp, c("chains", "class"), sep=' ') %>% 
      dplyr::rename(
        par = tail,
        val = scd
      ) %>% 
      mutate(par = paste("scd", par, sep='_'))
  )

# raw SCD profiles (probably for SI)
plot_scd = simdata_all %>% 
  filter(str_detect(par, "scd")) %>% 
  filter(press %in% c(0, 1000)) %>% 
  arrange(press) %>% 
  mutate(press = factor(press, levels=unique(press))) %>% 
  ggplot(
    aes(
      x = cnum,
      y = val,
      color = class,
      linetype = press
    )
  ) +
  facet_grid(cols = vars(chains), rows = vars(par)) +
  geom_line() +
  scale_color_manual(values = chroma_ppe) +
  guides(
    color = "none",
    linetype = "none"
  ) +
  theme_pubr() +
  labs(
    x = "Carbon #",
    y = "Order parameter SCD"
  )

# let's see it
plot_scd
ggsave(here("04-suppfigs", "20220630_raw_scd_1-1000bar.pdf"), width=8, height=6)

# plot delta SCD between two specified pressures
presslims = c(0, 1000)
plot_dscd = simdata_all %>% 
  filter(str_detect(par, "scd")) %>% 
  filter(press %in% presslims) %>% 
  group_by(comp, class, chains, par, cnum) %>% 
  arrange(press) %>% 
  summarise(val = val[[2]] - val[[1]]) %>% 
  ungroup() %>% 
  mutate(par = paste('d', par, sep='')) %>% 
  ggplot(
    aes(
      x = cnum,
      y = val,
      color = class
    )
  ) +
  facet_grid(cols = vars(chains), rows = vars(par)) +
  geom_line() +
  scale_color_manual(values = chroma_ppe) +
  guides(
    color = "none",
    linetype = "none"
  ) +
  theme_pubr() +
  labs(
    x = "Carbon #",
    y = str_glue("Delta SCD over ", presslims[[1]], "-", presslims[[2]], " bar")
  )

plot_dscd
ggsave(here("03-mainfigs", "synsim", "20220630_delta_scd_1-1000bar.pdf"), width=8, height=6)

# Now, how about "plateau" SCD vs. pressure?
# Defining "plateau" as the carbon with maximum SCD @ 1 bar
simdata_all %>% 
  filter(str_detect(par, "scd")) %>% 
  group_by(comp, class, chains, par) %>% 
  filter(press == 0) %>% 
  filter(val == max(val)) %>% 
  left_join(
    simdata_all %>% 
      filter(str_detect(par, "scd")),
    by = c("comp", "class", "chains", "par", "cnum")
  ) %>% 
  ggplot(
    aes(
      x = press.y,
      y = val.y,
      color = class,
      group = class
    )
  ) +
  facet_grid(rows = vars(par), cols = vars(chains), scales = "free_y", switch = 'y') +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  scale_color_manual(values = chroma_ppe) +
  theme_pubr() +
  labs(
    #title = "MD results summaries",
    x = "Simulation pressure",
    y = "Parameter value"
  ) +
  guides(color = "none")

#SCRATCH

# how about a raster plot of carbon vs pressure (SCD data dump)?  
grid_scd = crossing(
  cnum = seq(3, 20, 0.1),
  press = seq(0, max(scds$press), max(scds$press)/100)
)
  
scd_models = scds %>% 
  group_by(comp, tail) %>% 
  summarise(mod_loe = loess(scd ~ press * cnum, data = cur_data()) %>% list())

scd_scapes = grid_scd %>% 
  full_join(scd_models, by = character()) %>% 
  group_by(comp, tail) %>% 
  mutate(scd = predict(mod_loe[[1]], cur_data())) %>% 
  select(-mod_loe)

comps_sel = c("18:0/18:1 PE", "18:0/18:1 PPE")
scd_scapes %>% 
  #filter(tail == "sn2") %>% 
  #filter(comp %in% comps_sel) %>% 
  drop_na() %>% 
  ggplot(
    aes(
      x = press,
      y = cnum,
      z = scd,
      fill = scd
    )
  ) +
  facet_grid(
    rows = vars(comp),
    cols = vars(tail)
  ) +
  geom_raster() +
  #geom_point(data = scds) +
  scale_fill_viridis_c(option = "plasma") +
  scale_y_reverse() +
  theme_pubr()
ggsave(here("04-suppfigs", "20220609_scd_datadump.pdf"), width=6, height=12)
 
  