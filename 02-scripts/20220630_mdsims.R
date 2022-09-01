# Plotting for MD of low-complexity synthetic lipid systems
# 20220630 version includes SCD and can plot plateau value across pressure

library(here)
library(tidyverse)
library(broom)
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

# load SCD data
scds = read_tsv(here("01-rawdata", "20220609_Scd.tsv"))  %>% 
  pivot_longer(cols = c(sn1, sn2), names_to = "tail", values_to = "scd") %>% 
  separate(comp, c("chains", "class"), sep=' ', remove=FALSE) %>% 
  dplyr::rename(
    par = tail,
    val = scd
  ) %>% 
  mutate(par = paste("scd", par, sep='_'))

# get plateau SCD data
# Defining "plateau" as the carbon with maximum SCD @ 1 bar
scd_plat = scds %>% 
  group_by(comp, class, chains, par) %>% 
  filter(press_sim == 0) %>% 
  filter(val == max(val)) %>% 
  left_join(
    scds %>% 
      filter(str_detect(par, "scd")),
    by = c("comp", "class", "chains", "par", "cnum"),
    suffix = c(".x", '')
  ) %>% 
  select(val, press_sim)

# bind in plateau SCD data
simdata_all = simdata_summ %>% 
  bind_rows(scd_plat %>% dplyr::rename(mean = val)) %>% 
  replace_na(list("serr" = 0))

# FIG. 3 PANELS

# slope comparison tests: is class*press_sim significant within chains?
# the p-vals are extracted and printed in facet strips below for convenience
# appear to be significant for plateau SCD and for PUFA PPD
# (at least without serr weighting)
simdata_mods = simdata_all %>% 
  # replace SEMs of 0 with 1 for this computation
  # to avoid errors
  mutate(serr = ifelse(serr, serr, 1)) %>% 
  group_by(chains, par) %>% 
  summarize(
    mod_par = 
      ifelse(
        par != "dapl",
        lm(
          data = cur_data(),
          formula = mean ~ class*press_sim,
          weights = 1/(serr**2)
        ) %>% 
          # check out this less horrible new way to unpack model fits!
          broom::tidy() %>% 
          list(),
        lm(
          data = cur_data(),
          # lock to origin for % APL
          formula = mean ~ class*press_sim+0,
          weights = 1/(serr**2)
        ) %>% 
          # check out this less horrible new way to unpack model fits!
          broom::tidy() %>% 
          list()
      )
  ) %>% 
  rowwise() %>% 
  mutate(
    pval = mod_par %>% 
      filter(term == "classPE:press_sim") %>% 
      .$p.value
  )

library(BSDA)
# run pairwise t-tests
nmin = 600 # min number of independent sample blocks
simdata_ttests = simdata_all %>% 
  group_by(chains, par, press_sim) %>%
  arrange(class) %>% 
  summarise(
    ttest = tsum.test(
      mean.x = cur_data()$mean[[1]],
      s.x    = cur_data()$serr[[1]] * sqrt(nmin),
      n.x    = nmin,
      mean.y = cur_data()$mean[[2]],
      s.y    = cur_data()$serr[[2]] * sqrt(nmin),
      n.y    = nmin
    ) %>% tidy() %>% list()
  ) %>% 
  unnest_wider(ttest) %>% 
  group_by(chains, par) %>% 
  # CONSERVATIVE correction
  mutate(p.value = p.adjust(p.value, method = "bonferroni"))

# show the p-values
simdata_ttests %>% 
  filter(par == "dapl")

# turn it on end!
# no dioleoyl, for main text
simdata_all %>% 
  filter(
    #(ph == "la") &
      (chains != "18:1/18:1")
  ) %>% 
  left_join(
    simdata_mods %>% 
      select(chains, par, pval),
    by = c("chains", "par")
  ) %>% 
  arrange(par,desc(chains)) %>% 
  rowwise() %>% 
  mutate(
    chainpar = paste(chains, par, '\n', round(pval, 3)),
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
  scale_color_manual(values = chroma_ppe) +
  scale_y_reverse() +
  theme_pubr() +
  labs(
    title = "MD results summaries",
    y = "Simulation pressure",
    x = "Parameter value"
  ) +
  guides(color = "none")
ggsave(here("03-mainfigs", "synsim", "20220630_sims_1.pdf"), width=14, height=3)

# SCD profiles
# raw SCD profiles (probably for SI)
scds %>% 
  filter(press_sim %in% c(0, 1000)) %>% 
  arrange(press_sim) %>% 
  mutate(press_sim = press_sim %>% factor(., levels = unique(.))) %>% 
  ggplot(
    aes(
      x = cnum,
      y = val,
      color = class,
      linetype = press_sim
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
ggsave(here("04-suppfigs", "20220630_scd_profiles_1-1000bar.pdf"), width=4, height=3)

# plot delta SCD profiles between two specified pressures
presslims = c(0, 1000)
scds %>% 
  filter(str_detect(par, "scd")) %>% 
  filter(press_sim %in% presslims) %>% 
  group_by(comp, class, chains, par, cnum) %>% 
  arrange(press_sim) %>% 
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
ggsave(here("03-mainfigs", "synsim", "20220630_delta_scd_1-1000bar.pdf"), width=4, height=3)

# SCRATCH
# some really simple facets (from 20220510)

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


