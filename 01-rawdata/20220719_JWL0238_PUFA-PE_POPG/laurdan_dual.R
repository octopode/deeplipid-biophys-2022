### DUAL-WAVELENGTH LAURDAN DASHBOARD

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)

se = function(x){sd(x)/sqrt(length(x))}

## USER PARAMETERS
lop = 2 # how many leading readings to dump

# smoothing grid
grid = crossing(
  wl_ex = c(340, 410),
  T_act_mean = seq(0, 30, 0.1),
  P_act_mean = seq(0, 500, 2)
)

## smaller
#grid = crossing(
#  wl_ex = c(340, 410),
#  T_act_mean = seq(0, 30, 1),
#  P_act_mean = seq(0, 500, 20)
#)

datatypes = list(
  "P_act"     = col_number(),
  "P_set"     = col_number(),
  "T_act"     = col_number(),
  "T_ext"     = col_number(),
  "T_int"     = col_number(),
  "T_set"     = col_number(),
  "air"       = col_logical(),
  "clock"     = col_datetime(),
  "dewpt"     = col_number(),
  "dht_H"     = col_number(),
  "dht_T"     = col_number(),
  "intensity" = col_number(),
  "mlx_amb"   = col_number(),
  "mlx_inf"   = col_number(),
  "msg"       = col_character(),
  "n_read"    = col_integer(),
  "shutter"   = col_logical(),
  "slit_em"   = col_number(),
  "slit_ex"   = col_number(),
  "state"     = col_integer(),
  "time"      = col_number(),
  "vol"       = col_number(),
  "watch"     = col_number(),
  "wl_em"     = col_number(),
  "wl_ex"     = col_number()
)

# main function
load_visco = function(files){
  lapply(
    files,
    function(file){read_tsv(file, col_types = datatypes)}
  ) %>% 
    bind_rows()
}

# for the intensity and GP panels
raw2fluor = function(rawdata, lop = 2){
  rawdata %>% 
    filter(shutter) %>% 
    #separate(msg, into=c("key", "val"), sep='_')
    group_by(state, wl_ex, wl_em) %>% 
    arrange(watch) %>% 
    # drop the first reading(s)
    filter(row_number() > lop) %>% 
    mutate(n = n()) %>% 
    group_by(state, wl_ex, wl_em, P_set, T_set, n) %>% 
    summarize(across(everything(), list(mean = mean, se = se)))
}

# for the GP panels
fl2gp = function(fluordata){
  fluordata %>% 
    group_by(state, wl_ex) %>% 
    #group_by(across(c(-intensity_mean, -intensity_se, -wl_em))) %>%
    arrange(wl_em) %>% 
    mutate(across(contains("_act_mean"), mean)) %>% 
    group_by(across(contains("_act_mean")), state, wl_ex) %>% 
    summarize(
      gp = (first(intensity_mean)-last(intensity_mean))/sum(intensity_mean),
      clock_mean = mean(clock_mean)
    )
}

gp2scape = function(gpdata, grid, ...){
  # run LOESS smoothing
  gpdata %>% 
    group_by(wl_ex) %>% 
    summarize(mod = loess(
      data = cur_data(),
      formula = gp ~ T_act_mean * P_act_mean,
      control = ...
    ) %>% list()) %>% 
    right_join(
      grid %>% 
        # DO NOT EXTRAPOLATE
        filter(
          P_act_mean %>% between(min(gpdata$P_act_mean), max(gpdata$P_act_mean)) &
            T_act_mean %>% between(min(gpdata$T_act_mean), max(gpdata$T_act_mean))
          ),
      by="wl_ex"
    ) %>% 
    group_by(wl_ex) %>% 
    mutate(gp = safely(predict, otherwise = rep(NA, nrow(cur_data())), quiet=FALSE)(mod[[1]], newdata=cur_data()) %>% .$result) %>% 
    select(-mod)
}

gpfacet = function(smdata, gpdata){
  # plot GP contours overlaid with data points
  smdata %>% 
    ggplot(aes(x = T_act_mean, y = P_act_mean)) +
    facet_wrap(~wl_ex, ncol = 1) +
    geom_contour_filled(
      aes(z = gp, fill = stat(level_mid)),
      binwidth = 0.025#,
      #breaks = seq(-0.5, 0.5, 0.025)
    ) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1) +
    geom_point(
      data = gpdata,
      color = "black",
      alpha = 0.1,
      shape = 'o'
    ) +
    theme_pubr() +
    scale_y_reverse(limits = rev(range(smdata$P_act_mean))) +
    scale_x_continuous(limits = range(smdata$T_act_mean)) +
    labs(
      x = "Temperature (°C)",
      y = "Pressure (bar)",
      fill = "GP ratio"
    )
}

# SPECIFY FILES
#dir_data = here("20220308_JWL166P")
#dir_data = here("20220308_JWL174P")
#dir_data = here("20220309_JWL177P")
dir_data = here("20220719_JWL0238_PUFA-PE_POPG")
pat_data = "\\.tsv"

files = list.files(dir_data, pat_data, full.names = TRUE)

if(!exists("data_gp")){
  data_gp = tibble()
}

# load and parse
data_raw = files %>% 
  load_visco()

data_gp = data_raw %>% 
  raw2fluor() %>% 
  fl2gp() %>% 
  mutate(samp = dir_data %>% strsplit('_') %>% unlist() %>% .[[2]])# %>% 
  #bind_rows(data_gp)# %>% 
  #drop_na(clock_mean)

bwid = 0.025

## original plot
#data_gp %>% 
#  gp2scape(grid) %>% 
#  gpfacet(data_gp) %>% 
#  (function(x){x + ggtitle("JWL191P: Beroe")})(.)