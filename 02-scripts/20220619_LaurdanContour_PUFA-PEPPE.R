# refresh plot of Cubette data at a specified interval

library(tidyverse)
library(ggpubr)
library(gridExtra)

img_out = "~/Data/20220617_JWL0235_PUFA-PPE_POPG/JWL0235_dash.pdf"

dir_data = here("01-rawdata")
pat_data = ".*JWL023.*.tsv$" # filename pattern

# load files (once!)
files_data = list.files(path = dir_data, pattern = pat_data, full.names = T, recursive = T)

se = function(x){sd(x)/sqrt(length(x))}

# dark mode
bgcolor = "white"

lop = 2 # how many leading readings to dump

# GP contour bin width
bwid = 0.025

# smoothing grid
grid = crossing(
  wl_ex = c(340, 410),
  T_act_mean = seq(0, 30, 0.1),
  P_act_mean = seq(0, 626, 2)
)

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

## HELPER FUNCTIONS

# main function
load_visco = function(dir, pat){
  lapply(
    list.files(path = dir, pattern = pat, full.names = T, recursive = T),
    function(file){read_tsv(file, col_types = datatypes) %>% mutate(filename = basename(file))}
  ) %>% 
    bind_rows()
}

# for the intensity and GP panels
raw2fluor = function(rawdata, lop = 2){
  rawdata %>% 
    filter(shutter) %>% 
    #separate(msg, into=c("key", "val"), sep='_')
    group_by(samp, state, wl_ex, wl_em) %>% 
    arrange(watch) %>% 
    # drop the first reading(s)
    filter(row_number() > lop) %>% 
    mutate(n = n()) %>% 
    group_by(samp, state, wl_ex, wl_em, P_set, T_set, n) %>% 
    summarize(across(everything(), list(mean = mean, se = se)))
}

# for the GP panels
fl2gp = function(fluordata){
  fluordata %>% 
    group_by(samp, state, wl_ex) %>% 
    #group_by(across(c(-intensity_mean, -intensity_se, -wl_em))) %>%
    arrange(wl_em) %>% 
    mutate(across(contains("_act_mean"), mean)) %>% 
    group_by(across(contains("_act_mean")), samp, state, wl_ex) %>% 
    summarize(
      gp = (first(intensity_mean)-last(intensity_mean))/sum(intensity_mean),
      clock_mean = mean(clock_mean)
    )
}

gp2scape = function(gpdata, grid, ...){
  # run LOESS smoothing
  gpdata %>% 
    group_by(wl_ex, samp) %>% 
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
    group_by(wl_ex, samp) %>% 
    mutate(gp = safely(predict, otherwise = rep(NA, nrow(cur_data())), quiet=FALSE)(mod[[1]], newdata=cur_data()) %>% .$result) %>% 
    select(-mod)
}

# "Cubette state" plot
plot_pt = function(data_raw, bgcolor = "white"){
  data_raw %>%
    select(clock, T_act, P_act) %>%
    arrange(clock) %>%
    mutate(age = difftime(last(clock), clock, units = "sec") %>% as.numeric()) %>%
    ggplot(aes(x = T_act, y = P_act, color = age)) +
    geom_point(alpha = 0.017) +
    xlim(limits = range(grid$T_act_mean)) +
    scale_y_reverse(limits = rev(range(grid$P_act_mean))) +
    scale_color_gradient(low="red", high=ifelse(bgcolor == "black", "white", "black"), na.value = "white", limits = c(0, 600)) +
    theme_pubk() +
    #ifelse(bgcolor == "black", theme_pubk(), theme_pubr()) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank()
    ) +
    labs(
      x = "Experimental temperature (°C)",
      #y = "Experimental pressure (bar)",
      y = " ",
      title = "Cubette state"
    ) +
    guides(
      x = "none"
    )
}

# WonB theme for slide figures
theme_pubk = function(...){
  theme_pubr(...) +
    theme(
      # axis options
      axis.line  = element_line(color = "white"),
      axis.ticks = element_line(color = "white"),
      axis.text  = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      # legend
      legend.background = element_blank(),
      legend.key   = element_rect(color = NA,  fill = NA),
      legend.text  = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      # panel
      panel.background = element_blank(),
      #panel.border = element_blank()
      # facetting
      strip.background = element_rect(fill = NA, color = "white"),
      strip.text       = element_text(color = "white"),
      # plot options
      #plot.background = element_blank(),
      plot.background = element_rect(color = NA,  fill = "black"),
      plot.title = element_text(color = "white")
    )
}

## END HELPERS
## START LOOP

# load
data_raw = load_visco(dir_data, pat_data) %>% 
  # parse sample number
  separate(filename, sep='_', into = c("date", "samp", "dye", "fnum"))

# parse
data_gp = data_raw %>% 
  raw2fluor() %>% 
  fl2gp()# %>% 
  #mutate(samp = dir_data %>% strsplit('_') %>% unlist() %>% .[[2]])

## create state plot
#ptplot = data_raw %>% 
#  plot_pt(bgcolor = bgcolor)

# create contour plot
scape_gp = data_gp %>% 
  filter(wl_ex == 340) %>% 
  drop_na(gp) %>% 
  # run out loesses
  gp2scape(grid)

#print(scape_gp) #TEST
scape_gp %>% 
  # keep it reasonable to avoid too many contour lines
  filter((gp >= -1) & (gp <= 1)) %>% 
  # plot GP contours overlaid with data points
  ggplot(aes(x = T_act_mean, y = P_act_mean, z = gp)) +
  facet_grid(cols = vars(samp)) +
  geom_contour_filled(
    aes(fill = stat(level_mid)),
    binwidth = bwid
  ) +
  geom_contour(
    binwidth = bwid,
    color = "black"
  ) +
  scale_fill_distiller(
    palette = "YlGnBu", 
    direction = 1,
    guide = guide_colorbar(
      direction = "vertical"
    )
  ) +
  geom_point(
    data = data_gp,
    color = "black",
    alpha = 0.1,
    shape = 'o'
  ) +
  theme_pubr() +
  #ifelse(bgcolor == "black", theme_pubk(), theme_pubr()) +
  scale_y_reverse(limits = rev(range(grid$P_act_mean))) +
  scale_x_continuous(limits = range(grid$T_act_mean)) +
  labs(
    #x = element_blank(),
    x= "Experimental temperature (°C)",
    y = "Experimental pressure (bar)",
    fill = "C-laurdan\nGP"
  ) +
  #guides(x = "none") +
  theme(legend.position = "right") +
  ggtitle("C-laurdan GP")
ggsave(here("05-slides", "20220616_PUFA-PEPPE_claurdan.pdf"), width=12, height=4.75)
