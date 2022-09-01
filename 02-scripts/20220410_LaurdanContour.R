library(here)
library(tidyverse)
library(ggplot2)
library(metR)
library(ggpubr)

source(here("02-scripts", "20220502_hpsaxs_helpers.R"))

# presently this works thru symlinks
#NTS: copy data into repo prior to "deployment"
subscripts = list.files("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids", pattern = ".*laurdan_dual.R", recursive = TRUE)

# run all the individual sample scripts
# this should compile a data_gp tibble
data_gp = tibble()
subscripts %>% lapply(., source)

gp2scape = function(gpdata, grid){
  # run LOESS smoothing
  gpdata %>% 
    group_by(wl_ex, samp) %>% 
    summarize(mod = loess(
      data = cur_data(),
      formula = gp ~ T_act_mean * P_act_mean#,
      #control = loess.control(
      #  # allow extrapolation
      #  surface = "direct"
      #)
    ) %>% list()) %>% 
    right_join(grid, by="wl_ex") %>% 
    group_by(wl_ex, samp) %>% 
    mutate(gp = predict(mod[[1]], newdata=cur_data())) %>% 
    select(-mod)
}

# smoothing grid
grid = crossing(
  wl_ex = c(340, 410),
  T_act_mean = seq(0, 30, 0.1),
  P_act_mean = seq(0, 500, 2)
)

# identify and order the samples
data_gp_ord = data_gp %>% 
  mutate(
    samp = factor(
      samp,
      # set plotting order here
      levels = c(
        "JWL166P",
        "JWL191P",
        "JWL174P",
        "JWL177P"
      )
    ),
    samp = recode_factor(
      samp,
      # must match order above
      JWL166P = "Leucothea JWL166",
      JWL191P = "Beroe Sval. JWL191",
      JWL174P = "Cydippida sp. K JWL174",
      JWL177P = "Platyctene sp. T JWL177"
    )
  )

# calculate 410/340 ratio or diff
data_dgp = data_gp_ord %>% 
  mutate(state = ifelse(!(state%%2), state-1, state)) %>% # round state down to the nearest odd number
  group_by(samp, state) %>% 
  arrange() %>% 
  summarize(
    across(contains("mean"), mean),
    gp = last(gp) - first(gp)
  )

# define the contour breaks
bwid = 0.025

# from https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
sigfig <- function(vec, n=3){ 
  ### function to round values to N significant digits
  # input:   vec       vector of numeric
  #          n         integer is the required sigfig  
  # output:  outvec    vector of numeric rounded to N sigfig
  
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
  
}      # end of function   sigfig

data_gp_ord %>% 
  # run out loesses
  gp2scape(grid) %>% 
  # only one WL?
  filter(wl_ex == 340) %>% 
  # for now
  filter(!str_detect(samp, "JWL174")) %>% 
  # plot GP contours overlaid with data points
  ggplot(aes(x = T_act_mean, y = P_act_mean, z = gp)) +
  facet_grid(rows = vars(samp)) +
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
    data = data_gp_ord %>% 
    # for now
    filter(!str_detect(samp, "JWL174")) %>% 
    filter(wl_ex == 340),
    color = "black",
    alpha = 0.1,
    shape = 'o'
  ) +
  geom_text_contour(
    binwidth = bwid,
    #stroke = 0.2,
    label.placement = label_placement_fraction(0.5),
    skip = 0,
    check_overlap = TRUE,
    aes(label = sigfig(stat(level), 3))
  ) +
  theme_pubr() +
  scale_y_reverse(limits = rev(range(grid$P_act_mean))) +
  scale_x_continuous(limits = range(grid$T_act_mean)) +
  labs(
    #x = element_blank(),
    x= "Temperature (°C)",
    y = "Pressure (bar)",
    fill = "GP ratio"
  ) +
  #guides(x = "none") +
  theme(legend.position = "none") +
  ggtitle("C-laurdan GP landscapes for four ctenophores")

ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220410_ContourLaurdan.pdf", width = 4, height = 8)
