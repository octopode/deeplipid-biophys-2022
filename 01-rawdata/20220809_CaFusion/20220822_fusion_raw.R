library(here)
library(tidyverse)
library(tidymodels)
library(ggplot2)
library(ggpubr)

# locate data files
files_kindat = list.files(here(), pattern = "2022081[7-9].*.tsv")

# load them all
kindat_raw = files_kindat %>% 
  lapply(., function(x){read_tsv(x, col_types = cols(vol = col_double())) %>% mutate(filename = x)}) %>% 
  bind_rows() %>% 
  separate(filename, into = c("date", "id", "conc_ca", "press", "rep"), sep='_', remove = FALSE) %>% 
  mutate(press = press %>% str_extract("[0-9]+") %>% as.numeric())

# plot individually for QC and cropping
# takes awhile!
kindat_raw %>% 
  #filter(between(watch, 0, 1000)) %>% 
  ggplot(
    aes(
      x = watch,
      y = intensity,
      color = inj
    )
  ) +
  facet_wrap(~filename, ncol=1, scales = "free_x") +
  scale_x_continuous(breaks=function(lims){seq(round(lims[[1]]), lims[[2]], 30)}) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
ggsave(here("20220821_QC_traces_raw_full.pdf"), height = 40, width = 12)

# a function to do the math for my silly dilution scheme
dilfac = function(n_dils){
  v_init = 400
  v_curr = v_init + 65
  if(!n_dils){return(v_init/v_curr)} # 0-dil case
  v_curr = v_curr + 135
  v_curr = v_curr + (100*(n_dils-1))
  return(v_init/v_curr)
}

# do the dequench normalization
norm_facs = kindat_raw %>% 
  left_join(read_tsv(here("dequench_intervals.tsv")), by = "filename") %>% 
  #group_by(filename) %>% 
  rowwise() %>% 
  filter(between(watch, deq_beg, deq_end)) %>% 
  group_by(filename, n_dils) %>% 
  # adjust the dequenched intensity using the dil factor
  summarize(intensity_deq = mean(intensity)) %>% 
  mutate(intensity_deq = intensity_deq/dilfac(n_dils))

# window for baseline
min_win = c(60, 120)
# align injection times and normalize
kindat_alg = kindat_raw %>% 
  group_by(filename) %>% 
  mutate(
    baseline = cur_data() %>% 
      filter(between(watch, min_win[[1]], min_win[[2]])) %>% 
      .$intensity %>% 
      mean()
  ) %>% 
  # cut the garbage
  filter(
    !(filename %in%
        c(
          "20220818_SAPPE_1mMCa_200bar_c.tsv"
        )
    )
  ) %>% 
  # normalize by dequenched intensity
  left_join(norm_facs, by = "filename") %>% 
  mutate(intensity = (intensity-baseline)/(intensity_deq-baseline)) %>% 
  # align the time traces
  group_by(filename) %>% 
  filter(inj) %>% 
  mutate(watch = watch - min(watch)) %>% 
  # now try starting them all in the same place
  filter(watch >= 300) %>% 
  mutate(idx = row_number()) %>% 
  arrange(idx) %>% 
  group_by(id, press, idx) %>% 
  mutate(intens_mean = mean(intensity)) %>% 
  group_by(id, press, filename) %>% 
  # shift to mean start of group
  mutate(intensity = intensity-intensity[[1]])

chroma = c("SAPE" = "#BBBDBF", "SAPPE" = "#BAE4B3")

# set the actual kinetic window
kin_win = c(300, 900) # a nice 10 min!

slopes = kindat_alg %>% 
  group_by(filename, id, press, date) %>% 
  filter(between(watch, kin_win[[1]], kin_win[[2]])) %>% 
  summarize(
    linreg = cur_data() %>% 
      lm(intensity ~ watch, .) %>% 
      tidy() %>% 
      slice(2) %>% 
      list()
  ) %>% 
  unnest_wider(linreg) %>% 
  # fresh runs only
  filter(
    ((id == "SAPE") & (date == "20220818")) |
      ((id == "SAPPE") & (date == "20220817"))
  ) %>% 
  # only to 100 bar
  filter(press <= 100)

# plot those slopes!
slopes %>%
  group_by(id, press) %>% 
  mutate(
    mean_slope = mean(estimate),
    mean_slope = ifelse(row_number()>1, NA, mean_slope)
  ) %>% 
  ggplot(
    aes(
      x = id,
      label = filename,
      #color = date
      fill = id
    )
  ) +
  facet_wrap(~press, ncol=2) +
  geom_col(
    aes(y = mean_slope)
  ) +
  #geom_point(
  #  aes(y = estimate),
  #  color="black"
  #) +
  scale_fill_manual(values = chroma) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(
    y = "Fusion rate (fraction/s)"
  )
ggsave(here("20220822_normd_slopes.pdf"), height = 3, width = 3)

# plot actual traces
kindat_alg %>% 
  group_by(filename, id, press, date) %>% 
  filter(between(watch, kin_win[[1]], kin_win[[2]])) %>%
  mutate(watch = watch-min(watch)) %>% 
  # fresh runs only
  filter(
    ((id == "SAPE") & (date == "20220818")) |
      ((id == "SAPPE") & (date == "20220817"))
  ) %>% 
  # only to 100 bar
  filter(press <= 100) %>% 
  # rolling average to compress the data 5?-fold
  mutate(bin = as.integer(row_number()/5)) %>% 
  group_by(id, press, filename, bin) %>% 
  summarize(across(c(watch, intensity), mean)) %>% 
  # calc linregs to shift baseline
  mutate(
    linreg = cur_data() %>% 
      lm(intensity ~ watch, .) %>% 
      tidy() %>% 
      list()
  ) %>% 
  unnest(linreg) %>% 
  mutate(term = c("(Intercept)" = "xcept", "watch" = "slope")[term]) %>% 
  pivot_wider(id_cols = c(group_cols(), "watch", "intensity"), names_from = "term", values_from = "estimate") %>% 
  group_by(id, press) %>% 
  # only show the higher-slope replicate
  filter(slope == max(slope)) %>% 
  mutate(
    intensity = intensity - xcept#,
    #intensity = ifelse(!press, intensity + 0.01, intensity)
  ) %>% 
  ggplot(
    aes(
      x = watch,
      y = intensity,
      group = filename,
      color = id
    )
  ) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", se=FALSE) +
    scale_color_manual(values = chroma) +
    facet_wrap(~press, ncol = 2) +
    theme_pubr() +
    #lims(
    #  x = c(300, 900),
    #  y = c(-0.01, 0.02)
    #) +
    theme(legend.position = "none") +
    labs(
      x = "Time (s)",
      y = "Fraction fused"
    )
ggsave(here("20220822_normd_traces.pdf"), height = 3, width = 6)


#-----

# is there a significant slope within the high-P region?
kindat_all %>% 
  filter(str_detect(filename, "250bar")) %>%
  group_by(filename) %>%
  arrange(watch) %>% 
  filter(
    (P_act == 250) &
      (lead(P_act) == 250) &
      (lag(P_act) == 250)
  ) %>% 
  summarize(linreg = lm(intensity ~ watch, data = cur_data()) %>% tidy() %>% list()) %>% 
  unnest(linreg) %>% 
  rename(slope_hiP = estimate) %>% 
  filter(term == "watch")

thres_dequench = 300

kindat_all %>% 
  #filter(intensity > 20) %>% 
  filter(str_detect(filename, "20220810b")) %>%
  group_by(filename) %>% 
  # normalize pressure to max intensity
  mutate(P_act = max(intensity)*P_act/max(P_act)) %>% 
  ggplot(
    aes(
      x = watch
    )
  ) +
  geom_line(
    aes(y = P_act),
    color = "blue"
  ) +
  geom_line(
    aes(y = intensity),
    color = "black"
  ) +
  facet_wrap(~filename, ncol = 1, scales = "free_y") +
  theme_pubr() +
  lims(
    #x = c(0, 1500),
    y = c(20, 50)
  ) +
  labs(
    x = "Time (s)",
    y = "Intensity (RFU)"
  )
ggsave(here("20220814_SAPPE_100bar_15minEa.pdf"), width=6, height=4)

kindat_all %>% 
  #filter(intensity > 20) %>% 
  filter(str_detect(filename, "20220809a")) %>%
  group_by(filename) %>% 
  # normalize pressure to max intensity
  mutate(P_act = max(intensity)*P_act/max(P_act)) %>% 
  ggplot(
    aes(
      x = watch
    )
  ) +
  geom_line(
    aes(y = P_act),
    color = "blue"
  ) +
  geom_line(
    aes(y = intensity),
    color = "black"
  ) +
  facet_wrap(~filename, ncol = 1, scales = "free_y") +
  theme_pubr() +
  lims(
    x = c(250, 1000),
    y = c(0, 50)
  ) +
  labs(
    x = "Time (s)",
    y = "Intensity (RFU)"
  )
ggsave(here("20220814_SAPE_100bar_4minEa.pdf"), width=6, height=4)

kindat_all %>% 
  #filter(intensity > 20) %>% 
  filter(str_detect(filename, "longrun")) %>%
  group_by(filename) %>% 
  # normalize pressure to max intensity
  mutate(P_act = max(intensity)*P_act/max(P_act)) %>% 
  ggplot(
    aes(
      x = watch
    )
  ) +
  geom_line(
    aes(y = P_act),
    color = "blue"
  ) +
  geom_line(
    aes(y = intensity),
    color = "black"
  ) +
  facet_wrap(~filename, ncol = 1, scales = "free_y") +
  theme_pubr() +
  lims(
    #x = c(0, 1500),
    #y = c(20, 50)
  ) +
  labs(
    x = "Time (s)",
    y = "Intensity (RFU)"
  )
ggsave(here("20220814_SAPPE_longrun.pdf"), width=6, height=4)
