library(here)
library(tidyverse)
library(tidymodels)
library(ggplot2)
library(ggpubr)

# locate data files
files_kindat = list.files(here(), pattern = "*.tsv")

# load them all
kindat_all = files_kindat %>% 
  lapply(., function(x){read_tsv(x, col_types = cols(vol = col_double())) %>% mutate(filename = x)}) %>% 
  bind_rows()

# simple viz SAPPE 250 bar
kindat_all %>% 
  filter(str_detect(filename, "250bar")) %>%
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
  labs(
    x = "Time (s)",
    y = "Intensity (RFU)"
  )
ggsave(here("20220809_SAPPE_250bar_raw.pdf"), height = 4, width = 6)

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
