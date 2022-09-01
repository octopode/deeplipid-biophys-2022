library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pbapply)

dir_lsm_gp   = here("01-rawdata", "20220724_LSM_GP_final","G0.25", "Histograms")
files_lsm_gp = list.files(dir_lsm_gp, full.names =  TRUE)

# load ALL datafiles
data_lsm_gp = files_lsm_gp %>% 
  pblapply(
    .,
    function(file){
      fname = basename(file)
      file %>% 
        read_tsv() %>% 
        mutate(file = fname)
    }
  ) %>% 
  bind_rows() %>% 
  separate(file, into = c("sp", "temp_incu"), remove = FALSE)  %>% 
  mutate(
    sp        = sp %>% tolower() %>% c("boli" = "Boli_infu", "batho" = "Bath_fost")[.],
    temp_incu = temp_incu %>% str_extract("[0-9]+") %>% as.numeric()
  ) %>% 
  # discard the old normalization
  select(
    `GP values (GFactor-corrected)`,
    `Counts (Pixels)`,
    file,
    sp,
    temp_incu
  ) %>% 
  rename(
    gp = `GP values (GFactor-corrected)`,
    count = `Counts (Pixels)`
  )

data_lsm_norm = data_lsm_gp %>% 
  group_by(file, sp, temp_incu) %>% 
  # normalize
  mutate(count = count/sum(count))

# These are the image #s used in Fig. 1B
fig1set = c(
  29,
  69,
  73,
  83,
  101,
  104
)

data_fig1 = data_lsm_norm %>% 
  full_join(tibble(imgnum = fig1set), by = character()) %>% 
  filter(str_detect(file, imgnum %>% as.character() %>% paste(., "masked", sep='_')))

# Panels for Fig. 1B
data_fig1 %>% 
  arrange(temp_incu) %>% 
  mutate(temp_incu = temp_incu %>% factor(., levels = unique(.))) %>% 
  ggplot(
    aes(
      x = gp,
      y = count,
      color = temp_incu,
      group = file
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
  ) +
  lims(x = c(-1, 1))
ggsave(here("03-mainfigs", "depthdists", "20220726_clau_histos.pdf"), width=6, height=3)

# Barplot for Fig. S
data_gp_mean = data_lsm_norm %>% 
  group_by(sp, temp_incu, file) %>% 
  summarise(
    gp_mean = weighted.mean(gp, count)
  )

data_gp_mean %>% 
  ungroup() %>% 
  arrange(desc(sp)) %>% 
  mutate(sp = sp %>% factor(., levels=unique(.))) %>% 
  ggplot(
    aes(
      x = temp_incu %>% as.factor(),
      y = gp_mean,
      fill = temp_incu %>% as.character(),
      group = temp_incu
    )
  ) +
  facet_wrap(~sp, scales = "free") +
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(
    values = c(
      "2"  = "#1A2A6C",
      "10" = "#FDBB2C",
      "20" = "#B2201F"
    )
  ) +
  theme_pubr() +
  labs(
    x = "Incubation temperature (°C)",
    y = "Mean C-Laurdan GP"
  ) +
  theme(legend.position = "none")
ggsave(here("04-suppfigs", "20220726_lsm_gp_dists.pdf"), width=6, height=3)
