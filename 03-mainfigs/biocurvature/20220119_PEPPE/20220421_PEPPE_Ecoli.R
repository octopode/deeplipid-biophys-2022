library(tidyverse)
library(ggplot2)
library(ggpubr)

file_od600 = "/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/03-mainfigs/biocurvature/20220422_OD600-master.tsv"

col_width = 150
jit_width = 20

# lipid class colors
chroma_cl = c(
  brewer.pal(4, "Blues")[2:4], #PCs
  brewer.pal(4, "Greens")[2:4], #PEs
  brewer.pal(12, "Paired")[5:12] #PSs
) %>%
  setNames(
    c(
      "P-PC", "PC", "LPC",
      "P-PE", "PE", "LPE",
      "PS", "LPS",
      "DG",
      "PI",
      "SM",
      "Cer",
      "PG", "TG"
    )
  )

# load data and calc changes
data_ecoli = file_od600 %>% 
  read_tsv() %>% 
  mutate(
    delta_od = od_end - od_beg,
    doubling = log(od_end / od_beg, 2),
    strain = factor(strain, levels = c("PE", "PPE", "PC"))
  )

# get means
data_ecoli_summ = data_ecoli %>% 
  group_by(expt, strain, press) %>% 
  summarize(
    avg_delta_od = mean(delta_od),
    ser_delta_od = sd(delta_od)/sqrt(n()),
    avg_doubling = mean(doubling),
    ser_doubling = sd(doubling)/sqrt(n())
  )

# plot growth
data_ecoli_summ %>% 
  ggplot(aes(x = press, y = avg_delta_od, fill = strain)) +
  facet_grid(cols = vars(expt)) +
  geom_col(
    position = position_dodge(),
    width = col_width
  ) +
  geom_point(
    data = data_ecoli,
    aes(y = delta_od),
    position = position_jitterdodge(
      jitter.width = jit_width,
      dodge.width = col_width
    ),
    color = "black"#,
    #shape = 'o'
  ) +
  geom_errorbar(
    aes(
      group = strain,
      ymin = avg_delta_od - ser_delta_od,
      ymax = avg_delta_od + ser_delta_od
    ),
    position = position_dodge(width = col_width),
    width = 25
  ) +
  scale_x_continuous(breaks = c(0, 250, 500)) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(
    title = "PE vs PPE: delta OD",
    x = "Pressure (bar) for 24 h @ 37°C",
    y = "Change in OD600"
  )
 ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/biocurvature/20220119_PEPPE/20220421_PEPPE_delta_OD600.pdf", width = 4, height = 4)

# plot doublings
data_ecoli_summ %>% 
  ggplot(aes(x = press, y = avg_doubling, fill = strain)) +
  geom_col(
    position = position_dodge(),
    width = col_width
  ) +
  geom_point(
    data = data_ecoli,
    aes(y = doubling),
    position = position_jitterdodge(
      dodge.width = col_width
    ),
    color = "black",
    shape = 'o'
  ) +
  geom_errorbar(
    aes(
      group = strain,
      ymin = avg_doubling - ser_doubling,
      ymax = avg_doubling + ser_doubling
    ),
    position = position_dodge(width = col_width),
    width = 25
  ) +
  scale_x_continuous(breaks = c(0, 250, 500)) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(
    title = "PE vs PC: doublings by OD600",
    x = "Pressure (bar) for 24 h @ 37°C",
    y = "Doublings"
  )
ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/Ecoli/20220328_PEPC/20220329_PEPC_doublings.pdf", width = 6, height = 4)

## SURVIVAL
file_surv = "/Users/jwinnikoff/Documents/MBARI/Lipids/Ecoli/20220328_PEPC/20220329_PEPC_survival.tsv"

data_surv = file_surv %>% 
  read_tsv() %>% 
  mutate(
    # calculate cfu/µL; frogger deposits 3 µL
    surv = count * dil / 3,
    strain = factor(strain, levels = c("PE", "PC"))
  )

data_surv_summ = data_surv %>% 
  group_by(strain, press_set) %>% 
  summarize(
    avg_surv = mean(surv),
    ser_surv = sd(surv)/sqrt(n())
  )

# plot survival
data_surv_summ %>% 
  ggplot(aes(x = press_set, y = avg_surv, fill = strain)) +
  geom_col(
    position = position_dodge(),
    width = col_width
  ) +
  geom_point(
    data = data_surv,
    aes(y = surv),
    position = position_jitterdodge(
      dodge.width = col_width
    ),
    color = "black",
    shape = 'o'
  ) +
  geom_errorbar(
    aes(
      group = strain,
      ymin = avg_surv - ser_surv,
      ymax = avg_surv + ser_surv
    ),
    position = position_dodge(width = col_width),
    width = 25
  ) +
  scale_x_continuous(breaks = c(0, 250, 500)) +
  scale_y_log10() +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(
    title = "PE vs PC: survival post-decompression",
    x = "Pressure (bar) for 24 h @ 37°C",
    y = "Survival (cfu/µL)"
  )
ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/Ecoli/20220328_PEPC/20220329_PEPC_survival_raw.pdf", width = 6, height = 4)
