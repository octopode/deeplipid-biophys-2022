library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pbapply)
library(nsprcomp)
library(nnls)
library(lubridate)
library(magick)
library(abind)

source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_Nov2021/JWL157.R")
source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_May2021/JWL158.R")
source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_May2021/JWL161.R")
source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_May2021/JWL164.R")
source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_Nov2021/JWL167.R")
source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_Nov2021/JWL174.R")
source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_Nov2021/JWL177.R")
source("/Users/jwinnikoff/Documents/MBARI/SAXS/HPSAXS_7A_Nov2021/JWL191.R")

eid2sp = c(
  JWL157 = "Tetr_angu", # script done but the raw files have problems (bad config)
  JWL158 = "Leuc_pulc", # script done
  JWL161 = "Bero_cucu", # script done
  JWL164 = "Tjal_pink", # script done
  JWL167 = "Leuc_pulc", # script done
  JWL174 = "Cydi_blac", # script done
  JWL177 = "Tjal_pink", # script done
  JWL191 = "Bero_cucu"  # script done
)

# collate data from multiple samples
data_saxs_updn = bind_rows(
  data_saxs_157,
  data_saxs_158,
  data_saxs_161,
  data_saxs_164,
  data_saxs_167,
  data_saxs_174,
  data_saxs_177,
  data_saxs_191,
) %>% 
  filter(pdir %in% c("up", "dn"))# %>% 
  mutate(
    # identify samples
    samp = filename %>% str_split('_') %>% unlist() %>% .[[1]]#,
    # log-transform for interpolation
    #iq = log(iq)
  )# %>% 
  # downsample to q-res of 0.001
  #mutate(q = round(q, 3)) %>% 
  #group_by(across(!matches("iq"))) %>% 
  #group_by(temp, press, pdir, samp, q) %>% 
  #summarize(iq = mean(iq))

# add rows for sweep average
data_saxs_avg = data_saxs_updn %>% 
  group_by(temp, press, samp, q) %>% 
  summarise(iq = mean(iq)) %>% 
  mutate(pdir = "av")

# multisample version
data_working = data_saxs_updn %>% 
  filter(
    #samp == "JWL164",
    # cut out the hiP series 00 MPa data
    #!str_detect(filename, "hiP") # cut out the whole hiP series
    !(str_detect(filename, "hiP") & press < 1000) # cut out inconsistent low-pressure shots
  ) %>% 
  # crop beamstop
  filter(q >= 0.05) %>% 
  group_by(samp, press, temp, pdir) %>% 
  arrange(rep) %>% 
  # get ONLY last shot of each stop
  filter(rep == last(rep)) %>% 
  # normalize to area
  mutate(
    iq = iq - min(iq),
    iq = iq / sum(iq)
  ) %>% 
  filter(press != 800) %>% 
  # count replicates for weighting
  group_by(samp, q, press, temp) %>% 
  mutate(n = n()) %>% 
  # keep a sensible order
  arrange(samp, temp, press, q) %>% 
  group_by(samp) %>% 
  mutate(sp = eid2sp[samp]) %>% 
  # not sure why this mixup is happening!
  group_by(filename) %>% 
  filter(str_detect(filename, samp))

