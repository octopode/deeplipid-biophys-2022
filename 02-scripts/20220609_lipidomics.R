## NOTE: This script is a copy of 20220502_lipidomics with A. Armando's side-by-side correction factors
## applied to data collected in 2021.

library(here)
library(tidyverse)
library(ggpubr)
library(pbapply)
library(RColorBrewer)
library(spatstat)

source(here("02-scripts", "20220502_lipidomics_helpers.R"))

# this directory should contain Lipid Maps spreadsheets as TSVs
dir_data   = here("01-rawdata", "lipidmaps")
pat_data   = "tsv" # filename pattern; UNIX glob
files_data = list.files(path = dir_data, pattern = pat_data, full.names = T)

# this file contains corrections for erroneous response factors in 2021 data
file_corr = here("01-rawdata", "20220609_lcms_corrfacs.tsv")

# this file joins EID (JWL#) to species, SID (e.g. D1234-SS1), and environmental data for ALL samples
file_meta = here("01-rawdata", "20220509_metadata.tsv")

# this file contains intrinsic curvature data gleaned from Dymond 2021 and elsewhere
file_curv = here("01-rawdata", "c0_dymond_plus.tsv")

# set cutoffs for shallow and cold slices, analogous to Winnikoff et al. 2021
#cold = 5.0
#shal = 200
cold = 10
shal = 500

# collate all the Lipid Maps datafiles
lcmsdata_wide = files_data %>%
  pblapply(
    cl  = 1L,
    function(file_data){
      message(file_data)
      data_this_file <- file_data %>%
        # strip comment lines off the top
        read_tsv(skip = 7) %>%
        dplyr::rename(eid = `Sample ID`) %>%
        # replace "ND" with zero
        mutate_all(Vectorize(function(x){ifelse(x=="ND", 0, x)})) %>%
        # drop "invented" columns and rows
        select(-contains("...")) %>%
        mutate(across(!one_of("eid"), as.numeric)) %>% 
        drop_na()
      return(data_this_file)
    }
  ) %>%
  #plyr::join_all(by="eid", type="left") %>% 
  # more robust way to do the combined join/bind:
  bind_rows() %>% 
  group_by(eid) %>%
  summarise(across(everything(), function(x){c(na.omit(x), 0)[[1]]})) %>% 
  # strip the totals
  #rowwise() %>%
  filter(str_detect(eid, "JWL"))# %>%
  #group_by(eid)

# pivot data to long format and parse the LM nomenclature
lcmsdata_long = lcmsdata_wide %>%
  # `rab` stands for "relative/raw abundance"
  pivot_longer(cols = -one_of("eid"), names_to = "id", values_to = "rab") %>%
  #separate(id, into = c("class", "stuff")) %>% pull(id) %>% unique
  # see *lipidomics_helpers.R to understand/check outputs
  parse_lipidmaps_id() %>%
  # corrections
  group_by(id) %>%
  mutate(
    orig = eid %in% eids_orig,
    class_str = as.character(class)
  ) %>% 
  left_join(
    read_tsv(file_corr) %>% mutate(orig = TRUE),
    by = c("class_str", "orig")
  ) %>% 
  replace_na(list("corrfac" = 1)) %>% 
  mutate(rab = rab*corrfac) %>%
  # renormalize
  group_by(eid) %>%
  mutate(frac_molar = rab/sum(rab, na.rm = TRUE)) %>%
  # dealing with compound non-overlap
  replace_na(list("frac_molar" = 0)) %>% 
  group_by(id)

# load metadata
metadata = file_meta %>% read_tsv()

# join metadata
lcmsdata = lcmsdata_long %>%
  left_join(metadata, by = "eid") %>% 
  # create a factor allowing the samples to be ordered by species, temp, and depth
  group_by(sp) %>% 
  mutate(
    across(ends_with("_col"), mean, .names = "{.col}_spmean"),
    cold = temp_col_spmean <= cold,
    shal = depth_col_spmean <= shal
  ) %>% 
  ungroup() %>% 
  # this can be adjusted, enables precise systematic control of sample order in plots
  arrange(
    -shal, 
    cold,
    -temp_col_spmean,
    depth_col_spmean, 
    depth_col, 
    -temp_col
  ) %>% 
  mutate(
    sp_eid = factor(paste(sp, eid), levels = unique(paste(sp, eid))),
    sp = factor(sp, levels = unique(sp))
  ) %>% 
  # likewise make id a factor ordered by carbon count and unsaturation
  #NTS 20220503 can this be made more compact?
  mutate(
    id = factor(
      id,
      levels = lcmsdata_long %>%
        distinct(id, carbon, dbonds) %>%
        arrange(carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  )

# filter to phospholipids
lcmsdata_pl = lcmsdata %>%
  group_by(class) %>% 
  filter(str_detect(class, 'P') & !str_detect(class, 'Cer')) %>% 
  group_by(eid) %>%
  # renormalize
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # QC
  filter(
    !(eid %in% c("JWL0169")) & # this one is crap
      #tissue == "whole" & 
      !(sp %in% c("Ecol_wild", "Ecol_plpe"))
  )

lcmsdata_pl_wholectenos = lcmsdata_pl %>% 
  filter(tissue %in% c("whole", "body"))

lcmsdata_pl_origset = lcmsdata_pl %>% 
  filter(eid %in% eids_orig)

# load intrinsic curvature values
icurv = read_tsv(file_curv) %>% 
  mutate(
    c0 = mean,
    class = factor(
      class,
      levels = c(
        "PC", "P-PC", "LPC",
        "PE", "P-PE", "LPE",
        "PS", "LPS",
        "PI", "PG",
        "Cer", "SM",
        "DG", "TG"
      )
    )
  ) %>% 
  # set unknown errors to 0; they are minor overall contributors
  rowwise() %>% 
  mutate(
    lo = min(c0, lo, na.rm = TRUE),
    hi = max(c0, hi, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  select(class, c0, lo, hi)

#MAINTEXT plots

#NTS 20220503: maybe some of this can be decomposed into helpers
# get curvature contributions for each individual
curvdata = lcmsdata_pl_wholectenos %>%
  group_by(sp, eid, sp_eid, class, depth_col, temp_col) %>%
  summarize(frac_molar = sum(frac_molar)) %>%
  #group_by(sp, class) %>% 
  #summarize(frac_molar = mean(frac_molar)) %>%
  left_join(icurv, by = "class") %>% 
  # plci stands for PhosphoLipid Curvature Index
  mutate(plci = frac_molar * c0) %>% 
  drop_na() 

##plot Fig. 2D (curvature contributions)

# plot species means for a select few
curvdata %>% 
  filter(eid %in% eids_orig) %>% 
  group_by(sp, class, c0) %>% 
  summarize(plci = mean(plci)) %>%
  mutate(neg = plci<0) %>% 
  # wrangle for geom_tile
  group_by(sp) %>% 
  arrange(neg, -c0) %>% 
  mutate(
    # running total
    end = cumsum(plci),
    y = end - plci/2
  ) %>% 
  ggplot(
    aes(x = neg)
  ) +
  facet_wrap(~sp, nrow = 1, strip.position = "bottom") +
  geom_tile(
    aes(
      y = y,
      height = plci,
      fill = class
    ),
    width = 0.95,
    size = 0.05,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  # ind'lwise error bars
  geom_errorbar(
    data = lcmsdata_pl_wholectenos %>%
      filter(eid %in% eids_fig2d) %>% 
      group_by(sp, sp_eid, class) %>%
      summarize(frac_molar = sum(frac_molar)) %>%
      left_join(icurv, by = "class") %>% 
      mutate(plci = frac_molar * c0) %>% 
      drop_na() %>% 
      summarize(plci = sum(plci)) %>% 
      group_by(sp) %>% 
      summarize(
        plci_mean = mean(plci),
        plci_serr = sd(plci)/sqrt(n()),
        neg = TRUE
      ),
    aes(
      ymin = plci_mean - plci_serr,
      ymax = plci_mean + plci_serr
    ),
    width = 0.25
  ) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(color = "none") +
  labs(x = "Species", y = "PL curvature index", fill = "Lipid class", title = "Contributions to PL curvature index")
ggsave(here("03-mainfigs", "ctenolipids", "20220503_panelD.pdf"), width=6, height=3)

#SUPPINFO plots

#PLOT curvature estimates for PL classes
icurv %>% 
  # phospholipids only
  # drop unused levels, order by c0
  mutate(class = class %>% as.character()) %>% 
  filter(str_detect(class, 'P')) %>% 
  arrange(c0) %>% 
  mutate(class = class %>% factor(., levels = .)) %>% 
  ggplot(
    aes(
      x = class,
      y = c0,
      ymin = lo,
      ymax = hi,
      fill = class
    )
  ) +
  geom_col() +
  geom_errorbar(width = 0.1) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  guides(fill = "none") +
  labs(
    title = "Instrinsic curvature estimates for headgroup classes",
    y = "c0 (1/Å)"
  )
ggsave(here("04-suppfigs", "20220502_c0_estimates.pdf"), width=6, height=4)

#PLOT monster all-PLs plot
lcmsdata_pl_wholectenos %>%
  #filter(eid %in% eids_orig) %>% 
  group_by(sp_eid, class, annot) %>%
  arrange(id) %>% 
  ggplot(
    aes(
      x = sp_eid,
      y = frac_molar,
      fill = class
    )
  ) +
  # draw little hairlines in between compounds
  geom_col(width = 0.95, size = 0.05, color = "white") +
  geom_text(
    aes(
      label = ifelse(frac_molar >= 0.015, annot, ''),
      color = class
    ),
    size = 1.5,
    position = position_stack(vjust=0.5)
  ) +
  scale_fill_manual( values = chroma_cl) +
  scale_color_manual(values = chroma_tx) +
  theme_pubr() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "sample", y = "mole fraction", fill = "lipid class", title = "Phospholipid composition, all samples")
ggsave(here("04-suppfigs", "20220619_all_phospholipids.pdf"), width=18, height=7.5)

# HG classes averaged by species
lcmsdata_pl_wholectenos %>%
  group_by(sp, sp_eid, class) %>%
  summarize(frac_molar = sum(frac_molar)) %>% 
  # average
  group_by(sp, class) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # and renorm
  group_by(sp) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  arrange(class) %>% 
  ggplot(
    aes(
      x = sp,
      y = frac_molar,
      fill = class
    )
  ) +
  # draw little hairlines in between compounds
  geom_col(width = 0.95, size = 0.05, color = "black") +
  scale_fill_manual( values = chroma_cl) +
  theme_pubr() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(
    x = "Species", 
    y = "Mole fraction", 
    fill = "Headgroup", 
    title = "Phospholipid composition, species means"
  )
ggsave(here("04-suppfigs", "20220609_all_phospholipids_spmeans.pdf"), width=11, height=7.5)

## PCs and PEs plotted back-to-back (panel B)
lcmsdata_pl_wholectenos %>%
  # original species only
  filter(sp %in% c(
    "Boli_vitr",
    "Boli_infu",
    "Boli_arct",
    "Lamp_crue",
    "Tjal_pink"
  )) %>% 
  group_by(sp, sp_eid, class) %>%
  summarize(frac_molar = sum(frac_molar)) %>% 
  # average
  group_by(sp, class) %>% 
  summarize(frac_molar = mean(frac_molar)) %>% 
  # and renorm
  group_by(sp) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  arrange(class) %>% 
  # for PC/PE ratio
  group_by(sp, class) %>%
  filter(class %in% c("PC", "P-PC", "LPC", "PE", "P-PE", "LPE")) %>%
  mutate(frac_molar = ifelse(str_detect(class, "PC"), frac_molar, -1*frac_molar)) %>% 
  summarize(frac_molar = sum(frac_molar)) %>%
  # average indls of each species
  # comment to show individuals
  group_by(sp, class) %>%
  summarize(
    sem = sd(frac_molar)/sqrt(n()),
    frac_molar = mean(frac_molar)
  ) %>%
  ggplot(
    aes(
      x = sp,
      y = frac_molar,
      fill = class
    )
  ) +
  geom_col(
    width = 0.50,
    size = 0.05,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  # for displaying ratios
  scale_y_continuous(
    limits = c(-.75,.75), 
    breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),
    labels = abs(c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75))
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "sample", y = "mole fraction", fill = "lipid class", title = "Total PCs vs. PEs")
ggsave(here("03-mainfigs", "20220619_PCPE.pdf"), width=6, height=3)

#PLOT monster curvature plot
## might end up becoming Fig. S2B, if we can squeeze it under the above bar plot?
##NTS 20220503: best to do that stacking in Illustrator
curvdata %>% 
  mutate(neg = plci<0) %>% 
  # wrangle for geom_tile
  group_by(sp_eid) %>% 
  arrange(neg, -c0) %>% 
  mutate(
    # running total
    end = cumsum(plci),
    y = end - plci/2
  ) %>% #group_by(sp) %>% group_split() %>% .[[1]] %>% View
  ggplot(
    aes(x = neg)
  ) +
  facet_wrap(~sp_eid, nrow = 1, strip.position = "bottom") +
  geom_tile(
    aes(
      y = y,
      height = plci,
      fill = class
    ),
    width = 0.95,
    size = 0.05,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Species", y = "PL curvature index", fill = "Lipid class", title = "Contributions to PL curvature index")
ggsave(here("04-suppfigs", "20220506_all_curvatures.pdf"), width=10, height=7.5)

# acyl chain supp plots
# our favorite phospholipids
major_pls = c("PC", "P-PC", "PE", "P-PE", "PS")

#PLOT acyl carbon species means
lcmsdata_pl_wholectenos %>%
  # subset while maintaining class order
  arrange(class) %>% 
  filter(class %in% major_pls) %>%
  mutate(class = class %>% factor(., levels = unique(.))) %>% 
  # get total of each dbond count
  group_by(eid, sp, class, carbon) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  # renormalize within each class and sample
  group_by(eid, sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # now average samplewise
  group_by(sp, class, carbon) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # compute weighted median. ifelse() prevents overplotting
  group_by(sp, class) %>% 
  mutate(
    wtmed = ifelse(row_number() == 1, weighted.median(carbon, frac_molar, type=4), NA),
    wtavg = ifelse(row_number() == 1, sum(carbon * frac_molar), NA)
  ) %>%
  ggplot(aes(x=carbon, y=frac_molar, fill=class, alpha=-(carbon %% 2))) +
  facet_grid(rows = vars(sp), cols = vars(class)) +
  geom_col(size = 0.05, color = "white") +
  geom_vline(
    # It's an open question for now as to whether weighted mean or median is better; think about it.
    # As it is we're looking at the 50% weighted percentile; area to left and right of line are equal.
    aes(xintercept = wtmed),
    color = "black",
    alpha = 0.8,
    size = 1,
    linetype="dashed"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous() +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.4, 2)) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha="none", fill="none") +
  labs(
    title = "Acyl chain length distributions",
    x = "Total acyl carbons",
    y = "Mole fraction of lipid class"
  )
ggsave(here("04-suppfigs", "20220609_all_carbon.pdf"), width=10, height=7.5)

#PLOT acyl carbon species means, original 15
lcmsdata_pl_wholectenos %>%
  filter(eid %in% eids_orig) %>% 
  # subset while maintaining class order
  arrange(class) %>% 
  filter(class %in% major_pls) %>%
  mutate(class = class %>% factor(., levels = unique(.))) %>% 
  # get total of each dbond count
  group_by(eid, sp, class, carbon) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  # renormalize within each class and sample
  group_by(eid, sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # now average samplewise
  group_by(sp, class, carbon) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # compute weighted median. ifelse() prevents overplotting
  group_by(sp, class) %>% 
  mutate(
    wtmed = ifelse(row_number() == 1, weighted.median(carbon, frac_molar, type=4), NA),
    wtavg = ifelse(row_number() == 1, sum(carbon * frac_molar), NA)
  ) %>%
  ggplot(aes(x=carbon, y=frac_molar, fill=class, alpha=-(carbon %% 2))) +
  facet_grid(rows = vars(sp), cols = vars(class)) +
  geom_col(size = 0.05, color = "white") +
  geom_vline(
    # It's an open question for now as to whether weighted mean or median is better; think about it.
    # As it is we're looking at the 50% weighted percentile; area to left and right of line are equal.
    aes(xintercept = wtmed),
    color = "black",
    alpha = 0.8,
    size = 1,
    linetype="dashed"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous() +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.4, 2)) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha="none", fill="none") +
  labs(
    #title = "Acyl chain length distributions",
    x = "Total acyl carbons",
    y = "Mole fraction of lipid class"
  )
ggsave(here("04-suppfigs", "20220609_orig_carbon.pdf"), width=10, height=7.5)

#PLOT unsaturation species means
lcmsdata_pl_wholectenos %>%
  # subset while maintaining class order
  arrange(class) %>% 
  filter(class %in% major_pls) %>%
  mutate(class = class %>% factor(., levels = unique(.))) %>% 
  # get total of each dbond count
  group_by(eid, sp, class, dbonds) %>%
  replace_na(list("dbonds" = 0)) %>% # why needed?
  summarise(frac_molar = sum(frac_molar)) %>%
  # renormalize within each class and sample
  group_by(eid, sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # now average samplewise
  group_by(sp, class, dbonds) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # compute weighted median. ifelse() prevents overplotting
  group_by(sp, class) %>% 
  mutate(
    wtmed = ifelse(row_number() == 1, weighted.median(dbonds, frac_molar, type=4), NA),
    wtavg = ifelse(row_number() == 1, sum(dbonds * frac_molar), NA)
  ) %>%
  ggplot(aes(x=dbonds, y=frac_molar, fill=class)) +
  facet_grid(rows = vars(sp), cols = vars(class)) +
  geom_col(size = 0.05, color = "white") +
  geom_vline(
    # It's an open question for now as to whether weighted mean or median is better; think about it.
    # As it is we're looking at the 50% weighted percentile; area to left and right of line are equal.
    aes(xintercept = wtmed),
    color = "black",
    alpha = 0.8,
    size = 1,
    linetype="dashed"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous() +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.4, 2)) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha="none", fill="none") +
  labs(
    title = "Double bond distributions: both chains",
    x = "Number of double bonds",
    y = "Mole fraction of lipid class"
  )
ggsave(here("04-suppfigs", "20220609_all_dbonds.pdf"), width=10, height=7.5)

#PLOT unsaturation species means
lcmsdata_pl_wholectenos %>%
  filter(eid %in% eids_orig) %>% 
  # subset while maintaining class order
  arrange(class) %>% 
  filter(class %in% major_pls) %>%
  mutate(class = class %>% factor(., levels = unique(.))) %>% 
  # get total of each dbond count
  group_by(eid, sp, class, dbonds) %>%
  replace_na(list("dbonds" = 0)) %>% # why needed?
  summarise(frac_molar = sum(frac_molar)) %>%
  # renormalize within each class and sample
  group_by(eid, sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # now average samplewise
  group_by(sp, class, dbonds) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # compute weighted median. ifelse() prevents overplotting
  group_by(sp, class) %>% 
  mutate(
    wtmed = ifelse(row_number() == 1, weighted.median(dbonds, frac_molar, type=4), NA),
    wtavg = ifelse(row_number() == 1, sum(dbonds * frac_molar), NA)
  ) %>%
  ggplot(aes(x=dbonds, y=frac_molar, fill=class)) +
  facet_grid(rows = vars(sp), cols = vars(class)) +
  geom_col(size = 0.05, color = "white") +
  geom_vline(
    # It's an open question for now as to whether weighted mean or median is better; think about it.
    # As it is we're looking at the 50% weighted percentile; area to left and right of line are equal.
    aes(xintercept = wtmed),
    color = "black",
    alpha = 0.8,
    size = 1,
    linetype="dashed"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous() +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.4, 2)) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha="none", fill="none") +
  labs(
    #title = "Double bond distributions: both chains",
    x = "Number of double bonds",
    y = "Mole fraction of lipid class"
  )
ggsave(here("04-suppfigs", "20220524_orig_dbonds.pdf"), width=10, height=7.5)

# Tentacle vs. body comparison
lcmsdata_pl %>% 
  filter(sp == "Tjal_pink") %>% 
  group_by(sp_eid, class, annot) %>%
  arrange(id) %>% 
  mutate(depth_col = round(depth_col)) %>% 
  unite("sid", sid, depth_col, eid) %>% 
  ungroup() %>% 
  ggplot(
    aes(
      x = sid,
      y = frac_molar,
      fill = class
    )
  ) +
  facet_wrap(~tissue, scales = "free_x") +
  # draw little hairlines in between compounds
  geom_col(width = 0.95, size = 0.05, color = "white") +
  geom_text(
    aes(
      label = ifelse(frac_molar >= 0.015, annot, ''),
      color = class
    ),
    size = 1.5,
    position = position_stack(vjust=0.5)
  ) +
  scale_fill_manual( values = chroma_cl) +
  scale_color_manual(values = chroma_tx) +
  theme_pubr() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "Sample", y = "Mole fraction", fill = "Lipid class", title = "Phospholipid composition, Tjalfiella")
ggsave(here("04-suppfigs", "20220517_Tjalfie-tentacles.pdf"), width=10, height=7.5)

## E. coli STUFF

lcmsdata_ecoli = lcmsdata %>%
  group_by(class) %>% 
  filter(str_detect(class, 'P') & !str_detect(class, 'Cer')) %>% 
  # spurious PPC
  filter(class != "P-PC") %>% 
  group_by(eid) %>%
  # renormalize
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # QC
  filter(sp %in% c("Ecol_wild", "Ecol_plpe"))

# E. coli full lipidomics plot
lcmsdata_ecoli %>%
  # make a new x-axis var
  group_by(sp, depth_col) %>% 
  arrange(eid) %>% 
  mutate(rep = eid %>% as.factor() %>% as.integer()) %>% 
  group_by(sp_eid, class, annot) %>%
  arrange(id) %>% 
  ggplot(
    aes(
      x = rep,
      y = frac_molar,
      fill = class
    )
  ) +
  facet_grid(rows = vars(depth_col), cols = vars(sp), scales = "free") +
  # draw little hairlines in between compounds
  geom_col(width = 0.95, size = 0.05, color = "white") +
  geom_text(
    aes(
      label = ifelse(frac_molar >= 0.015, annot, ''),
      color = class
    ),
    size = 1.5,
    position = position_stack(vjust=0.5)
  ) +
  scale_fill_manual( values = chroma_cl) +
  scale_color_manual(values = chroma_tx) +
  theme_pubr() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "Replicate", y = "Mole fraction", fill = "Lipid class", title = "E. coli PL composition")
ggsave(here("04-suppfigs", "20220517_Ecoli_all_PLs.pdf"), width=8.5, height=7.5)

curvdata_ecoli = lcmsdata_ecoli %>%
  group_by(sp, eid, sp_eid, class, depth_col, temp_col) %>%
  summarize(frac_molar = sum(frac_molar)) %>%
  #group_by(sp, class) %>% 
  #summarize(frac_molar = mean(frac_molar)) %>%
  left_join(icurv, by = "class") %>% 
  # plci stands for PhosphoLipid Curvature Index
  mutate(plci = frac_molar * c0) %>% 
  drop_na() 

curvdata_ecoli %>% 
  # make a new x-axis var
  group_by(sp, depth_col) %>% 
  #arrange(eid) %>% 
  #mutate(rep = eid %>% as.factor() %>% as.integer()) %>% 
  arrange(sp, eid) %>% 
  mutate(rep = paste(sp, eid %>% as.factor() %>% as.integer())) %>% 
  mutate(rep = factor(rep, levels = unique(rep))) %>% 
  mutate(neg = plci<0) %>% 
  # wrangle for geom_tile
  group_by(sp_eid) %>% 
  arrange(neg, -c0) %>% 
  mutate(
    # running total
    end = cumsum(plci),
    y = end - plci/2
  ) %>% #group_by(sp) %>% group_split() %>% .[[1]] %>% View
  ggplot(
    aes(x = neg)
  ) +
  facet_grid(rows = vars(depth_col), cols = vars(rep), switch = 'x') +
  geom_tile(
    aes(
      y = y,
      height = plci,
      fill = class
    ),
    width = 0.95,
    size = 0.05,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Replicate", y = "PL curvature index", fill = "Lipid class", title = "E. coli: contributions to PL curvature index")
ggsave(here("04-suppfigs", "20220506_Ecoli_curvatures.pdf"), width=8.5, height=7.5)

## CTENO METADATA TABLE
lcmsdata_pl_wholectenos %>% 
  group_by(
    eid,
    sid,
    sp,
    depth_col,
    temp_col
  ) %>% 
  summarise() %>% 
  arrange(sp, -temp_col, depth_col) %>% 
  write_tsv(here("04-suppfigs", "20220525_metadata_supptable.tsv"))

## TEXT CLAIMS

# max PPE
lcmsdata_pl %>% filter(class == "P-PE") %>% group_by(eid) %>% summarize(frac_molar = sum(frac_molar)) %>% arrange(-frac_molar)
# Di-PUFAs
lcmsdata_pl_wholectenos %>% filter(depth_col >= 1000) %>% filter(dbonds > 7) %>% group_by(eid) %>% summarize(frac_molar = sum(frac_molar)) %>% .$frac_molar %>% mean()

# how many lipids detected?
lcmsdata_pl_wholectenos %>% filter(frac_molar > 0) %>% drop_na(frac_molar) %>% mutate(id = as.character(id)) %>% group_by(id) %>% summarise()
# how many samples?
lcmsdata_pl_wholectenos$eid %>% unique() %>% length()

# spmeans summary table (wide format)
lcmsdata_pl_wholectenos %>%
  group_by(sp, sp_eid, class) %>%
  summarize(frac_molar = sum(frac_molar)) %>% 
  # average
  group_by(sp, class) %>% 
  summarize(
    n = n(),
    mean_mf = mean(frac_molar),
    serr_mf = serr(frac_molar)
  ) %>% 
  # and renorm
  group_by(sp, n) %>% 
  mutate(
    adj = sum(mean_mf),
    mean_mf = mean_mf/adj,
    serr_mf = serr_mf/adj
  ) %>% 
  select(-adj) %>% 
  pivot_longer(cols = contains("mf"), names_to = "var", values_to = "val") %>% 
  arrange(class, var) %>% 
  # widen
  unite(var, c(class, var)) %>% 
  pivot_wider(names_from = var, values_from = val) %>% 
  write_tsv(here("04-suppfigs", "20220619_spmeans.tsv"))
