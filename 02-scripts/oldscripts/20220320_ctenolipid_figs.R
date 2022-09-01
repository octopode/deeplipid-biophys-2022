library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(spatstat)

dir_data <- "/Users/jwinnikoff/Documents/MBARI/Lipids/20201122_hplc/"
pat_data <- "tsv" # filename pattern; UNIX glob
# for now, this file links EIDs to SIDs
file_meta <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/00-metadata/old/20201022_lipid_metadata.tsv"
dir_envi <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/00-metadata/"

## helpers
class2tails <- Vectorize(function(class){
  if(str_detect(class, 'L')){
    return(1L)
  }
  if(str_detect(class, 'T')){
    return(3L)
  }else{
    return(2L)
  }
})

files_data <- list.files(path = dir_data, pattern = pat_data, full.names = T)

hplcdata_wide <- files_data %>%
  lapply(
    function(file_data){
      data_this_file <- file_data %>%
        read_tsv(skip = 7) %>%
        #mutate(file = file_data %>% basename()) %>%
        dplyr::rename(eid = `Sample ID`) %>%
        # replace "ND" with zero
        mutate_all(Vectorize(function(x){ifelse(x=="ND", 0, x)})) %>%
        # drop empty columns and rows
        select(-X1, -X3) %>%
        drop_na()
      return(data_this_file)
    }
  ) %>%
  plyr::join_all(by="eid", type="left") %>%
  # strip the totals
  rowwise() %>%
  filter(str_detect(eid, "JWL")) %>%
  group_by(eid) %>%
  mutate_all(as.numeric)

# reshape data and parse compound ids
hplcdata_long <- hplcdata_wide %>%
  pivot_longer(cols = -one_of("eid"), names_to = "id", values_to = "rab") %>%
  # parse LIPID MAPS analyte nomenclature
  group_by(id) %>%
  mutate(
    # parse headgroup-level structural data
    struc  = id %>%
      str_split('[ _:/]'),
    # headgroup class
    class =  unlist(struc)[[1]] %>%
      factor(
        levels = c(
          "P-PC", "PC", "LPC",
          "P-PE", "PE", "LPE",
          "PS", "LPS",
          "PI", "PG",
          "Cer", "SM",
          "DG", "TG"
        )
      ),
    # tails
    tails = class2tails(class),
    # total # hydroxy groups
    hydrox = unlist(struc)[[2]] %>%
      str_extract("[nmdt]") %>%
      c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
    # total carbons
    carbon = unlist(struc)[[2]] %>%
      str_remove("[nmdt]") %>%
      as.integer(),
    # total double bonds
    dbonds  = unlist(struc)[[3]] %>%
      as.integer(),
    rt      = unlist(struc)[[4]] %>%
      as.numeric(),
    # fatty acid-level structural data
    hydrsn1 = ifelse(
      "|" %in% unlist(struc),
      unlist(struc)[[6]] %>%
        str_extract("[nmdt]") %>%
        c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
      NA
      ),
    carbsn1 = ifelse(
      "|" %in% unlist(struc),
      unlist(struc)[[6]] %>%
        str_remove("[nmdt]") %>%
        as.numeric(),
      NA
      ),
    dbonsn1 = ifelse(
      "|" %in% unlist(struc),
      unlist(struc)[[7]] %>%
        as.integer(),
      NA
      ),
    hydrsn2 = ifelse(
      "|" %in% unlist(struc),
      unlist(struc)[[8]] %>%
        str_extract("[nmdt]") %>%
        c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
      NA
    ),
    carbsn2 = ifelse(
      "|" %in% unlist(struc),
      unlist(struc)[[8]] %>%
        str_remove("[nmdt]") %>%
        as.numeric(),
      NA
    ),
    dbonsn2 = ifelse(
      "|" %in% unlist(struc),
      unlist(struc)[[9]] %>%
        as.integer(),
      NA
    ),
    hydrsn3 = ifelse(
      ("|" %in% unlist(struc)) & (tails == 3),
      unlist(struc)[[10]] %>%
        str_extract("[nmdt]") %>%
        c('n'=0L, 'm'=1L, 'd'=2L, 't'=3L)[.],
      NA
    ),
    carbsn3 = ifelse(
      ("|" %in% unlist(struc)) & (tails == 3),
      unlist(struc)[[10]] %>%
        str_remove("[nmdt]") %>%
        as.numeric(),
      NA
    ),
    dbonsn3 = ifelse(
      ("|" %in% unlist(struc)) & (tails == 3),
      unlist(struc)[[11]] %>%
        as.integer(),
      NA
    ),
    # sn1/2 positions resolved?
    stereo = str_detect(id, '/') %>% unlist(),
    # make an annotation string that is as specific as possible
    annot = ifelse(
      '|' %in% unlist(struc),
      # if there is acyl chain info
      id %>%
        str_split(' \\| ') %>%
        unlist() %>%
        .[[2]],
      # if there's no acyl chain info
      paste(unlist(struc)[[2]], unlist(struc)[[3]], sep=':')
    )
  ) %>%
  # correct for DAG, TAG stds being 10X of PLs
  group_by(id) %>%
  mutate(rab = ifelse(class %in% c("DG", "TG"), 10*rab, rab)) %>%
  # normalize!
  group_by(eid) %>%
  mutate(frac_molar = rab/sum(rab)) %>%
  group_by(id)

## load and join metadata
# read VARS-format files
envi_data <- list.files(path = dir_envi, full.names = T, pattern = "20201025.*.tsv$") %>%
  lapply(
    function(file_data){
      data_this_file <- file_data %>%
        read_tsv() %>%
        mutate(file = file_data %>% basename())
      return(data_this_file)
    }
  ) %>%
  do.call(rbind, .) %>%
  dplyr::rename(
    sid = samplenum,
    spfull = sp
  ) %>%
  mutate(sid = str_to_upper(sid)) %>%
  right_join(read_tsv(file_meta), by="sid") %>%
  group_by(eid, sp, sid) %>%
  summarize(
    depth_col = first(depth),
    temp_col = first(temp),
    do_col = first(o2_ppt),
    tissue = first(tissue)
  )

# export environmental data
envi_data %>% write_tsv("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/01-rawdata/20220502_metadata.tsv")

# join environmental data
hplcdata <- hplcdata_long %>%
  left_join(envi_data, by = "eid") %>%
  # make the an ordered factor running hi->lo temp, then shallow->deep
  mutate(sp = factor(sp, levels = c("Boli_vitr", "Boli_infu", "Boli_arct", "Lamp_crue", "Tjal_pink"))) %>%
  arrange(sp, depth_col, -temp_col) %>%
  mutate(sp_eid = factor(paste(sp, eid), levels = unique(paste(sp, eid))))

# remove acylglycerols and renorm
hplcdata_membrane <- hplcdata %>%
  group_by(class) %>%
  filter(!(class %in% c("DG", "TG"))) %>%
  group_by(eid) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar))

# exploratory plotting

# lipid class colors
chroma_cl <- c(
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

# species colors
chroma_sp <- RColorBrewer::brewer.pal(9, "Set1") %>%
  setNames(c(
    "Lamp_crue",
    "Boli_arct",
    "Boli_infu",
    "Boli_vitr",
    "Leuc_pulc",
    "Bath_fost",
    "Pleu_bach",
    "Tjal_pink",
    "other"))

# break generator for acyl chain lengths
even_breaks <- function(x, n=5){
  x = sort(round(x))
  if((x[[1]] %% 2) > 0){x[[1]] <- x[[1]]+1}
  spacing <- 2
  breaks <- seq(min(x), max(x), spacing)
  while(length(breaks) > n){
    spacing <- spacing + 2
    breaks <- seq(min(x), max(x), spacing)
  }
  breaks
}

## PCs and PEs plotted back-to-back (panel B)
hplcdata_membrane %>%
  # order the compounds by carbon count
  ungroup() %>%
  mutate(
    id = factor(
      id,
      levels = hplcdata %>%
        arrange(carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  ) %>%
  group_by(sp_eid, sp, class) %>%
  # for PC/PE ratio
  filter(class %in% c("PC", "P-PC", "LPC", "PE", "P-PE", "LPE")) %>%
  mutate(frac_molar = ifelse(str_detect(class, "PC"), frac_molar, -1*frac_molar)) %>% 
  summarize(frac_molar = sum(frac_molar)) %>%
  # average indls of each species
  # comment to show individuals
  group_by(sp, class) %>%
  summarize(
    sem = sd(frac_molar)/sqrt(n()),
    frac_molar = mean(frac_molar),
    sp_eid = first(sp)
  ) %>%
  ggplot(
    aes(
      x = sp_eid,
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
  scale_y_continuous(limits=c(-.5,.5), labels=c(0,.25,.5,.75,1)) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "sample", y = "mole fraction", fill = "lipid class", title = "Total PCs vs. PEs")

ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220323_PCPE.pdf", width=6, height=3)

## Membrane HG composition (NOT USED as of 20220321)
hplcdata_membrane %>%
  # filter classes
  group_by(sp, sp_eid, eid) %>% 
  filter(!(class %in% c("Cer", "SM"))) %>% 
  # renorm
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # order the compounds by carbon count
  ungroup() %>%
  mutate(
    id = factor(
      id,
      levels = hplcdata %>%
        arrange(carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  ) %>%
  group_by(sp_eid, sp, class) %>%
  summarize(frac_molar = sum(frac_molar)) %>%
  # average indls of each species
  # comment to show individuals
  group_by(sp, class) %>%
  summarize(
    sem = sd(frac_molar)/sqrt(n()),
    frac_molar = mean(frac_molar),
    sp_eid = first(sp)
  ) %>%
  ggplot(
    aes(
      x = sp_eid,
      y = frac_molar,
      fill = class
    )
  ) +
  geom_col(
    width = 0.95,
    size = 0.05,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "Species", y = "Mole fraction", fill = "Lipid class", title = "Membrane lipid headgroup composition")

ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220321_membraneHGs_noSL.pdf", width=6, height=4)

## PL curvature index (panel C)
icurv = read_tsv("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/c0_dymond_plus.tsv") %>% 
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

# let's see it! (for SUPP)
icurv %>% 
  filter(class != "DG") %>% 
  drop_na() %>% 
  arrange(c0) %>% 
  mutate(
    class = class %>% 
      as.character() %>% 
      factor(., levels = .)
  ) %>% 
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
    y = "c0"
  )

ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220323_c0_estimates.pdf", width=6, height=4)

## PLCI contributions plot (panel C)
hplcdata_membrane %>%
  group_by(sp, sp_eid, class) %>%
  summarize(frac_molar = sum(frac_molar)) %>%
  group_by(sp, class) %>% 
  summarize(frac_molar = mean(frac_molar)) %>%
  left_join(icurv, by = "class") %>% 
  mutate(plci = frac_molar * c0) %>% 
  drop_na() %>% 
  ## pare down legend, keeping order
  ## doesn't actually affect scale_manual() :(
  #arrange(class) %>% 
  #mutate(
  #  class = class %>% 
  #    as.character() %>% 
  #    factor(., levels = unique(.))
  #) %>% 
  mutate(neg = plci<0) %>% 
  # wrangle for geom_tile
  group_by(sp) %>% 
  arrange(neg, -c0) %>% 
  mutate(
    # running total
    end = cumsum(plci),
    y = end - plci/2
  ) %>% #group_by(sp) %>% group_split() %>% .[[1]] %>% View
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
    data = hplcdata_membrane %>%
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
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(color = "none") +
  labs(x = "Species", y = "PL curvature index", fill = "Lipid class", title = "Contributions to PL curvature index")

ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220323_plci_contribs.pdf", width=6, height=3)

hplcdata_membrane %>%
  group_by(sp_eid, depth_col, sp, class) %>%
  # sum lipids classwise
  summarize(frac_molar = sum(frac_molar)) %>%
  left_join(icurv, by = "class") %>% 
  mutate(
    plci = frac_molar * c0,
    plci_hi = plci + frac_molar * hi,
    plci_lo = plci - frac_molar * lo
  ) %>% 
  # sum contributions from classes
  group_by(sp_eid, sp, depth_col) %>% 
  # IGNORING unknown curvatures, which happen to be for minor components
  summarize(across(contains("plci"), function(x){sum(x, na.rm = TRUE)})) %>% 
  ggplot(
    aes(
      x = depth_col,
      y = plci,
      ymin = plci_lo,
      ymax = plci_hi,
      color = sp
    )
  ) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  scale_color_manual(values = chroma_sp) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  scale_x_log10() +
  labs(x = "depth", y = "PL curvature index", color = "species", title = "PL curvature index vs. depth")
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220320_plci_indl_scatter.pdf", width=6, height=4)

# PLCI scatterplots (NOT USED as of 20220324)
hplcdata_membrane %>%
  group_by(sp_eid, depth_col, sp, class) %>%
  # sum lipids classwise
  summarize(frac_molar = sum(frac_molar)) %>%
  left_join(icurv, by = "class") %>% 
  mutate(plci = frac_molar * c0) %>% 
  # sum contributions from classes
  group_by(sp_eid, sp, depth_col) %>% 
  # IGNORING unknown curvatures, which happen to be for minor components
  summarize(plci = sum(plci, na.rm = TRUE)) %>% 
  ggplot(
    aes(
      x = sp,
      y = plci,
      color = sp
    )
  ) +
  geom_point() +
  scale_color_manual(values = chroma_sp) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "species", y = "PL curvature index", color = "species", title = "PL curvature index vs. species")
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220320_plci_sp_scatter.pdf", width=6, height=4)

# P-PE plot with error bars (NOT USED as of 20220320)
hplcdata_membrane %>%
  # order the compounds by carbon count
  ungroup() %>%
  mutate(
    id = factor(
      id,
      levels = hplcdata %>%
        arrange(carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  ) %>%
  group_by(sp_eid, sp, class) %>%
  # for PC/PE ratio
  filter(class %in% c("P-PE")) %>%
  summarize(frac_molar = sum(frac_molar)) %>%
  # comment to show individuals
  group_by(sp, class) %>%
  summarize(
    sem = sd(frac_molar)/sqrt(n()),
    frac_molar = mean(frac_molar),
    sp_eid = first(sp)
  ) %>%
  ggplot(
    aes(
      x = sp_eid,
      y = frac_molar,
      fill = class
    )
  ) +
  geom_col(
    width = 0.50,
    size = 0.05,
    color = "white"
  ) +
  geom_errorbar(
    aes(
      ymin = frac_molar - sem,
      ymax = frac_molar + sem
    ),
    width = 0.125#,
    #color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  scale_y_continuous(limits=c(0,0.35)) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "sample", y = "mole fraction", fill = "lipid class", title = "Mole fraction P-PE by species")
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220320_P-PE_spmean.pdf", width=6, height=4)

## Acyl chain figs (for SUPP)
class_subset = c("PC", "P-PC", "PE", "P-PE", "PS")

# dbonsn2 by class and species
hplcdata %>%
  # subset
  filter(class %in% class_subset) %>%
  # get total of each chain length
  group_by(eid, sp, class, dbonsn2) %>%
  replace_na(list("dbonsn2" = 0)) %>% # why needed?
  summarise(frac_molar = sum(frac_molar)) %>%
  # renormalize within each class and sample
  group_by(eid, sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # now average samplewise
  group_by(sp, class, dbonsn2) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # compute weighted median. ifelse() prevents overplotting
  group_by(sp, class) %>% 
  mutate(wtmed = ifelse(row_number() == 1, weighted.median(dbonsn2, frac_molar, type=4), NA)) %>%
  # stretch x-axis to 12 (AFTER computing median, bc adds an observation!!)
  full_join(
    crossing(
      sp = hplcdata$sp %>% unique(),
      class = class_subset,
      dbonsn2 = 12,
      frac_molar = 0
    )
  ) %>% 
  ggplot(aes(x=dbonsn2, y=frac_molar, fill=class)) +
  facet_grid(rows = vars(sp), cols = vars(class)) +
  geom_col(size = 0.05, color = "white") +
  geom_vline(
    aes(xintercept = wtmed),
    color = "black",
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
    title = "Sn-2 double bond distributions",
    x = "Number of double bonds",
    y = "Mole fraction of lipid class"
  )
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220324_dbonsn2.pdf", width=6, height=6)
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220324_dbonsn2_wide.pdf", width=12, height=6)

# dbonsn1 by class and species
hplcdata %>%
  # subset
  filter(class %in% class_subset) %>%
  # get total of each chain length
  group_by(eid, sp, class, dbonsn1) %>%
  replace_na(list("dbonsn1" = 0)) %>% # why needed?
  summarise(frac_molar = sum(frac_molar)) %>%
  # renormalize within each class and sample
  group_by(eid, sp, class) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # now average samplewise
  group_by(sp, class, dbonsn1) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # compute weighted median. ifelse() prevents overplotting
  group_by(sp, class) %>% 
  mutate(wtmed = ifelse(row_number() == 1, weighted.median(dbonsn1, frac_molar, type=4), NA)) %>%
  # stretch x-axis to 12 (AFTER computing median, bc adds an observation!!)
  full_join(
    crossing(
      sp = hplcdata$sp %>% unique(),
      class = class_subset,
      dbonsn1 = 12,
      frac_molar = 0
    )
  ) %>% 
  ggplot(aes(x=dbonsn1, y=frac_molar, fill=class)) +
  facet_grid(rows = vars(sp), cols = vars(class)) +
  geom_col(size = 0.05, color = "white") +
  geom_vline(
    aes(xintercept = wtmed),
    color = "black",
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
    title = "Sn-1 double bond distributions",
    x = "Number of double bonds",
    y = "Mole fraction of lipid class"
  )
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220324_dbonsn1.pdf", width=6, height=6)
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220324_dbonsn1_wide.pdf", width=12, height=6)

# total dbonds by class and species
hplcdata %>%
  # subset
  filter(class %in% class_subset) %>%
  mutate(class = factor(class, levels = class_subset)) %>% 
  # get total of each chain length
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
  mutate(wtmed = ifelse(row_number() == 1, weighted.median(dbonds, frac_molar, type=4), NA)) %>%
  ggplot(aes(x=dbonds, y=frac_molar, fill=class)) +
  facet_grid(rows = vars(sp), cols = vars(class)) +
  geom_col(size = 0.05, color = "white") +
  geom_vline(
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
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220324_dbonds.pdf", width=6, height=6)
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220324_dbonds_wide.pdf", width=11, height=8.5)

# PUFA content (Fig. 2D)
hplcdata %>%
  # subset
  filter(class %in% c(class_subset, c("LPC", "LPE"))) %>%
  select(eid, id, sp, class, frac_molar, dbonsn1, dbonsn2) %>%
  pivot_longer(c(dbonsn1, dbonsn2), names_to = "chain", values_to = "dbonds") %>% 
  group_by(eid, sp) %>%
  replace_na(list("dbonds" = 0)) %>% # why needed?
  # renormalize
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  filter(dbonds == 1) %>% # PUFAs only
  summarise(frac_molar = sum(frac_molar)) %>%
  group_by(sp) %>% 
  summarise(
    mean_frac_molar = mean(frac_molar),
    serr_frac_molar = sd(frac_molar)/sqrt(n())
  ) %>% 
  ggplot(
    aes(
      x = sp,
      y = mean_frac_molar,
      ymin = mean_frac_molar - serr_frac_molar,
      ymax = mean_frac_molar + serr_frac_molar
    )
  ) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "species", y = "Mole fraction", color = "species", title = "MUFA content vs. species")

# Showcase different desat metrics
hplcdata %>%
  # subset
  filter(class %in% c("PC", "P-PC", "PE", "P-PE", "LPC", "LPE")) %>%
  select(eid, id, sp, class, frac_molar, dbonsn1, dbonsn2) %>%
  pivot_longer(c(dbonsn1, dbonsn2), names_to = "chain", values_to = "dbonds") %>% 
  group_by(eid, sp) %>%
  replace_na(list("dbonds" = 0)) %>% # why needed?
  # renormalize
  mutate(
    frac_molar = frac_molar/sum(frac_molar),
    fm_ufa  = ifelse(dbonds > 0, frac_molar, 0),
    fm_mufa = ifelse(dbonds == 1, frac_molar, 0),
    fm_pufa = ifelse(dbonds > 1, frac_molar, 0),
    pc = ifelse(str_detect(class, "PC"), 1, -1)
  ) %>% 
  group_by(sp, eid, pc) %>% 
  summarise(across(contains("fm"), sum)) %>%
  pivot_longer(contains("fm"), names_to = "unsatcl", values_to = "frac_molar") %>% 
  group_by(sp, unsatcl, pc) %>% 
  mutate(
    frac_molar = frac_molar * pc,
    mean_frac_molar = mean(frac_molar),
    mean_frac_molar = ifelse(row_number() == 1, mean_frac_molar, NA),
    serr_frac_molar = sd(frac_molar*pc)/sqrt(n()),
    serr_frac_molar = ifelse(row_number() == 1, serr_frac_molar, NA)
  ) %>% 
  ggplot(
    aes(
      x = sp,
      y = frac_molar,
      ymin = mean_frac_molar - serr_frac_molar,
      ymax = mean_frac_molar + serr_frac_molar#,
      #color = pc
    )
  ) +
  facet_grid(cols = vars(unsatcl)) +
  geom_hline(yintercept = 0) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_errorbar(width = 0.2) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "species", y = "Mole fraction", color = "species", title = "UFA content vs. species")

# PUFA only
hplcdata %>%
  # subset
  filter(class %in% c("PC", "P-PC", "PE", "P-PE", "LPC", "LPE")) %>%
  select(eid, id, sp, class, frac_molar, dbonsn1, dbonsn2) %>%
  pivot_longer(c(dbonsn1, dbonsn2), names_to = "chain", values_to = "dbonds") %>% 
  group_by(eid, sp) %>%
  replace_na(list("dbonds" = 0)) %>% # why needed?
  # renormalize
  mutate(
    frac_molar = frac_molar/sum(frac_molar),
    fm_ufa  = ifelse(dbonds > 0, frac_molar, 0),
    fm_mufa = ifelse(dbonds == 1, frac_molar, 0),
    fm_pufa = ifelse(dbonds > 1, frac_molar, 0),
    pc = ifelse(str_detect(class, "PC"), 1, -1)
  ) %>% 
  group_by(sp, eid, pc) %>% 
  summarise(across(contains("fm"), sum)) %>%
  pivot_longer(contains("fm"), names_to = "unsatcl", values_to = "frac_molar") %>% 
  group_by(sp, unsatcl, pc) %>% 
  mutate(
    frac_molar = frac_molar * pc,
    mean_frac_molar = mean(frac_molar),
    mean_frac_molar = ifelse(row_number() == 1, mean_frac_molar, NA),
    serr_frac_molar = sd(frac_molar*pc)/sqrt(n()),
    serr_frac_molar = ifelse(row_number() == 1, serr_frac_molar, NA)
  ) %>% 
  #filter(unsatcl == "fm_ufa") %>% 
  ggplot(
    aes(
      x = sp,
      y = frac_molar,
      ymin = mean_frac_molar - serr_frac_molar,
      ymax = mean_frac_molar + serr_frac_molar#,
      #color = pc
    )
  ) +
  facet_grid(cols = vars(unsatcl)) +
  geom_hline(yintercept = 0) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_errorbar(width = 0.2) +
  theme_pubr() +
  lims(y = c(-0.52, 0.5)) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "species", y = "Mole fraction", color = "species", title = "UFA content vs. species")
ggsave("/Users/jwinnikoff/Documents/MBARI/deeplipid-biophys-2022/ctenolipids/20220414_ufa_layer.pdf", width=6, height=4)
