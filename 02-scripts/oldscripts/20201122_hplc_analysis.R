library(tidyverse)
library(ggpubr)
library(RColorBrewer)

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
        select(-1, -3) %>%
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
          "PC", "P-PC", "LPC",
          "PE", "P-PE", "LPE",
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
  group_by(eid, sp) %>%
  summarize(
    depth_col = depth[[1]],
    temp_col =  temp[[1]],
    do_col =    o2_ppt[[1]],
    tissue =    tissue[[1]]
  )

# join environmental data
hplcdata <- hplcdata_long %>%
  left_join(envi_data, by = "eid") %>%
  # make the an ordered factor running hi->lo temp, then shallow->deep
  mutate(sp = factor(sp, levels = c("Boli_vitr", "Boli_infu", "Boli_arct", "Lamp_crue", "Tjal_pink"))) %>%
  arrange(sp, depth_col, -temp_col) %>%
  mutate(sp_eid = factor(paste(sp, eid), levels = unique(paste(sp, eid))))

# show where the critters came from
hplcdata %>% group_by(eid, sp, depth_col, temp_col) %>% summarise() %>% arrange(sp, eid) %>% 
  group_by(sp) %>% 
  summarise(
    depth_min = min(depth_col),
    depth_max = max(depth_col),
    press_min = min(depth_col)/10,
    press_max = max(depth_col)/10,
    temp_min  = min(temp_col),
    temp_max  = max(temp_col)
  )

# remove acylglycerols and renorm
hplcdata_membrane <- hplcdata %>%
  group_by(class) %>%
  filter(!(class %in% c("DG", "TG"))) %>%
  group_by(eid) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar))

# exploratory plotting

# WonB theme
theme_pubblack <- function(...){
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
      plot.background = element_blank(),
      plot.title = element_text(color = "white")
    )
}

# lipid class colors
chroma_cl <- c(
  brewer.pal(4, "Blues")[2:4], #PCs
  brewer.pal(4, "Greens")[2:4], #PEs
  brewer.pal(12, "Paired")[5:12] #PSs
) %>%
  setNames(
    c(
      "PC", "P-PC", "LPC",
      "PE", "P-PE", "LPE",
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

# plot classes by sample
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
  # filter a subset of classes
  #filter(class %in% c("LPC", "LPE", "LPS")) %>%
  filter(class == "P-PE") %>%
  # for PC/PE ratio
  #filter(class %in% c("PC", "P-PC", "LPC", "PE", "P-PE", "LPE")) %>%
  #mutate(frac_molar = ifelse(class %in% c("PC", "P-PC", "LPC"), frac_molar, -1*frac_molar)) %>%
  summarize(frac_molar = sum(frac_molar)) %>%
  # to combine classes
  #group_by(sp_eid) %>%
  #mutate(frac_molar = sum(frac_molar)) %>%
  # average indls of each species
  # comment to show individuals
  group_by(sp, class) %>%
  summarize(
    sem = sd(frac_molar)/sqrt(n()),
    frac_molar = mean(frac_molar),
    sp_eid = sp[[1]]
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
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = frac_molar - sem,
      ymax = frac_molar + sem
    ),
    width = 0.125,
    color = "white"
  ) +
  scale_fill_manual(values = chroma_cl) +
  # for displaying ratios
  #scale_y_continuous(limits=c(-.5,.5), labels=c(0,.25,.5,.75,1)) +
  scale_y_continuous(limits=c(0,0.35)) +
  # regular scale
  #scale_y_continuous() +
  theme_pubblack() +
  #theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = FALSE) +
  labs(x = "sample", y = "mole fraction", fill = "lipid class", title = "lipid class by species, polar only")

ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/20201122_hplc/20210316_PPE_spmean.pdf", width=12, height=6)

# now, divide individual compounds and label the major ones
hplcdata_membrane %>%
  # order the compounds by carbon count
  ungroup() %>%
  mutate(
    id = factor(
      id,
      levels = hplcdata %>%
        distinct(id, carbon, dbonds) %>%
        arrange(carbon, dbonds) %>%
        pull(id) %>%
        unique()
    )
  ) %>%
  # average indls of each species
  # comment to show individuals
  group_by(sp, class, id, annot) %>%
  summarize(
    sp_eid = sp[[1]],
    frac_molar = mean(frac_molar)
  ) %>%
  group_by(sp_eid, class, annot) %>%
  ggplot(
    aes(
      x = sp_eid,
      y = frac_molar,
      fill = class
    )
  ) +
  geom_col(width = 0.95, size = 0.05, color = "black") +
  geom_text(
    aes(
      label = ifelse(frac_molar >= 0.015, annot, ''),
      color = class
    ),
    size = 1.5,
    position = position_stack(vjust=0.5)
  ) +
  scale_fill_manual(values = chroma_cl) +
  scale_color_manual(values = c(
    "PC"="black", "P-PC"="black", "LPC"="white",
    "PE"="black", "P-PE"="black", "LPE"="white",
    "PS"="black", "LPS"="white",
    "DG"="black",
    "PI"="black",
    "SM"="black",
    "Cer"="white",
    "PG"="black", "TG"="black"
  )) +
  theme_pubblack() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = FALSE) +
  labs(x = "sample", y = "mole fraction", fill = "lipid class", title = "lipid class by sample")

ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/20201122_hplc/20201130_lipidclass_spmean_cpdsplit.pdf", width=12, height=6)

# chain length by class and species
hplcdata %>%
  group_by(eid, sp, class, carbon) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  group_by(sp, class, carbon) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # no PGs, DAGs, TAGs, SMs
  filter(!(class %in% c("PG", "DG", "TG", "SM"))) %>%
  ggplot(aes(x=carbon, y=frac_molar, fill=sp, alpha=-(carbon %% 2))) +
    facet_grid(rows = vars(sp), cols = vars(class), scales = "free") +
    geom_col() +
    scale_x_continuous(breaks = even_breaks) +
    scale_fill_manual(values = chroma_sp) +
    scale_alpha(range = c(0.4, 2)) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    guides(alpha=FALSE) +
    ggtitle("chain length by class and species")

# dbonds by class and species
hplcdata %>%
  group_by(eid, sp, class, dbonds) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  group_by(sp, class, dbonds) %>%
  summarize(frac_molar = mean(frac_molar)) %>%
  # no PGs, DAGs, TAGs, SMs
  filter(!(class %in% c("PG", "DG", "TG", "SM"))) %>%
  ggplot(aes(x=dbonds, y=frac_molar, fill=sp)) +
  facet_grid(rows = vars(sp), cols = vars(class), scales = "free") +
  geom_col() +
  scale_x_continuous(breaks = even_breaks) +
  scale_fill_manual(values = chroma_sp) +
  theme_pubr() +
  ggtitle("double bonds by class and species")

class_subset <- c("DG", "TG")

## as above, but normalized by class
# chain length by class and species
hplcdata %>%
  # subset
  #group_by(class) %>%
  #filter(class %in% class_subset) %>%
  # get total of each chain length
  group_by(eid, sp, class, carbon) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  group_by(sp, class, carbon) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # normalize
  group_by(sp, class) %>%
  mutate(prop_molar = frac_molar/max(frac_molar)) %>%
  ggplot(aes(x=carbon, y=prop_molar, fill=class, alpha=-(carbon %% 2))) +
  facet_grid(rows = vars(sp), cols = vars(class), scales = "free") +
  geom_col(size = 0.05, color = "black") +
  geom_vline(
    data = hplcdata %>%
      # subset
      group_by(class) %>%
      #filter(class %in% class_subset) %>%
      group_by(eid, sp, class, carbon) %>%
      summarise(frac_molar = sum(frac_molar)) %>%
      group_by(sp, class, carbon) %>%
      summarise(frac_molar = mean(frac_molar)) %>%
      # normalize
      group_by(sp, class) %>%
      mutate(frac_molar = frac_molar/sum(frac_molar)) %>%
      summarise(
        # weighted mean
        #carbon = sum(frac_molar*carbon, na.rm=TRUE)
        carbon = spatstat::weighted.median(carbon, frac_molar)
      ),
    aes(xintercept = carbon),
    color = "white",
    alpha = 0.8,
    linetype="dotted"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.6, 2)) +
  theme_pubblack() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha=FALSE, fill=guide_legend(nrow=1)) +
  ggtitle("chain length by class and species")

ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/20201122_hplc/20201130_chainlen_medlines.pdf", width=16, height=8.5)

# the above, with compounds split and labeled
# chain length by class and species
hplcdata %>%
  # subset
  group_by(class) %>%
  filter(class %in% class_subset) %>%
  group_by(sp, class, carbon, id, annot) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # normalize
  group_by(sp, class, carbon) %>%
  mutate(tot_frac_molar = sum(frac_molar)) %>%
  group_by(sp, class) %>%
  mutate(prop_molar = frac_molar/max(tot_frac_molar)) %>%
  ggplot(aes(x=carbon, y=prop_molar, fill=class, alpha=-(carbon %% 2))) +
  facet_grid(rows = vars(sp), cols = vars(class), scales = "free") +
  geom_col(size = 0.05, color = "black") +
  # compound labels
  geom_text(
    aes(label = ifelse(prop_molar >= 0.25, annot, '')),
    color = "white",
    size = 1.5,
    position = position_stack(vjust=0.5),
    angle = 90,
    alpha = 1
  ) +
  # median lines
  geom_vline(
    data = hplcdata %>%
      # subset
      group_by(class) %>%
      filter(class %in% class_subset) %>%
      group_by(eid, sp, class, carbon) %>%
      summarise(frac_molar = sum(frac_molar)) %>%
      group_by(sp, class, carbon) %>%
      summarise(frac_molar = mean(frac_molar)) %>%
      # normalize
      group_by(sp, class) %>%
      mutate(frac_molar = frac_molar/sum(frac_molar)) %>%
      summarise(
        # weighted mean
        #carbon = sum(frac_molar*carbon, na.rm=TRUE)
        carbon = spatstat::weighted.median(carbon, frac_molar)
      ),
    aes(xintercept = carbon),
    color = "white",
    alpha = 0.8,
    linetype="dotted"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.4, 2)) +
  theme_pubblack() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha=FALSE) +
  ggtitle("chain length by class and species")

ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/20201122_hplc/20201129_chainlen_CerSM_cpdsplit_classcolor.pdf", width=12, height=8.5)

class_subset <- c("PC","PE","PS")

# dbonds by class and species
hplcdata %>%
  # subset
  #group_by(class) %>%
  filter(class %in% class_subset) %>%
  # get total of each chain length
  group_by(eid, sp, class, dbonds) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  group_by(sp, class, dbonds) %>%
  mutate(dbonds = ifelse(is.na(dbonds), 0 , dbonds)) %>% 
  summarise(frac_molar = mean(frac_molar)) %>%
  # normalize
  group_by(sp, class) %>%
  mutate(prop_molar = frac_molar/sum(frac_molar)) %>%
  ggplot(aes(x=dbonds, y=prop_molar, fill=class)) +
  facet_grid(rows = vars(sp), cols = vars(class), scales = "free_y") +
  geom_col(size = 0.05, color = "black") +
  geom_vline(
    data = hplcdata %>%
      # subset
      group_by(class) %>%
      filter(class %in% class_subset) %>%
      group_by(eid, sp, class, dbonds) %>%
      summarise(frac_molar = sum(frac_molar)) %>%
      group_by(sp, class, dbonds) %>%
      summarise(frac_molar = mean(frac_molar)) %>%
      # normalize
      group_by(sp, class) %>%
      mutate(frac_molar = frac_molar/sum(frac_molar)) %>%
      summarise(
        # weighted mean
        #carbon = sum(frac_molar*carbon, na.rm=TRUE)
        dbonds = spatstat::weighted.median(dbonds, frac_molar)
      ),
    aes(xintercept = dbonds),
    color = "white",
    alpha = 0.8,
    linetype="dotted"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.4, 2)) +
  theme_pubblack() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha=FALSE, fill=guide_legend(nrow=1)) +
  ggtitle("double bonds by class and species")

ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/20201122_hplc/20210726_dbonds_medlines.pdf", width=12, height=8.5)

class_subset <- c("PC","PE","PS")

# dbonsn2 by class and species
hplcdata %>%
  # subset
  #group_by(class) %>%
  filter(class %in% class_subset) %>%
  # get total of each chain length
  group_by(eid, sp, class, dbonsn2) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  mutate(dbonsn2 = ifelse(is.na(dbonsn2), 0 , dbonsn2)) %>% 
  group_by(sp, class, dbonsn2) %>%
  summarise(frac_molar = mean(frac_molar)) %>%
  # normalize
  group_by(sp, class) %>%
  mutate(prop_molar = frac_molar/sum(frac_molar)) %>%
  ggplot(aes(x=dbonsn2, y=prop_molar, fill=class)) +
  facet_grid(rows = vars(sp), cols = vars(class), scales = "free_y") +
  geom_col(size = 0.05, color = "black") +
  geom_vline(
    data = hplcdata %>%
      # subset
      group_by(class) %>%
      filter(class %in% class_subset) %>%
      group_by(eid, sp, class, dbonsn2) %>%
      summarise(frac_molar = sum(frac_molar)) %>%
      group_by(sp, class, dbonsn2) %>%
      summarise(frac_molar = mean(frac_molar)) %>%
      # normalize
      group_by(sp, class) %>%
      mutate(frac_molar = frac_molar/sum(frac_molar)) %>%
      summarise(
        # weighted mean
        #carbon = sum(frac_molar*carbon, na.rm=TRUE)
        dbonsn2 = spatstat::weighted.median(dbonsn2, frac_molar)
      ),
    aes(xintercept = dbonsn2),
    color = "white",
    alpha = 0.8,
    linetype="dotted"
  ) +
  scale_x_continuous(breaks = even_breaks) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = chroma_cl) +
  scale_alpha(range = c(0.4, 2)) +
  theme_pubblack() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(alpha=FALSE, fill=guide_legend(nrow=1)) +
  ggtitle("double bonds by class and species, sn2 only")
ggsave("/Users/jwinnikoff/Documents/MBARI/Lipids/20201122_hplc/20210723_dbonsn2_medlines.pdf", width=12, height=8.5)

# what about ceramide hydroxylation?
hplcdata %>%
  filter(class=="Cer") %>%
  group_by(id, hydrox, sp) %>%
  summarize(frac_molar = mean(frac_molar)) %>%
  group_by(hydrox, sp) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = hydrox, y = frac_molar, fill=sp)) +
  facet_wrap(~sp, ncol = 1) + geom_col() +
  theme_pubr() +
  scale_fill_manual(values = chroma_sp) +
  ggtitle("ceramide acyl hydroxylation by species")

# what are the most abundant chemical species?
hplcdata %>%
  group_by(eid) %>%
  arrange(desc(frac_molar)) %>%
  filter(frac_molar >= frac_molar[[5]]) %>%
  arrange(eid, -frac_molar) %>% View()

## evaluate the FAME conclusions
# plot average # carbons
hplcdata %>%
  # get average # carbons
  group_by(eid, sp, temp_col, depth_col) %>%
  summarise(carbon = sum(carbon * frac_molar)) %>%
  mutate(shal = ifelse(depth_col <= 200, 1.0, 0.4)) %>%
  ggplot(aes(x = temp_col, y = carbon, color = sp, alpha = shal)) +
  geom_point() +
  scale_color_manual(values = chroma_sp) +
  scale_alpha(range=c(0.4, 1)) +
  theme_pubr() +
  theme(legend.position = c(0.6, 0.7)) +
  guides(alpha=FALSE) +
  labs(
    x = "collection temperature (˚C)",
    y = "mean chain length",
    color = "species",
    title = "Chain length vs. temperature"
  )
# this is not what the FAME data had!

# plot average # double bonds
hplcdata %>%
  group_by(eid, sp, temp_col, depth_col) %>%
  summarise(dbonds = sum(dbonds * frac_molar)) %>%
  mutate(cold = ifelse(temp_col <= 7.5, 1, 0.4)) %>%
  ggplot(aes(x = depth_col, y = dbonds, color = sp, alpha = cold)) +
  geom_point() +
  scale_color_manual(values = chroma_sp) +
  scale_alpha(range=c(0.4, 1)) +
  theme_pubr() +
  theme(legend.position = c(0.6, 0.4)) +
  guides(alpha=FALSE) +
  labs(
    x = "collection depth (m)",
    y = "mean # double bonds (DBI)",
    color = "species",
    title = "DBI vs. depth"
  )

## which fatty acids are driving the trends?
# extract fatty acid data in long format: 1 entry/FA
# NTS 20201124: Itay says not important
fadata <- bind_rows(
  mutate(across(where(str_detect("sn1"))))
)

# did the injected masses affect sample complexity?
masses <- read_tsv(paste(dir_data, "lipidmaps_sample_masses.tsv", sep='/'))

complexity_v_mass <- hplcdata %>% 
  filter(frac_molar > 0) %>% 
  group_by(eid, sp) %>% 
  summarize(n_cpds = length(unique(id))) %>% 
  left_join(masses, by = "eid")

complexity_v_mass %>% 
  ggplot(aes(x = mass, y = n_cpds, color = sp)) +
  geom_point() +
  theme_pubr() +
  labs(x = "sample mass (µg)", y = "# compounds detected")

# 20220121 PE-O

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
  group_by(sp_eid, sp, class, depth_col, temp_col) %>%
  filter(sp %in% c("Boli_arct", "Lamp_crue", "Tjal_pink")) %>% # deep only
  # filter a subset of classes
  filter(class == "P-PE") %>%
  summarize(frac_molar = sum(frac_molar)) %>%
  ggplot(
    aes(
      x = depth_col,
      y = frac_molar,
      color = sp
    )
  ) +
  geom_smooth(
    method="lm",
    color="black"
  ) +
  geom_point(
    size=3
  ) +
  scale_color_manual(values = chroma_sp) +
  # regular scale
  scale_x_log10() +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  #guides(color = FALSE) +
  labs(x = "Depth (m)", y = "mole fraction", fill = "lipid class", title = "P-PE vs. depth")

data_annotated = data_long %>% 
  mutate(
    sn1_len = name2sn1l(cmpd),
    sn1_dbl = name2sn1d(cmpd),
    
  )

## 20220318 PL curvature index
icurv = tibble

