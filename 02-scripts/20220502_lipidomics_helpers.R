# HELPER FUNCTIONS for lipidomics analyses and plotting

library(here)
library(tidyverse)
library(phytools)
library(ggpubr)
library(RColorBrewer)

## PARSING HELPERS

serr = function(x){sd(x)/sqrt(length(x))}

# How many acyl chains does passed class of lipids have?
class2tails = Vectorize(function(class){
  if(str_detect(class, 'L|AC|CE')){
    return(1L)
  }
  if(str_detect(class, 'T')){
    return(3L)
  }else{
    return(2L)
  }
})

# Take a long-format tibble with an `id` column containing lipid maps compound names,
# break that column out into several informative columns
#NTS 20220502: May need to be updated for changes in nomenclature between 2020, 2022.
parse_lipidmaps_id = function(longdf){
  longdf %>% 
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
            "Cer", "Cer1P", "SM",
            "DG", "TG",
            "AC", "CE"
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
    )
}

## PLOTTING HELPERS

# full 14-headgroup color map
chroma_cl = c(
  brewer.pal(4, "Blues")[2:4], #PCs
  brewer.pal(4, "Greens")[2:4], #PEs
  brewer.pal(12, "Paired")[5:12], #PSs,
  "grey50", "grey25"
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
      "PG", "TG",
      "AC", "CE"
    )
  )

# a special little palette for when we want ester PE to be grey
chroma_ppe = chroma_cl
chroma_ppe[["PE"]] = "#BBBDBF"

# colors for labeling compounds on a bar plot
chroma_tx = c(
  "PC"="black", "P-PC"="black", "LPC"="white",
  "PE"="black", "P-PE"="black", "LPE"="white",
  "PS"="black", "LPS"="white",
  "DG"="black",
  "PI"="black",
  "SM"="black",
  "Cer"="white",
  "PG"="black", "TG"="black",
  "AC"="white", "CE"="white"
)

# species colors
# 20220511 new ones for the species in this paper
chroma_sp = c(
  "Leuc_pulc" = "#FF7F00", # orange
  "Bero_cucu" = "#3A7DB4",
  "Boli_vitr" = "#984EA3", # the original purple
  "Boli_infu" = "#109292", # muted teal used in 2021 paper
  "Lamp_crue" = "#E41A1C",
  "Bath_fost" = "#A65628", # brick red
  "Llyr_deep" = "#FCC06F", # peach gold
  "Llyr_bent" = "#FCC06F", 
  "Cydi_blac" = "#000000", # straight up black
  "Tjal_pink" = "#F781BF",
  #"other"     = "#999999"
  "Mert_angu" = "#8A8B8A" # 50% grey bc related to black cyd
)

spp_figS2E = c(
  "Mert_angu",
  "Cydi_blac",
  "Boli_arct",
  "Lamp_crue",
  "Bath_fost"#,
  #  "Tjal_pink"
)

##NTS 20220502: I should amend this to match final figs in the FAME paper
#chroma_sp = RColorBrewer::brewer.pal(9, "Set1") %>%
#  setNames(c(
#    "Lamp_crue",
#    "Bero_cucu",
#    "Boli_infu",
#    "Mert_angu",
#    "Llyr_deep",
#    "other",
#    "Bath_fost",
#    "Tjal_pink",
#    "Cydi_blac")) %>% 
#  c(., c("Llyr_bent" = "#FF7F00"))

## these are from the FAME dataset
#chroma_sp = RColorBrewer::brewer.pal(9, "Set1") %>%
#  setNames(c(
#    "Lamp_crue",
#    "Boli_arct",
#    "Boli_infu",
#    "Boli_vitr",
#    "Leuc_pulc",
#    "Bath_fost",
#    "Pleu_bach",
#    "Tjal_pink",
#    "other"))

# break generator for acyl chain lengths
even_breaks = function(x, n=5){
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

## PGLS helpers
# uniting the sovereign realms of phytools and tidyverse
# refactored, tidied-up PGLS routines from the JEB paper

# match phenotypic data to a phylogeny
# "samp" and "sp" are the default column names mapping sample to species
# return namedlist containing tibble and phylo
# tip2sp is the function that converts tip names to species
match_tree = function(pheno, phylo, tip2sp=identity){
  
  # species in the pheno dataset
  sp = unique(pheno$sp)
  # warn about any species not in the tree
  sp_missing = sp[which(!(sp %in% tip2sp(phylo$tip.label)))]
  if(length(sp_missing) > 0){warning(paste(paste(sp_missing, collapse=", "), "not found in tree!", sep=' '))}
  
  # remove those from the table
  pheno_ok = pheno %>% filter(sp %in% phylo$tip.label)
  
  # prune taxa not in pheno
  phylo_pruned = phylo %>% keep.tip(phylo$tip.label[which(phylo$tip.label %in% pheno_ok$sp)])
  
  # add conspecifics as polytomies, use for loop because node matching recurses
  samps = pheno_ok %>% rowwise() %>% group_split()
  phylo_indl = phylo_pruned
  for(row in samps){
    phylo_indl = phylo_indl %>%
      bind.tip(tip.label = row$sp, where = match(row$sp, phylo_indl$tip.label), edge.length = 0)
  }
  
  # drop leaves that do not match samps verbatim
  phylo_ok = phylo_indl %>% 
    # keep.tip automatically drops redundantly named leaves
    keep.tip(.$tip.label[which(.$tip.label %in% pheno_ok$sp)])
  
  return(list("pheno" = pheno_ok, "phylo" = phylo_ok))
}

# this function takes a phylo object and a table of phenotypic data with column "sp"
# It consolidates duplicate tip labels, then makes polytomies to match the number of indl's/sp in pheno.
expand_tree = function(pheno, phylo, tip2sp=function(x){str_remove_all(x, "[0-9]")}){
  # construct branch length table
  branches = cbind(phylo$edge, phylo$edge.length) %>%
    as_tibble() %>%
    set_names(c("mrca", "node", "dist")) %>%
    mutate(
      node = as.integer(node),
      mrca = as.integer(mrca)
    )
  
  # get duplicated taxa and their intraspecific MRCAs
  leaves = tibble(sp = phylo$tip.label) %>%
    mutate(node = row_number())
  
  leaves_dup = leaves %>%
    group_by(sp) %>%
    mutate(count = n()) %>%
    filter(count > 1) %>%
    mutate(mrca = phylo %>% getMRCA(node) %>% as.integer())
  
  # duplicate taxa and their average distance from MRCA
  taxa_dup = leaves_dup %>%
    rowwise() %>%
    mutate(dist = tibble(node, mrca) %>% left_join(branches, by = c("node", "mrca")) %>% .$dist) %>%
    group_by(sp, mrca) %>%
    quietly(summarize)(
      node = min(node),
      dist = mean(dist, na.rm = TRUE)
    ) %>% .$result %>% 
    ungroup()
  
  # leaves to drop (all but the first, per sp)
  leaves_drop = leaves_dup %>%
    anti_join(leaves_dup %>% summarize(node = first(node)), by = c("sp", "node"))
  
  # set first conspecific branch to average length
  branches_adj = bind_rows(anti_join(branches, taxa_dup, by = c("node", "mrca")), taxa_dup)
  
  # replace conspecific clades with a single tip
  phylo_pruned = phylo %>%
    compute.brlen(
      phylo$edge %>%
        as_tibble() %>%
        set_names(c("mrca", "node")) %>%
        # put in order
        left_join(branches_adj, by = c("mrca", "node")) %>%
        pull(dist)
    ) %>%
    drop.tip(leaves_drop %>% pull(node)) %>%
    # finally, prune taxa not in phenodata
    drop.tip(phylo$tip.label %>% .[which(!(. %in% pheno$sp))])
  
  # Supplement pheno table for output
  pheno_indl = pheno %>% 
    group_by(sp) %>% 
    mutate(indl = paste(sp, row_number(), sep=''))
  
  # Start generating a list of the to-be tip names (indls)
  leaves_add = pheno %>% 
    summarize() %>% 
    group_by(sp) %>% 
    summarise(indl = paste(sp, seq(n()), sep='') %>% list()) %>% 
    unnest(indl) %>% 
    group_by(sp) %>% 
    arrange(indl) %>% 
    filter(row_number() > 1)
  
  phylo_indl = phylo_pruned
  # concat "1" to all the existing tip labels
  phylo_indl$tip.label = paste(phylo_indl$tip.label, '1', sep='')
  
  # iteratively add leaves to tree
  edge_polytomy = 0
  for(row in leaves_add %>% rowwise() %>% group_split()){
    phylo_indl = phylo_indl %>%
      bind.tip(tip.label = row$indl, where = match(row$sp, tip2sp(phylo_indl$tip.label)), edge.length = edge_polytomy)
  }
  # extra char on the end of the tip names; just cut it!
  phylo_indl$tip.label = str_sub(phylo_indl$tip.label, end=-2)
  
  # number the individuals in the tree_rooted and the table
  phylo_indl$tip.label = tibble(sp = phylo_indl$tip.label) %>%
    group_by(sp) %>%
    mutate(indl = paste(sp, row_number(), sep="")) %>%
    pull(indl)
  
  return(list("pheno" = pheno_indl, "phylo" = phylo_indl))
}

# A PGLS wrapper that takes a tibble with sp column, so you don't have to worry that the data table and
# tree tips are in the same order. You do, however, have to have a column in the tibble that corresponds
# 1:1 with tips in the tree.
pgls_tbl = function(pheno, phylo, form, startval = 0, gamma=1, tipcol="indl", corfunc=ape::corMartins, ...){
  
  # make the individual names rownames in the data
  pheno_rownames = pheno %>%
    column_to_rownames(var = tipcol)
  
  # rescale branches; this matters for OU models
  phylo$edge.length = phylo$edge.length * gamma
  
  # do.call is needed to preserve the formula call in the object
  fit = do.call(
    nlme::gls,
    list(
      data = pheno_rownames,
      model = as.formula(form),
      correlation = corfunc(startval, phy = phylo, ...)
    )
  )
  
  return(fit)
}

# helper function to properly decode the tTable
tTable_rename = function(ttablist){
  if(length(ttablist)){
    tt = as.data.frame(ttablist) %>% 
      t() %>% 
      as_tibble() %>% 
      dplyr::rename(name = V1) %>% 
      mutate(col = (row_number()-1)%%4 + 1) %>% 
      group_by(col) %>% 
      mutate(row = row_number()) %>% 
      left_join(
        tibble(
          col = seq(4),
          par = c("val", "se", "tval", "pval")
        ),
        by = "col"
      ) %>% 
      mutate(name = paste(par, "_var", row, sep = ''))
    return(ttablist %>% setNames(tt$name))
  }else{
    return(NULL)
  }
}

# Take a gls object, return namedlist of its pertinent contents
unpack_gls = function(gls){
  summlist = gls %>% 
    summary() %>% 
    unlist()
  c(
    summlist %>% 
      .[!str_detect(names(.), "fitted|residuals|call|tTable")],
    summlist %>% 
      .[str_detect(names(.), "tTable")] %>% 
      tTable_rename()
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
      plot.background = element_blank(),
      plot.title = element_text(color = "white")
    )
}

# the original 3x5 dataset
eids_orig = c(
  "JWL0021",
  "JWL0025",
  "JWL0028",
  "JWL0034",
  "JWL0039",
  "JWL0042",
  "JWL0043",
  "JWL0044",
  "JWL0058",
  "JWL0062",
  "JWL0070",
  "JWL0099",
  "JWL0100",
  "JWL0131",
  "JWL0148"
)

major_pls = c("PC", "P-PC", "PE", "P-PE", "PS")
