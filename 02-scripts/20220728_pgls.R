# PHYLOGENETIC AND ORDINARY REGRESSIONS OF LIPIDOME FEATURES AGAINST DEPTH

library(here)
library(tidyverse)
library(scales)
library(ggpubr)
library(RColorBrewer)
library(ape)
library(AICcmodavg) # provides predictSE.gls

source(here("02-scripts", "20220502_lipidomics_helpers.R"))
# Load up all the data for the rest of Fig. 2
# this should furnish lcmsdata_pl, icurv, curvdata
source(here("02-scripts", "20220609_lipidomics.R"))

# load the phylogeny
file_tree = here("03-mainfigs", "ctenolipids", "20220506_iq.tree")
phylodata = ape::read.tree(file_tree)
plot(phylodata) # have a quick look

# prep the phenotypic data in long form
# this can handle multiple variables
# moving forward, there are just "pred" and "resp" variables
# all else is handled by dplyr grouping
phenodata = bind_rows(
    # composition data
    lcmsdata_pl_wholectenos %>% 
      filter(class %in% major_pls) %>% 
      group_by(eid, sp, class, depth_col, temp_col) %>% 
      summarise(frac_molar = sum(frac_molar)) %>% 
      pivot_longer(frac_molar, names_to = "resp_var", values_to = "resp"),
    # total plasmalogen
    lcmsdata_pl_wholectenos %>% 
      filter(class %in% c("P-PE", "P-PC")) %>% 
      group_by(eid, sp, depth_col, temp_col) %>% 
      summarise(frac_molar = sum(frac_molar)) %>% 
      mutate(class = "plasm") %>% 
      pivot_longer(frac_molar, names_to = "resp_var", values_to = "resp"),
    # curvature data
    curvdata %>% 
      group_by(eid, sp, depth_col, temp_col) %>% 
      summarise(plci = sum(plci)) %>% 
      pivot_longer(plci, names_to = "resp_var", values_to = "resp")
  ) %>% 
  group_by(class, resp_var) %>% 
  mutate(resp_var = paste(resp_var, c(na.omit(as.character(class)), '', sep='')[[1]])) %>%
  # extra predictors could be added here
  pivot_longer(depth_col, names_to = "pred_var", values_to = "pred") %>% 
  group_by(sp) %>% 
  filter(mean(temp_col) <= cold) %>% # subset!
  #mutate(pred = log10(pred)) %>%  # transform
  # since we don't have an Arctic Boli txome
  mutate(sp = ifelse(sp == "Boli_arct", "Boli_infu", as.character(sp)) %>% factor(levels = levels(sp)))
  
# Reduce the phylo- and phenodata to their taxonomic intersection
# give warnings if species in phenodata are not found in phylogeny
phenophylo = match_tree(phenodata, phylodata, tip2sp = function(x){str_remove_all(x, "[0-9]")})
pheno = phenophylo$pheno
phylo = phenophylo$phylo

# root it?
# this seems to be the conservative choice - LBA affects PGLS!
message("midpoint-rooting phylogeny")
phylo = phylo %>% midpoint.root()

# ultrametrize?
message("making phylogeny ultrametric")
phylo = phylo %>% chronos()

# Does the tree of species look right?
plot(phylo)

# Make a little pruned tree for Fig. S2E
phylo %>% 
  drop.tip(which(!(.$tip.label %in% c(spp_figS2E, "Boli_infu")))) %>% 
  write.tree(file = here("04-suppfigs", "20220814_small.tre"))

# Expand the tree so it has a tip for every individual (distinguished by trailing numbers)
# the pheno input grouping is essential!
phenophylo_indl = expand_tree(pheno %>% group_by(sp, eid), phylo)
pheno_indl = phenophylo_indl$pheno
phylo_indl = phenophylo_indl$phylo

# Does the tree of individuals look right?
plot(phylo_indl)

# now we are ready to run regressions!
pgls_mods = pheno_indl %>% 
  # for the OU model, we are gonna try a range of different gamma values, i.e. branch length multipliers
  full_join(tibble(gamma = lapply(seq(0, 20, 1), function(x){10**x}) %>% unlist()), by = character()) %>% 
  group_by(pred_var, resp_var, class, gamma) %>% 
  summarise(
    mod_ou = safely(pgls_tbl, NULL)(
      cur_data(), 
      phylo_indl, 
      "resp ~ pred", # if we want to apply a transform to resp or pred, do it here
      startval = 1,
      gamma = gamma[[1]],
      tipcol = "indl", # a default anyway
      corfunc=ape::corMartins
    ) %>% .$result %>% list()
  )

# extract all the params nicely
pgls_results = pgls_mods %>% 
  rowwise() %>% 
  mutate(
    par_ou = mod_ou %>% 
      unpack_gls() %>% 
      list()
  ) %>% 
  unnest_wider(par_ou)

# find the best fits
params_best = pgls_results %>%
  drop_na(`logLik`) %>% 
  group_by(pred_var, resp_var, class) %>% 
  dplyr::rename(alpha = modelStruct.corStruct) %>% 
  # can adjust the fit selection criteria here
  arrange(-`logLik`, alpha) %>% 
  #filter(gamma == 1) %>% 
  slice(1) %>% 
  select(pred_var, resp_var, alpha, gamma, `logLik`, contains("coefficient"), mod_ou)

# calculate fits frames for each model across the domain
domain_fit = seq(0, 4000, 10)
fits_best = params_best %>% 
  rowwise() %>% 
  mutate(
    pred = list(domain_fit),
    resp = predictSE.gls(mod_ou, newdata = list("pred" = unlist(pred))) %>% list()
  ) %>% 
  select(-contains("mod")) %>% 
  ungroup() %>% 
  unnest_wider(resp) %>% 
  unnest(c(pred, fit, se.fit)) %>% 
  dplyr::rename(
    resp = fit,
    resp_serr = se.fit
  )
  
# plot them fits!
subs_resp = c(
  "frac_molar plasm",
  "frac_molar P-PE", 
  "frac_molar PE"#, 
  #"plci "
)
phenodata %>% 
  filter(resp_var %in% subs_resp) %>% 
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  facet_grid(cols = vars(pred_var), rows = vars(resp_var), scales = "free_y") +
  # underlay precalculated fits, robust to transformation unlike abline
  geom_ribbon(
    data = fits_best %>% filter(resp_var %in% subs_resp),
    aes(
      ymin = resp - resp_serr,
      ymax = resp + resp_serr
    ),
    #color = "black", # needed to override missing mapping
    fill = "black",
    alpha = 0.2
  ) +
  geom_line(
    data = fits_best %>% filter(resp_var %in% subs_resp),
    color = "black"
  ) +
  geom_point(
    aes(fill = sp),
    shape = 21,
    size = 2,
    color = "white"
  ) +
  theme_pubr() +
  scale_fill_manual(values = chroma_sp) +
  #coord_trans(x = scales::exp_trans(10)) +
  ggtitle("PGLS: [PPE] and PL curvature vs. depth")
ggsave(here("03-mainfigs", "ctenolipids", "20220609_pgls.pdf"), width=6, height=8)

# ordinary least-squares plot
phenodata %>% 
  #filter(resp_var %in% subs_resp) %>% 
  ggplot(
    aes(
      x = pred,
      y = resp
    )
  ) +
  facet_grid(
    cols = vars(pred_var), 
    rows = vars(resp_var), 
    scales = "free_y"
  ) +
  # underlay precalculated fits, robust to transformation unlike abline
  #geom_ribbon(
  #  data = fits_best %>% filter(resp_var %in% subs_resp),
  #  aes(
  #    ymin = resp - resp_serr,
  #    ymax = resp + resp_serr
  #  ),
  #  #color = "black", # needed to override missing mapping
  #  fill = "black",
  #  alpha = 0.2
  #) +
  #geom_line(
  #  data = fits_best %>% filter(resp_var %in% subs_resp),
  #  color = "black"
  #) +
  geom_smooth(
    method = "lm",
    #alpha = 0.2,
    color = "black"
  ) +
  geom_point(
    aes(fill = sp),
    shape = 21,
    size = 2,
    color = "white"
  ) +
  theme_pubr() +
  scale_fill_manual(values = chroma_sp) +
  #coord_trans(x = scales::exp_trans(10)) +
  guides(fill  = "none") +
  #ggtitle("OLS: [PPE] and [PE]\nvs. depth")
  labs(
    x = "Depth (m)",
    y = "Mole fraction",
    title = "Major PLs vs. depth"
  )
ggsave(here("04-suppfigs", "20220728_ols_all-major-pl.pdf"), width=3, height=8)

# rn need to check significance with e.g.:
params_best %>% .$mod_ou %>% lapply(., summary) %>% setNames(params_best$resp_var)
# [PPE] is significant, PLCI is not

#NTS 20220506: something is terribly wrong with the ttable unpacker; I should scrap it
# with REML, the `coefficients` columns are ust fine!