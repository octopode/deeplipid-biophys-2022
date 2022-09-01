# SLIDE FIGURES FOR VARIOUS TALKS

library(here)
library(tidyverse)
library(ggpubr)
library(pbapply)

source(here("02-scripts", "20220502_lipidomics_helpers.R"))

## LCMS barplots
#source(here("02-scripts", "20220609_lipidomics.R"))

# correction factors
lcmsdata_long %>% 
  group_by(class, corrfac) %>% 
  summarise() %>% 
  filter(corrfac != 1) %>% 
  ungroup() %>% 
  #mutate(corrfac = log(corrfac, 2)) %>% 
  mutate(corrfac = ifelse(class %in% c("DG", "TG"), corrfac/10, corrfac)) %>% 
  arrange(-corrfac) %>% 
  mutate(class = class %>% factor(., levels = .)) %>% 
  ggplot(
    aes(
      x = class,
      y = corrfac,
      fill = class
    )
  ) +
  # draw little hairlines in between compounds
  geom_col(width = 0.95, size = 0.05, color = "black") +
  geom_hline(yintercept = 1) +
  scale_fill_manual( values = chroma_cl) +
  scale_color_manual(values = chroma_tx) +
  theme_pubk() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(
    x = "Lipid class", 
    y = "Fold difference in response factor"
  )
ggsave(here("05-slides", "20220711_corrfacs.pdf"), width=12, height=5.5)

# old vs new
# Tentacle vs. body comparison
lcmsdata_pl %>% 
  filter(sp == "Tjal_pink") %>% 
  group_by(sp_eid, class, annot) %>%
  # UNCORRECT
  mutate(frac_molar = frac_molar/corrfac) %>% 
  group_by(sp_eid) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
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
  geom_col(width = 0.95, size = 0.05, color = "black") +
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
  theme_pubk() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none", fill = "none") +
  labs(x = "Sample", y = "Mole fraction", fill = "Lipid class")
ggsave(here("05-slides", "20220711_Tjalfie-Old.pdf"), width=10, height=7.5)

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
  geom_col(width = 0.95, size = 0.05, color = "black") +
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
  theme_pubk() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(x = "Replicate", y = "Mole fraction", fill = "Lipid class", title = "E. coli PL composition")
ggsave(here("05-slides", "20220711_Ecoli_all_PLs.pdf"), width=8.5, height=7.5)

#All PLs all ind'ls
lcmsdata_pl %>%
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
  geom_col(width = 0.95, size = 0.05, color = "black") +
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
  theme_pubk() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(
    x = "Sample", 
    y = "Mole fraction", 
    fill = "Headgroup"#, 
    #title = "Phospholipid composition, all samples"
  )
ggsave(here("05-slides", "20220711_all_phospholipids_indls.pdf"), width=12, height=5.5)

# HG classes averaged by species
lcmsdata_pl_origset %>%
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
  theme_pubk() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(
    x = "Species", 
    y = "Mole fraction", 
    fill = "Headgroup"#, 
    #title = "Phospholipid composition, all samples"
  )
ggsave(here("05-slides", "20220609_all_phospholipids_spmeans.pdf"), width=12, height=4.75)

# HG classes averaged by species
lcmsdata_pl_origset %>%
  # sum within headgroups
  group_by(sp, sp_eid, class) %>%
  summarize(frac_molar = sum(frac_molar)) %>% 
  # average
  group_by(sp, class) %>% 
  summarize(
    serr_frac_molar = sd(frac_molar)/sqrt(n()),
    frac_molar = mean(frac_molar)
  ) %>% 
  # and renorm
  group_by(sp) %>% 
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>% 
  # get out PPE
  #filter(class %in% c("P-PE", "P-PC")) %>% 
  filter(class == "P-PE") %>% 
  arrange(class) %>% 
  ggplot(
    aes(
      x = sp,
      y = frac_molar,
      fill = class,
      ymin = frac_molar - serr_frac_molar,
      ymax = frac_molar + serr_frac_molar,
    )
  ) +
  # draw little hairlines in between compounds
  geom_col(width = 0.95, size = 0.05, color = "black") +
  geom_errorbar(width = 0.125, color = "grey50") +
  scale_fill_manual( values = chroma_cl) +
  theme_pubk() +
  theme(
    legend.position = "right",
    #legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.5)
  ) +
  guides(color = "none") +
  labs(
    x = "Species", 
    y = "Mole fraction", 
    fill = "Headgroup"#, 
    #title = "Phospholipid composition, all samples"
  )
ggsave(here("05-slides", "20220609_PPE_spmeans.pdf"), width=12, height=4.75)

## PGLS scatterplot
# always check dated references!
# only rerun the source script initially, then go interactive with it
#source(here("02-scripts", "20220504_pgls.R"))
# plot them fits!
subs_resp = c("frac_molar P-PE")
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
    fill = "white",
    alpha = 0.2
  ) +
  geom_line(
    data = fits_best %>% filter(resp_var %in% subs_resp),
    color = "white"
  ) +
  geom_point(
    aes(fill = sp),
    shape = 21,
    size = 3,
    color = "white"
  ) +
  theme_pubk() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = chroma_sp) +
  guides(fill = "none")
ggsave(here("05-slides", "20220609_pgls.pdf"), width=3.5, height=3.5)

# how about a nice tree for the PGLS?
pdf(here("05-slides", "20220511_timetree.pdf"), width=6, height=4)
plot(phylo)
dev.off()

# it would be nice to rotate some tips in FigTree!
phylo %>% write.tree(file = here("05-slides", "20220511_timetree.tre"))

## OLS scatterplot
# always check dated references!
# only rerun the source script initially, then go interactive with it
#source(here("02-scripts", "20220504_pgls.R"))
# plot them fits!
subs_resp = c("frac_molar P-PE")
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
  geom_smooth(
    method = "lm",
    fill = "white",
    alpha = 0.2,
    color = "white"
  ) +
  geom_point(
    aes(fill = sp),
    shape = 21,
    size = 3,
    color = "white"
  ) +
  theme_pubk() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = chroma_sp) +
  guides(fill = "none") +
  labs(
    x = "Collection depth (m)",
    y = "Mole fraction PPE"
  )
ggsave(here("05-slides", "20220609_ols.pdf"), width=3.5, height=3.5)

## Synthlipid phase diagrams!
#source(here("02-scripts", "20220510_hpsaxs_synthlipids.R"))
# master plot
data_saxs_fits %>% 
  ggplot(
    aes(
      x = press, 
      y = frac, 
      fill = ph
    )
  ) +
  facet_grid(cols = vars(comp)) +
  geom_area() + # the models, stacked by default
  geom_point(
    data = data_saxs_frac %>% 
      # avoid redundant plotting
      filter(
        ((logreg == 1) & (ph == "lb")) |
          ((logreg == 2) & (ph == "la"))
      ) %>% 
      # recode "bottom" points
      mutate(pdir = ifelse(press == 2000, "bt", pdir)),
    aes(shape = pdir),
    #position = position_dodge(100),
    fill = "black",
    color= "white"
  ) +
  theme_pubk() +
  # orientation control
  coord_flip() +
  scale_x_reverse() +
  scale_y_continuous(breaks = c(0, 1)) +
  # show bidirectional data
  scale_fill_manual(values = chroma_ph) +
  scale_shape_manual(values = c("dn" = 24, "up" = 25, "bt" = 22)) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(
    panel.spacing = unit(0.25, "in")
  ) +
  guides(
    fill  = "none",
    shape = "none"
  ) +
  labs(
    x = "Experimental pressure (bar)",
    y = "Intensity ratio"
  )
ggsave(here("05-slides", "20220511_hpsaxs_synlipids.pdf"), width=5, height=4)

# Some nice raw profiles to show what's going on?

## MD SIMULATIONS!
#source(here("02-scripts", "20220510_mdsims.R"))
barwid = 100
simdata_summ %>% 
  filter(
    (press_sim <= 1000) |
      (chains == "18:0/20:4")
  ) %>% 
  # kbc0 should be dropped in place of ka when we get it
  filter(par %in% c("apl", "dpp", "ka")) %>% 
  filter(chains %in% c("18:0/18:1", "18:0/20:4")) %>% 
  filter((par != "ka") | (chains == "18:0/18:1")) %>% 
  group_by(par, press_sim) %>% 
  #arrange(desc(chains)) %>% 
  arrange(chains) %>% 
  mutate(chains = chains %>% factor(., levels = unique(.))) %>% 
  #mutate(comp = str_glue(comp, ' ', temp_sim, "°C")) %>% 
  ggplot(
    aes(
      x = press_sim,
      y = mean,
      ymin = mean - serr,
      ymax = mean + serr,
      color = class,
      group = class
    )
  ) +
  facet_grid(rows = vars(par), cols = vars(chains), scales = "free_y") +
  #geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.5) +
  geom_errorbar(width = barwid, position = position_dodge(width = 0.75*barwid)) +
  geom_point(position = position_dodge(width = 0.75*barwid)) +
  scale_color_manual(values = chroma_ppe) +
  theme_pubk() +
  labs(
    #title = "MD results summaries",
    x = "Simulation pressure",
    y = "Parameter value"
  ) +
  guides(color = "none")
ggsave(here("05-slides", "20220511_mdsims.pdf"), width=6.5, height=4)
