# Generate a Cubette animation in the visual language of the 2022 biophys paper
library(here)
library(pbapply)
library(gridExtra)

# load up all the raw data for sample JWL191P
source(here("20220719_JWL0238_PUFA-PE_POPG", "laurdan_dual.R"))

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
      #plot.background = element_blank(),
      plot.background = element_rect(color = NA,  fill = "black"),
      plot.title = element_text(color = "white")
    )
}

# "Cubette state" plot
plot_pt = function(data_raw, bgcolor = "white"){
  data_raw %>%
  select(clock, T_act, P_act) %>%
  arrange(clock) %>%
  mutate(age = difftime(last(clock), clock, units = "sec") %>% as.numeric()) %>%
  ggplot(aes(x = T_act, y = P_act, color = age)) +
  geom_point(alpha = 0.017) +
  xlim(limits = range(grid$T_act_mean)) +
  scale_y_reverse(limits = rev(range(grid$P_act_mean))) +
  scale_color_gradient(low="red", high=ifelse(bgcolor == "black", "white", "black"), na.value = "black", limits = c(0, 600)) +
  theme_pubr() +
  #ifelse(bgcolor == "black", theme_pubk(), theme_pubr()) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  labs(
    x = "Experimental temperature (°C)",
    #y = "Experimental pressure (bar)",
    y = " ",
    title = "Cubette state"
  ) +
  guides(
    x = "none"
  )
}

# 20220604 from 20220410_LaurdanContour.R
plot_gp = function(data_gp, bgcolor = "white"){
  scape_gp = data_gp %>% 
    # run out loesses
    #gp2scape(grid, surface = "direct")
    gp2scape(grid)
  #print(scape_gp) #TEST
  scape_gp %>% 
    # keep it reasonable to avoid too many contour lines
    filter((gp >= -1) & (gp <= 1)) %>% 
    # plot GP contours overlaid with data points
    ggplot(aes(x = T_act_mean, y = P_act_mean, z = gp)) +
    geom_contour_filled(
      aes(fill = stat(level_mid)),
      binwidth = bwid
    ) +
    geom_contour(
      binwidth = bwid,
      color = "black"
    ) +
    scale_fill_distiller(
      palette = "YlGnBu", 
      direction = 1,
      guide = guide_colorbar(
        direction = "vertical"
      )
    ) +
    geom_point(
      data = data_gp,
      color = "black",
      alpha = 0.1,
      shape = 'o'
    ) +
    theme_pubr() +
    #ifelse(bgcolor == "black", theme_pubk(), theme_pubr()) +
    scale_y_reverse(limits = rev(range(grid$P_act_mean))) +
    scale_x_continuous(limits = range(grid$T_act_mean)) +
    labs(
      #x = element_blank(),
      x= "Experimental temperature (°C)",
      y = "Experimental pressure (bar)",
      fill = "C-laurdan GP"
    ) +
    #guides(x = "none") +
    theme(legend.position = "none") +
    ggtitle("C-laurdan GP")
}

# how many frames?
nframes = 100
interval = diff(data_raw$clock %>% range()) / nframes
bgcolor = "white"

seq(nframes) %>% 
  .[48] %>% 
  #.[3:length(.)] %>% 
  pblapply(
    cl = 1L,
    function(framenum){
      img_out = here("20220719_JWL0238_PUFA-PE_POPG", "JWL0238_reanim", paste(framenum, "pdf", sep='.'))
      begtime = data_raw$clock[[1]]
      endtime = begtime + (framenum * interval)
      print(endtime)
      # generate state plot
      ptplot = data_raw %>% 
        filter(clock <= endtime) %>% 
        plot_pt(bgcolor = bgcolor)
      # and contour plot
      #print(data_gp%>% 
      #        filter(wl_ex == 340) %>% 
      #        filter(clock_mean <= endtime))
      gpplot = data_gp %>% 
        filter(wl_ex == 340) %>% 
        filter(clock_mean <= endtime) %>% 
        # just in case it's empty
        safely(function(data){plot_gp(data, bgcolor = bgcolor)}, otherwise = ggplot(), quiet = FALSE)() %>% 
        .$result
      message(img_out)
      image = arrangeGrob(ptplot, gpplot, ncol=1, heights=c(1,2))
      ggsave(file = img_out, plot = image, width = 5.2, height = 6.4)
    }
  )
