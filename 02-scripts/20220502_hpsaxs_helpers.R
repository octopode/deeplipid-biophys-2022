# HELPER FUNCTIONS for HP-SAXS analyses

library(tidyverse)
library(ggplot2)
library(lubridate)
library(rjson)
library(pbapply)
library(ggpubr)
library(calecopal)

# load .dat file output by BIOXTAS-RAW, return profile dataframe and metadata JSON object
load_saxsdat = function(file){
  # read profile data
  scatdat = file %>% 
    read_table(
      comment = '#', 
      col_names = c("q", "iq", "err"),
      col_types = list(col_double(), col_double(), col_double())#,
    ) %>% 
    mutate_all(as.numeric)
  # read header data
  header = file %>% 
    read_file() %>% 
    sub(".*### HEADER:", '', .) %>% 
    gsub('#', '', .) %>% 
    rjson::fromJSON()
  return(list("data" = scatdat, "header" = header))
}

# Load HP-SAXS data, parsing P, T, etc into columns
# can run parallel if desired
load_hp_saxsdat = function(files, clust = 1L){
  files %>% 
    pblapply(
    cl = clust,
    function(file){
      # parse metadata from filename
      metadata = basename(file) %>% 
        strsplit('_') %>% .[[1]] %>% 
        .[1:5] %>% 
        setNames(c("samp", "press", "pdir", "temp", "rep")) %>% 
        bind_rows()
      # load [meta]data from file
      datahead = quietly(load_saxsdat)(file)$result
      # bind all metadata and return
      datahead$data %>% 
        bind_cols(metadata) %>% 
        bind_cols(tibble(
          datetime = datahead$header$counters$date,
          filename = datahead$header$filename
        ))
    }
  ) %>% 
    bind_rows() %>%  
    # clean up press and temp columns
    mutate(across(c(press, temp), ~as.numeric(gsub("[^0-9.-]", "", .x)))) %>% 
    # convert MPa to bar
    mutate(
      press = press*10,
      datetime = parse_date_time(datetime, orders = "%a %B %d %H:%M:%S %Y")
    )
}

# muted palette for phases
chroma_ph = cal_palette("bigsur")[order(c(3,1,2))] %>% 
  setNames(c("hii", "lb", "la"))

# from https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
sigfig <- function(vec, n=3){ 
  ### function to round values to N significant digits
  # input:   vec       vector of numeric
  #          n         integer is the required sigfig  
  # output:  outvec    vector of numeric rounded to N sigfig
  
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
  
}      # end of function   sigfig