library(tidyverse)
library(data.table)
roh <- fread("output/ROH/roh_cM.hom") %>%
      rename(id = IID) %>% 
      rename(cM = KB) %>%
      mutate(cM = cM/1e3) 
max(roh$cM)
mean(roh$cM)

roh %>% 
      group_by(id) %>% 
      summarise(prop_roh = sum(cM)/ 3146) %>% 
      summarise(mean(prop_roh), range(prop_roh))
