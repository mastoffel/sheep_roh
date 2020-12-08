library(tidyverse)
library(data.table)
load("data/fitness_roh.RData") 

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

juv_survival <- fitness_data %>% 
   filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
   filter(age == 0)

juv_survival %>% 
   summarise(across(starts_with("froh_short"), list(mean = mean, sd = sd, min = min, max = max)))

             