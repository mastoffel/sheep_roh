# garlic bed output to data frame with froh_* 

library(tidyverse)
library(data.table)

# get genome size from assembly
autosomal_genome_size <- read_delim("../sheep/data/sheep_genome/chromosome_info_oar31.txt", delim = "\t") %>% 
      rename(size_BP = Length,
             CHR = Part) %>% 
      mutate(size_KB = size_BP / 1000)%>% 
      .[2:27, ] %>% 
      summarise(sum_KB = sum(size_KB)) %>% 
      as.numeric()

# check ROH
roh_garlic_ids <- fread("output/ROH_garlic/roh.roh.bed", fill = TRUE) %>% 
      filter(V1 == "track") %>% 
      select(V3) %>% 
      as_tibble()
roh_garlic <- fread("output/ROH_garlic/roh.roh.bed", fill = TRUE) %>% 
      filter(V1 != "track") %>% 
      dplyr::select(V1) %>% 
      separate(V1, sep = "\t", into = c("chr", "roh_start",
                                        "roh_end", "size_class",
                                        "roh_length", "p1", "p2", "p3", "RGB")) %>% 
      as_tibble()

num_roh_p_id <- fread("output/ROH_garlic/roh.roh.bed", fill = TRUE) %>% 
      select(V1) 
num_roh <- which(num_roh_p_id$V1 == "track")
num_roh_p_id <- c(1, num_roh[2:length(num_roh)] - 1:(length(num_roh)-1))

roh_garlic$ids <- NA_integer_
roh_garlic[num_roh_p_id, "ids"] <- roh_garlic_ids$V3

# all roh
roh <- roh_garlic %>% 
      fill(ids, .direction = "down") %>% 
      mutate(KB = as.numeric(roh_length)/1000)

# froh by size class
froh <- roh %>% 
            group_by(ids, size_class) %>% 
            summarise(sum_KB = sum(KB)) %>% 
            mutate(froh = sum_KB/autosomal_genome_size) %>% 
            pivot_wider(id_cols = ids, names_from = size_class, values_from = froh) %>% 
            rename(froh_shortg = A, froh_mediumg = B, froh_longg = C)

froh_all <- roh %>% 
                  group_by(ids) %>% 
                  summarise(sum_KB = sum(KB)) %>% 
                  mutate(frohg = sum_KB/autosomal_genome_size)

froh_full <- froh_all %>% 
                  left_join(froh) %>% 
                  rename(id = ids) %>% 
                  mutate(id = as.factor(id))


load("../sheep_ID/data/survival_mods_data.RData")

rohs <- fitness_data %>% 
      filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
      select(id, contains("froh")) %>% 
      group_by(id) %>% 
      summarise(across(.cols = contains("froh"), mean, na.rm = TRUE)) %>% 
      left_join(froh_full) 

kde <- read_delim("output/ROH_garlic/roh.200SNPs.kde", " ", col_names = F) %>% 
      plot()

plot(rohs$froh_short, rohs$froh_shortg)
mean(rohs$froh_all)
mean(rohs$frohg)
