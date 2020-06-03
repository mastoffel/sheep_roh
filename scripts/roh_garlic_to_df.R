library(tidyverse)
library(data.table)

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
roh <- roh_garlic %>% 
      fill(ids, .direction = "down") %>% 
      mutate(KB = as.numeric(roh_length)/1000)

froh <- roh %>% 
      group_by(ids, size_class) %>% 
      summarise(sum_KB = sum(KB)) %>% 
      mutate(froh = sum_KB/autosomal_genome_size)

ggplot(froh, aes(size_class, froh)) +
      geom_boxplot() +
      geom_point()

froh %>% 
      select(-sum_KB) %>% 
      pivot_wider(names_from = size_class, values_from = froh) %>% 
      ggpairs(columns = 2:4)


