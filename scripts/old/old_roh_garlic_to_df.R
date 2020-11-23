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

ggplot(roh, aes(KB)) +
   geom_histogram(bins = 100) +
   scale_x_continuous(limits = c(1, 5000))
   #facet_wrap(~size_class, scales = "free_x")

# froh by size class
froh <- roh %>% 
            group_by(ids, size_class) %>% 
            summarise(sum_KB = sum(KB)) %>% 
            mutate(froh = sum_KB/autosomal_genome_size) %>% 
            pivot_wider(id_cols = ids, names_from = size_class, values_from = froh) %>% 
            #dplyr::rename(froh_shortg = A, froh_longg = B)
            dplyr::rename(froh_shortg = A, froh_mediumg = B, froh_longg = C)

froh_all <- roh %>% 
                  group_by(ids) %>% 
                  summarise(sum_KB = sum(KB)) %>% 
                  mutate(froh_allg = sum_KB/autosomal_genome_size)

froh_full <- froh_all %>% 
                  left_join(froh) %>% 
                  rename(id = ids) %>% 
                  mutate(id = as.factor(id))

cor(froh_full$froh_shortg, froh_full$froh_longg)
hist(froh_full$froh_longg)

load("../sheep_ID/data/survival_mods_data.RData")

froh_complete <- fitness_data %>% 
      filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
      select(id, contains("froh")) %>% 
      group_by(id) %>% 
      summarise(across(.cols = contains("froh"), mean, na.rm = TRUE)) %>% 
      left_join(froh_full) 

kde <- read_delim("output/ROH_garlic/roh.100SNPs.kde", " ", col_names = F) %>% 
      plot()

hist(as.numeric(roh$roh_length)/1000)

# make dataset for modeling
fitroh_data <- fitness_data %>% 
                  left_join(froh_full %>% mutate(id = as.character(id)))

save(fitroh_data, file = "data/fitroh_data.RData")
write_delim(roh, "output/ROH_garlic/roh_df.txt")
