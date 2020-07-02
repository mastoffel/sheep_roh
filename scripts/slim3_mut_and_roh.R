# slim output to roh

library(data.table)
library(tidyverse)
library(vcfR)

# cutoff for short/long ROH in KB
roh_cutoff <- 3000

# use vcf output to call ROH
system(paste0("/usr/local/bin/plink --vcf slim_sim/sheep_recap.vcf --sheep --out output/roh ",
              # "--keep output/ROH/ids_surv.txt ",
              "--homozyg --homozyg-window-snp 20 --homozyg-snp 20 --homozyg-kb 1200 ",
              "--homozyg-gap 100 --homozyg-density 100 --homozyg-window-missing 0 ",
              "--homozyg-het 0 ",
              "--homozyg-window-het 0"))

file_path <- "output/roh.hom"
roh <- fread(file_path) %>% mutate(IID = str_replace(IID, "indv", ""))
hist(roh$KB)

froh <- roh %>% 
      group_by(IID) %>% 
      summarise(froh = sum(KB)/100000) 
hist(froh$froh, breaks = 100) 


# load mutations per individual
muts <- read_delim("slim_sim/mutperind2.txt", " ") %>% 
         group_by(pedigree_id, pos) %>%  
         summarise(pos = mean(pos), s = mean(s), copies = n()) %>% 
         rename(ind_id = pedigree_id)

mut_in_roh <- function(ind) {
   muts_ind <- muts %>% filter(ind_id == ind)
   roh_ind <- roh %>% filter(IID == ind)
   
   mut_roh <- map_dbl(1:nrow(muts_ind), function(x) {
      in_roh <- (roh_ind$POS1 <= muts_ind$pos[x]) & (roh_ind$POS2 >= muts_ind$pos[x])
      out <- ifelse(any(in_roh), roh_ind$KB[in_roh], NA)
   })
   mut_roh 
}

# if mutation is in roh, add roh length
muts$roh_kb <- map(unique(muts$ind_id), mut_in_roh) %>% 
            unlist()

muts <- muts %>% 
         mutate(roh_class = case_when(
            roh_kb > roh_cutoff ~ "long",
            roh_kb < roh_cutoff ~ "short",
            is.na(roh_kb) ~ "outside_roh"
         ) ) %>% 
        mutate(mut_in_roh = ifelse(is.na(roh_kb), FALSE, TRUE))

roh_ind <- roh %>% 
   mutate(roh_class = case_when(
      KB > roh_cutoff ~ "long",
      KB < roh_cutoff ~ "short"
   )) %>% 
   group_by(IID, roh_class) %>% 
   summarise(sum_kb = sum(KB)) %>% 
   pivot_wider(names_from = roh_class, values_from = sum_kb) %>% 
   setNames(c("ind_id", "long", "short")) %>% 
   replace_na(list(long = 0, short = 0)) %>% 
   mutate(outside_roh = 1e5 - (long + short)) %>% 
   mutate(ind_id = str_replace(ind_id, "indv", "")) %>% 
   pivot_longer(names_to = "roh_class", values_to = "roh_class_genome_cov", cols = long:outside_roh) 


mut_and_roh <- muts %>% 
                      mutate(ind_id = as.character(ind_id)) %>% 
                      left_join(roh_ind) %>% 
                      group_by(ind_id, roh_class) %>% 
                      add_count(roh_class) %>% 
                      mutate(s = ifelse(copies == 2, s, s * 0.1)) %>% 
                      summarise(sum_s = sum(s),
                                mean_s = mean(s),
                                var_s = var(s),
                                num_mut = mean(n)) %>% 
                      left_join(roh_ind) %>% 
                      mutate(s_per_MB = (sum_s / roh_class_genome_cov) * 1000,
                             num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000) %>% 
                      mutate(roh_class = factor(roh_class, levels = c("long", "short", "outside_roh")))
                             
source("../sheep_ID/theme_simple.R")           
ggplot(mut_and_roh, aes(roh_class, s_per_MB)) + 
   geom_boxplot() +
  # scale_y_continuous(limits = c(-0.01, 0)) +
  # geom_jitter(size = 2, alpha = 0.1) +
   theme_simple()
