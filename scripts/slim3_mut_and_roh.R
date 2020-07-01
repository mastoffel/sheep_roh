# slim output to roh

library(data.table)
library(tidyverse)
library(vcfR)

# use vcf output to call ROH
system(paste0("/usr/local/bin/plink --vcf slim_sim/sheep_recap.vcf --sheep --out output/roh ",
              # "--keep output/ROH/ids_surv.txt ",
              "--homozyg --homozyg-window-snp 30 --homozyg-snp 30 --homozyg-kb 600 ",
              "--homozyg-gap 100 --homozyg-density 100 --homozyg-window-missing 0 ",
              "--homozyg-het 0 ",
              "--homozyg-window-het 0"))

sheep_vcf <- read.vcfR("slim_sim/sheep.vcf")
gt <- extract.gt(sheep_vcf, element = "GT", as.numeric = TRUE)
gt
chromoqc(sheep_vcf)
summary(sheep_vcf@gt)

file_path <- "output/roh.hom"
roh <- fread(file_path)
hist(roh$KB)

froh <- roh %>% 
      group_by(IID) %>% 
      summarise(froh = sum(KB)/100000) 
hist(froh$froh, breaks = 100) 

# load mutation details
# del_muts <- read_delim("slim_sim/mut.txt", " ") %>% 
#                rename(mut_id = id)

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
            roh_kb > 3000 ~ "long",
            roh_kb < 3000 ~ "short",
            is.na(roh_kb) ~ "outside_roh"
         ) ) %>% 
        mutate(mut_in_roh = ifelse(is.na(roh_kb), FALSE, TRUE))

roh_ind <- roh %>% 
   mutate(roh_class = case_when(
      KB > 3000 ~ "long",
      KB < 3000 ~ "short"
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
                      summarise(sum_s = sum(s),
                                mean_s = mean(s),
                                var_s = var(s),
                                num_mut = mean(n)) %>% 
                      left_join(roh_ind) %>% 
                      mutate(s_per_MB = (sum_s / roh_class_genome_cov) * 1000,
                             num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000) %>% 
                      mutate(roh_class = factor(roh_class, levels = c("long", "short", "outside_roh")))
                             
                             
ggplot(mut_and_roh, aes(roh_class, s_per_MB)) + 
   geom_boxplot() +
   geom_jitter()
