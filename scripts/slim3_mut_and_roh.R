# combine ROH and mutations from slim

combine_mut_roh <- function(run_name, roh_cutoff = 5000) {
   
   # cutoff for short/long ROH in KB
   #roh_cutoff <- 5000
   
   # get roh
   file_path_roh <- paste0("slim_sim/sims/roh/", run_name, ".hom")
   roh <- fread(file_path_roh) %>% mutate(IID = str_replace(IID, "indv", ""))
   
   # froh <- roh %>% 
   #       group_by(IID) %>% 
   #       summarise(froh = sum(KB)/100000) 
   # hist(froh$froh, breaks = 100) 
   
   # load mutations per individual
   run_name_mut <- str_replace(run_name, "sheep", "mutperind")
   file_path_mut <- paste0("slim_sim/sims/muts/", run_name_mut, ".txt")
   
   muts <- read_delim(file_path_mut, " ") %>% 
            dplyr::group_by(pedigree_id, pos) %>%  
            dplyr::summarise(pos = mean(pos), s = mean(s), copies = n()) %>% 
            dplyr::rename(ind_id = pedigree_id)
   
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
   muts$roh_kb <- map(unique(muts$ind_id), mut_in_roh) %>% unlist()
   
   # calculate categories: long, short and outside of roh
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
      setNames(c("ind_id", "short", "long")) %>% 
      replace_na(list(long = 0, short = 0)) %>% 
      mutate(outside_roh = 1e5 - (long + short)) %>% 
      mutate(ind_id = str_replace(ind_id, "indv", "")) %>% 
      pivot_longer(names_to = "roh_class", values_to = "roh_class_genome_cov", cols = short:outside_roh) 

   
   mut_and_roh <- muts %>% 
      mutate(ind_id = as.character(ind_id)) %>% 
      group_by(ind_id, roh_class) %>% 
      add_count(roh_class) %>%
      mutate(s = ifelse(copies == 2, s, s * 0.1)) %>% 
      summarise(sum_s = sum(s, na.rm = TRUE),
                num_mut = mean(n, na.rm = TRUE)#,
                #roh_class_genome_cov = sample(roh_class_genome_cov, 1)
      ) %>% 
      full_join(roh_ind, by = c("ind_id", "roh_class")) %>% 
      arrange(ind_id) %>% 
      mutate(s_per_MB = (sum_s / roh_class_genome_cov) * 1000,
             num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000) %>% 
      mutate(roh_class = factor(roh_class, levels = c("long", "short", "outside_roh"))) %>% 
      mutate(seed = str_replace(run_name, "sheep_", ""))
   
   
   
# source("../sheep_ID/theme_simple.R")
# ggplot(mut_and_roh, aes(roh_class, s_per_MB)) +
#    geom_boxplot(width = 0.5) +
#   # scale_y_continuous(limits = c(-0.01, 0)) +
#    geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
#    theme_simple()

}
