# combine ROH and mutations from slim

combine_mut_roh <- function(run_name, out_path, roh_cutoff_small = 1221, roh_cutoff_long = 4885) {
      
      # cutoff for short/long ROH in KB
      #roh_cutoff <- 5000
      
      # get roh
      file_path_roh <- paste0(out_path, "/roh/", run_name, ".hom")
      roh <- fread(file_path_roh) %>% 
            mutate(IID = str_replace(IID, "indv", "")) %>% 
            as_tibble()
      
      # load mutations per individual
      run_name_mut <- str_replace(run_name, "sheep", "mutperind")
      file_path_mut <- paste0(out_path, "/muts/", run_name_mut, ".txt")
      
      # same mutation id, same position
      # collapse each mutation into one row and add variable for homo/heterozygosity
      muts <- read_delim(file_path_mut, " ") %>% 
            dplyr::group_by(pedigree_id, mut_id) %>%  
            dplyr::summarise(pos = first(pos), s = first(s), copies = n()) %>% 
            dplyr::rename(ind_id = pedigree_id)
      
      # checks whether a mutation is part of an ROH and gives the ROH length if so, else NA
      mut_in_roh <- function(ind) {
            
            muts_ind <- muts %>% filter(ind_id == ind)
            roh_ind <- roh %>% filter(IID == ind)
            mut_roh <- map_dbl(1:nrow(muts_ind), function(x) {
                  # if the mutation is within a ROH, mark its row index with a logical
                  in_roh <- (roh_ind$POS1 <= muts_ind$pos[x]) & (roh_ind$POS2 >= muts_ind$pos[x])
                  # if the mutation is within a ROH, add ROH length
                  out <- ifelse(any(in_roh), roh_ind$KB[in_roh], NA)
            })
            #tibble(muts_ind, roh_kb = mut_roh) 
            muts_ind$roh_kb <- mut_roh
            muts_ind
      }
      
      # if mutation is in roh, add roh length
      muts_roh <- map_dfr(unique(muts$ind_id), mut_in_roh)
      
      # calculate categories: long, short and outside of roh
      muts_roh <- muts_roh %>% 
            mutate(roh_class = case_when(
                  roh_kb >= roh_cutoff_long ~ "long",
                  (roh_kb >= roh_cutoff_small) & (roh_kb < roh_cutoff_long) ~ "medium",
                  roh_kb < roh_cutoff_small ~ "short",
                  is.na(roh_kb) ~ "outside_roh"
            ) ) #%>% 
           # mutate(mut_in_roh = ifelse(is.na(roh_kb), FALSE, TRUE))
      
      # proportion of the genome in long and short ROH and outside of ROH
      roh_ind <- roh %>% 
            mutate(roh_class = case_when(
                  KB >= roh_cutoff_long ~ "long",
                  (KB >= roh_cutoff_small) & (KB < roh_cutoff_long) ~ "medium",
                  KB < roh_cutoff_small ~ "short"
            )) %>% 
            group_by(IID, roh_class) %>% 
            summarise(sum_kb = sum(KB)) %>% 
            pivot_wider(names_from = roh_class, values_from = sum_kb) %>% 
            rename(ind_id = IID) %>% 
            replace_na(list(long = 0, medium = 0, short = 0)) %>% 
            mutate(outside_roh = 1e5 - long - medium - short) %>% 
           # mutate(ind_id = str_replace(ind_id, "indv", "")) %>% 
            pivot_longer(names_to = "roh_class", values_to = "roh_class_genome_cov",
                         cols = -ind_id) # cols = long:outside_roh 
      
      
      muts_roh_full <- muts_roh %>% 
            mutate(ind_id = as.character(ind_id)) %>%
            group_by(ind_id, roh_class) %>% 
            # filter only homozygous sites
            filter(copies == 2) %>% 
            # how many deleterious homozygotes?
            add_count(roh_class) %>%
            # count selection coefficients only in homozygous state here
            #mutate(s = ifelse(copies == 2, s, s * 0.1)) %>% 
            summarise(sum_s = sum(s, na.rm = TRUE),
                      mean_s = mean(s, na.rm = TRUE),
                      num_mut = mean(n, na.rm = TRUE)#,
                      #roh_class_genome_cov = sample(roh_class_genome_cov, 1)
            ) %>% 
            full_join(roh_ind, by = c("ind_id", "roh_class")) %>% 
            arrange(ind_id) %>% 
            mutate(s_sum_per_MB = (sum_s / roh_class_genome_cov) * 1000,
                   num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000) %>% 
            mutate(roh_class = factor(roh_class, levels = c("long", "medium", "short", "outside_roh"))) %>% 
            mutate(seed = str_replace(run_name, "sheep_", ""))
      
      muts_roh_full
      
      
}
