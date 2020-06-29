library(data.table)
library(tidyverse)
system(paste0("/usr/local/bin/plink --vcf slim_sim/sheep_recap.vcf --sheep --out output/roh ",
              # "--keep output/ROH/ids_surv.txt ",
              "--homozyg --homozyg-window-snp 20 --homozyg-snp 20 --homozyg-kb 700 ",
              "--homozyg-gap 100 --homozyg-density 100 --homozyg-window-missing 0 ",
              "--homozyg-het 0 ",
              "--homozyg-window-het 0"))

file_path <- "output/roh.hom"
roh <- fread(file_path)
hist(roh$KB)

froh <- roh %>% 
      group_by(IID) %>% 
      summarise(froh = sum(KB)/100000) 

hist(froh$froh, breaks = 100) 

# load mutations
del_muts <- read_delim("slim_sim/mut.txt", " ") %>% 
                  filter(s < 0)

count_del <- function(ind) {
   del_mut <- del_muts[ind, ]
   out <- data.frame(del = rep(NA, length(roh$POS1)),
                     s = rep(NA, length(roh$POS1)))
   pos <- del_mut$pos
   out$del <- (pos >= roh$POS1) & (pos <= roh$POS2)
   out[out$del, "s"] <- del_mut$s
   out
}

all_del <- map(1:nrow(del_muts), count_del) %>% 
            map("del") %>% 
            bind_cols() %>% 
            as.matrix() %>% 
            rowSums()

all_s <- map(1:nrow(del_muts), count_del) %>% 
   map("s") %>% 
   bind_cols() %>% 
   as.matrix() %>% 
   rowMeans(na.rm = TRUE)

all_del <- map(del_muts$pos, count_del) %>% 
            bind_cols() %>% 
            as.matrix() %>% 
            rowSums()

roh <- roh %>% 
            mutate(num_del = all_del,
                   s = all_s)
#autosomal_genome_size <- 2600
calc_froh_classes <- function(roh_crit, roh_lengths) {
      
      roh_filt <- dplyr::case_when(
            roh_crit == "short"  ~ expr(KB < 1200),
            roh_crit == "medium" ~ expr((KB > 1200)&(KB < 3000)),
            roh_crit == "long"   ~ expr(KB > 3000),
            roh_crit == "all" ~ expr(KB > 0)
      )
      
      roh_lengths %>%
            dplyr::group_by(IID) %>%
            #filter({{ roh_filt }}) %>% 
            filter(!!roh_filt) %>% 
            dplyr::summarise(KBSUM = sum(KB), del_vars = sum(num_del), s_del = sum(s)) %>% 
            mutate(FROH = KBSUM / 100000) %>% 
            dplyr::select(IID, FROH, del_vars, s_del) %>% 
            mutate(class = paste0("froh_", roh_crit)) %>% 
            rename(id = IID)
           # rename(ID = IID, !! paste0("FROH_", roh_crit) := FROH)
      
}

ROH_classes <- c("short", "medium", "long", "all")
froh <- purrr::map(ROH_classes, calc_froh_classes, roh) %>% 
     # purrr::reduce(left_join, by = "ID") %>% 
      bind_rows() %>% 
      mutate(del_prop = del_vars * (FROH * 1000))
      #replace_na(list(FROH_long = 0, FROH_medium = 0, FROH_short = 0)) %>% 

ggplot(froh, aes(del_vars, s_del)) + 
      geom_point() +
      facet_wrap(~class)
      
ggplot(froh, aes(class, s_del)) + 
   geom_boxplot()
