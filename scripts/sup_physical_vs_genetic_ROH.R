library(tidyverse)
library(data.table)
library(patchwork)
source("scripts/theme_simple.R")
load("data/fitness_roh.RData") 

# call physical map ROH
system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar31_17052020 ",
              "--sheep --out output/ROH/roh_bp ",
              "--homozyg --homozyg-window-snp 25 --homozyg-snp 25 --homozyg-kb 390 ",
              "--homozyg-gap 250 --homozyg-density 100 --homozyg-window-missing 2 ",
              "--homozyg-het 2 ",
              "--homozyg-window-het 2"))

# ROH in cM
froh_cM <- fitness_data %>% 
      filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
      filter(age == 0) %>% 
      select(c(id, starts_with("froh"))) %>% 
      rename_with(~paste0(., "_cM"), starts_with("froh"))

# autosomal genome size
chr_data <- read_delim("data/chromosome_info_oar31.txt", delim = "\t") %>% 
      rename(size_BP = Length,
             CHR = Part) %>% 
      mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
      .[2:27, ] %>% 
      summarise(sum_KB = sum(size_KB)) %>% 
      as.numeric()

roh_bp <- fread("output/ROH/roh_bp.hom") %>%
      rename(id = IID) %>% 
      select(id, KB) %>% 
      mutate(mb = KB/1000)

# define ROH length classes
calc_froh_classes <- function(roh_crit, roh_lengths) {
      
      roh_filt <- dplyr::case_when(
            roh_crit == "short"  ~ expr(mb < 1.5625),
            roh_crit == "medium" ~ expr((mb >= 1.5625) & (mb < 12.5)),
            roh_crit == "long"   ~ expr(mb >= 12.5),
            roh_crit == "all" ~ expr(mb > 0)
      )
      
      roh_lengths %>%
            dplyr::group_by(id) %>%
            #filter({{ roh_filt }}) %>% 
            filter(!!roh_filt) %>% 
            dplyr::summarise(sum_mb = sum(mb)) %>% 
            mutate(froh = sum_mb/ (autosomal_genome_size/1000)) %>% 
            dplyr::select(id, froh) %>% 
            rename(!! paste0("froh_", roh_crit) := froh)
      
}

# proportion of ROH length classes in each genome. Individuals which
# do not have long ROH have 0 for this class.
ROH_classes <- c("short","medium", "long", "all")
froh_bp <- purrr::map(ROH_classes, calc_froh_classes, roh_bp) %>% 
      purrr::reduce(left_join, by = "id") %>% 
      replace_na(list(froh_long = 0)) %>% 
      #rename(id = ID) %>% 
      mutate(id = as.character(id)) %>% 
      rename_with(~paste0(., "_Mb"), starts_with("froh"))

froh <- froh_cM %>% 
         left_join(froh_bp)

# make plot

library(ggpubr)

cols <- viridis::viridis(3)

point_alpha <- 0.3
point_size <- 2

p1 <- ggplot(froh, aes(froh_all_cM, froh_all_Mb)) +
   geom_point(shape = 21, alpha = point_alpha, fill = "black",
              size = point_size) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
   xlab(expression(F[ROH]~(genetic))) +
   ylab(expression(F[ROH]~(physical))) +
   stat_cor(method = "pearson", aes(label = ..r.label..))
p1
p2 <- ggplot(froh, aes(froh_long_cM, froh_long_Mb)) +
   geom_point(shape = 21, alpha = point_alpha, fill = cols[1],
              size = point_size) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
   scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
   xlab(expression(F[ROH[long]]~(genetic))) +
   ylab(expression(F[ROH[long]]~(physical))) +
   stat_cor(method = "pearson", aes(label = ..r.label..))

p3 <- ggplot(froh, aes(froh_medium_cM, froh_medium_Mb)) +
   geom_point(shape = 21, alpha = point_alpha, fill = cols[2],
              size = point_size)+
   theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
   xlab(expression(F[ROH[medium]]~(genetic))) +
   ylab(expression(F[ROH[medium]]~(physical))) +
   stat_cor(method = "pearson", aes(label = ..r.label..))

p4 <- ggplot(froh, aes(froh_short_cM, froh_short_Mb)) +
   geom_point(shape = 21, alpha = point_alpha, fill = cols[3],
              size = point_size)+
   theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
   scale_x_continuous(breaks = c(0.07, 0.1, 0.13)) +
   xlab(expression(F[ROH[short]]~(genetic))) +
   ylab(expression(F[ROH[short]]~(physical))) +
   stat_cor(method = "pearson", aes(label = ..r.label..))
p4

p_final <- (p1 + p2) / (p3 + p4) +
   plot_annotation(tag_levels = "a") &
   theme(axis.title = element_text(margin=margin(t=-20, r=-20)),
         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
         plot.tag = element_text(face = 'bold'))

ggsave(plot = p_final, filename = "figs/cor_physical_genetic_roh.jpg",
       width = 6.6, height = 6)
         