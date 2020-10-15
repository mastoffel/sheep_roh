# call ROH
library(tidyverse)
library(data.table)
library(snpStats)
library(GGally)
library(gghalves)
library(furrr)
source("../sheep_ID/theme_simple.R")
library(janitor)

# Chr lengths
chr_data <- read_delim("../sheep/data/sheep_genome/chromosome_info_oar31.txt", delim = "\t") %>% 
      rename(size_BP = Length,
             CHR = Part) %>% 
      mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
      .[2:27, ] %>% 
      summarise(sum_KB = sum(size_KB)) %>% 
      as.numeric()

# load sample qc
sample_qc <- read_delim("output/sample_qc.txt", delim = " ") %>% 
               mutate(id = as.character(id))

# fitness data
# annual measures of traits and fitness
fitness_path <- "../sheep/data/1_Annual_Fitness_Measures_April_20190501.txt"
annual_fitness <- read_delim(fitness_path, delim = "\t") %>% clean_names()
names(annual_fitness)

# prepare, order pedigree
sheep_ped <- read_delim("../sheep/data/SNP_chip/20190711_Soay_Pedigree.txt", 
                        delim = "\t",
                        col_types = "ccc") %>%
      as.data.frame() %>%
      MasterBayes::orderPed() 

# roh
roh <- fread("output/ROH/roh_cM.hom") %>%
      rename(id = IID) %>% 
      rename(cM = KB) %>%
      mutate(cM = cM/1e3) 

# linkage map
lmap <- read_delim("../sheep_ID/data/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
      rename(snp = SNP.Name,
             chr = Chr) %>% 
      select(chr, snp, cMPosition, cMPosition.Female, cMPosition.Male) %>% 
      filter(chr %in% 1:26)

full_length_map <- lmap %>% group_by(chr) %>% 
      summarise(max_per_chr = max(cMPosition)) %>% 
      summarise(genetic_map_length = sum(max_per_chr)) %>% 
      .[[1]]

# lmap %>% group_by(chr) %>% 
#       summarise(max_per_chr = max(cMPosition.Male)) %>% 
#       summarise(genetic_map_length = sum(max_per_chr))
# ROH classes ------------------------------------------------------------------

# expected ROH length using cM/Mb from Johnston et al (2016)
length_dist <- data.frame(g = c(1, 2,2^2, 2^3, 2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13)) %>%
      mutate(ROH_length_cM = 100 / (2*g))

IBD_df <- roh %>%
      mutate(length_cM = cM) %>%
      mutate(class = case_when(length_cM >= 25 ~ 1,
                               length_cM < 25 &      length_cM >= 12.5 ~ 2,
                               length_cM < 12.5 &    length_cM >= 6.25 ~ 4,
                               length_cM < 6.25 &    length_cM >=  3.125 ~ 8,
                               length_cM <  3.125 &  length_cM >=1.5625 ~ 16,
                               length_cM < 1.5625 &  length_cM >= 0.78125 ~ 32,
                               length_cM < 0.78125 ~ 64)) %>% # 0.610683047
      mutate(length_class = case_when(
            class == 1 ~ "> 25 (2g)",
            class == 2 ~ "25-12.5 (2-4g)",
            class == 4 ~ "12.5-6.25 (4-8g)",
            class == 8 ~ "6.25-3.13 (8-16g)",
            class == 16 ~ "3.13-1.56 (16-32g)",
            class == 32 ~ "1.56-0.78 (32-64g)",
            class == 64 ~ "0.78-0.39 (64-128g)"
           # class == 64 ~ ">0.6-1.2 (>32-64g)"
            # class == 128 ~ "0.6-0.3 (128G)"
      )) %>% 
      mutate(length_class = fct_reorder(length_class, class)) %>% 
      mutate(id = as.character(id)) %>% 
      group_by(id, class, length_class) %>%
      dplyr::summarise(prop_IBD = sum(length_cM) / full_length_map,
                       sum_cM = sum(length_cM)) #%>% 

IBD_df_with_0 <- IBD_df %>% 
      # add missing length classes as 0
      ungroup() %>% 
      tidyr::complete(length_class, nesting(id)) %>% 
      mutate(class = ifelse(is.na(class), length_class, class)) %>% 
      mutate(sum_cM = ifelse(is.na(sum_cM), 0, sum_cM)) %>% 
      mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD))

IBD_df_with_0 %>% arrange(id)

library(viridis)
col_pal <- viridis(7)
p_roh_dist <- IBD_df_with_0 %>% 
      #mutate(prop_IBD = prop_IBD / 10) %>% 
      ggplot(aes(length_class, prop_IBD, fill = length_class)) +
      geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
                      transformation_params = list(height = 0, width = 1.3, seed = 1)) +
      geom_half_boxplot(side = "r", outlier.color = NA,
                        width = 0.6, lwd = 0.3, color = "black",
                        alpha = 0.8) +
      theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
      ylab("% genome") +
      scale_y_continuous(labels = c("0", "5", "10"), breaks = c(0, 0.05, 0.1)) +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
      theme(legend.position = "none",
            #axis.ticks.x = element_blank(),
            axis.title=element_text(size = rel(1.1)), 
            axis.text = element_text(color = "black")) + 
      xlab("ROH classes in cM") 
p_roh_dist 
ggsave("figs/roh_len_dist_cM.jpg", p_roh_dist, height = 3, width = 8)
quantile(roh$cM, probs = c(0.25,0.5,0.75))

# define ROH length classes
calc_froh_classes <- function(roh_crit, roh_lengths) {
      
      roh_filt <- dplyr::case_when(
            roh_crit == "short"  ~ expr(cM < 1.56),
            roh_crit == "medium" ~ expr((cM >= 1.56) & (cM < 6.25)),
            roh_crit == "long"   ~ expr(cM >= 6.25),
            roh_crit == "all" ~ expr(cM > 0)
      )
      
      roh_lengths %>%
            dplyr::group_by(id) %>%
            #filter({{ roh_filt }}) %>% 
            filter(!!roh_filt) %>% 
            dplyr::summarise(sum_cM = sum(cM)) %>% 
            mutate(froh = sum_cM/ full_length_map) %>% 
            dplyr::select(id, froh) %>% 
            rename(!! paste0("froh_", roh_crit) := froh)
      
}

# proportion of ROH length classes in each genome. Individuals which
# do not have long ROH have 0 for this class.
ROH_classes <- c("short","medium", "long", "all")
froh <- purrr::map(ROH_classes, calc_froh_classes, roh) %>% 
      purrr::reduce(left_join, by = "id") %>% 
      replace_na(list(froh_long = 0)) %>% 
      #rename(id = ID) %>% 
      mutate(id = as.character(id))

#plot(froh$froh_all, homs_df$nsnps)

# fitness data
fitness_data <- annual_fitness %>% 
      mutate(id = as.character(id)) %>% 
      left_join(froh, by = "id") %>% 
      left_join(sample_qc)


# make data.frame for analysis
fitness_data <- fitness_data %>% 
      clean_names() %>% 
      rename(birth_year = birthyear,
             mum_id = mother) %>% 
      dplyr::select(id, survival, sheep_year, age, birth_year,
                    sex, mum_id, twin, froh_all, froh_long, froh_medium, froh_short, 
                    starts_with("offspring"), everything()) %>% 
      # for sex, check what Cas is
      filter(sex %in% c("F", "M")) %>% 
      clean_names() %>% 
      # only take individuals with very high genotyping call (imputation) rate.
      filter(call_rate > 0.99) %>% 
      # na rows should be discarded
      mutate_at(c("id", "sheep_year", "birth_year", "sex",
                  "mum_id", "twin"), as.factor)


save(fitness_data, file = "data/fitness_roh.RData") # formerly fitness_roh_df
save(sheep_ped, file = "data/sheep_ped.RData") # ordered ped

