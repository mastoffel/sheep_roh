library(tidyverse)
library(ggplot2)
library(gghalves)
source("../sheep_ID/theme_simple.R")

mut_df <- read_delim("output/qm_slim/par_combs_popsize1_1000_popsize2_200.txt", " ") %>% 
      mutate(popsize = 1000) 


# mut_p <- mut_df %>% 
#          filter(copies == 2) %>% 
#          filter(roh_class != "outside_roh") %>% 
#          group_by(id, roh_class) %>% 
#          add_count(roh_class) %>%
#          mutate(sum_s = sum(s, na.rm = TRUE),
#                 num_mut = mean(n, na.rm = TRUE)) %>% 
#          mutate(s_sum_per_MB = (sum_s / roh_class_genome_cov) * 1000,
#                 num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000) %>% 
#          group_by(mut1_dom_coeff, mut1_gam_mean, seed, roh_class) %>% 
#          summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
#                    num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE),
#                    mean_origin = mean(originG, na.rm = TRUE),
#                    mean_freq = mean(mut_freq, na.rm = TRUE)) %>% 
#          mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short"))))


mut_all <- mut_df %>% 
      #sample_frac(0.001) %>% 
      # filter homozygous sites
      filter(copies == 2) %>% 
      group_by(id, roh_class) %>% 
      add_count(roh_class) %>%
      summarise(sum_s = sum(s, na.rm = TRUE),
                num_mut = mean(n, na.rm = TRUE),
                mean_freq = mean(mut_freq, na.rm = TRUE),
                mean_originG = mean(originG),
                roh_class_genome_cov = first(roh_class_genome_cov),
                mut1_gam_mean = first(mut1_gam_mean),
                mut1_dom_coeff = first(mut1_dom_coeff),
                seed = first(seed)) %>% 
      mutate(s_sum_per_MB = (sum_s / roh_class_genome_cov) * 1000,
              num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000) #%>% 

mut_p <- mut_all %>% 
      filter(roh_class != "outside_roh") %>% 
      group_by(seed, mut1_dom_coeff, mut1_gam_mean, roh_class) %>% 
      summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
                num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE),
                mean_origin = mean(mean_originG, na.rm = TRUE),
                mean_freq = mean(mean_freq, na.rm = TRUE)) %>% 
      mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short")))) %>% 
      ungroup()

p1 <- mut_p %>% 
      ggplot(aes(roh_class, s_sum_per_MB, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.8, stroke = 0.1, size =2,
                      transformation_params = list(height = 0, width = 1.3, seed = 1)) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.8, lwd = 0.5, color = "black",
                        alpha = 0.8) +
      #scale_y_continuous(limits = c(-0.02, 0)) +
      scale_fill_viridis_d("ROH length class", direction = -1,
                           guide = guide_legend(reverse = TRUE),
                           labels=rev(c("long (>6.25cM)", "medium (>1.56cM & <6.25cM)", "short (<1.56cM)"))) +
      ylab("Selection coefficient per cM") +
      xlab("ROH length class") +
      facet_wrap(mut1_dom_coeff ~  mut1_gam_mean ) +
      scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
      #  ggtitle("Sim: weakly del\nROH cutoff 4900KB, \nmut1. dist: -0.03, 2, dom coeff 0.1 \nmut2. dist: -0.2, 3, dom coeff 0.01") +
      # geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
      theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
      theme(
            # panel.grid.major = element_blank(),
            #panel.grid.minor = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            #strip.text.y.right = element_text(angle = 0),
            # legend.position = "none",
            axis.text = element_text(color = "black"),
            legend.position = "top"
      ) +
      #ggtitle("Ancestral N = 1000, current N = 200")
      coord_flip()
p1
