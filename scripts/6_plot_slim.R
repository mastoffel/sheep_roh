library(tidyverse)
library(ggplot2)
library(gghalves)
source("../sheep_ID/theme_simple.R")
mut_df1 <- read_delim("output/qm_slim/par_combs_popsize1_1000_popsize2_200.txt", " ") %>% 
            mutate(popsize = 1000)

# mut_df1 %>% 
#    group_by(seed) %>% 
#    tally() -> df
# any(df$n != 800)
# 
# 
# mut_df2 <- read_delim("output/qm_slim/par_combs_popsize1_1000_popsize2_1000.txt", " ") %>% 
#    mutate(popsize = 10001000)


mut_df2 <- read_delim("output/qm_slim/par_combs_popsize1_5000_popsize2_200.txt", " ")%>% 
   mutate(popsize = 5000)
mut_df3 <- read_delim("output/qm_slim/par_combs_popsize1_10000_popsize2_200.txt", " ")%>% 
   mutate(popsize = 10000)

mut_all <- bind_rows(mut_df1, mut_df2, mut_df3) 

# filter out one parameter combination
mut_p <- mut_all %>% 
      #filter(roh_class != "outside_roh") %>% 
      group_by(popsize, mut1_dom_coeff, mut1_gam_mean, seed,  roh_class) %>% 
      summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
                num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE)) %>% 
      mutate(roh_class = factor(roh_class, levels = c("long", "medium","short")))
      #mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short"))))

# % deleteriousness lost
mut_p %>% 
   group_by(mut1_dom_coeff, roh_class) %>% 
   summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
             num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE)) 

p1 <- mut_p %>% 
      #filter(popsize == 1000) %>% 
      rename(selection = mut1_gam_mean,
             dominance = mut1_dom_coeff) %>% 
      ggplot(aes(roh_class, s_sum_per_MB, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.8, stroke = 0.1, size =2,
                      transformation_params = list(height = 0, width = 1.3, seed = 1)) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "black",
                        alpha = 0.8) +
      #scale_y_continuous(limits = c(-0.02, 0)) +
      scale_fill_viridis_d(direction = 1) +
      ylab("Selection coefficient per cM") +
      xlab("ROH length class") +
      scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
      facet_grid(popsize + dominance ~ selection, 
                 labeller = label_both, #scales = "free",
                 switch = "y") +
      #  ggtitle("Sim: weakly del\nROH cutoff 4900KB, \nmut1. dist: -0.03, 2, dom coeff 0.1 \nmut2. dist: -0.2, 3, dom coeff 0.01") +
      # geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
      theme_simple(grid_lines = TRUE, axis_lines = TRUE,  base_size = 12) +
      theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
           # axis.title.y = element_blank(),
           strip.text.y.right = element_text(angle = 0),
            legend.position = "none",
            axis.text = element_text(color = "black")
      ) +
      #ggtitle("Ancestral N = 1000, current N = 200")
      coord_flip()
p1
#ggsave("figs/sim_s_h_ne1000.jpg", p1, height = 6, width = 7)


mut_p <- mut_df2 %>% 
   filter(roh_class != "outside_roh") %>% 
   group_by(mut1_dom_coeff, mut1_gam_mean, seed,  roh_class) %>% 
   summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
             num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE)) %>% 
   mutate(roh_class = factor(roh_class, levels = c("long", "medium","short")))
#mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short"))))

p2 <- mut_p %>% 
   rename(selection = mut1_gam_mean,
          dominance = mut1_dom_coeff) %>% 
   ggplot(aes(roh_class, s_sum_per_MB, fill = roh_class)) +
   geom_half_point(side = "r", shape = 21, alpha = 0.8, stroke = 0.1, size =2,
                   transformation_params = list(height = 0, width = 1.3, seed = 1)) +
   geom_half_boxplot(side = "l", outlier.color = NA,
                     width = 0.5, lwd = 0.5, color = "black",
                     alpha = 0.8) +
   #scale_y_continuous(limits = c(-0.02, 0)) +
   scale_fill_viridis_d(direction = 1) +
   ylab("Selection coefficient per cM") +
   # xlab("ROH length class") +
   scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
   facet_grid(dominance ~ selection, labeller = label_both, scales = "free") +
   #  ggtitle("Sim: weakly del\nROH cutoff 4900KB, \nmut1. dist: -0.03, 2, dom coeff 0.1 \nmut2. dist: -0.2, 3, dom coeff 0.01") +
   # geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
   theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = "black")
   ) #+
#coord_flip()

library(patchwork)
p1 + p2
