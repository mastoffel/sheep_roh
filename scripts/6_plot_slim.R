library(tidyverse)
library(ggplot2)
library(gghalves)
source("../sheep_ID/theme_simple.R")
mut_df1 <- read_delim("output/qm_slim/slim1000200/out/par_combs_popsize1_1000_popsize2_200.txt", " ") %>% 
            mutate(popsize = 1000)
mut_df2 <- read_delim("output/qm_slim/slim5000200/out/par_combs_popsize1_5000_popsize2_200.txt", " ")%>% 
   mutate(popsize = 5000)
mut_df3 <- read_delim("output/qm_slim/slim10000200/out/par_combs_popsize1_10000_popsize2_200.txt", " ")%>% 
   mutate(popsize = 10000)

mut_all <- bind_rows(mut_df1, mut_df2, mut_df3) 

# filter out one parameter combination
mut_p <- mut_all %>% 
      filter(roh_class != "outside_roh") %>% 
      group_by(popsize, mut1_dom_coeff, mut1_gam_mean, seed,  roh_class) %>% 
      summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
                num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE)) %>% 
      mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short"))))
      #mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short"))))

# % deleteriousness lost
mut_p %>% 
   group_by(mut1_dom_coeff, roh_class) %>% 
   summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
             num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE)) 

p1 <- mut_p %>% 
      filter(!is.na(roh_class)) %>% 
      mutate(popsize = case_when(
         popsize == 1000 ~ "1k",
         popsize == 5000 ~ "5k",
         popsize == 10000 ~ "10k",
      )) %>% 
      mutate(popsize = factor(popsize, levels = c("1k", "5k", "10k"))) %>% 
      rename(s = mut1_gam_mean,
             h = mut1_dom_coeff,
             Nanc = popsize) %>% 
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
      scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
      scale_y_continuous(breaks = c(0, -0.005, -0.01),
                         labels = c("0", "-0.005", "-0.01")) +
      facet_grid(h + Nanc ~ s, 
                 labeller = label_both, #scales = "free",
                 switch = "y") +
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
ggsave("figs/sim_s_h_ne_all.jpg", p1, height = 8, width = 7)



# % differences

mut_p %>% 
   group_by(roh_class) %>% 
   summarise(mean(s_sum_per_MB))

# % decrease long to medium
1 - (-0.00249 /   -0.00339)
1 - ( -0.00136 /   -0.00339)

mut_p %>% 
   group_by(roh_class, mut1_gam_mean) %>% 
   summarise(mean(s_sum_per_MB))

1- (-0.00133 /  -0.00402)
1- (-0.00140 /  -0.00264)

mut_p %>% 
   group_by(roh_class, mut1_dom_coeff) %>% 
   summarise(mean(s_sum_per_MB))

1- (-0.00194 /   -0.00477)
1- (-0.000800 /   -0.00195 )

mut_p <- mut_df2 %>% 
   filter(roh_class != "outside_roh") %>% 
   group_by(mut1_dom_coeff, mut1_gam_mean, seed,  roh_class) %>% 
   summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
             num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE)) %>% 
   mutate(roh_class = factor(roh_class, levels = c("long", "medium","short")))
#mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short"))))

p2 <- mut_p %>% 
   filter(popsize = 5000) %>% 
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
