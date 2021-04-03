library(tidyverse)
library(ggplot2)
library(gghalves)
library(data.table)
library(patchwork)
source("../sheep_ID/theme_simple.R")
library(vroom)
library(colorspace)

mut_df <- fread("output/qm_slim/slim1000200_bot_7030del/par_combs_popsize1_1000_popsize2_200.txt")

mut_all <- mut_df %>% 
      filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05)# %>% 
      #group_by(seed) %>% 
     # sample_frac(0.01)


# get mean fitness at each locus (mutation position)
mean_fit <- mut_all %>% 
      #sample_n(1000) %>% 
      mutate(s2 = case_when(
         copies == 2 ~ s,
         copies == 1 ~ s * mut1_dom_coeff
      )) %>% 
      group_by(seed, pos) %>% 
      mutate(mean_fit = sum(s2, na.rm = TRUE)/200) %>%
      mutate(inb_load = s2/mean_fit) %>%
      #filter(copies == 2) %>%
      filter(roh_class != "outside_roh") %>%
      group_by(seed, roh_class) %>%
      summarise(mean_load = mean(inb_load, na.rm = TRUE)) %>%
      mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short")))) %>%
      ungroup()

mean_fit <- mut_all %>% 
   #sample_n(1000) %>% 
   mutate(s2 = case_when(
      copies == 2 ~ s,
      copies == 1 ~ s * mut1_dom_coeff
   )) %>% 
   group_by(seed, pos) %>% 
   # not quite correct
   #mutate(mean_fit = (1+s2) * mut_freq^2 + 2*mut_freq*(1-mut_freq)*(1-mut1_dom_coeff*s2) + (1-mut_freq)^2) %>%
   #mutate(inb_load = (1+s2)/mean_fit) %>%
   mutate(mean_fit = (s2) * mut_freq^2 + 2*mut_freq*(1-mut_freq)*(mut1_dom_coeff*s2) + (1-mut_freq)^2) %>%
   mutate(inb_load = (s2)/mean_fit) %>%
   #filter(copies == 2) %>%
   filter(roh_class != "outside_roh") %>%
   group_by(seed, roh_class) %>%
   summarise(mean_load = mean(inb_load, na.rm = TRUE)) %>%
   mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short")))) %>%
   ungroup()

# %>% 
#       group_by(seed, pos) %>% 
#       rowwise() %>% 
#       mutate(mean_fit = sample(c(s, 0), size = 2, 
#                         replace = TRUE, prob = c(mut_freq, 1-mut_freq)))
#    
#     
#       # 400 alleles

# estimate B 
B <- mut_all %>% 
   group_by(seed, pos) %>% 
   add_count() %>% 
   mutate(pi = n / 400) %>% 
   mutate(B = s * (1 - 2 * mut1_dom_coeff) * pi * (1-pi)) %>% 
   filter(copies == 2) %>% 
   filter(roh_class != "outside_roh") %>% 
   group_by(seed, roh_class) %>% 
   summarise(mean_load = mean(B, na.rm = TRUE)) %>% 
   mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short")))) %>% 
   ungroup()

ggplot(mean_fit, aes(roh_class, mean_load, fill = roh_class)) +
      #geom_line(mapping = aes(group = seed), size = 0.2, alpha = 0.1,
      #        color = "#4C566A") +
      geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, 
                      size = 2,  width = 0.5) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "black",
                        alpha = 0.8, notch = FALSE) +
      scale_fill_viridis_d("ROH length class", direction = -1, option = "D",
                          guide = guide_legend(reverse = TRUE),
                          labels=rev(c("long (>12.5cM | < 4g)",
                                       "medium (1.56cM-12.5cM | 4-32g)",
                                       "short (0.39cM-1.56cM | 32-128g)"))) +
      theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
      ylab("inbreeding load") +
      theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            #axis.title.x = element_blank(),
            # panel.spacing = unit(2, "lines"),
            axis.text = element_text(color = "black")#,
           # legend.position = "none"
      ) +
      coord_flip()
ggsave("figs/inb_load.jpg", width = 7, height = 3)
