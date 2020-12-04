library(tidyverse)
library(ggplot2)
library(gghalves)
library(data.table)
library(patchwork)
source("../sheep_ID/theme_simple.R")
library(vroom)
mut_df <- read_delim("output/qm_slim/slim1000200_bot_3e9/par_combs_popsize1_1000_popsize2_200.txt", " ",
                     n_max = 1000) %>% 
      mutate(popsize = 1000) 

#mut_df <- vroom("output/qm_slim/slim1000200_bot_3e9/par_combs_popsize1_1000_popsize2_200.txt")

mut_df <- fread("output/qm_slim/slim1000200_bot_7030del/par_combs_popsize1_1000_popsize2_200.txt")
mut_all <- mut_df %>% 
      #sample_frac(0.001) %>% 
      # filter homozygous sites
      filter(copies == 2) %>%
      #filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
      #filter(s < -0.005) %>% 
      group_by(id, roh_class) %>% 
      #add_count(id, roh_class, name = "num_mut0") %>%
      summarise(sum_s = sum(s, na.rm = TRUE),
                #num_mut = mean(num_mut0, na.rm = TRUE),
                num_mut = n(),
                mean_freq = mean(mut_freq, na.rm = TRUE),
                median_freq = median(mut_freq, na.rm = TRUE),
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

# combine 
mut_sub <- mut_p %>%
   filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
   mutate(mean_origin = 11000 - mean_origin)
 
make_plot <- function(ss, axis_title, df) {
   ss <- ensym(ss)
   p <- ggplot(df, aes(roh_class, !!ss, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.3, stroke = 0.1, 
                      size = 2,  width = 0.5) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "black",
                        alpha = 0.8, notch = FALSE) +
      scale_fill_viridis_d("ROH length class", direction = -1, option = "D",
                         guide = guide_legend(reverse = TRUE),
                         labels=rev(c("long (>12.5cM | < 4g)", 
                                      "medium (>1.56cM | 4-32g)", 
                                      "short (>0.39cM | 32-128g)"))) +
     # scale_fill_manual("ROH length class", values = c("#D9B08C","#116466","#2C3531"),
     #                      guide = guide_legend(reverse = TRUE),
     #                      labels=rev(c("long (>12.5cM | < 4g)", "medium (>1.56cM | < 32g)", "short (>0.39cM | <128g)"))) +
      #scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
      theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
      ylab(axis_title) +
      theme(
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(), 
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         #axis.title.x = element_blank(),
        # panel.spacing = unit(2, "lines"),
         axis.text = element_text(color = "black"),
         legend.position = "right"
      ) +
      coord_flip()
   p
}

p1 <- make_plot(s_sum_per_MB,"Selection coefficient per cM", mut_sub)  
p2 <- make_plot(mean_freq,"Mutation frequency", mut_sub)  
p3 <- make_plot(num_mut_per_MB, "Number of mutations per cM", mut_sub)
p4 <- make_plot(mean_origin,"Mutation age in generations", mut_sub)
p_final <- p1 + p2 + p3 + p4 + 
   plot_layout(guides = 'collect') #+
  # plot_annotation(tag_levels = 'a')
p_final

ggsave(plot = p_final, filename = "figs/Fig2_gen4_32_7030.jpg", width = 7, height = 4.5)





ggplot(mut_sub, aes(roh_class, s_sum_per_MB, fill = roh_class)) +
   geom_line(mapping = aes(group = seed), size = 0.2, alpha = 0.1) +
   geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, size = 2) +
   geom_half_boxplot(side = "l", outlier.color = NA,
                     width = 0.5, lwd = 0.5, color = "black",
                     alpha = 0.8, notch = TRUE) +
   scale_fill_viridis_d("ROH length class", direction = -1,
                        guide = guide_legend(reverse = TRUE),
                        labels=rev(c("long (>6.25cM)", "medium (1.56cM - 6.25cM)", "short (<1.56cM)"))) +
   #scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
  # ylab(axis_title) +
   theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      #axis.title.x = element_blank(),
      # panel.spacing = unit(2, "lines"),
      axis.text = element_text(color = "black"),
      legend.position = "right"
   ) +
   coord_flip()












library(ggbeeswarm)
mut_sub_summary <- mut_sub %>% group_by(roh_class) %>%  summarise(s_sum_per_MB = mean(s_sum_per_MB))
ggplot(mut_sub, aes(roh_class, s_sum_per_MB, color = roh_class)) +
   geom_point() +
   geom_quasirandom(width = 0.4, method = "frowney") + 
  # geom_beeswarm() +
  # scale_y_continuous(limits = c(-0.005, 0)) +
   geom_crossbar(data=mut_sub_summary, aes(ymin = s_sum_per_MB, ymax = s_sum_per_MB),
                 size=0.5,col="grey", width = 0.4) +
   scale_color_viridis_d("ROH length class", direction = -1,
                        guide = guide_legend(reverse = TRUE),
                        labels=rev(c("long (>12.5cM)", "medium (1.56cM - 12.5cM)", "short (<1.56cM)"))) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
   ylab("Selection coefficient per cM") +
   theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(color = "black"),
      legend.position = "right"
   ) +
   coord_flip()







# no grouping before plotting
mut_plot <- mut_df %>% 
               filter(roh_class != "outside_roh") %>% 
               filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
               filter(copies == 2) %>% 
             #  sample_frac(0.0001) %>% 
               group_by(id, roh_class) %>% 
               summarise(s_per_cM = (sum(s)/roh_class_genome_cov) * 1000,
                         num_mut_per_cM = (sum(n()) / roh_class_genome_cov)*1000)

ggplot(mut_plot, aes(roh_class, s_per_cM, fill = roh_class)) +
   geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, size = 2,
                     transformation_params = list(height = 0, width = 1.3, seed = 1)) +
   geom_half_boxplot(side = "l", outlier.color = NA,
                     width = 0.5, lwd = 0.5, color = "black",
                     alpha = 0.8) +
  # scale_y_continuous(limits = c(-0.03, 0)) +
   #facet_wrap(mut1_gam_mean ~ mut1_gam_shape) +
   scale_fill_viridis_d("ROH length class", direction = -1,
                        guide = guide_legend(reverse = TRUE),
                        labels=rev(c("long (>6.25cM)", "medium (1.56cM - 6.25cM)", "short (<1.56cM)"))) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
   #ylab(axis_title) +
   theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(color = "black"),
      legend.position = "right"
   ) +
   coord_flip()



# different grouping
mut_all <- mut_df %>% 
   sample_frac(1) %>% 
   # filter homozygous sites
   filter(copies == 2) %>% 
   filter(roh_class != "outside_roh") %>% 
   filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
   group_by(id, roh_class) %>%  # mut1_dom_coeff, mut1_gam_mean, 
   #add_count(id, roh_class, name = "num_mut0") %>%
   summarise(sum_s = sum(s, na.rm = TRUE),
             #num_mut = mean(num_mut0, na.rm = TRUE),
             num_mut = n(),
             mean_freq = mean(mut_freq, na.rm = TRUE),
             median_freq = median(mut_freq, na.rm = TRUE),
             mean_origin = 11000 - mean(originG),
             roh_class_genome_cov = first(roh_class_genome_cov),
             mut1_gam_mean = first(mut1_gam_mean),
             mut1_dom_coeff = first(mut1_dom_coeff),
             seed = first(seed)) %>% 
   mutate(s_sum_per_MB = (sum_s / roh_class_genome_cov) * 1000,
          num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000) %>% 
   mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short")))) 

make_plot <- function(ss, axis_title, df) {
   ss <- ensym(ss)
   p <- ggplot(df, aes(roh_class, !!ss, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, size = 2,
                      transformation_params = list(height = 0, width = 1.3, seed = 1)) +
      # geom_half_boxplot(side = "l", outlier.color = NA,
      #                   width = 0.5, lwd = 0.5, color = "black",
      #                   alpha = 0.8) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "black",
                        alpha = 0.8) +
      #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
      scale_fill_viridis_d("ROH length class", direction = -1,
                           guide = guide_legend(reverse = TRUE),
                           labels=rev(c("long (>6.25cM)", "medium (1.56cM - 6.25cM)", "short (<1.56cM)"))) +
      #scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
      theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
      ylab(axis_title) +
      theme(
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(), 
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         #axis.title.x = element_blank(),
         # panel.spacing = unit(2, "lines"),
         axis.text = element_text(color = "black"),
         legend.position = "right"
      ) +
      coord_flip()
   p
}

p1 <- make_plot(s_sum_per_MB,"Selection coefficient per cM", mut_all)  
p2 <- make_plot(mean_freq,"Mutation frequency", mut_all)  
p3 <- make_plot(num_mut_per_MB, "Number of mutations per cM", mut_all)
p4 <- make_plot(mean_origin,"Mutation age in generations", mut_all)
p_final <- p1 + p2 + p3 + p4 + 
   plot_layout(guides = 'collect') +
   plot_annotation(tag_levels = 'a')
p_final


ggplot(mut_all, aes(roh_class, s_sum_per_MB, fill = roh_class)) +
   geom_boxplot(outlier.shape = NA) +
 #  geom_point(alpha = 0.05) +
  # scale_y_continuous(limits = c(0, -0.05)) +
   #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
   scale_fill_viridis_d("ROH length class", direction = -1,
                        guide = guide_legend(reverse = TRUE),
                        labels=rev(c("long (>6.25cM)", "medium (1.56cM - 6.25cM)", "short (<1.56cM)"))) +
   #scale_x_discrete(labels = c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short]))) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
   scale_y_discrete(limits = c(-0.02, 0)) +
   theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      #axis.title.x = element_blank(),
      # panel.spacing = unit(2, "lines"),
      axis.text = element_text(color = "black"),
      legend.position = "right"
   ) 
