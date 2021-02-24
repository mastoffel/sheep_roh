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
     # geom_line(mapping = aes(group = seed), size = 0.2, alpha = 0.1,
      #          color = "#4C566A") +
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
         legend.position = "none"
      ) +
      coord_flip()
   p
}

p1 <- make_plot(s_sum_per_MB,"Selection coefficient per cM", mut_sub)  
p2 <- make_plot(mean_freq,"Mean mutation frequency", mut_sub)  
p3 <- make_plot(num_mut_per_MB, "Number of mutations per cM", mut_sub)
#p4 <- make_plot(mean_origin,"Mutation age in generations", mut_sub)
p_final <- p1 + p3 + p2 + 
   plot_annotation(tag_levels = 'a') +
   plot_layout(guides = 'collect') & theme(legend.position = 'bottom') 
 
p_final

ggsave(plot = p_final, filename = "figs/Fig2_gen4_32_7030.jpg", width = 8, height = 3)


# check 
mut_sub %>% 
   group_by(roh_class) %>% 
   summarise(mean = across(c(s_sum_per_MB, num_mut_per_MB), mean))




p <- ggplot(mut_sub, aes(roh_class, s_sum_per_MB, fill = roh_class)) +
   # geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, 
   #                 size = 2,  width = 0.5) +
   geom_beeswarm() +
   # geom_half_boxplot(#mapping = aes(fill = lighten(fill, .5)),
   #                                 side = "l", outlier.color = NA,
   #                   width = 0.5, lwd = 0.5, color = "black",
   #                   alpha = 0.8, notch = FALSE) +
   scale_fill_viridis_d("ROH length class", direction = -1, option = "D",
                        guide = guide_legend(reverse = TRUE),
                        labels=rev(c("long (>12.5cM | < 4g)", 
                                     "medium (>1.56cM | 4-32g)", 
                                     "short (>0.39cM | 32-128g)"))) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
   #ylab(axis_title) +
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




# how much is purged?

# selection coefficient
mut_sub %>% 
   group_by(roh_class) %>% 
   summarise(mean(s_sum_per_MB))

mut_sub %>% 
   group_by(roh_class) %>% 
   summarise(mean(num_mut_per_MB))

# Supplementary plot with all s and h
mut_sub2 <- mut_p %>%
   #filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
   mutate(mean_origin = 11000 - mean_origin) %>%
   rename(s = mut1_gam_mean,
          h = mut1_dom_coeff)
  # pivot_longer(cols = s_sum_per_MB:mean_freq, names_to = "feature")
library(scales)
make_plot2 <- function(ss, axis_title, df) {
   ss <- ensym(ss)
   ggplot(df, aes(roh_class, !!ss, fill = roh_class)) +
      #geom_boxploth() +
      geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, size = 2) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "black",
                        alpha = 0.8, notch = FALSE) +
      facet_wrap(s ~ h, scales = "free_x",
                 labeller = "label_both") +
      scale_fill_viridis_d("ROH length class", direction = -1, option = "D",
                           guide = guide_legend(reverse = TRUE),
                           labels=rev(c("long (>12.5cM | < 4g)", 
                                        "medium (>1.56cM | 4-32g)", 
                                        "short (>0.39cM | 32-128g)"))) +
      scale_y_continuous(breaks = pretty_breaks(n = 4)) +
      theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
      ylab(axis_title) +
      theme(
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(), 
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         #axis.title.x = element_blank(),
         axis.text.x = element_text(size = 8),
         panel.spacing = unit(2, "lines"),
         axis.text = element_text(color = "black"),
         legend.position = "right"
      ) +
      coord_flip()
}

p1 <- make_plot2(s_sum_per_MB, "Selection coefficient per cM", mut_sub2)
p2 <- make_plot2(num_mut_per_MB, "Number of mutations per cM", mut_sub2)
p3 <- make_plot2(mean_freq, "Mutation frequency", mut_sub2)
p4 <- make_plot2(mean_origin, "Mutation age in generations", mut_sub2)

ggsave("figs/SupSelCoef.jpg",plot = p1, width = 8.5, height = 6)
ggsave("figs/MutFreq.jpg",plot = p3, width = 8.5, height = 6)
ggsave("figs/NumMut.jpg",plot = p2, width = 8.5, height = 6)
ggsave("figs/MutOrigin.jpg",plot = p4, width = 8.5, height = 6)




library(ggforce)
library(ggbeeswarm)
make_plot <- function(ss, axis_title, df) {
   ss <- ensym(ss)
   p <- ggplot(df, aes(roh_class, !!ss, fill = roh_class)) +
      # geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1,
      #                 size = 2,  width = 0.5) +
      geom_beeswarm(beeswarmArgs=list(side=1)) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "black",
                        alpha = 0.8, notch = FALSE) +
      #geom_sina() +
      scale_color_viridis_d("ROH length class", direction = -1, option = "D",
                           guide = guide_legend(reverse = TRUE),
                           labels=rev(c("long (>12.5cM | < 4g)", 
                                        "medium (>1.56cM | 4-32g)", 
                                        "short (>0.39cM | 32-128g)"))) +
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
   plot_layout(guides = 'collect') +
   plot_annotation(tag_levels = 'a')
p_final

ggplot(mut_sub, aes(roh_class, s_sum_per_MB, fill = roh_class)) +
   # geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1,
   #                 size = 2,  width = 0.5) +
   geom_half_dotplot() + 
   geom_half_boxplot(side = "l", outlier.color = NA,
                     width = 0.5, lwd = 0.5, color = "black",
                     alpha = 0.8, notch = FALSE) 








# individual scale
mut_p <- mut_all %>% 
   filter(roh_class != "outside_roh") %>% 
   mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short"))))

# combine 
mut_sub <- mut_p %>%
   filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
   mutate(mean_origin = 11000 - mean_originG)



make_plot <- function(ss, axis_title, df) {
   ss <- ensym(ss)
   p <- ggplot(df, aes(roh_class, !!ss, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, 
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
   plot_layout(guides = 'collect') +
   plot_annotation(tag_levels = 'a')
p_final





# only highly deleterious muts

mut_all <- mut_df %>% 
   #sample_frac(0.001) %>% 
   # filter homozygous sites
   filter(copies == 2) %>%
   #filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
   filter(s > -0.1) %>% 
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
          num_mut_per_MB = (num_mut / roh_class_genome_cov) * 1000)

mut_p <- mut_all %>% 
   filter(roh_class != "outside_roh") %>% 
   group_by(seed, mut1_dom_coeff, mut1_gam_mean, roh_class) %>% 
   summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE),
             num_mut_per_MB = mean(num_mut_per_MB, na.rm = TRUE),
             mean_origin = mean(mean_originG, na.rm = TRUE),
             mean_freq = mean(mean_freq, na.rm = TRUE)) %>% 
   mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium","short")))) %>% 
   ungroup()

mut_sub <- mut_p %>%
   filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05) %>% 
   mutate(mean_origin = 11000 - mean_origin)

make_plot <- function(ss, axis_title, df) {
   ss <- ensym(ss)
   p <- ggplot(df, aes(roh_class, !!ss, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.5, stroke = 0.1, 
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
   plot_layout(guides = 'collect') +
   plot_annotation(tag_levels = 'a')
p_final

ggsave(plot = p_final, filename = "figs/Fig2_gen4_32_7030_s_under_01.jpg", width = 7, height = 4.5)
