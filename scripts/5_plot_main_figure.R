library(tidyverse)
library(broom.mixed)
source("../sheep_ID/theme_simple.R")
library(gghalves)
library(janitor)
library(sjPlot)
library(patchwork)

# first year survival model
mod <- readRDS("output/juv_survival_model.RDS")
summary(mod)

# CIs
mod_sum <- tidy(mod, conf.int = TRUE) #  conf.method = "boot", nsim = 1000

mod_sum_filt <- mod_sum %>% 
      clean_names() %>% 
      filter(str_detect(term, "froh")) %>% 
      mutate(across(.cols = c("estimate", "conf_low", "conf_high"), exp)) %>% 
      mutate(term = fct_rev(factor(term)))


p1 <- ggplot(mod_sum_filt, aes(estimate, term, xmax = conf_high, xmin = conf_low, color = term, fill = term)) +
      geom_vline(xintercept = 1, linetype='dashed', size = 0.3) +  #colour =  "#4c566a",   # "#4c566a"  "#eceff4"
      geom_errorbarh(alpha = 1, height = 0, size = 1) +
      geom_point(size = 3.5, shape = 21, #col = "#4c566a", #fill = "#eceff4", # "grey69"
                 alpha = 1, stroke = 0.3) + 
      scale_y_discrete(labels = rev(c(expression(F[ROH[long]]),expression(F[ROH[medium]]), expression(F[ROH[short]])))) +
     # ggtitle("Empirical analysis") +
      theme_simple(axis_lines = TRUE, base_size = 12) +
      scale_fill_viridis_d(direction = -1) +
      scale_color_viridis_d(direction = -1) +
      theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            axis.text = element_text(color = "black", size = 12),
            axis.title.x = element_text(size = 12)
      ) +
      xlab("Odds ratio and 95% CI")
p1


# load simulated data
# 6 parameter combinations,100 runs with 200 individuals each
mut_df <- read_delim("output/qm_slim/slim5000200/out/par_combs_popsize1_5000_popsize2_200.txt", " ")
# filter out one parameter combination
mut_p <- mut_df %>% 
      filter(mut1_gam_mean == -0.03, mut1_dom_coeff == 0.05, roh_class != "outside_roh") %>% 
      group_by(mut1_dom_coeff, mut1_gam_mean, seed,  roh_class) %>% 
      summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE)) %>% 
      mutate(roh_class = factor(roh_class, levels = rev(c("long", "medium",
                                                      "short"))))

p2 <- mut_p %>% 
      rename(selection = mut1_gam_mean,
             dominance = mut1_dom_coeff) %>% 
      ggplot(aes(roh_class, s_sum_per_MB, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.6, stroke = 0.3, size =2, color = "#2E3440",
                      transformation_params = list(height = 0, width = 1.3, seed = 1)) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "#2E3440",
                        alpha = 0.8) +
      scale_fill_viridis_d(direction = -1) +
      ylab("Selection coefficient per cM") +
      scale_x_discrete(labels = rev(c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short])))) +
      theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
      #ggtitle("Simulation") +
      theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none",
            axis.title.x = element_text(size = 12),
            axis.text = element_text(color = "black", size = 12)
      ) +
      coord_flip()
p2


p_final <- p2 + p1 + plot_annotation(tag_levels = "a")
p_final

ggsave("figs/Fig1.jpg", p_final, width = 7, height = 3)
