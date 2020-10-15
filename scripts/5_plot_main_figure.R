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
      # when only showing roh related effs
      #scale_x_log10(breaks = c(0.1, 0.3, 1, 2, 5), limits = c(0.1, 5), labels = c( "0.1", "0.3", "1", "2", "5")) +
     # scale_y_discrete(labels = rev(c(expression(F[ROH]), expression(F[ROH]~'*'~Age), expression(F[ROH]~'*'~Lamb)))) +
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
            axis.text = element_text(color = "black")
      ) +
      xlab("Odds ratio and 95% CI")
p1


# load simulated data
# 6 parameter combinations,100 runs with 200 individuals each
mut_df <- read_delim("output/eddie_slim/par_combs_popsize1_1000.txt", " ")
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
      geom_half_point(side = "r", shape = 21, alpha = 0.8, stroke = 0.1, size =2,
                      transformation_params = list(height = 0, width = 1.3, seed = 1)) +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "black",
                        alpha = 0.8) +
      #scale_y_continuous(limits = c(-0.02, 0)) +
      scale_fill_viridis_d(direction = -1) +
      ylab("Selection coefficient per cM") +
     # xlab("ROH length class") +
      scale_x_discrete(labels = rev(c(expression(ROH[long]),expression(ROH[medium]), expression(ROH[short])))) +
      #facet_grid(dominance ~ selection, labeller = label_both) +
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
      ) +
      coord_flip()
p2


p_final <- p2 + p1 + plot_annotation(tag_levels = "A")
p_final

ggsave("figs/Fig1.jpg", p_final, width = 6, height = 2.6)
