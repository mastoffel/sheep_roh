library(tidyverse)

# 6 parameter combinations,100 runs with 200 individuals each
mut_df <- read_delim("output/eddie_slim/par_combs_popsize1_1000.txt", " ")

mut_df %>% 
      group_by(mut1_dom_coeff) %>% 
      tally()

mut_p <- mut_df %>% 
      group_by(mut1_dom_coeff, mut1_gam_mean, seed,  roh_class) %>% 
      summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE)) %>% 
      mutate(roh_class = factor(roh_class, levels = c("long", "medium",
                                                      "short", "outside_roh")))

# mut_p <- mut_df %>% 
#    group_by(seed, roh_class) %>% 
#    summarise(s_sum_per_MB = mean(s_sum_per_MB, na.rm = TRUE))
# 
source("../sheep_ID/theme_simple.R")
library(gghalves)
# 
p1 <- mut_p %>% 
      rename(selection = mut1_gam_mean,
             dominance = mut1_dom_coeff) %>% 
   ggplot(aes(roh_class, s_sum_per_MB, fill = roh_class)) +
   geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
                   transformation_params = list(height = 0, width = 1.3, seed = 1)) +
   geom_half_boxplot(side = "r", outlier.color = NA,
                     width = 0.6, lwd = 0.3, color = "black",
                     alpha = 0.8) +
   #scale_y_continuous(limits = c(-0.02, 0)) +
   scale_fill_viridis_d() +
   ylab("selection coefficient per cM") +
   facet_grid(dominance ~ selection, labeller = label_both) +
 #  ggtitle("Sim: weakly del\nROH cutoff 4900KB, \nmut1. dist: -0.03, 2, dom coeff 0.1 \nmut2. dist: -0.2, 3, dom coeff 0.01") +
   # geom_jitter(size = 2, alpha = 0.3, width = 0.2) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE)
p1
ggsave("figs/del_roh_popsize1000.jpg", plot = p1, width = 10, height = 4)
# 
# mut_df
# 
# library(lme4)
# 
# mod <- lmer(s_sum_per_MB ~ roh_class + (1|seed), data = mut_df)
# summary(mod)
# tidy(mod, conf.int = TRUE)
# library(broom.mixed)
# 
# tidy(mod, conf.int = TRUE)
