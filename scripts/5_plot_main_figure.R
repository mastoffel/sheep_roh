library(tidyverse)
library(broom.mixed)
source("../sheep_ID/theme_simple.R")
library(gghalves)
library(janitor)
library(sjPlot)
library(patchwork)
library(tidybayes)
library(bayesplot)
library(brms)
load("data/fitness_roh.RData") 
library(gt)
# first year survival model
mod <- readRDS("output/juv_survival_model_nonstd_brm.RDS")
summary(mod)

juv_survival <- fitness_data %>% 
   filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
   mutate(age_std = as.numeric(scale(age)),
          age_std2 = age_std^2,
          froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
          froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
          froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
          froh_all_std = as.numeric(scale(froh_all)),
          froh_short_std = as.numeric(scale(froh_short)),
          froh_long_std = as.numeric(scale(froh_long)),
          # hom_std = as.numeric(scale(hom1_all))) %>% 
          froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
          froh_medium_std = scale(froh_medium)) %>% 
   filter(age == 0)

# distribution 
p_dist <- juv_survival %>% 
   pivot_longer(cols = froh_long:froh_short, names_to = "roh_class", 
                values_to = "froh") %>% 
   mutate(roh_class = factor(roh_class, levels = rev(c("froh_long", "froh_medium", "froh_short")))) %>% 
   select(roh_class, froh) %>% 
   mutate(froh = 100*froh) %>% 
   ggplot(aes(roh_class, froh, fill = roh_class)) +
      geom_half_point(side = "r", shape = 21, alpha = 0.1, stroke = 0.3, size = 2, color = "#2E3440") +
      geom_half_boxplot(side = "l", outlier.color = NA,
                        width = 0.5, lwd = 0.5, color = "#2E3440",
                        alpha = 0.8) +
      scale_fill_viridis_d(direction = -1) +
      #ylab("Selection coefficient per cM") +
      scale_x_discrete(labels = rev(c(expression(paste(F[ROH[long]], "(>12.5cM | <4g)")),
                                      expression(paste(F[ROH[medium]], "(1.56 - 12.5cM | 4-32g)")),
                                      expression(paste(F[ROH[short]] ~ "(>1.56g | >32g)"))))) +
      theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
      ylab("% genome") + 
      #ggtitle("Simulation") +
      theme(
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         legend.position = "none",
         axis.title.x = element_text(size = 12),
         axis.text = element_text(color = "black", size = 12)
      ) +
      coord_flip()
p_dist


# 
# brm model
post <- posterior_samples(mod, pars = c("froh_long" ,"froh_medium", "froh_short"))

p_emp <- post %>% 
   pivot_longer(cols = starts_with("b"), names_to = "froh", values_to = "beta") %>% 
   mutate(beta = exp(beta),
          froh =  fct_rev(factor(froh))) %>% 
   ggplot(aes(beta, froh, color = froh)) +
   geom_vline(xintercept = 1, linetype='dashed', size = 0.3) +
   stat_slab(alpha = 0.6, slab_type = "pdf", height = 0.7, size = 0.4,  fill = "#ECEFF4") +
   stat_pointinterval(.width = c(.95), point_size = 3) +
   scale_x_continuous(breaks = c(.8, 1, 1.2, 1.4)) + 
   scale_y_discrete(labels = rev(c(expression(atop(F[ROH[long]], "(>12.5cM | <4g)")),
                                   expression(atop(F[ROH[medium]], "(1.56-12.5cM | 4-32g)")),
                                   expression(atop(F[ROH[short]], "(<1.56cM | >32g)"))))) +
   # ggtitle("Empirical analysis") +
   theme_simple(axis_lines = TRUE, base_size = 12) +
   #scale_fill_viridis_d(direction = -1) +
   scale_color_viridis_d(direction = -1) +
   theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(hjust=0.5),
      legend.position = "none",
      axis.text = element_text(color = "black", size = 12),
      axis.title.x = element_text(size = 12)
   ) +
   xlab("Odds-ratio and 95% CI\n(First-year survival)")


library(patchwork)
p_final <- p_dist + p_emp + plot_annotation(tag_levels = "a") +
   plot_layout(widths = c(1.5, 2)) +
   theme(plot.tag.position = c(0.36, 0.94),
         plot.tag = element_text( hjust = 0, vjust = 0))
p_final

ggsave("figs/Fig1_bayes.jpg", p_final, width = 7, height = 3)


# figure of differences
p_sup <- post %>% 
   mutate(d_long_medium = b_froh_long - b_froh_medium,
          d_long_short = b_froh_long - b_froh_short,
          d_medium_short = b_froh_medium - b_froh_short) %>% 
   pivot_longer(cols = starts_with("d"), names_to = "froh", values_to = "beta") %>% 
   mutate(froh =  fct_rev(factor(froh))) %>% 
   ggplot(aes(beta, froh, fill = froh)) +
   geom_vline(xintercept = 0, linetype='dashed', size = 0.3 ) +
   stat_slab(alpha = 0.6, slab_type = "pdf", height = 0.7, size = 0.4, color = "#3B4252", fill = "#ECEFF4") +
   stat_pointinterval(.width = c(.95), point_size = 3, color = "#3B4252") +
   #scale_x_continuous(limits = c(0.73, 1.47)) + 
   scale_y_discrete(labels = rev(c(expression(beta[F[ROHlong]] - beta[F[ROHmedium]]),
                                   expression(beta[F[ROHlong]] - beta[F[ROHshort]]),
                                   expression(beta[F[ROHmedium]] - beta[F[ROHshort]])))) +
   # ggtitle("Empirical analysis") +
   theme_simple(axis_lines = TRUE, base_size = 12) +
   theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(hjust=0.5),
      legend.position = "none",
      axis.text = element_text(color = "black", size = 12),
      axis.title.x = element_text(size = 12)
   ) +
   xlab(expression(Delta~log~odds~ratio~and~'95%'~CI)) # "log odds-ratio and 95% CI"

ggsave("figs/Sup_fig1_bayes_diff.jpg", p_sup, width = 4, height = 2.7)


# TABLE
pred_names = c(Intercept = "Intercept", froh_long = "FROHlong", 
               froh_medium = "FROHmedium", froh_short = "FROHshort",
               sexM = "Sex [0=Male, 1=Female]", twin1 = "Twin [0=singleton, 1=twin]", 
               birth_year = "Birth year", 
               mum_id = "Mother ID")

tab_model(mod, pred.labels = pred_names,
          show.r2 = FALSE, show.icc = FALSE)


tidy_mod <- tidy(mod) %>% 
               select(term:conf.high) %>% 
               mutate(term = c("Intercept", 
                               "F<sub>ROHlong</sub>",
                               "F<sub>ROHmedium</sub>",
                               "F<sub>ROHshort</sub>",
                               "Sex",
                               "Twin",
                               "Birth year",
                               "Mother ID")) %>% 
               mutate(add_info = c("", rep("continuous", 3), 
                                   "categorical (0=male, 1=female)",
                                   "categorical (0=singleton, 1=twin)",
                                   "n = 1118",
                                   "n = 39")) %>% 
               mutate(effect = c("", rep("Population level/fixed effects", 5),
                                 rep("Group level/random effects (standard deviation)", 2))) %>% 
               select(7, 1:6) %>% 
               mutate_if(is.numeric, function(x) paste0(round(x, 3)," (",round(exp(x), 3), ")")) %>% 
               setNames(c("effect", "Term", "Post.Mean", "Std.Error", "CI (2.5%)", "CI (97.5%)", "Info"))
                                    
tidy_mod %>% gt(
   rowname_col = "term",
   groupname_col = "effect") %>% 
   tab_style(
      style = cell_text( weight = "bold"),
      locations = cells_column_labels(columns = TRUE)
   ) %>% 
   fmt_markdown(columns = TRUE) %>% 
   gtsave("JS_model_table.png", path = "tables/")









# Model and CIs
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
      # scale_y_discrete(labels = rev(c(expression(F[ROH[long]]~"">12.5~cM),
      #                                 expression(F[ROH[medium]]~1.56 - 12.5~cM),
      #                                 expression(F[ROH[short]]~""<1.56~cM)))) +
      # scale_y_discrete(labels = rev(c(expression(paste(F[ROH[long]], "(>12.5cM | <4g)")),
      #                             expression(paste(F[ROH[medium]], "(1.56 - 12.5cM | 4-32g)")),
      #                             expression(paste(F[ROH[short]] ~ "(>1.56g | >32g)"))))) +
      scale_y_discrete(labels = rev(c(expression(atop(F[ROH[long]], "(>12.5cM | <4g)")),
                                   expression(atop(F[ROH[medium]], "(1.56-12.5cM | 4-32g)")),
                                   expression(atop(F[ROH[short]], "(<1.56cM | >32g)"))))) +
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
            axis.text.y = element_text(hjust=0.5),
            legend.position = "none",
            axis.text = element_text(color = "black", size = 12),
            axis.title.x = element_text(size = 12)
      ) +
      xlab("Odds-ratio and 95% CI")
p1
ggsave("figs/Fig1_simple.jpg", p1, width = 6, height = 3)

library(patchwork)
p_final <- p_dist + p1 + plot_annotation(tag_levels = "a") +
   theme(plot.tag.position = c(0.36, 0.94),
         plot.tag = element_text( hjust = 0, vjust = 0))
p_final

ggsave("figs/Fig1.jpg", p_final, width = 7, height = 3)



# brm model
post <- posterior_samples(mod, pars = c("froh_long" ,"froh_medium", "froh_short"))
post %>% 
   pivot_longer(cols = starts_with("b"), names_to = "froh", values_to = "beta") %>% 
   mutate(beta = exp(beta),
          froh =  fct_rev(factor(froh))) %>% 
   ggplot(aes(beta, froh, fill = froh)) +
   geom_vline(xintercept = 1, linetype='dashed', size = 0.3) +
   stat_slab(alpha = 0.6, slab_type = "pdf", height = 0.5) +
   stat_pointinterval(.width = c(.50,.90)) +
   scale_x_continuous(limits = c(0.73, 1.47)) + 
   scale_y_discrete(labels = rev(c(expression(atop(F[ROH[long]], "(>12.5cM | <4g)")),
                                   expression(atop(F[ROH[medium]], "(1.56-12.5cM | 4-32g)")),
                                   expression(atop(F[ROH[short]], "(<1.56cM | >32g)"))))) +
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
      axis.text.y = element_text(hjust=0.5),
      legend.position = "none",
      axis.text = element_text(color = "black", size = 12),
      axis.title.x = element_text(size = 12)
   ) +
   xlab("Odds-ratio and 95% CI")

post %>% 
   mutate(d_long_medium = b_froh_long - b_froh_medium,
          d_long_short = b_froh_long - b_froh_short,
          d_medium_short = b_froh_medium - b_froh_short) %>% 
   pivot_longer(cols = starts_with("d"), names_to = "froh", values_to = "beta") %>% 
   mutate(beta = exp(beta),
          froh =  fct_rev(factor(froh))) %>% 
   ggplot(aes(beta, froh, fill = froh)) +
   geom_vline(xintercept = 1, linetype='dashed', size = 0.3) +
   stat_slab(alpha = 0.6, slab_type = "pdf", height = 0.5) +
   stat_pointinterval(.width = c(.50,.90)) +
   scale_x_continuous(limits = c(0.73, 1.47)) + 
   scale_y_discrete(labels = rev(c(expression(atop(F[ROH[long]], "(>12.5cM | <4g)")),
                                   expression(atop(F[ROH[medium]], "(1.56-12.5cM | 4-32g)")),
                                   expression(atop(F[ROH[short]], "(<1.56cM | >32g)"))))) +
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
      axis.text.y = element_text(hjust=0.5),
      legend.position = "none",
      axis.text = element_text(color = "black", size = 12),
      axis.title.x = element_text(size = 12)
   ) +
   xlab("Odds-ratio and 95% CI")












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
