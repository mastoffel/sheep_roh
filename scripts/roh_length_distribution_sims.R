library(data.table)
library(tidyverse)
source("../sheep_ID/theme_simple.R")
library(gghalves)
# roh lenght classes in simulations

# combinations
pop_size1 <- 1000
pop_size2 <- 200
time1 <- pop_size1 * 10
time2 <- time1 + 1000
mut1_dom_coeff <- c(0, 0.05, 0.2)
mut1_gam_mean <- c(-0.01, -0.03, -0.05)
mut1_gam_shape <- 0.2
genome_size <- 1e8
mut_rate_del <- 7e-9
recomb_rate <- 1e-8

params <- expand.grid(pop_size1, pop_size2, mut1_dom_coeff, mut1_gam_mean, mut1_gam_shape,
                      genome_size, mut_rate_del, recomb_rate) %>%
   setNames(c("pop_size1", "pop_size2", "mut1_dom_coeff", "mut1_gam_mean", "mut1_gam_shape",
              "genome_size", "mut_rate_del", "recomb_rate")) %>%
   as_tibble() %>%
   mutate(time1 = 10000,        #pop_size1 * 10,
          time2 = time1 + 1000)

# try 10 simulation with only weakly deleterious alleles
num_sim_per_parset <- 50
set.seed(123)
seeds <- sample(1:1e5, num_sim_per_parset * nrow(params))
# replicate each parameter set num_sim_per_parset times
params_sim <- params[rep(1:nrow(params), each =num_sim_per_parset), ] %>%
   mutate(seed = seeds)

# subset sims
seed_sub <- params_sim %>% 
   filter(mut1_gam_mean == -0.03,
          mut1_dom_coeff == 0.05) %>% 
   select(seed)


roh_paths <- list.files("output/qm_slim/slim1000200_bot_7030del/roh", 
                        pattern = ".hom$", full.names = TRUE)

roh_paths %>% 
   str_split("_") %>% 
   map_chr(~.[5]) %>% 
   str_split("\\.") %>% 
   map_chr(~.[1]) -> seeds

seed_ind <- as.numeric(seeds) %in% seed_sub$seed

roh_paths_sub <- roh_paths[seed_ind]

roh <- map_dfr(roh_paths_sub, function(x) as_tibble(fread(x)), .id = "run")
mean(roh$KB)

IBD_df <- roh %>%
      mutate(length_Mb = KB / 1000) %>%
      mutate(length_class = case_when(length_Mb >= 12.5 ~ "roh_long",
                               length_Mb > 1.5625 &  length_Mb < 12.5 ~ "roh_medium",
                               length_Mb <= 1.5625 ~ "roh_short")) %>% # 0.610683047
      mutate(id = as.character(IID)) %>% 
      group_by(run, id, length_class) %>%
      dplyr::summarise(prop_IBD = sum(length_Mb) / 100,
                       sum_Mb = sum(length_Mb)) #%>% 

IBD_df_with_0 <- IBD_df %>% 
      # add missing length classes as 0
      ungroup() %>% 
      tidyr::complete(length_class, nesting(run, id)) %>% 
      mutate(sum_Mb = ifelse(is.na(sum_Mb), 0, sum_Mb)) %>% 
      mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD))

IBD_df_with_0 %>% arrange(id)

p_dist_sim <- IBD_df_with_0 %>% 
   ungroup() %>% 
   rename(froh = prop_IBD) %>% 
   mutate(froh = 100*froh) %>% 
   mutate(length_class = fct_rev(length_class)) %>% 
   #mutate(prop_IBD = prop_IBD / 10) %>% 
   ggplot(aes(length_class, froh, fill = length_class)) +
   geom_half_point(side = "r", shape = 21, alpha = 0.1, stroke = 0.3, size = 2, color = "#2E3440") +
   geom_half_boxplot(side = "l", outlier.color = NA,
                     width = 0.5, lwd = 0.5, color = "#2E3440",
                     alpha = 0.8) +
   scale_fill_viridis_d(direction = -1, guide = guide_legend(reverse = TRUE)) +
   #ylab("Selection coefficient per cM") +
   scale_x_discrete(labels = rev(c(expression(paste(ROH[long], "(>12.5cM | <4g)")),
                                   expression(paste(ROH[medium], "(1.56 - 12.5cM | 4-32g)")),
                                   expression(paste(ROH[short] ~ "(>1.56g | >32g)"))))) +
   theme_simple(grid_lines = FALSE, axis_lines = TRUE,  base_size = 12) +
   ylab("% simulated 100Mb genome") + 
   #ggtitle("Simulation") +
   theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #axis.line.y = element_blank(),
      #axis.ticks.y = element_blank(),
      #axis.title.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text = element_text(color = "black", size = 12)
   ) +
   coord_flip()
#p_dist_sim

ggsave("figs/roh_dist_sim.jpg", p_dist_sim, width = 5.5, height = 3)





ibd %>% 
   ungroup() %>% 
   #mutate(prop_IBD = prop_IBD / 10) %>% 
   ggplot(aes(three_classes, prop_IBD, fill = three_classes)) +
   geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
                   transformation_params = list(height = 0, width = 1.3, seed = 1)) +
   geom_half_boxplot(side = "r", outlier.color = NA,
                     width = 0.6, lwd = 0.3, color = "black",
                     alpha = 0.8) +
   theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
   ylab("% genome") +
   # scale_y_continuous(labels = c("0", "5", "10"), breaks = c(0, 0.05, 0.1)) +
   scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
   scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
   theme(legend.position = "none",
         #axis.ticks.x = element_blank(),
         axis.title=element_text(size = rel(1.1)), 
         axis.text = element_text(color = "black")) + 
   xlab("ROH classes in cM") 





library(viridis)
col_pal <- viridis(7)
p_roh_dist <- IBD_df_with_0 %>% 
   #mutate(prop_IBD = prop_IBD / 10) %>% 
   ggplot(aes(length_class, prop_IBD, fill = length_class)) +
   geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
                   transformation_params = list(height = 0, width = 1.3, seed = 1)) +
   geom_half_boxplot(side = "r", outlier.color = NA,
                     width = 0.6, lwd = 0.3, color = "black",
                     alpha = 0.8) +
   theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
   ylab("% genome") +
   # scale_y_continuous(labels = c("0", "5", "10"), breaks = c(0, 0.05, 0.1)) +
   scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
   scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
   theme(legend.position = "none",
         #axis.ticks.x = element_blank(),
         axis.title=element_text(size = rel(1.1)), 
         axis.text = element_text(color = "black")) + 
   xlab("ROH classes in cM") 
p_roh_dist 
