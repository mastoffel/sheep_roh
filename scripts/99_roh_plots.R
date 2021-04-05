library(data.table)
library(tidyverse)
library(ggplot2)
library(viridis)
# roh
roh <- fread("output/ROH/roh_cM.hom") %>%
      rename(id = IID) %>% 
      rename(cM = KB) %>%
      mutate(cM = cM/1e3) 
max(roh$cM)
mean(roh$cM)

# linkage map
lmap <- read_delim("../sheep_ID/data/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
      rename(snp = SNP.Name,
             chr = Chr) %>% 
      select(chr, snp, cMPosition, cMPosition.Female, cMPosition.Male) %>% 
      filter(chr %in% 1:26)

full_length_map <- lmap %>% group_by(chr) %>% 
      summarise(max_per_chr = max(cMPosition)) %>% 
      summarise(genetic_map_length = sum(max_per_chr)) %>% 
      .[[1]]


# ROH classes ------------------------------------------------------------------

# expected ROH length using cM/Mb from Johnston et al (2016)
length_dist <- data.frame(g = c(1, 2,2^2, 2^3, 2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13)) %>%
      mutate(ROH_length_cM = 100 / (2*g))

IBD_df <- roh %>%
      mutate(length_cM = cM) %>%
      mutate(class = case_when(length_cM >= 25 ~ 1,
                               length_cM < 25 &      length_cM >= 12.5 ~ 2,
                               length_cM < 12.5 &    length_cM >= 6.25 ~ 4,
                               length_cM < 6.25 &    length_cM >=  3.125 ~ 8,
                               length_cM <  3.125 &  length_cM >=1.5625 ~ 16,
                               length_cM < 1.5625 &  length_cM >= 0.78125 ~ 32,
                               length_cM < 0.78125 ~ 64)) %>% # 0.610683047
      mutate(length_class = case_when(
            class == 1 ~ "> 25 (2g)",
            class == 2 ~ "25-12.5 (2-4g)",
            class == 4 ~ "12.5-6.25 (4-8g)",
            class == 8 ~ "6.25-3.13 (8-16g)",
            class == 16 ~ "3.13-1.56 (16-32g)",
            class == 32 ~ "1.56-0.78 (32-64g)",
            class == 64 ~ "0.78-0.39 (64-128g)"
            # class == 64 ~ ">0.6-1.2 (>32-64g)"
            # class == 128 ~ "0.6-0.3 (128G)"
      )) %>% 
      mutate(length_class = fct_reorder(length_class, class)) %>% 
      mutate(id = as.character(id)) %>% 
      group_by(id, class, length_class) %>%
      dplyr::summarise(prop_IBD = sum(length_cM) / full_length_map,
                       sum_cM = sum(length_cM)) #%>% 

IBD_df_with_0 <- IBD_df %>% 
      # add missing length classes as 0
      ungroup() %>% 
      tidyr::complete(length_class, nesting(id)) %>% 
      mutate(class = ifelse(is.na(class), length_class, class)) %>% 
      mutate(sum_cM = ifelse(is.na(sum_cM), 0, sum_cM)) %>% 
      mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD))

IBD_df_with_0 %>% arrange(id)


col_pal <- viridis(7)
p_roh_dist <- IBD_df_with_0 %>% 
      #mutate(prop_IBD = prop_IBD / 10) %>% 
      ggplot(aes(length_class, prop_IBD, fill = length_class)) +
      geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2) +
      geom_half_boxplot(side = "r", outlier.color = NA,
                        width = 0.6, lwd = 0.3, color = "black",
                        alpha = 0.8) +
      theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
      ylab("% genome") +
      scale_y_continuous(labels = c("0", "5", "10"), breaks = c(0, 0.05, 0.1)) +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
      theme(legend.position = "none",
            #axis.ticks.x = element_blank(),
            axis.title=element_text(size = rel(1.1)), 
            axis.text = element_text(color = "black")) + 
      xlab("ROH classes in cM") 
p_roh_dist 



# roh for some individuals
all_roh <- roh %>% 
      group_by(id) %>% 
      summarise(sum_roh = sum(cM)) %>% 
      ungroup() %>% 
      arrange(desc(sum_roh))

longest_roh <- all_roh %>% 
      top_n(3)
shortest_roh <- all_roh %>% 
      top_n(-3)
num_ind <- 6

extreme_roh <- rbind(longest_roh, shortest_roh)

df <- roh %>%
      mutate(POS1 = POS1 / 1e+6,
             POS2 = POS2 / 1e+6)

df <- df %>% filter(id %in% extreme_roh$id) %>% 
      mutate(id = factor(id, levels = extreme_roh$id))

yax <- data.frame(id = fct_inorder(levels(df$id))) %>%
      mutate(yax = seq(from = 2,
                       to = 2*length(unique(df$id)),
                       by = 2)) 

df <- left_join(df, yax, by = "id")

shade <- df %>%
      group_by(CHR) %>%
      summarise(min = min(POS1), max = max(POS2)) %>%
      mutate(min = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
                                   CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
                                   CHR == 22 | CHR == 24 | CHR == 26 ~ 0,
                             TRUE ~ min)) %>%
      mutate(max = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
                                   CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
                                   CHR == 22 | CHR == 24 | CHR == 26 ~ 0,
                             TRUE ~ max))

col <- c("#4c566a", "#d8dee9")
col <- c("#1E3231", "#9aadbf")
chr_names <- as.character(1:26)
names(chr_names) <- as.character(1:26)
chr_names[c(11, 13, 15, 17, 19, 21, 23, 25)] <- ""

df %>% 
      #filter(cM > 1) %>% 
      filter(CHR %in% 3) %>% 
      mutate(length_class = case_when(
         cM > 12.5 ~ "long",
         cM < 12.5 & cM > 1.56 ~ "medium",
         TRUE ~ "short"
      )) %>% 
      ggplot() +
      #geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
      #          alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
      geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.4) +
      geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.4, fill = length_class), 
                   col = "grey", size = 0.2, alpha = 1) +   #  fill = "#E5E9F0"
      scale_fill_viridis_d() +
      #scale_fill_manual(values = rep(col, 18)) + 
     # scale_color_manual(values = rep(col, 18)) +
      scale_y_reverse(expand = c(0, 0)) +
      theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
     # facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
       #          labeller = as_labeller(chr_names)) +
      theme(#strip.placement = 'outside',
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.spacing = unit(0, "lines"),
            plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
            axis.line.x = element_blank(),
            legend.position="none",
           # axis.title.x = element_text(margin=margin(t=0)),
            axis.title.y = element_text(margin=margin(r=0)),
            axis.text.y = element_blank(),
            axis.line.y = element_blank()) +
      coord_cartesian(clip = 'off') +
      xlab("ROH along chromosome 3") +
      ylab("Individuals") -> ROH_per_ind
ROH_per_ind

ggsave("figs/roh_per_individual_color_talk.jpg",ROH_per_ind, width = 9, height = 2.5)
ggsave("figs/roh_per_individual_color_talk.jpg",ROH_per_ind, width = 5, height = 2.5)
