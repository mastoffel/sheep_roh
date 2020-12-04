library(ggplot2)
source("../sheep_ID/theme_simple.R")
dat <- data.frame(gamma = -rgamma(20000, 0.2, rate = 6.666))

p1 <- ggplot(dat, aes(gamma)) +
      geom_histogram(bins = 200) +
      scale_x_continuous(limits = c(-0.05, 0)) +
      scale_y_continuous(limits = c(0,750)) +
      xlab("selection coefficient") +
      theme_simple(grid_lines = FALSE) +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank())
ggsave(filename = "figs/gamma_dist_003_02.jpg", p1, width = 3, height = 2.3)
