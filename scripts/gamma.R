library(ggplot2)
source("../sheep_ID/theme_simple.R")
library(simstudy)
gammaGetShapeRate(0.05, 5)
gammaGetShapeRate(0.03, 5)
gammaGetShapeRate(0.01, 5)
gammaGetShapeRate(0.001, 5)

hist(-rgamma(10000, 0.2, rate = 6.666), breaks = 1000, xlim = c(0,-0.2),
     main = "Gamma distribution with mean -0.03 and shape 0.2",
     xlab = "Selection coefficient s")

muts <- rgamma(10000, 0.2, rate = 6.666) * 200
muts <- rgamma(100000, 0.2, rate = 200) * 200
mean(muts)
median(muts)

sum(muts < 1)/10000
sum(muts >= 1 & muts <= 10)/10000
sum(muts >= 10)/10000

dat <- data.frame(gamma = -rgamma(20000, 0.2, rate = 6.666),
                  main = "Gamma")

p1 <- ggplot(dat, aes(gamma)) +
      geom_histogram(bins = 200) +
      scale_x_continuous(limits = c(-0.05, 0)) +
      scale_y_continuous(limits = c(0,750)) +
      xlab("selection coefficient") +
      theme_simple(grid_lines = FALSE) +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank())
ggsave(filename = "figs/gamma_dist_003_02.jpg", p1, width = 3, height = 2.3)

