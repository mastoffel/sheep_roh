# this is a simulation to see whther mixed models can estimate effects of 
# correlated predictors.
library(broom.mixed)
library(tidyverse)
# froh correlation
# cor_lm = long, medium
# cor_ls = long, short
# cor_ms = medium, short 

# froh slopes: b_long, b_medium, b_short

sim_froh_mod <- function(iter, cor_lm = -0.15, cor_ls = -0.45, cor_ms = -0.08,
                         b_long = -0.4, b_medium = -0.2, b_short = 0) {
      
      b1 <- b_long
      b2 <- b_medium
      b3 <- b_short
      
      n_years <- 20
      sheep_per_year <- 100
      n_data <- n_years * sheep_per_year
      b0 <- 0.5
      sds <- 2
      sd <- 1
      
      #set.seed(15)
      birth_year <- rep(LETTERS[1:n_years], each = sheep_per_year)
      year_eff <- rep(rnorm(n_years, 0, sds), each = sheep_per_year)
      res <- rnorm(n_data, 0, sd)
      
      # calculated from var(juv_survival$froh_long*100) and mean()
      sigma_long <- 3.58
      sigma_medium <- 1.00
      sigma_short <- 0.288
      mu_long <- 3.67
      mu_medium <- 12
      mu_short <- 8.3
      
      cor_mat <- matrix(c(sigma_long, cor_lm, cor_ls, cor_lm, sigma_medium, cor_ms, cor_ls, cor_ms, sigma_short), ncol = 3)
      vars <- MASS::mvrnorm(n_data, c(mu_long, mu_medium, mu_short), cor_mat)
      froh_long <- vars[, 1]
      froh_medium <- vars[, 2]
      froh_short <- vars[, 3]
      
      # make response variable
      resp <- b0 + b1 * froh_long + b2 * froh_medium + b3 * froh_short + year_eff + res
      resp_bin <- plogis(resp)
      resp_bin <- ifelse(resp_bin >= 0.5, 1, 0)
      
      df <- data.frame(resp, resp_bin, froh_long, froh_medium, froh_short, birth_year) %>% 
            mutate(froh_long_std = (froh_long - mean(froh_long))/sd(froh_long),
                   froh_medium_std = (froh_medium - mean(froh_medium))/sd(froh_medium),
                   froh_short_std = (froh_short - mean(froh_short))/sd(froh_short))
      
      # model gaussian
      mod1 <- lmer(resp ~ froh_long + froh_medium + froh_short + (1|birth_year),
                   data = df)
      mod1_tidy <- tidy(mod1)
      
      mod2 <- lmer(resp ~ froh_long_std + froh_medium_std + froh_short_std + (1|birth_year), 
                   data = df)
      mod2_tidy <- tidy(mod2)
      
      list(mod1_tidy, mod2_tidy)
      # model bin
      # mod <- glmer(resp_bin ~ froh_long + froh_medium + froh_short + (1|birth_year),
      #              family = "binomial", data = df)
      # tidy(mod)
}

all_mods <- map(1:100, sim_froh_mod)

all_mods %>% 
      map(1) %>% 
      bind_rows() %>% 
      filter(str_detect(term, "froh")) %>% 
      ggplot(aes(estimate, fill = term)) +
      geom_histogram(bins = 50) #+
      #facet_wrap(~term, scale = "free")

# mod2 <- glmer(resp_bin ~ froh_long + froh_medium + froh_short + (1|birth_year),
#               family = "binomial", data = df)
# summary(mod2)
