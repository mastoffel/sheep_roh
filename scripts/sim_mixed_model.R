library(lme4)
library(broom.mixed)

# froh correlation
# cor_lm = long, medium
# cor_ls = long, short
# cor_ms = medium, short 

# froh slopes: b_long, b_medium, b_short

sim_froh_mod <- function(iter, cor_lm = -0.3, cor_ls = -0.45, cor_ms = -0.03,
                         b_long = -0.4, b_medium = -0.2, b_short = -0.1) {
      
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
      
      cor_mat <- matrix(c(1, cor_lm, cor_ls, cor_lm, 1, cor_ms, cor_ls, cor_ms, 1), ncol = 3)
      vars <- MASS::mvrnorm(n_data,c(0,0,0), cor_mat)
      froh_long <- vars[, 1]
      froh_medium <- vars[, 2]
      froh_short <- vars[, 3]
      
      # make response variable
      resp <- b0 + b1 * froh_long + b2 * froh_medium + b3 * froh_short + year_eff + res
      resp_bin <- plogis(resp)
      resp_bin <- ifelse(resp_bin >= 0.5, 1, 0)
      
      df <- data.frame(resp, resp_bin, froh_long, froh_medium, froh_short, birth_year)
      
      # model gaussian
      mod <- lmer(resp ~ froh_long + froh_medium + froh_short + (1|birth_year), data = df)
      tidy(mod)
      
      mod <- lmer(resp ~ scale(froh_long) + scale(froh_medium) + scale(froh_short) + (1|birth_year), data = df)
      tidy(mod)
      
      # model bin
      # mod <- glmer(resp_bin ~ froh_long + froh_medium + froh_short + (1|birth_year),
      #              family = "binomial", data = df)
      # tidy(mod)
}

all_mods <- map_dfr(1:100, sim_froh_mod, cor_lm = -0.3, cor_ls = -0.5, cor_ms = -0.03,
                                         b_long = -0.4, b_medium = -0.2, b_short = 0, .id = "sim")

all_mods %>% 
      filter(term %in% c("froh_long", "froh_medium", "froh_short")) %>% 
      ggplot(aes(estimate)) +
      geom_histogram() +
      facet_wrap(~term, scale = "free")

# mod2 <- glmer(resp_bin ~ froh_long + froh_medium + froh_short + (1|birth_year),
#               family = "binomial", data = df)
# summary(mod2)
