library(tidyverse)
library(lme4)
library(GGally)
library(partR2)
library(broom.mixed)
library(performance)
library(brms)
library(ggeffects)
library(GGally)
library(sjPlot)
library(optiSel)
load("data/fitness_roh.RData") 
#load("data/sheep_ped.RData")
ped <- sheep_ped

# check for froh across years
# fitness_data %>% 
#       select(starts_with("froh"), birth_year) %>% 
#       pivot_longer(cols = c(starts_with("froh"))) %>% 
#       mutate(birth_year = as.numeric(as.character(birth_year))) %>% 
#       ggplot(aes(birth_year, value)) +
#       geom_point() +
#       geom_smooth(method = "lm") +
#       facet_wrap(~name, scales = "free_y")

# survival data preprocessing
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

#ggpairs(juv_survival, columns = 9:12)

# time saver function for modeling
nlopt <- function(par, fn, lower, upper, control) {
      .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper,
                                opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                            maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
      list(par = res$solution,
           fval = res$objective,
           conv = if (res$status > 0) 0 else res$status,
           message = res$message
      )
}
# froh_long + froh_short + 


# no std:
juv_survival2 <- juv_survival %>% 
   mutate(froh_long = froh_long * 100,
          froh_medium = froh_medium*100,
          froh_short = froh_short *100,
          froh_all = (froh_all * 100))

m2 <- glmer(survival ~ froh_long + froh_medium + froh_short + sex + twin + (1|mum_id) + (1|birth_year), 
            family = binomial, data = juv_survival2,
            control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
mod_tidy <- tidy(m2, conf.int = TRUE)
mod_tidy

saveRDS(m2, "output/juv_survival_model_nonstd.RDS")

# jarrod suggestion
# m3 <- glmer(survival ~ froh_all + mean_cM + sex + twin + (1|mum_id) + (1|birth_year), 
#             family = binomial, data = juv_survival2,
#             #control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
#             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
# mod_tidy <- tidy(m3, conf.int = TRUE)
# mod_tidy

# brms
juv_survival2 <- juv_survival %>% 
   mutate(froh_long = froh_long * 100,
          froh_medium = froh_medium*100,
          froh_short = froh_short *100,
          froh_all = (froh_all * 100))
brm_fit <- brm(survival ~ froh_long + froh_medium + froh_short  + sex + twin + (1|mum_id) + (1|birth_year), 
               family = bernoulli(), data = juv_survival2, cores = 4, iter = 10000,
               set_prior("normal(0,5)", class = "b"))

saveRDS(brm_fit, "output/juv_survival_model_nonstd_brm.RDS")
brm_fit <- readRDS("output/juv_survival_model_nonstd_brm.RDS")
prior_summary(brm_fit)
summary(brm_fit)
plot(brm_fit, ask = FALSE)
loo(brm_fit)
pp_check(brm_fit, nsamples = 100 )
summary(mod_brm)
conditional_effectnsamples = s(mod_brm)


# suggestion jarrod
brm_fit2 <- brm(survival ~ froh_all + ROH_len_cM + sex + twin + (1|mum_id) + (1|birth_year), 
               family = bernoulli(), data = juv_survival2, cores = 4, iter = 10000,
               set_prior("normal(0,5)", class = "b"))
summary(brm_fit2)
saveRDS(brm_fit2, "output/juv_survival_model_nonstd_brm_roh_all_length.RDS")
brm_fit2 <- readRDS("output/juv_survival_model_nonstd_brm_roh_all_length.RDS")



# trying out things here
# 
# plot(mod_brm)
# ggpredict(mod_brm, c("froh_long_std [all]")) 
# brms::conditional_effects(mod_brm)
# me <- brms::conditional_effects(mod_brm, "froh_long_std")
# summary(mod_brm)
# plot_model(mod_brm, terms = c("froh_long", "froh_medium", "froh_short"))
# 
# library(tidybayes)
# post1 <- posterior_samples(mod_brm, pars = c("froh_long" ,"froh_medium", "froh_short"))
# post <- spread_draws(mod_brm, )
# post1 %>% 
#    ggplot(aes)
# 
# library(bayesplot)
# mcmc_areas(post1, prob = 0.9)
# 
# library(brmstools)
# forest(mod_brm, pars =  c("froh_long" ,"froh_medium", "froh_short"))
# 
# library(ggridges)
# 
# mod_brm_summary <- post1 %>% 
#    pivot_longer(cols = starts_with("b"), names_to = "froh", values_to = "beta") %>% 
#    group_by(froh) %>% 
#    mutate(beta = exp(beta)) %>% 
#    mean_qi(.width = c(.50, .90))
# 
# post1 %>% 
#    pivot_longer(cols = starts_with("b"), names_to = "froh", values_to = "beta") %>% 
#    mutate(beta = exp(beta)) %>% 
#    ggplot(aes(beta, froh, fill = froh)) +
#    geom_vline(xintercept = 1, color = "black", size = 1) +
#    stat_halfeye() +
#    scale_fill_viridis_d() +
#    theme_minimal()
# 
# 
