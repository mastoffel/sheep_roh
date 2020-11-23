library(tidyverse)
library(lme4)
library(GGally)
library(partR2)
library(broom.mixed)
library(performance)
load("data/fitness_roh.RData") 
load("data/sheep_ped.RData")
ped <- sheep_ped

# check for froh across years
fitness_data %>% 
      select(starts_with("froh"), birth_year) %>% 
      pivot_longer(cols = c(starts_with("froh"))) %>% 
      mutate(birth_year = as.numeric(as.character(birth_year))) %>% 
      ggplot(aes(birth_year, value)) +
      geom_point() +
      geom_smooth(method = "lm") +
      facet_wrap(~name, scales = "free_y")

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
mod <- glmer(survival ~ froh_long_std + froh_medium_std + froh_short_std  + sex + twin + (1|mum_id) + (1|birth_year), 
             family = binomial, data = juv_survival,
             #control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
#mod_out <- tidy(mod, conf.int = TRUE, conf.method = "boot", nsim = 1000)
tidy(mod, conf.int = TRUE)
summary(mod)
library(sjPlot)
tab_model(mod,show.r2 = FALSE, show.icc = FALSE, auto.label = TRUE, transform = NULL)




# model weight
juv_weight <- fitness_data %>% 
   filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id, weight, hindleg), ~ !is.na(.)) %>% 
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

mod <- lmer(weight ~ froh_long_std + froh_medium_std + froh_short_std  + sex + twin + (1|mum_id) + (1|birth_year), 
            data = juv_weight)

tidy(mod, conf.int = TRUE)
summary(mod)
library(sjPlot)
tab_model(mod,show.r2 = FALSE, show.icc = FALSE, auto.label = TRUE, transform = NULL)



