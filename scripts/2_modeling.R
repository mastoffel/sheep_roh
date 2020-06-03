library(tidyverse)
library(lme4)

load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")
ped <- sheep_ped
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# survival data preprocessing
juv_survival <- fitness_data %>% 
      # filter na rows
      #filter_at(vars(survival, froh_all, birth_year, sheep_year, mum_id), ~ !is.na(.)) %>% 
      filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
      mutate(age_std = as.numeric(scale(age)),
             age_std2 = age_std^2,
             froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
             froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
             froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
             froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
             froh_all_std = as.numeric(scale(froh_all)),
             froh_shorter_cent = (froh_short + froh_medium) - mean((froh_short + froh_medium)),
             froh_short_std = as.numeric(scale(froh_short)),
             froh_medium_std = as.numeric(scale(froh_medium)),
             froh_long_std = as.numeric(scale(froh_long)),
             lamb = as.factor(ifelse(age == 0, 1, 0))) %>% 
      filter(age == 0)

cor(juv_survival$hom, juv_survival$froh_all)
sd(juv_survival$hom) 
sd(juv_survival$froh_all)
# lme4 -------------------------------------------------------------------------
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


mod <- glmer(survival ~ froh_long_cent + froh_shorter_cent + twin + sex + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
             family = binomial, data = juv_survival,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod)

mod_std <- glmer(survival ~ froh_long_std + froh_medium_std + froh_short_std + twin + sex + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
             family = binomial, data = juv_survival,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod_std)


ad_survival <- fitness_data %>% 
      # filter na rows
      #filter_at(vars(survival, froh_all, birth_year, sheep_year, mum_id), ~ !is.na(.)) %>% 
      filter_at(vars(survival, froh_all, sheep_year, mum_id), ~ !is.na(.)) %>% 
      mutate(age_std = as.numeric(scale(age)),
             age_std2 = age_std^2,
             froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
             froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
             froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
             froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
             froh_short2 = (froh_short + froh_medium) - mean(froh_short + froh_medium),
             froh_all_std = as.numeric(scale(froh_all)),
             froh_short_std = as.numeric(scale(froh_short)),
             froh_medium_std = as.numeric(scale(froh_medium)),
             froh_long_std = as.numeric(scale(froh_long)),
             lamb = as.factor(ifelse(age == 0, 1, 0))) %>% 
      filter(age > 0)

mod_ad <- glmer(survival ~ froh_long_cent + froh_short2 + twin + sex + age_std + age_std2 + (1|id) + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
             family = binomial, data = ad_survival,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod_ad)

mod_ad <- glmer(survival ~ froh_long_cent + froh_medium_cent + twin + sex  + (1|sheep_year) + (1|mum_id), 
                family = binomial, data = ad_survival,
                control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

summary(mod_ad)

mod_ad_std <- glmer(survival ~ froh_long_std + froh_medium_std + froh_short_std + twin + sex + age_std + age_std2 + (1|id) + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
                family = binomial, data = ad_survival,
                control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
options(scipen = 999)
summary(mod_ad_std)




ggpairs(juv_survival[c("froh_short_std", "froh_medium_std", "froh_long_std")])

mod <- glmer(survival ~  froh_short_cent + froh_medium_cent + froh_long_cent +  twin + sex + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
             family = binomial, data = juv_survival,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod)

# try with reproductive success
rep_suc <-  fitness_data %>% 
      filter_at(vars(froh_all, sheep_year, mum_id, offspr_surv), ~ !is.na(.)) %>% 
      mutate(age_std = as.numeric(scale(age)),
             age_std2 = age_std^2,
             froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
             froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
             froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
             froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
             lamb = as.factor(ifelse(age == 0, 1, 0))) %>% 
      filter(sex == "M") %>% 
      mutate(olre = 1:nrow(.)) 
hist(rep_suc$offspr_surv)

mod2 <- glmer(offspr_surv ~ froh_long_cent + froh_medium_cent + froh_short_cent + twin + age + age^2 + (1|sheep_year) + (1|olre) + (1|mum_id) + (1|birth_year), 
             family = poisson, data = rep_suc)
summary(mod2)

# try with body weight 
weight_df <-  fitness_data %>% 
      filter_at(vars(froh_all, sheep_year, mum_id, weight), ~ !is.na(.)) %>% 
      mutate(age_std = as.numeric(scale(age)),
             age_std2 = age_std^2,
             froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
             froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
             froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
             froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
             lamb = as.factor(ifelse(age == 0, 1, 0)),
             rel_weight = weight/(hindleg/10)) %>% 
      filter(age == 0)

mod3 <- lmer(rel_weight ~ froh_long_cent + froh_medium_cent + froh_short_cent + twin + sex + (1|sheep_year) + (1|mum_id), 
              data = weight_df)
summary(mod3)
