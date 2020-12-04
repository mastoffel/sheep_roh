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
load("data/fitness_roh2.RData") 
load("data/sheep_ped.RData")
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

ggpairs(juv_survival, columns = 9:12)

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

# lme4
m1 <- glmer(survival ~ froh_long_std + froh_medium_std + froh_short_std  + sex + twin + (1|mum_id) + (1|birth_year), 
             family = binomial, data = juv_survival,
             #control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
#mod_out <- tidy(mod, conf.int = TRUE, conf.method = "boot", nsim = 1000)
tidy(m1, conf.int = TRUE)
saveRDS(m1, "output/juv_survival_model_std.RDS")
summary(m1)
tab_model(m1,show.r2 = FALSE, show.icc = FALSE, auto.label = TRUE, transform = NULL)


# no std:
juv_survival2 <- juv_survival %>% 
   mutate(froh_long = (froh_long-mean(froh_long)) * 100,
          froh_medium = (froh_medium-mean(froh_medium))*100,
          froh_short = (froh_short - mean(froh_short)) *100,
          froh_all = (froh_all * 100))
m2 <- glmer(survival ~ froh_long + froh_medium + froh_short + sex + twin + (1|mum_id) + (1|birth_year), 
            family = binomial, data = juv_survival2,
            #control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
            control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
#mod_out <- tidy(mod, conf.int = TRUE, conf.method = "boot", nsim = 1000)
mod_tidy <- tidy(m2, conf.int = TRUE)
# mod_tidy %>% 
#    filter(str_detect(term, "froh")) %>% 
#    #mutate(across(.cols = c("estimate", "conf.low", "conf.high"), exp)) %>% 
#    ggplot(aes(estimate, term, xmax = conf.high, xmin = conf.low, color = term, fill = term)) +
#    geom_vline(xintercept = 1, linetype='dashed', size = 0.3) +  #colour =  "#4c566a",   # "#4c566a"  "#eceff4"
#    geom_errorbarh(alpha = 1, height = 0, size = 1) +
#    geom_point(size = 3.5, shape = 21, #col = "#4c566a", #fill = "#eceff4", # "grey69"
#               alpha = 1, stroke = 0.3)

saveRDS(m2, "output/juv_survival_model_nonstd.RDS")


m3 <- glmer(survival ~ froh_all + sex + twin + (1|mum_id) + (1|birth_year), 
            family = binomial, data = juv_survival,
            #control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
            control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(m3, conf.int = TRUE)




# brm
mod_brm <- brm(survival ~ froh_long + froh_medium + froh_short  + sex + twin + (1|mum_id) + (1|birth_year), 
               family = bernoulli(), data = juv_survival2, cores = 4, iter = 5000)
summary(mod_brm)
conditional_effects(mod_brm)

plot(mod_brm)
ggpredict(mod_brm, c("froh_long_std [all]")) 
brms::conditional_effects(mod_brm)
me <- brms::conditional_effects(mod_brm, "froh_long_std")
summary(mod_brm)
plot_model(mod_brm, terms = c("froh_long", "froh_medium", "froh_short"))

library(tidybayes)
post1 <- posterior_samples(mod_brm, pars = c("froh_long" ,"froh_medium", "froh_short"))
post <- spread_draws(mod_brm, )
post1 %>% 
   ggplot(aes)

library(bayesplot)
mcmc_areas(post1, prob = 0.9)

library(brmstools)
forest(mod_brm, pars =  c("froh_long" ,"froh_medium", "froh_short"))

library(ggridges)

mod_brm_summary <- post1 %>% 
   pivot_longer(cols = starts_with("b"), names_to = "froh", values_to = "beta") %>% 
   group_by(froh) %>% 
   mutate(beta = exp(beta)) %>% 
   mean_qi(.width = c(.50, .90))

post1 %>% 
   pivot_longer(cols = starts_with("b"), names_to = "froh", values_to = "beta") %>% 
   mutate(beta = exp(beta)) %>% 
   ggplot(aes(beta, froh, fill = froh)) +
   geom_vline(xintercept = 1, color = "black", size = 1) +
   stat_halfeye() +
   scale_fill_viridis_d() +
   theme_minimal()


