library(tidyverse)
library(lme4)
library(GGally)
library(partR2)
library(broom.mixed)
load("data/roh_mods.RData") 
load("data/sheep_ped.RData")
ped <- sheep_ped

# general fitness_data checks
fitness_data %>% 
   select(starts_with("froh")) %>% 
   ggpairs()

# check for froh across years
fitness_data %>% 
   select(starts_with("froh"), birth_year) %>% 
   pivot_longer(cols = c(starts_with("froh"))) %>% 
   mutate(birth_year = as.numeric(as.character(birth_year))) %>% 
   ggplot(aes(birth_year, value)) +
   geom_point() +
   geom_smooth(method = "lm") +
   facet_wrap(~name, scales = "free_y")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# survival data preprocessing
juv_survival <- fitness_data %>% 
      filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
      mutate(age_std = as.numeric(scale(age)),
             age_std2 = age_std^2,
             froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
             froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
             froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
             froh_all_std = scale(froh_all),
             froh_short_std = scale(froh_short),
             froh_long_std = scale(froh_long),
             froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
             froh_medium_std = scale(froh_medium)) %>% 
             # hom1_all_cent = hom1_all - mean(hom1_all, na.rm = TRUE),
             # hom2_out_cent = hom2_out - mean(hom2_out, na.rm = TRUE),
             # hom1_all_std = scale(hom1_all),
             # hom2_out_std = scale(hom2_out)) %>% 
      filter(age == 0)# %>% 
      #filter(sex == "F")

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
# froh_long + froh_short + 
mod <- glmer(survival ~ froh_long + froh_medium + froh_short + sex + twin + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
             family = binomial, data = juv_survival,
             #control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod)
VarCorr(mod)



ad_survival_m <- fitness_data %>% 
   filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
   filter(sex == "F") %>% 
   mutate(age_std = as.numeric(scale(age)),
          age_std2 = age_std^2,
          froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
          froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
          froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
          froh_all_std = scale(froh_all),
          froh_short_std = scale(froh_short),
          froh_long_std = scale(froh_long),
          froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
          froh_medium_std = scale(froh_medium),
          lamb = as.factor(ifelse(age == 0, 1, 0))) 
     # filter(age > 0) %>% 
 

mod_ad <- glmer(survival ~ froh_long + froh_medium + froh_short + twin + sex + age_std + age_std2 + (1|id) + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
             family = binomial, data = ad_survival,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod_ad)

mod_ad <- glmer(survival ~ froh_short + twin + age_std + age_std2 + (1|id) + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
                family = binomial, data = ad_survival_m,
                control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod_ad)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Body Weight~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# survival data preprocessing
juv_weight <- fitness_data %>% 
   filter_at(vars(hindleg, weight, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
   mutate(age_std = as.numeric(scale(age)),
          age_std2 = age_std^2,
          froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
          froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
          froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
          froh_all_std = scale(froh_all),
          froh_short_std = scale(froh_short),
          froh_long_std = scale(froh_long),
          froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
          froh_medium_std = scale(froh_medium),
          hindweight = (weight/hindleg) * 10) %>% 
   # hom1_all_cent = hom1_all - mean(hom1_all, na.rm = TRUE),
   # hom2_out_cent = hom2_out - mean(hom2_out, na.rm = TRUE),
   # hom1_all_std = scale(hom1_all),
   # hom2_out_std = scale(hom2_out)) %>% 
   filter(age == 0)

library(INLA)
library(AnimalINLA)
sheep_ped_inla <- ped %>% 
   as_tibble() %>% 
   dplyr::rename(id = ID,
                 mother = MOTHER,
                 father = FATHER) %>% 
   mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "F", "888")) %>% 
   mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "M", "999")) %>% 
   mutate(father = ifelse(is.na(father), 0, father),
          mother = ifelse(is.na(mother), 0, mother)) %>% 
   mutate_if(is.character, list(as.numeric)) %>% 
   as.data.frame() 

# compute Ainverse and map
comp_inv <- AnimalINLA::compute.Ainverse(sheep_ped_inla)

# make Cmatrix and map
ainv <- comp_inv$Ainverse
ainv_map <- comp_inv$map
Cmatrix <- sparseMatrix(i=ainv[,1],j=ainv[,2],x=ainv[,3])
juv_weight$IndexA <-  match(juv_weight$id, ainv_map[, 1])
juv_weight$IndexA2 <- juv_weight$IndexA

prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
# model
formula_weight <- as.formula(paste('hindweight ~ froh_long_std + froh_medium_std + froh_short_std + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                # 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                  'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
#control.family1 = list(control.link=list(model="logit"))
mod_inla <- inla(formula=formula_weight, family="gaussian",
                 data=juv_weight, 
                # control.family = control.family1,
                 control.compute = list(dic = TRUE, config=TRUE))
summary(mod_inla)



# survival data preprocessing
ad_weight <- fitness_data %>% 
   filter_at(vars(weight, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
   mutate(age_std = as.numeric(scale(age)),
          age_std2 = age_std^2,
          froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
          froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
          froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
          froh_all_std = scale(froh_all),
          froh_short_std = scale(froh_short),
          froh_long_std = scale(froh_long),
          froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
          froh_medium_std = scale(froh_medium)) %>% 
   # hom1_all_cent = hom1_all - mean(hom1_all, na.rm = TRUE),
   # hom2_out_cent = hom2_out - mean(hom2_out, na.rm = TRUE),
   # hom1_all_std = scale(hom1_all),
   # hom2_out_std = scale(hom2_out)) %>% 
   filter(age > 3)

sheep_ped_inla <- ped %>% 
   as_tibble() %>% 
   dplyr::rename(id = ID,
                 mother = MOTHER,
                 father = FATHER) %>% 
   mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "F", "888")) %>% 
   mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "M", "999")) %>% 
   mutate(father = ifelse(is.na(father), 0, father),
          mother = ifelse(is.na(mother), 0, mother)) %>% 
   mutate_if(is.character, list(as.numeric)) %>% 
   as.data.frame() 

# compute Ainverse and map
comp_inv <- AnimalINLA::compute.Ainverse(sheep_ped_inla)

# make Cmatrix and map
ainv <- comp_inv$Ainverse
ainv_map <- comp_inv$map
Cmatrix <- sparseMatrix(i=ainv[,1],j=ainv[,2],x=ainv[,3])
ad_weight$IndexA <-  match(ad_weight$id, ainv_map[, 1])
ad_weight$IndexA2 <- ad_weight$IndexA

prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
# model
formula_weight <- as.formula(paste('weight ~ froh_long + froh_medium + froh_short + sex + twin + 1', 
                                   'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                   'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                   'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                   'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                   'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
#control.family1 = list(control.link=list(model="logit"))
mod_inla <- inla(formula=formula_weight, family="gaussian",
                 data=ad_weight, 
                 # control.family = control.family1,
                 control.compute = list(dic = TRUE, config=TRUE))
summary(mod_inla)










#~~~~~~~~~~~~~~~~~~~~~~~~~~~~by age ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
dat_survival <- fitness_data %>% 
   filter_at(vars(survival, froh_all, sheep_year, birth_year, mum_id), ~ !is.na(.)) %>% 
   mutate(age_std = as.numeric(scale(age)),
          age_std2 = age_std^2,
          froh_all_cent =    froh_all - mean(froh_all, na.rm = TRUE),
          froh_short_cent =  froh_short - mean(froh_short, na.rm = TRUE),
          froh_long_cent =   froh_long - mean(froh_long, na.rm = TRUE),
          froh_all_std = scale(froh_all),
          froh_short_std = scale(froh_short),
          froh_long_std = scale(froh_long),
          froh_medium_cent = froh_medium - mean(froh_medium, na.rm = TRUE),
          froh_medium_std = scale(froh_medium)) 


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

b_froh_per_age <- function(ageclass) {
   dat_survival_age <- dat_survival %>% filter(age == ageclass) #%>% filter(sex == sexclass)
   # froh_long + froh_short + 
   mod <- glmer(survival ~ froh_long + froh_medium + froh_short + sex + twin + (1|sheep_year) + (1|mum_id) + (1|birth_year), 
                family = binomial, data = dat_survival_age,
                control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
  tidy(mod)
}

#mods_f <- map_dfr(0:2, b_froh_per_age, "F", .id = "age") %>% mutate(sex = "F")
#mods_m <- map_dfr(0:2, b_froh_per_age, "M", .id = "age") %>% mutate(sex = "M")
#mods <- rbind(mods_f, mods_m)

mods <- map_dfr(0:4, b_froh_per_age, .id = "age") 

mods %>% 
   filter(term %in% c("froh_long", "froh_medium", "froh_short")) %>% 
   mutate(age = fct_inorder(age)) %>% 
   ggplot(aes(age, estimate, color = term)) +
      geom_point(position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymax = estimate + std.error, ymin = estimate - std.error), 
                 width = 0.2, position = position_dodge(width = 0.5))


















