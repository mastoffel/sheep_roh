library(tidyverse)
library(data.table)
library(furrr)
library(gtable)

### prepare LD data with cM ####################################################

# make ped file to replace with genetic distances
system(paste0("/usr/local/bin/plink --bfile data/Plates_1to87_QC3 ",
              "--sheep --out data/Plates_1to87_QC3 ",
              "--recode --tab --chr 1-26"))

# replace physical with genetic distances in map file
map_file <- read_delim("data/Plates_1to87_QC3.map", delim = "\t",
                       col_names = FALSE) %>% 
      setNames(c("chr", "snp", "cM", "bp"))
head(map_file)

# linkage map from Johnston et al. (2020)
lmap <- read_delim("../sheep_ID/data/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
      rename(snp = SNP.Name,
             chr = Chr) %>% 
      select(chr, snp, cMPosition, cMPosition.Female, cMPosition.Male)

map_cM <- map_file %>% 
      left_join(lmap[c("chr", "snp", "cMPosition")]) %>% 
      mutate(org_map = ifelse(!is.na(cMPosition), TRUE, FALSE),
             index = 1:nrow(.))

intpol_cM <- function(row, snp_map) {
      # current SNP
      row_cur <- snp_map[row, ]
      
      # return row if SNP is on genetic map
      if (!is.na(row_cur$cMPosition)) return(row_cur)
      
      # take current chromosome
      snp_map_chr <- snp_map %>% filter(chr == row_cur$chr)
      row_chr <- which(snp_map_chr$snp == row_cur$snp)
      
      # start of chromosome
      if (row_chr < min(which(snp_map_chr$org_map))) {
            row_cur$cMPosition <- 0
            return(row_cur)
      }
      
      # if not start, lmap position before
      row_before <- snp_map_chr %>% 
            filter((index < row) & org_map) %>% 
            filter(index == max(index)) 
      
      # end of chromosome
      if (row_chr > max(which(snp_map_chr$org_map))) {
            row_cur$cMPosition <- row_before$cMPosition
            return(row_cur)
      }
      
      # if not end, lmap position after
      row_next <- snp_map_chr %>% 
            filter((index > row) & org_map) %>% 
            filter(index == min(index)) 
      
      # if flanking SNPs have the same cM position
      if (row_before$cMPosition == row_next$cMPosition) {
            row_cur$cMPosition <- row_before$cMPosition
            return(row_cur)
      }
      
      # interpolate if nothing of the above is true
      diff_bp <- row_next$bp - row_before$bp
      diff_cM <- row_next$cMPosition - row_before$cMPosition
      
      row_cur$cMPosition <- row_before$cMPosition + (row_cur$bp - row_before$bp)/diff_bp * diff_cM
      return(row_cur)
}


# takes a while
plan(multisession, workers = 6)
map_cM_sex_av <- future_map(1:nrow(map_cM), intpol_cM, map_cM) %>% 
      bind_rows()

map_cM_sex_av_int <- map_cM_sex_av %>% 
      mutate(bp = cMPosition * 1e6) %>% 
      dplyr::select(-cMPosition, -org_map, -index)

write_delim(map_cM_sex_av_int, file = "data/Plates_1to87_QC3_cM.map", 
            delim = "\t", col_names = FALSE)

# convert to bed
system(paste0("/usr/local/bin/plink --ped data/Plates_1to87_QC3.ped ",
              "--map data/Plates_1to87_QC3_cM.map ",
              "--sheep --out data/Plates_1to87_QC3_cM ",
              "--make-bed"))


### call LD ROH with cM ########################################################
system(paste0("/usr/local/bin/plink --bfile data/Plates_1to87_QC3_cM ",
              "--sheep --out output/ROH/roh_ld_cM ",
              "--homozyg --homozyg-window-snp 25 --homozyg-snp 25 --homozyg-kb 1500 ",
              "--homozyg-gap 500 --homozyg-density 200 --homozyg-window-missing 2 ",
              "--homozyg-het 1 ",
              "--homozyg-window-het 1"))

### calc FROH for LD ###########################################################

# linkage map length
lmap <- read_delim("../sheep_ID/data/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
      rename(snp = SNP.Name,
             chr = Chr) %>% 
      select(chr, snp, cMPosition, cMPosition.Female, cMPosition.Male) %>% 
      filter(chr %in% 1:26)

full_length_map <- lmap %>% group_by(chr) %>% 
      summarise(max_per_chr = max(cMPosition)) %>% 
      summarise(genetic_map_length = sum(max_per_chr)) %>% 
      .[[1]]

# calculate FROH for LD chip
calc_froh <- function(plink_roh) {
      roh <- fread(plink_roh)
      out <- roh %>% 
            rename(id = IID) %>% 
            rename(cM = KB) %>%
            mutate(cM = cM/1e3) %>% 
            group_by(id) %>% 
            summarise(froh = sum(cM)/full_length_map) %>% 
            mutate(id = as.character(id))
      out
}

froh_ld <- calc_froh("output/ROH/roh_ld_cM.hom")

# let's keep the SNP density constant
# SNPs imputed: 419281

# modeling
load("data/fitness_roh.RData") 
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

# add to fitness data
juv_surv_mod <- juv_survival %>% 
            left_join(froh_ld) %>% 
            mutate(froh_imp = froh_all * 100,
                   froh_ld = froh * 100)


# modeling lme4
library(lme4)
library(broom.mixed)
library(brms)
library(sjPlot)

brm_fit_imp <- brm(survival ~ froh_imp + sex + twin + (1|mum_id) + (1|birth_year), 
               family = bernoulli(), data = juv_surv_mod, cores = 4, iter = 10000,
               set_prior("normal(0,5)", class = "b"))
saveRDS(brm_fit_imp, "output/juv_survival_imp_cM.RDS")

brm_fit_ld <- brm(survival ~ froh_ld + sex + twin + (1|mum_id) + (1|birth_year), 
                   family = bernoulli(), data = juv_surv_mod, cores = 4, iter = 10000,
                   set_prior("normal(0,5)", class = "b"))
saveRDS(brm_fit_ld, "output/juv_survival_ld_cM.RDS")

tab_model(brm_fit_imp, brm_fit_ld)

# make tables for supplementary material
# TABLE

# make nice table
get_tidy_mod <- function(mod) {
      tidy(mod) %>% 
            select(term:conf.high) %>% 
            mutate(term = c("Intercept", 
                            "F<sub>ROH</sub>",
                            "Sex",
                            "Twin",
                            "Birth year",
                            "Mother ID")) %>% 
            mutate(add_info = c("", rep("continuous", 1), 
                                "categorical (0=male, 1=female)",
                                "categorical (0=singleton, 1=twin)",
                                "n = 1118",
                                "n = 39")) %>% 
            mutate(effect = c("", rep("Population level/fixed effects", 3),
                              rep("Group level/random effects (standard deviation)", 2))) %>% 
            select(7, 1:6) %>% 
            mutate_if(is.numeric, function(x) paste0(round(x, 3)," (",round(exp(x), 3), ")")) %>% 
            select(-std.error) %>% 
            setNames(c("effect", "Term", "Post.Mean", "CI (2.5%)", "CI (97.5%)", "Info"))
}

# imputed data
tidy_mod_imp <- get_tidy_mod(brm_fit_imp)
# ld data
tidy_mod_ld <- get_tidy_mod(brm_fit_ld)

tidy_mod_imp %>% gt(
      rowname_col = "term",
      groupname_col = "effect") %>% 
      tab_style(
            style = cell_text( weight = "bold"),
            locations = cells_column_labels(columns = TRUE)
      ) %>% 
      fmt_markdown(columns = TRUE) %>% 
      gtsave("JS_model_table_froh_all_imputed.png", path = "tables/")

tidy_mod_ld %>% gt(
      rowname_col = "term",
      groupname_col = "effect") %>% 
      tab_style(
            style = cell_text( weight = "bold"),
            locations = cells_column_labels(columns = TRUE)
      ) %>% 
      fmt_markdown(columns = TRUE) %>% 
      gtsave("JS_model_table_froh_all_ld.png", path = "tables/")

# lme4
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

m1 <- glmer(survival ~ froh_imp + sex + twin + (1|mum_id) + (1|birth_year), 
            family = binomial, data = juv_surv_mod,
            control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
#mod_out <- tidy(mod, conf.int = TRUE, conf.method = "boot", nsim = 1000)
mod_tidy <- tidy(m1, conf.int = TRUE)
mod_tidy

m2 <- glmer(survival ~ froh_ld + sex + twin + (1|mum_id) + (1|birth_year), 
            family = binomial, data = juv_surv_mod,
            control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
#mod_out <- tidy(mod, conf.int = TRUE, conf.method = "boot", nsim = 1000)
mod_tidy <- tidy(m2, conf.int = TRUE)
mod_tidy

tab_model(m1, m2)







# 
# # SNP density imputed, SNPs per Mb
# snp_dens_imp <- 419281/(autosomal_genome_size/1000)
# # SNP density per minimum ROH length
# snp_dens_imp * 0.39 
# # SNP density 50K, SNPs per Mb
# snp_dens_ld <- 39368/(autosomal_genome_size/1000)
# # minimum ROH length for ~ 25 SNPs
# 25/snp_dens_ld
# 
# system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar31_17052020 ",
#               "--sheep --out output/ROH/roh_imp_bp ",
#               "--homozyg --homozyg-window-snp 25 --homozyg-snp 25 --homozyg-kb 390 ",
#               "--homozyg-gap 250 --homozyg-density 100 --homozyg-window-missing 2 ",
#               "--homozyg-het 2 --chr 1-26 ",
#               "--homozyg-window-het 2"))
# 
# system(paste0("/usr/local/bin/plink --bfile ../sheep/data/SNP_chip/oar31_mapping/Plates_1to87_QC3 ",
#               "--sheep --out output/ROH/roh_ld_bp ",
#               "--homozyg --homozyg-window-snp 25 --homozyg-snp 25 --homozyg-kb 1500 ",
#               "--homozyg-gap 500 --homozyg-density 300 --homozyg-window-missing 2 ",
#               "--homozyg-het 2 ",
#               "--homozyg-window-het 2"))
