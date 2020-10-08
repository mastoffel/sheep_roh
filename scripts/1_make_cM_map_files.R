# call ROH
library(tidyverse)
library(data.table)
library(snpStats)
library(GGally)
library(gghalves)
library(furrr)
source("../sheep_ID/theme_simple.R")

# make ped file to replace with genetic distances
system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar31_17052020 ",
              "--sheep --out data/sheep_geno_imputed_oar31_17052020 ",
              "--recode --tab"))

# replace physical with genetic distances in map file
map_file <- read_delim("data/sheep_geno_imputed_oar31_17052020.map", delim = "\t",
                       col_names = FALSE) %>% 
            setNames(c("chr", "snp", "cM", "bp"))
head(map_file)

# linkage map from Johnston et al. (2020)
lmap <- read_delim("../sheep_ID/data/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
      rename(snp = SNP.Name,
             chr = Chr) %>% 
      select(chr, snp, cMPosition, cMPosition.Female, cMPosition.Male)

# join and interpolate cM positions (sex-averaged)
map_cM_sex_av <- map_file %>% 
                  left_join(lmap[c("chr", "snp", "cMPosition")]) %>% 
                  fill(cMPosition) %>% 
                  # give 0 at the beginning and make it cMPosition otherwise
                  mutate(cM = ifelse(is.na(cMPosition), 0, cMPosition)) %>% 
                  # trick plink by giving cM position to bp
                  mutate(bp = cM * 1e6) %>% 
                  select(-cMPosition)
write_delim(map_cM_sex_av, path = "data/sheep_geno_imputed_oar31_17052020_cM.map", 
                              delim = "\t", col_names = FALSE)
                  
# join and interpolate cM positions (female)
map_cM_female <- map_file %>% 
      left_join(lmap[c("chr", "snp", "cMPosition.Female")]) %>% 
      rename(cMPosition = cMPosition.Female) %>% 
      fill(cMPosition) %>% 
      # give 0 at the beginning and make it cMPosition otherwise
      mutate(cM = ifelse(is.na(cMPosition), 0, cMPosition)) %>% 
      # trick plink by giving cM position to bp
      mutate(bp = cM * 1e6) %>% 
      select(-cMPosition)

write_delim(map_cM_female, path = "data/sheep_geno_imputed_oar31_17052020_cM_female.map", 
            delim = "\t", col_names = FALSE)

# join and interpolate cM positions (male)
map_cM_male <- map_file %>% 
      left_join(lmap[c("chr", "snp", "cMPosition.Male")]) %>% 
      rename(cMPosition = cMPosition.Male) %>% 
      fill(cMPosition) %>% 
      # give 0 at the beginning and make it cMPosition otherwise
      mutate(cM = ifelse(is.na(cMPosition), 0, cMPosition)) %>% 
      # trick plink by giving cM position to bp
      mutate(bp = cM * 1e6) %>% 
      select(-cMPosition)

write_delim(map_cM_male, path = "data/sheep_geno_imputed_oar31_17052020_cM_male.map", 
            delim = "\t", col_names = FALSE)


# convert to bed
system(paste0("/usr/local/bin/plink --ped data/sheep_geno_imputed_oar31_17052020.ped ",
              "--map data/sheep_geno_imputed_oar31_17052020_cM.map ",
              "--sheep --out data/sheep_geno_imputed_oar31_17052020_cM ",
              "--make-bed"))

system(paste0("/usr/local/bin/plink --ped data/sheep_geno_imputed_oar31_17052020.ped ",
              "--map data/sheep_geno_imputed_oar31_17052020_cM_female.map ",
              "--sheep --out data/sheep_geno_imputed_oar31_17052020_cM_female ",
              "--make-bed"))

system(paste0("/usr/local/bin/plink --ped data/sheep_geno_imputed_oar31_17052020.ped ",
              "--map data/sheep_geno_imputed_oar31_17052020_cM_male.map ",
              "--sheep --out data/sheep_geno_imputed_oar31_17052020_cM_male ",
              "--make-bed"))
