# call ROH
library(tidyverse)
library(data.table)
library(snpStats)
library(GGally)
library(gghalves)
library(furrr)
library(qtl)
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

## create interpolated map, where constant recombination rate is assumed
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

# 
write_delim(map_cM_sex_av, file = "data/sheep_geno_imputed_oar31_17052020_cM_full.txt", 
            delim = "\t", col_names = FALSE)

# complete map file
map_cM_sex_av_int <- map_cM_sex_av %>% 
   mutate(bp = cMPosition * 1e6) %>% 
   dplyr::select(-cMPosition, -org_map, -index)

write_delim(map_cM_sex_av_int, file = "data/sheep_geno_imputed_oar31_17052020_cM.map", 
            delim = "\t", col_names = FALSE)

print(map_cM_sex_av, n = 100)

# compare to Jon's interpolation \ looks arlight
map_js <- read_delim("~/Downloads/Oar3.1_Interpolated.txt", delim = "\t") %>% 
            setNames(c("chr", "snp", "bp", "cM"))
map_cM_sex_av %>% 
   select(-bp,-cM) %>% 
   left_join(map_js) %>% 
   ggplot(aes(cMPosition, cM)) +
      geom_point() +
      facet_wrap(~chr)


# convert to bed
system(paste0("/usr/local/bin/plink --ped data/sheep_geno_imputed_oar31_17052020.ped ",
              "--map data/sheep_geno_imputed_oar31_17052020_cM.map ",
              "--sheep --out data/sheep_geno_imputed_oar31_17052020_cM ",
              "--make-bed"))



# old interpolation

# # join and interpolate cM positions (sex-averaged)
# map_cM_sex_av <- map_file %>% 
#    left_join(lmap[c("chr", "snp", "cMPosition")]) %>% 
#    fill(cMPosition) %>% 
#    # give 0 at the beginning and make it cMPosition otherwise
#    mutate(cM = ifelse(is.na(cMPosition), 0, cMPosition)) %>% 
#    # trick plink by giving cM position to bp
#    mutate(bp = cM * 1e6) %>% 
#    select(-cMPosition)
# 
# write_delim(map_cM_sex_av, path = "data/sheep_geno_imputed_oar31_17052020_cM.map", 
#             delim = "\t", col_names = FALSE)
# 
# # join and interpolate cM positions (female)
# map_cM_female <- map_file %>% 
#    left_join(lmap[c("chr", "snp", "cMPosition.Female")]) %>% 
#    rename(cMPosition = cMPosition.Female) %>% 
#    fill(cMPosition) %>% 
#    # give 0 at the beginning and make it cMPosition otherwise
#    mutate(cM = ifelse(is.na(cMPosition), 0, cMPosition)) %>% 
#    # trick plink by giving cM position to bp
#    mutate(bp = cM * 1e6) %>% 
#    select(-cMPosition)
# 
# write_delim(map_cM_female, path = "data/sheep_geno_imputed_oar31_17052020_cM_female.map", 
#             delim = "\t", col_names = FALSE)
# 
# # join and interpolate cM positions (male)
# map_cM_male <- map_file %>% 
#    left_join(lmap[c("chr", "snp", "cMPosition.Male")]) %>% 
#    rename(cMPosition = cMPosition.Male) %>% 
#    fill(cMPosition) %>% 
#    # give 0 at the beginning and make it cMPosition otherwise
#    mutate(cM = ifelse(is.na(cMPosition), 0, cMPosition)) %>% 
#    # trick plink by giving cM position to bp
#    mutate(bp = cM * 1e6) %>% 
#    select(-cMPosition)
# 
# write_delim(map_cM_male, path = "data/sheep_geno_imputed_oar31_17052020_cM_male.map", 
#             delim = "\t", col_names = FALSE)
# 
# 
# 
# system(paste0("/usr/local/bin/plink --ped data/sheep_geno_imputed_oar31_17052020.ped ",
#               "--map data/sheep_geno_imputed_oar31_17052020_cM_female.map ",
#               "--sheep --out data/sheep_geno_imputed_oar31_17052020_cM_female ",
#               "--make-bed"))
# 
# system(paste0("/usr/local/bin/plink --ped data/sheep_geno_imputed_oar31_17052020.ped ",
#               "--map data/sheep_geno_imputed_oar31_17052020_cM_male.map ",
#               "--sheep --out data/sheep_geno_imputed_oar31_17052020_cM_male ",
#               "--make-bed"))

