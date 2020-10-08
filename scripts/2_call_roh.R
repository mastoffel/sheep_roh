# call ROH
library(tidyverse)
library(data.table)
library(snpStats)
library(GGally)
library(gghalves)
library(furrr)
library(janitor)
source("../sheep_ID/theme_simple.R")

system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar31_17052020_cM ",
              "--sheep --out output/ROH/roh_cM ",
              "--homozyg --homozyg-window-snp 30 --homozyg-snp 30 --homozyg-kb 390 ",
              "--homozyg-gap 100 --homozyg-density 100 --homozyg-window-missing 2 ",
              "--homozyg-het 1 ",
              "--homozyg-window-het 1"))

file_path <- "output/ROH/roh_cM.hom"
roh <- fread(file_path)
range(roh$KB)
hist(roh$KB, breaks = 1000)
#roh %>% group_by(IID) %>% summarise(prop_ibd = sum(KB / 3146000)) %>% summarise(mean(prop_ibd))

# individual qc  ---------------------------------------------------------------
sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")

# load genotypes
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
sheep_geno <- full_sample$genotypes
sample_qc <- row.summary(sheep_geno) %>% as_tibble(rownames = "id")

write_delim(sample_qc %>% clean_names(), "output/sample_qc.txt")

# calculate homozygosity in the rest of the genome -----------------------------

# # filter roh and genotypes for individuals with call rate > 0.99
# high_call_rate_id <- sample_qc %>% filter(Call.rate >= 0.99) %>% .$id
# 
# # genotypes
# sheep_geno <- sheep_geno[rownames(sheep_geno) %in% high_call_rate_id, ]
# # roh
# roh <- roh %>% filter(IID %in% high_call_rate_id)
# 
# # get snp names
# snp_names <- read_delim("../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020.bim", "\t",
#                         col_names = FALSE) %>% .$X2
# snp_index <- 1:length(snp_names)
# snps_df <- data.frame(SNP1 = snp_names, SNP1_index = snp_index,
#                       stringsAsFactors = FALSE)
# 
# # find SNPs in ROH for each individual
# roh_snps <- roh %>%
#       as_tibble() %>%
#       #sample_frac(0.00001) %>%
#       left_join(snps_df, by = "SNP1") %>%
#       mutate(SNP2_index = SNP1_index + NSNP - 1) %>%
#       mutate(snp_indices = paste0(SNP1_index, ":", SNP2_index)) %>%
#       mutate(snp_indices_full = purrr::map(snp_indices, function(x) eval(parse(text = x)))) %>%
#       group_by(IID) %>%
#       summarise(snps_in_roh = list(unlist(c(snp_indices_full))))
# 
# 
# double_check <- unlist(apply(roh_snps,1, function(x) length(x[2][[1]])))
# hist(double_check)
# 
# 
# # sheep_geno <- as(full_sample$genotypes, Class = "numeric")
# # sheep_geno <- sheep_geno[rownames(full_sample$genotypes) %in% as.character(roh_snps$IID), ]
# 
# # calculate maf
# snps_stats <- col.summary(full_sample$genotypes) %>% as_tibble(rownames = "snp")
# 
# calc_hom <- function(id, sheep_geno, roh_snps) {
#       
#       id_index <- which(rownames(sheep_geno) == id)
#       roh_snps_sub <- roh_snps[roh_snps$IID == id, ]$snps_in_roh[[1]]
#       num_snps_roh <- length(roh_snps_sub)
#       
#       # convert to numeric matrix and make it a vector
#       genos <- as(sheep_geno[id_index, ], Class = "numeric")[1, ]
#       # how many non-na genotypes overall?
#       num_snps_all <- sum(!is.na(genos))
#       
#       # give roh snps NA so that they aren't counted in het calc
#       genos_not_roh <- genos
#       genos_not_roh[roh_snps_sub] <- NA 
#       
#       # calculate moments estimate of inbreeding coefficient (Wright, 1948)
#       # get expected homozygosity
#       # snp_stats_sub <- snps_stats[!is.na(genos_not_roh), ]
#       # mafs <- snp_stats_sub$MAF
#       # # # exp number of homozygous snps
#       # exp_hom <- sum(1 - (2*mafs*(1-mafs)))
#       # obs_hom <- sum(genos_not_roh != 1, na.rm = TRUE)
#       # 
#       # # see clark et al (2019), equation 4 in methods
#       # o_dash_hom <-   obs_hom - num_snps_roh
#       # e_dash_hom <- ((num_snps_all - num_snps_roh)/num_snps_all) * exp_hom
#       # n_dash <- num_snps_all - num_snps_roh
#       # f_snp_outside_roh <- (o_dash_hom - e_dash_hom) / (n_dash - e_dash_hom)
#       
#       # calculate number of SNPs outside ROH
#       num_snps_non_roh <- sum(!is.na(genos_not_roh))
#       num_snps_non_roh_hom <- sum(genos_not_roh!=1, na.rm = TRUE)
#       out <- data.frame(num_snps_non_roh_hom, num_snps_non_roh, num_snps_all)
# }
# 
# homs <- map(roh_snps$IID, calc_hom, sheep_geno, roh_snps)
# homs_df <- bind_rows(homs) %>% 
#       mutate(hom1_all = num_snps_non_roh_hom/num_snps_all,
#              hom2_out = num_snps_non_roh_hom/num_snps_non_roh)
# 
# saveRDS(homs_df, file = "output/homs_cM_df.rds")
# #homs_df <- readRDS("output/homs_df.rds")
# homs_df$id <- as.character(roh_snps$IID)
# 
# # combine homozygosity outside ROH and qc
# hom_and_qc <- homs_df %>%
#       left_join(sample_qc) %>%
#       dplyr::select(id, everything()) %>%
#       rename(call_rate = Call.rate, het_all = Heterozygosity) %>% 
#       #setNames(c("id", "f_snp", "hom_out_roh", "nsnps_hom", "call_rate", "certain", "het_all")) %>%
#       select(-Certain.calls)
# 
# saveRDS(hom_and_qc, file = "output/homs_and_qc_cM_df.rds")

# system(paste0("/usr/local/bin/plink --bfile ../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020 ",
#               "--sheep --out output/ROH/roh ",
#               "--homozyg --homozyg-window-snp 30 --homozyg-snp 30 --homozyg-kb 600 ",
#               "--homozyg-gap 300 --homozyg-density 50 --homozyg-window-missing 2 ",
#               "--homozyg-het 2 ",
#               "--homozyg-window-het 1"))
# tests ------------------------------------------------------------------------
roh_kB <- fread("output/ROH/roh.hom")
roh_cM <- fread("output/ROH/roh_cM.hom")

ind_roh_kB <- roh_kB %>% 
                  group_by(IID) %>% 
                  summarise(sum(KB))
ind_roh_cM <- roh_cM %>% 
      group_by(IID) %>% 
      summarise(sum(KB))

# compare cM to kB ROH using Froh
plot(ind_roh_cM$`sum(KB)`, ind_roh_kB$`sum(KB)`)


