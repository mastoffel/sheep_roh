# call ROH
library(tidyverse)
library(data.table)
library(snpStats)
library(GGally)
library(gghalves)
source("../sheep_ID/theme_simple.R")
# Chr lengths
chr_data <- read_delim("../sheep/data/sheep_genome/chromosome_info_oar31.txt", delim = "\t") %>% 
   rename(size_BP = Length,
          CHR = Part) %>% 
   mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
   .[2:27, ] %>% 
   summarise(sum_KB = sum(size_KB)) %>% 
   as.numeric()

# snps per 100Kb ~ 17
(417000/autosomal_genome_size) * 100

# 150 Kb ~ 256 generations ago
system(paste0("/usr/local/bin/plink --bfile ../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020 ",
              "--sheep --out output/ROH/roh ",
              "--homozyg --homozyg-window-snp 30 --homozyg-snp 30 --homozyg-kb 600 ",
              "--homozyg-gap 300 --homozyg-density 50 --homozyg-window-missing 2 ",
              "--homozyg-het 2 ",
              "--homozyg-window-het 2"))

file_path <- "output/ROH/roh.hom"
roh <- fread(file_path)
range(roh$KB)
hist(roh$KB, breaks = 1000)

# # calculate homozygosity in the rest of the genome -----------------------------
# sheep_plink_name <- "../../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_01052020"
# # read merged plink data
# sheep_bed <- paste0(sheep_plink_name, ".bed")
# sheep_bim <- paste0(sheep_plink_name, ".bim")
# sheep_fam <- paste0(sheep_plink_name, ".fam")
# #full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
# #summary(full_sample$genotypes)
# 
# snp_names <- read_delim("../../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_01052020.bim", "\t",
#                          col_names = FALSE) %>% .$X2
# snp_index <- 1:length(snp_names)
# snps_df <- data.frame(SNP1 = snp_names, SNP1_index = snp_index, 
#                       stringsAsFactors = FALSE)
# 
# # find SNPs in ROH for each individual
# roh_snps <- roh %>% 
#    #sample_frac(0.00001) %>% 
#    left_join(snps_df, by = "SNP1") %>% 
#    mutate(SNP2_index = SNP1_index + NSNP - 1) %>% 
#    mutate(snp_indices = paste0(SNP1_index, ":", SNP2_index)) %>% 
#    as_tibble() %>% 
#    mutate(snp_indices_full = purrr::map(snp_indices, function(x) eval(parse(text = x)))) %>% 
#    group_by(IID) %>% 
#    summarise(snps_in_roh = list(unlist(c(snp_indices_full))))
# 
# 
# double_check <- unlist(apply(roh_snps,1, function(x) length(x[2][[1]])))
# hist(double_check)
# 
# # load genotypes
# full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
# sheep_geno <- full_sample$genotypes
# sample_qc <- row.summary(sheep_geno) 
# sample_qc <- sample_qc %>% as_tibble(rownames = "id")
# 
# #sheep_geno <- as(full_sample$genotypes, Class = "numeric")
# # sheep_geno <- sheep_geno[rownames(full_sample$genotypes) %in% as.character(roh_snps$IID), ]
# 
# calc_hom <- function(id, sheep_geno, roh_snps) {
#    id_index <- which(rownames(sheep_geno) == id)
#    roh_snps_sub <- roh_snps[roh_snps$IID == id, ]$snps_in_roh[[1]]
#    genos <- as.numeric(as(sheep_geno[id_index, ], Class = "numeric"))
#    genos[roh_snps_sub] <- NA # give roh snps NA so that they aren't counted in het calc
#    # calc multilocus homozygosity
#    hom <- sum(genos!=1, na.rm = TRUE)/sum(!is.na(genos))
#    nsnps <- sum(!is.na(genos))
#    out <- data.frame(hom, nsnps)
# }
# 
# homs <- map(roh_snps$IID, calc_hom, sheep_geno, roh_snps)
# homs_df <- bind_rows(homs)
# #saveRDS(homs_df, file = "output/homs_df.rds")
# homs_df <- readRDS("output/homs_df.rds")
# homs_df$id <- as.character(roh_snps$IID)
# 
# # combine homozygosity outside ROH and qc
# hom_and_qc <- homs_df %>% 
#          left_join(sample_qc) %>% 
#          dplyr::select(id, everything()) %>% 
#          setNames(c("id", "hom_out_roh", "nsnps_hom", "call_rate", "certain", "het_all")) %>% 
#          select(-certain)
# 
# plot(test$hom, test$Call.rate)


# ROH classes ------------------------------------------------------------------

# expected ROH length using cM/Mb from Johnston et al (2016)
length_dist <- data.frame(g = c(1, 2,2^2, 2^3, 2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13)) %>%
   mutate(ROH_length_cM = 100 / (2*g)) %>% 
   mutate(ROH_length_Mb = ROH_length_cM * 0.7816743)

prop_IBD_df <- roh %>%
      mutate(length_Mb = KB/1000) %>%
      mutate(class = case_when(length_Mb >= 39.083715000 ~ 1,
                               length_Mb < 39.083715000 & length_Mb >= 19.541857500 ~ 2,
                               length_Mb < 19.541857500 & length_Mb >= 9.770928750 ~ 4,
                               #length_Mb < 9.770928750 & length_Mb >= 6.513952500 ~ 6,
                               length_Mb < 9.770928750& length_Mb >= 4.885464375 ~ 8,
                               # length_Mb < 4.885464375 & length_Mb >= 3.908371500 ~ 10,
                               length_Mb < 4.885464375 & length_Mb >= 2.442732188 ~ 16,
                               length_Mb < 2.442732188 & length_Mb >= 1.221366094 ~ 32,
                               length_Mb < 1.221366094 & length_Mb >= 0.6 ~ 64)) %>% # 0.610683047
      mutate(length_class = case_when(
            class == 1 ~ ">39 (1G)",
            class == 2 ~ ">19.5-39 (>1-2G)",
            class == 4 ~ ">9.8-19.5 (>2-4G)",
            # class == 6 ~ "6.5-9.7 (6G)",
            class == 8 ~ ">4.9-9.8 (>4-8G)",
            # class == 10 ~ "3.9-4.9 (10G",
            class == 16 ~ ">2.4-4.9 (>8-16G)",
            class == 32 ~ ">1.2-2.4 (>16-32G)",
            class == 64 ~ ">0.6-1.2 (>32-64G)"
           # class == 128 ~ "0.6-0.3 (128G)"
      )) %>% 
      mutate(length_class = fct_reorder(length_class, class)) %>% 
      mutate(IID = as.character(IID)) %>% 
      group_by(IID, class, length_class) %>%
      dplyr::summarise(prop_IBD = sum(length_Mb / (autosomal_genome_size/1000))) #%>% 

# add IBD of non-ROH snps if wanted
#  bind_rows(homs) 

prop_IBD_df_with_0 <- prop_IBD_df %>% 
      # add missing length classes as 0
      ungroup() %>% 
      tidyr::complete(length_class, nesting(IID)) %>% 
      mutate(class = ifelse(is.na(class), length_class, class)) %>% 
      mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD))

prop_IBD_df_with_0 %>% arrange(IID)


lowerFn <- function(data, mapping, method = "lm", ...) {
      p <- ggplot(data = data, mapping = mapping) +
            geom_point(colour = "blue") +
            geom_smooth(method = method, color = "red", ...)
      p
}
prop_IBD_df_with_0 %>% 
      select(-class) %>% 
      pivot_wider(names_from = length_class, values_from = prop_IBD) %>% 
      ggpairs(columns = 2:ncol(.),
              lower = list(continuous = wrap(lowerFn, method = "lm"))) +
      theme_simple()
ggsave("figs/froh_corr_mat.jpg", height = 11, width = 14)


# ROH class distribution across individuals
library(viridis)
col_pal <- viridis(8)
prop_IBD_df_with_0 %>% 
      mutate(prop_IBD = prop_IBD * 100) %>% 
      ggplot(aes(length_class, prop_IBD, fill = length_class)) +
      geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
                      transformation_params = list(height = 0, width = 1.3, seed = 1)) +
      geom_half_boxplot(side = "r", outlier.color = NA,
                        width = 0.6, lwd = 0.3, color = "black",
                        alpha = 0.8) +
     # theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
      ylab("% genome") +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
      theme(legend.position = "none",
            #axis.ticks.x = element_blank(),
            axis.title=element_text(size = rel(1.1)), 
            axis.text = element_text(color = "black")) + 
      xlab("ROH classes in Mb") 

quantile(roh$KB, probs = c(0.25,0.5,0.75))

# define ROH length classes
calc_froh_classes <- function(roh_crit, roh_lengths) {
   
   roh_filt <- dplyr::case_when(
      roh_crit == "short"  ~ expr(KB < 1221),
      roh_crit == "medium" ~ expr((KB >= 1221)&(KB < 4885)),
      roh_crit == "long"   ~ expr(KB >= 4885),
      roh_crit == "all" ~ expr(KB > 0)
   )
   
   roh_lengths %>%
      dplyr::group_by(IID) %>%
      #filter({{ roh_filt }}) %>% 
      filter(!!roh_filt) %>% 
      dplyr::summarise(KBSUM = sum(KB)) %>% 
      mutate(FROH = KBSUM / autosomal_genome_size) %>% 
      dplyr::select(IID, FROH) %>% 
      rename(ID = IID, !! paste0("froh_", roh_crit) := FROH)
   
}

# proportion of ROH length classes in each genome. Individuals which
# do not have long ROH have 0 for this class.
ROH_classes <- c("short", "medium", "long", "all")
froh <- purrr::map(ROH_classes, calc_froh_classes, roh) %>% 
   purrr::reduce(left_join, by = "ID") %>% 
   replace_na(list(FROH_long = 0)) %>% 
   rename(id = ID) %>% 
   mutate(id = as.character(id))

# join with hom
froh <- froh %>% 
         left_join(hom_and_qc)

froh %>% filter(call_rate > 0.95) %>% 
   ggplot(aes(froh_long, hom_out_roh)) +
   geom_point()
#plot(froh$froh_all, homs_df$nsnps)

# fitness data
# annual measures of traits and fitness
fitness_path <- "../sheep/data/1_Annual_Fitness_Measures_April_20190501.txt"
annual_fitness <- read_delim(fitness_path, delim = "\t")
names(annual_fitness)

# prepare, order pedigree
sheep_ped <- read_delim("../sheep/data/SNP_chip/20190711_Soay_Pedigree.txt", 
                        delim = "\t",
                        col_types = "ccc") %>%
   as.data.frame() %>%
   MasterBayes::orderPed() 

fitness_data <- annual_fitness %>% 
   left_join(froh, by = "ID")

fitness_data <- fitness_data %>% 
   dplyr::rename(birth_year = BIRTHYEAR,
                 sheep_year = SheepYear,
                 age = Age,
                 id = ID,
                 twin = TWIN,
                 sex = SEX,
                 mum_id = MOTHER,
                 froh_short =  froh_short,
                 froh_medium = froh_medium,
                 froh_long =   froh_long,
                 froh_all =    froh_all,
                 hom = hom,
                 nsnps_hom = nsnps,
                 # froh_not_roh = hom,
                 survival = Survival,
                 comment = Comment,
                 seen_in_rut = SeenInRut,
                 dad_id = FATHER,
                 weight = Weight, 
                 birth_weight = BIRTHWT,
                 cap_month = CapMonth,
                 hindleg = Hindleg,
                 foreleg = Foreleg,
                 rut_measure = RutMeasure,
                 horn_len = HornLen,
                 horn_circ = HornCirc,
                 bol_circ = BolCirc,
                 bol_len = BolLen,
                 horn = Horn,
                 freq = Freq,
                 offspr_born = OffspringBorn,
                 offspr_surv = OffspringSurvived,
                 offspr_surv_aug = AugSurvOffspring,
                 offspr_surv_oct = OctSurvOffspring,
                 offspr_surv_wint = OverWinterOffspring) %>% 
   dplyr::select(id, survival, comment, sheep_year, age, birth_year,
                 sex, mum_id, twin, froh_all, froh_long, froh_medium, 
                 froh_all, everything()) %>% 
   # for sex, check what Cas is
   filter(sex %in% c("F", "M")) %>% 
   # some individuals arent imputed well and should be discarded 
   #filter(!(id %in% IDs_lots_missing$id)) %>% 
   # na rows should be discarded
   mutate_at(c("id", "sheep_year", "birth_year", "sex",
               "mum_id", "twin"), as.factor)


save(fitness_data, file = "data/survival_mods_data.RData") # formerly fitness_roh_df
save(sheep_ped, file = "data/sheep_ped.RData") # ordered ped

