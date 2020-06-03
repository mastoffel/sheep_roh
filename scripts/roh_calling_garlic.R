# roh calling with GARLIC
library(tidyverse)
read_delim("data/hd_ids.txt", " ", col_names = FALSE) %>% 
      mutate(X1 = 1) %>% write_delim("data/hd_ids.txt")

chr_data <- read_delim("../sheep/data/sheep_genome/chromosome_info_oar31.txt", delim = "\t") %>% 
      rename(size_BP = Length,
             CHR = Part) %>% 
      mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
      .[2:27, ] %>% 
      summarise(sum_KB = sum(size_KB)) %>% 
      as.numeric()

# generate transposed files
system(paste0("/usr/local/bin/plink --bfile ../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020 ",
              "--sheep --out data/sheep_geno --recode --keep data/hd_ids.txt --transpose "))

# fake centromere file
data.frame(1:26, 1, 2) %>% write_delim("data/dummy_cent", col_names = FALSE)

system(paste0("/usr/local/bin/garlic-master/bin/osx/garlic --tped data/sheep_geno.tped ",
              "--tfam data/sheep_geno.tfam --centromere data/dummy_cent ",
              "--error 0.005 --winsize-multi 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 --out output/ROH_garlic/test"))

all_kde <- purrr::map(c(seq(10, 170, 10)), function(x) {
      read_delim(paste0("output/ROH_garlic/test.", x, "SNPs.kde"), " ", col_names = FALSE)
      }) %>% 
      map(as.data.frame) %>% 
      bind_rows(.id = "winsize") %>% 
      as_tibble()
            
ggplot(all_kde, aes(X1, X2)) +
      geom_point() +
      facet_wrap(~winsize, scales = "free")

# lets choose window size of 60 SNPs
system(paste0("/usr/local/bin/garlic-master/bin/osx/garlic --tped data/sheep_geno.tped ",
              "--tfam data/sheep_geno.tfam --centromere data/dummy_cent ",
              "--error 0.005 --winsize 60 --out output/ROH_garlic/roh"))

# check ROH
roh_garlic_ids <- fread("output/ROH_garlic/roh.roh.bed", fill = TRUE) %>% 
                  filter(V1 == "track") %>% 
                  select(V3) %>% 
                  as_tibble()
roh_garlic <- fread("output/ROH_garlic/roh.roh.bed", fill = TRUE) %>% 
                  filter(V1 != "track") %>% 
                  dplyr::select(V1) %>% 
                  separate(V1, sep = "\t", into = c("chr", "roh_start",
                                                    "roh_end", "size_class",
                                                    "roh_length", "p1", "p2", "p3", "RGB")) %>% 
                  as_tibble()

num_roh_p_id <- fread("output/ROH_garlic/roh.roh.bed", fill = TRUE) %>% 
                  select(V1) 
num_roh <- which(num_roh_p_id$V1 == "track")
num_roh_p_id <- c(1, num_roh[2:length(num_roh)] - 1:(length(num_roh)-1))

roh_garlic$ids <- NA_integer_
roh_garlic[num_roh_p_id, "ids"] <- roh_garlic_ids$V3
roh <- roh_garlic %>% 
      fill(ids, .direction = "down") %>% 
      mutate(KB = as.numeric(roh_length)/1000)

froh <- roh %>% 
            group_by(ids, size_class) %>% 
            summarise(sum_KB = sum(KB)) %>% 
            mutate(froh = sum_KB/autosomal_genome_size)

ggplot(froh, aes(size_class, froh)) +
      geom_boxplot() +
      geom_point()

froh %>% 
      select(-sum_KB) %>% 
      pivot_wider(names_from = size_class, values_from = froh) %>% 
      ggpairs(columns = 2:4)


# now with all individuals
# generate transposed files
system(paste0("/usr/local/bin/plink --bfile ../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020 ",
              "--sheep --out data/sheep_geno_full --recode --transpose --exclude ../sheep_id/data/oar_imp_low_call_snps095.txt "))

# lets choose window size of 60 SNPs
system(paste0("/usr/local/bin/garlic-master/bin/osx/garlic --tped data/sheep_geno_full.tped ",
              "--tfam data/sheep_geno_full.tfam --centromere data/dummy_cent ",
              "--error 0.005 --winsize 60 --out output/ROH_garlic/roh_full"))





                  
