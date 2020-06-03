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





                  
