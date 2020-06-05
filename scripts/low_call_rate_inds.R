library(snpStats)
library(tidyverse)

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

ind_sum <- row.summary(full_sample$genotypes) %>% 
            as_tibble(rownames = "id")

# write ids with call rates lower than 98% to file
ind_sum %>% filter(Call.rate < 0.98) %>% .$id %>% cbind(1,.) %>% as.data.frame() %>% write_delim("data/low_call_rate_inds.txt", col_names = FALSE)
# send to eddie 
system("scp data/low_call_rate_inds.txt mstoffel@eddie.ecdf.ed.ac.uk:/exports/csce/eddie/biology/groups/pemberton/martin/sheep_roh/data/")
