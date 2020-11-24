# currently, pyslim sometimes outputs weird REF and ALT positions in the vcf file
# these might be the mutations from SLiM
# When one of these is the first line in the VCF gt table, then REF is the empty
# string and ALT is a large number. plink doesn't like this and will find no ROH
# This means that some of the simulations can't be processed further. 

library(data.table)
library(tidyverse)
library(furrr)
source("scripts/slim3alt_mut_and_roh.R")

out_path <- "output/qm_slim/slim10000200"
file_names <- list.files(paste0(out_path, "/vcfs"), full.names = TRUE)

correct_vcf <- function(vcf_file) {
      # read meta data
      meta <- readr::read_delim(vcf_file, delim = "\t", n_max = 5,
                         col_names = FALSE)[[1]]
      # read genotypes, data.table does this correctly
      # i.e. jumping over the first 5 lines
      gt <- data.table::fread(vcf_file, skip = 5)
      # correct ref and alt columns
      gt$REF <- 0
      gt$ALT <- 1
      # write out meta
      readr::write_lines(meta, path = vcf_file)
      # append genotypes to make file complete
      fwrite(gt, file = vcf_file, append = TRUE,
             col.names = TRUE, sep = "\t")
      
      # now rerun ROH calling
      run_name <- stringr::str_split(vcf_file, "/") %>% 
                        unlist() %>% 
                        .[length(.)] %>% 
                        str_remove(".vcf")
            
      system(paste0("plink --vcf ", vcf_file, " ",  # /usr/local/bin/plink
                    "--sheep --out ", out_path, "/roh/", run_name, " ",
                    "--homozyg --homozyg-window-snp 30 --homozyg-snp 30 --homozyg-kb 390 ",
                    "--homozyg-gap 100 --homozyg-density 100 --homozyg-window-missing 0 ",
                    "--homozyg-het 0 ",
                    "--homozyg-window-het 0"))
      
}

# replace vcfs and re-run ROH calc
plan(multiprocess, workers = 6)
future_map(file_names, correct_vcf)

# extract run names
run_names <- stringr::str_split(file_names, "/") %>% 
      map_chr(5) %>% 
      str_remove(".vcf")

# make a safe combine function
combine_safe <- safely(combine_mut_roh)

#plan(multiprocess, workers = 4)
# combine mutations and roh data and calculate length classes
out <- future_map(run_names, combine_safe, out_path, roh_cutoff_small = 1560, roh_cutoff_long = 6250)

# params
mut1_dom_coeff <- c(0, 0.05, 0.2)
mut1_gam_mean <- c(-0.01, -0.03, -0.05)
mut1_gam_shape <- 0.2
genome_size <- 1e8
mut_rate_del <- 1e-9
recomb_rate <- 1e-8
pop_size1 <- 10000
pop_size2 <- 200
params <- expand.grid(pop_size1, pop_size2, mut1_dom_coeff, mut1_gam_mean, mut1_gam_shape,
                      genome_size, mut_rate_del, recomb_rate) %>% 
      setNames(c("pop_size1", "pop_size2", "mut1_dom_coeff", "mut1_gam_mean", "mut1_gam_shape",
                 "genome_size", "mut_rate_del", "recomb_rate")) %>% 
      as_tibble() %>% 
      mutate(time1 = 10000,        #pop_size1 * 10,
             time2 = time1 + 1000)

# try 10 simulation with only weakly deleterious alleles
num_sim_per_parset <- 100
set.seed(123)
seeds <- sample(1:1e5, num_sim_per_parset * nrow(params))
# replicate each parameter set num_sim_per_parset times
params_sim <- params[rep(1:nrow(params), each =num_sim_per_parset), ] %>% 
      mutate(seed = seeds)

# extract only non-error runs and combine
mut_df <- out %>% 
      purrr::transpose() %>% 
      .[[1]] %>% 
      bind_rows() %>% 
      mutate(id = paste0(ind_id, seed),
             seed = as.numeric(seed)) %>% 
      dplyr::arrange(id) %>% 
      left_join(params_sim, by = "seed")
# complete(nesting(id, roh_class)) %>% 
# arrange(id)

#saveRDS(mut_df, file = "slim_sim/sims/out/sims_weakly_03.RData")

#dir.create("slim_sim/sims/out", recursive = TRUE, showWarnings = TRUE)
write_delim(mut_df, paste0(out_path, "/out/par_combs_popsize1_", pop_size1,
                           "_popsize2_", pop_size2, ".txt"))
