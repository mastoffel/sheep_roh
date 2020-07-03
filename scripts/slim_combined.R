library(data.table)
library(tidyverse)
source("scripts/slim3_mut_and_roh.R")
library(furrr)
# run slim simulation

slim_roh <- function(seed) {
      
      run_name <- paste0("sheep_", seed)
      
      # run slim
      system(paste("slim -time -Memhist -seed", seed, "scripts/slim1_simulation.slim"))
      
      # recapitation and overlay of neutral mutations
      system(paste("python3.8 scripts/slim2_overlay_mut.py", run_name))
      
      # call ROH
      # use vcf output to call ROH
      system(paste0("/usr/local/bin/plink --vcf slim_sim/sims/vcfs/", run_name, ".vcf ",
                    "--sheep --out slim_sim/sims/roh/", run_name, " ",
                    "--homozyg --homozyg-window-snp 20 --homozyg-snp 20 --homozyg-kb 1200 ",
                    "--homozyg-gap 100 --homozyg-density 100 --homozyg-window-missing 0 ",
                    "--homozyg-het 0 ",
                    "--homozyg-window-het 0"))
      
      out <- combine_mut_roh(run_name = run_name) 
}


seeds <- sample(1:1e5, 8)
plan(multiprocess, workers = 4)
out <- future_map(seeds, slim_roh)
save(out, file = "slim_sim/sims/out_test.RData")



