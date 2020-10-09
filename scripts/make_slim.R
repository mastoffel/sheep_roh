library(glue)
library(tidyverse)

# makes a slim simulation file
# original params genome_size = 1e8, pop_size1 = 5000, pop_size2 = 200, 
#  time1 = 10000, time2 = 11000,  mut_rate_del = 1e-9, recomb_rate = 1.27e-8,
#    mut1_gam_shape = 0.2,
make_slim <- function(genome_size = NULL, pop_size1 = NULL, pop_size2 = NULL, 
                      time1 = NULL, time2 = NULL,
                      mut_rate_del = NULL, recomb_rate = NULL,
                      mut1_dom_coeff = NULL, mut1_gam_mean = NULL, 
                      mut1_gam_shape = NULL,
                      mut2_dom_coeff = NULL, mut2_gam_mean = NULL, 
                      mut2_gam_shape = NULL, mut2_rel_freq = NULL,
                      mut3_dom_coeff = NULL, mut3_gam_mean = NULL, 
                      mut3_gam_shape = NULL, mut3_rel_freq = NULL,
                      out_dir = "slim_sim/sims",
                      seed = NULL) {
      
      # output directory
      slim_code_out <- paste0(out_dir, "/slim_code")
      dir.create(slim_code_out, recursive = TRUE, showWarnings = TRUE)
      # output for mutations
      muts_out <- paste0(out_dir, "/muts")
      dir.create(muts_out, recursive = TRUE, showWarnings = TRUE)
      # output for trees
      trees_out <- paste0(out_dir, "/trees")
      dir.create(trees_out, recursive = TRUE, showWarnings = TRUE)
      # # output for vcfs
      # vcfs_out <- paste0(out_dir, "/vcfs")
      # dir.create(vcfs_out, recursive = TRUE, showWarnings = TRUE)
      # 
      # length starts from 0
      genome_size <- genome_size - 1
      
      # needs seed for naming
      if (is.null(seed)) stop("need seed for naming file")
      
      # construct mutation bit -------------------------------------------------------
      # up to three mutation types on top of neutral mut possible
      mut2 <- mut3 <- NULL
      
      mut1 <- glue('
                   initializeMutationType("m1", {mut1_dom_coeff}, "g", \\
                   {mut1_gam_mean}, {mut1_gam_shape});
                   ')
      
      # define mutation type 2
      if (!is.null(mut2_dom_coeff)) {
            if (is.null(mut2_gam_mean) | is.null(mut2_gam_shape) | 
                is.null(mut2_rel_freq)) {
                  stop("define mean, shape and frequency for mutation 2")
            }
            mut2 <- glue('
                         initializeMutationType("m2", {mut2_dom_coeff}, "g", \\
                         {mut2_gam_mean}, {mut2_gam_shape});
                         ')
      }
      
      # define mutation type 3
      if (!is.null(mut3_dom_coeff)) {
            if (is.null(mut3_gam_mean) | is.null(mut3_gam_shape) | 
                is.null(mut3_rel_freq)) {
                  stop("define mean, shape and frequency for mutation 3")
            }
            mut3 <- glue('
                         initializeMutationType("m3", {mut3_dom_coeff}, "g", \\
                         {mut3_gam_mean}, {mut3_gam_shape});
                         ')
      }
      
      if (is.null(mut2)) {
            muts <- glue('
                         {mut1}
                         initializeGenomicElementType("g1", c(m1), c(1.0));
                         ')
      } 
      
      if (!is.null(mut2) & is.null(mut3)) {
            muts <-  glue('
                          {mut1}
                          {mut2}
                          initializeGenomicElementType("g1", c(m1, m2),\\
                          c(1.0, {mut2_rel_freq}));
                         ')
      }
      
      if(!is.null(mut2) & !is.null(mut3)) {
            muts <-  glue('
                          {mut1}
                          {mut2}
                          {mut3}
                          initializeGenomicElementType("g1", c(m1, m2, m3), \\
                          c(1.0, {mut2_rel_freq}, {mut3_rel_freq}));
                         ')
      }
      
      # end construct mutation bit ---------------------------------------------------
      
      
      slim_file <- glue('
      
      // set up a simple neutral simulation
      initialize() {{
      	initializeSLiMOptions(keepPedigrees=T);
      	initializeTreeSeq();
      	initializeMutationRate({mut_rate_del});
      	
      	// mutation types and genomic elelement type
            {muts}
      	
      	// uniform chromosome of length 100 Mb with uniform recombination
      	initializeGenomicElement(g1, 0, {genome_size});
      	initializeRecombinationRate({recomb_rate});
      }}
           
      
      // create a population of pop_size individuals
      1 {{
      	defineConstant("simID", getSeed());
      	sim.addSubpop("p1", {pop_size1});
      }}
      
      {time1} {{ p1.setSubpopulationSize({pop_size2}); }}
      
      {time2} late() {{
      
      	// tree sequence recording for 
      	sim.treeSeqOutput("{trees_out}/sheep_"+simID+".trees");
      	
      	// mutations per genome
      	genomes = sim.subpopulations.genomes;
      	out1 = string(length(sim.subpopulations.genomes.mutations));
      	out1_index = 0;
      	
      	for (genome_index in seqAlong(genomes)) {{
      		genome = genomes[genome_index];
      		pedigree_id = genome.individual.pedigreeID;
      		muts = genome.mutations;
      	
      		for (mut_index in seqAlong(muts)){{
      			out1[out1_index]=paste(c(genome_index, asString(pedigree_id),\\
muts[mut_index].id, muts[mut_index].position, muts[mut_index].selectionCoeff));
      			out1_index = out1_index + 1;
      		}}
      		header1=paste(c("genome_id", "pedigree_id", "mut_id", "pos", "s"));
      		writeFile("{muts_out}/mutperind_"+simID+".txt", \\
paste(c(header1, out1), sep="\\n"));
      	}}
      	
      }}
      
      ')
      
      write_lines(slim_file, glue("{slim_code_out}/sheep_{seed}.slim"))
      
}




