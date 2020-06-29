library(glue)
library(tidyverse)


out_dir <- "slim_sim"

make_slim <- function(chr_length = 1e7, pop_size1 = 5000, pop_size2 = 200, out_dir = "slim_sim"){
      
# length starts from 0
chr_length <- chr_length - 1


slim_file <- glue('

// set up a simple neutral simulation
initialize() {{
	initializeMutationRate(1e-8);
	
	// m1 mutation type: neutral
	initializeMutationType("m1", 0.5, "f", 0.0);
	// weakly deleterious mutations
	initializeMutationType("m2", 0.1, "g", -0.03, 0.2);
	// g1 genomic element type: uses m1 for all mutations
	initializeGenomicElementType("g1", c(m1, m2), c(1.0, 0.2));
	
	// uniform chromosome of length 100 Mb with uniform recombination
	initializeGenomicElement(g1, 0, {chr_length});
	initializeRecombinationRate(1.27 * 1e-8);
}}
     
// create a population of 1000 individuals
1 {{
	sim.addSubpop("p1", 5000);
}}

10000 {{ p1.setSubpopulationSize(200); }}

11000 late() {{
	g = p1.sampleIndividuals(200).genomes;
	g.outputVCF(filePath="{out_dir}/sheep.vcf");
	
	out = string(length(sim.mutations.position));
	for (i in 0:(length(sim.mutations.position)-1)){{
                  out[i]=paste(c(asString(sim.mutations.position[i]),
                        sim.mutations.mutationType[i].id,
                        sim.mutations.mutationType[i].dominanceCoeff,
                        sim.mutations.selectionCoeff[i]));
	}}
	header=paste(c("pos","type","h","s"));

	writeFile("{out_dir}/mut.txt",paste(c(header,out),"\\n"));

}}

')

write_lines(slim_file, glue("{out_dir}/test.slim"))
      
}

make_slim(chr_length = 1e8, out_dir = "slim_sim")

system("slim -m -l 2 -t slim_sim/test.slim")


