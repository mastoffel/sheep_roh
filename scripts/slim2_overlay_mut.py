# Keywords: Python, tree-sequence recording, tree sequence recording

from sys import argv
import inspect
import msprime, pyslim, gzip
import numpy as np
#import matplotlib.pyplot as plt
import os

# get args
script, run_name, pop_size = argv
infile = "slim_sim/sims/trees/" + run_name + ".trees"
outfile = "slim_sim/sims/vcfs/" + run_name + ".vcf"

# tree sequence
ts = pyslim.load(infile) # "slim_sim/sheep.trees"

# alive = ts.individuals_alive_at(0)
# print(f"There are {len(alive)} individuals alive from the final generation.")
ts.individual(201).id

# why was each individuals retained in tree sequence?
# indiv_types = {"first_gen" : 0,
#                "remembered" : 0,
#                "alive" : 0}
# for ind in mutated.individuals():
#    if ind.flags & pyslim.INDIVIDUAL_FIRST_GEN:
#       indiv_types['first_gen'] += 1
#    if ind.flags & pyslim.INDIVIDUAL_REMEMBERED:
#       indiv_types['remembered'] += 1
#    if ind.flags & pyslim.INDIVIDUAL_ALIVE:
#       indiv_types['alive'] += 1
# 
# for k in indiv_types:
#    print(f"Number of individuals that are {k}: {indiv_types[k]}")
   
   
# is recapitation necessary?
# recapitation adds diversity present in the initial generation;
# will it make a difference? In fact, most segments of the genome have
# already coalesced, i.e. 13379 out of 15471 have one root
sum([t.num_roots == 1 for t in ts.trees()])
sum([t.num_roots > 0 for t in ts.trees()])

# recapitate anyway.
recap = ts.recapitate(recombination_rate=1e-8, Ne=pop_size, random_seed=np.random.randint(1,100000000))

# verify that it worked
# orig_max_roots = max(t.num_roots for t in ts.trees())
# recap_max_roots = max(t.num_roots for t in recap.trees())
# print(f"Before recapitation, the max number of roots was {orig_max_roots}, "
#       f"and after recapitation, it was {recap_max_roots}.")


# some individuals might not be simplified away because they have nodes that
# are required to describe the genealogies of the sample
# keep_indivs = np.random.choice(recap.individuals_alive_at(0), 200, replace=False)
# keep_nodes = []
# for i in keep_indivs:
#    keep_nodes.extend(recap.individual(i).nodes)
# recap_simp = recap.simplify(keep_nodes)

# add mutations. Also wrap in SlimTreeSequence so that pyslim can still work with it.
mutated = pyslim.SlimTreeSequence(msprime.mutate(recap, rate=1e-8, 
                                  random_seed=np.random.randint(1,100000000), 
                                  keep=True))
# print(f"The tree sequence now has {mutated.num_mutations} mutations, "
#       f"and mean pairwise nucleotide diversity is {mutated.diversity()}.")

# only keep individuals that are alive today
# mutated=mutated.simplify()

n_dip_indv = int(mutated.num_samples / 2)
ind_names=[0] * n_dip_indv
for i in range(n_dip_indv):
      ind_names[i] = mutated.individual(i).metadata["pedigree_id"]
      
#indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
indv_names = [f"tsk_{str(i)}indv" for i in ind_names]

alive = mutated.individuals_alive_at(0).tolist()

with open(outfile, "w") as vcf_file:
    mutated.write_vcf(vcf_file, individuals = alive, individual_names=indv_names)
    








  
