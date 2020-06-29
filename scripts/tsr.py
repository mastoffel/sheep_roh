# Keywords: Python, tree-sequence recording, tree sequence recording

from sys import argv
import msprime, pyslim, gzip
import numpy as np
import matplotlib.pyplot as plt

ts = pyslim.load("slim_sim/sheep.trees").simplify()

# is recapitation necessary?
# recapitation adds diversity present in the initial generation;
# will it make a difference? In fact, most segments of the genome have
# already coalesced, i.e. 13379 out of 15471 have one root
sum([t.num_roots == 1 for t in ts.trees()])
sum([t.num_roots > 0 for t in ts.trees()])

# recapitate anyway.
recap = ts.recapitate(recombination_rate=1.27e-8, Ne=5000, random_seed=1)

# verify that it worked
orig_max_roots = max(t.num_roots for t in ts.trees())
recap_max_roots = max(t.num_roots for t in recap.trees())
print(f"Before recapitation, the max number of roots was {orig_max_roots}, "
      f"and after recapitation, it was {recap_max_roots}.")

# add mutations. Also wrap in SlimTreeSequence so that pyslim can still work with it.
mutated = pyslim.SlimTreeSequence(msprime.mutate(recap, rate=1e-8, random_seed=np.random.randint(1,100000000), keep=True))

alive = mutated.individuals_alive_at(0)

print(f"The tree sequence now has {mutated.num_mutations} mutations, "
      f"and mean pairwise nucleotide diversity is {mutated.diversity()}.")

n_dip_indv = int(mutated.num_samples / 2)
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]

with open("slim_sim/sheep_recap.vcf", "w") as vcf_file:
    mutated.write_vcf(vcf_file,  individual_names=indv_names)
    

