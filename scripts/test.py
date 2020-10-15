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

print(run_name)
print(int(pop_size) * 3)
