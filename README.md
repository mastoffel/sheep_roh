### Mutation load decreases with haplotype age in wild Soay sheep
Stoffel, M.A, Johnston, S.E., Pilkington, J.G., Pemberton, J.M  
*bioRxiv* [https://doi.org/10.1101/2021.03.04.433700](https://doi.org/10.1101/2021.03.04.433700)  

### Overview
This repository contains the analysis code for our paper under scripts/

**1_make_cM_map_files:** Converts physical position PLINK .map file to cM positions, and uses interpolation to infer cM position for SNPs which are not part of the linkage map.

**2_call_roh:** ROH calling and quality control.

**3_roh_length_classes:** Divides ROH into length classes, calculates different Froh based on these classes and combines ROH and fitness data.

**4_modeling:** Mixed models to estimate inbreeding depression per ROH length class.

**5_plot_main_figure:** Formatting and plotting of the modeling results and Figure 1.

**6_plot_simulations:** Plotting SLiM simulations. 