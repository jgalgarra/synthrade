# k-decomposition and analysis of mutualistic networks

Authors: Javier Garcia-Algarra/ Juan Manuel Pastor (UPM, Spain)


## Description

This repository contains code for the kcore decomposition, analysis and visualization of mutualistic networks.

### Prerequisites

R 3.2 or newer
Python 3.0 or newer
git bash installed

Intall kcorebip package with `devtools` package:


`install_github("jgalgarra/kcorebip")`

### Reproducibility

- Clone the repository with `git clone https://github.com/jgalgarra/kcore_robustness.git`
- Move to `kcore_robustness` directory and set it as R working directory (RStudio use is recommended)
- The `data` directory contains the 89 network interaction matrices downloaded from the web of life site
- Run `testing-all.R`. The output file `results/datos_analisis.RData` stores the k-magnitudes in `.csv` format
- Run `kdegree_calc_store_results.R` to get the `results/datos_analisis_condegs.RData` (results + network degrees and correlation degree kdegree)
- Run `network_k_parameters.R` to create individual k-magnitude files in `analysis_indiv_extended` directory

- `destruction_first_algorithm.R`. Algorithm to find the number of primary extinctions of any guild to destroy half the giant component according to different indexes. Read the R file documentation header for detailed instructions

- Go to the python directory and run `extictions_compare_new.py` . This task may be slow, do not stop it

Go back to RStudio
- Run `best_1stalg.R` and `best_2ndalg` to have a fast count of comparative performances. 
- Run `paint_destructions_1stalg_network.R` . The directory `graphs/FIRST` contains the plots of individual network performance
- Run `paint_destructions_2ndalg_network.R` . The directory `graphs/pyhton` contains the plots of individual network performance for the two outcomes ofthe second extinction algorithm
- Run `paint_extinctions_1stalg_results.R` and `paint_extinctions_2ndtalg_results.R` to create the four comparative figures in the directory `graphs`
- Run `paint_degree_distribution_steps.R` to create the comparative figure of degree vs. kdegree of network M_PL_001
- Run `paint_2ndalg_areas.R` to build the destruction AUC for the second algorithm under `graphs/AREAS`.
- Run `create_007_destroy_2nd_alg` to create the stages of destruction of network M_PL_007
- Run `bipartite_graphs_007.R` to create the bipartite plots of destruction of networkl M_PL_007