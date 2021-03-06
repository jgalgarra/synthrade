# Synthrade: Synthetic model simulator of World Trade Network

Authors: Javier Garcia-Algarra (U-TAD, Spain)
         Mary Luz Mouronte-Lopez (UFV, Spain)
         Javier Galeano (UPM, Spain)


## Description

This repository contains the code of Synthrade simulator and input data.

### Prerequisites

R 3.2 or newer
Python 3.0 or newer
git bash installed


### Reproducibility

- Download or Clone the repository
- Run scripts from command line with Rscript

#### Getting original data

[CAUTION: This step takes a long time to run. Original .tsv files are huge and we do not provide it in this repository. We recommend to skip this preprocessing step]

- Download `.tsv` files from OEC [`https://atlas.media.mit.edu/en/`] and store them in raw_data directory
- Run python scripts `CalcTotalYears_sitc.py` & `CalcTotalYears_hs07.py` Results: File `TODOSYYYY.csv`, where `YYYY` is the year. `TODOSYYYY.csv` files may be found at `data/raw_data`

#### Creation of bipatrite trade matrices

- Run script `create_raw_trade_matrix.R init_year end_year`. Creates files `RedAdyComYYYY.txt` inside `data/ directory`. Each file stores the global ammount of trade among exporters (rows) and importers (columns), for one year.
- Run script `filter_matrixes.R init_year end_year`. Results: `RedAdyComYYYY_FILT.txt` files.

#### Run the synthetic model

- Run script `synthrade.R init_year end_year num_experiments` . This procedure is CPU intensive, try to run it in batch mode. Results: `RedAdyComYYYY_FILT_W_N.txt` files, where `YYYY` is the year and `N` the number of experiment.

#### Compute goodness of fit
- Run script `compute_gof_distributions.R init_year end_year`. Results are stored in `results`

#### Plots
- Run `Rscript paint_lillies_test.R`
- Run `Rscript paint_KS_plots.R`
- Run `Rscript paint_matrix.R init_year end_year`
- Run `Rscript paint_density_plots.R init_year end_year FS`
- Run `Rscript compute_strength_degree_slopes.R init_year end_year`
