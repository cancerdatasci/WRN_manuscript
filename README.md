# WRN MANUSCRIPT

Analysis scripts used in the paper "WRN Helicase is a Synthetic Lethal Target in Microsatellite Unstable Cancers" by Chan and Shibue et al. More info [here](https://depmap.org/WRN/)

## Data

An rds file containing a list of datasets used in the main analysis (**generate_figs.Rmd**) is available as a [Figshare repo](https://figshare.com/articles/DepMap_Datasets_for_WRN_manuscript/7712756).

## What's what.

* **process_CCLE_MSI_data.R**: This script takes the table of MSI calls from the CCLE Phase 2 manuscript (*ref*) and generates some normalized MSI indel stats used in **make_cell_line_info.R**

* **make_cell_line_info.R**: This script generates the table of cell line features (omics features and dependency scores) used to generate most analysis. This table is saved to **WRN_final_cell_line_table.csv**

* **WRN_helpers.R**: This script is sourced in **generate_figs.Rmd** and has some simpler stats and plotting helper functions.

* **generate_figs.Rmd**: This script generates all figures analyzing the CCLE/Achilles/GDSC data. It takes as input the cell line table mentioned above (Supplementary Table 1), and also takes as input a set of cell line feature matrices from DepMap 18Q4 data release (available [here](https://figshare.com/articles/DepMap_Datasets_for_WRN_manuscript/7712756)). In the analysis script these are pulled from our internal data server into a data list object "**dat**". For convenience, we also provide this list saved as an rds file on the Figshare repo.

* **WRN_diff_express.Rmd**: This script runs differenital expression analyses on the RNAseq data and generates figures.

* **in_vivo_KM12_analysis.R**: Short script to run linear mixed model analysis on in vivo xenograft data.

* **intensity_calcs.Rmd**: Runs linear model contrast tests on IF intensity values across cells.

* **WRN_stats_calcs.Rmd**: Script that runs some other small-scale t-test and ANOVA stats.