# ELISA Data Analysis Pipeline

An automated R pipeline for processing raw plate reader data, quantifying target protein concentrations via 4-parameter logistic (4PL) regression, and generating publication-ready quality control figures.

Dummy experimental data has been added into the `data/` folder for testing.

## Project Context
Built for multi-variable clinical studies, this pipeline evaluates plasma enzyme levels across disease cohorts. It integrates high-throughput 3D organoid metadata, leveraging custom color-coding, shape-based contamination tracking, and patient labeling to generate richly stratified visualizations for downstream applications.

## Overview
This script is designed for quantitative enzyme-linked immunosorbent assays (ELISA). It handles background subtraction, standard curve interpolation using the `drc` package, replicate-level tracking, and clinical metadata integration. 

## Key Features
* **Rigorous QC:** Performs blank-adjustment prior to replicate averaging to minimize background noise.
* **Non-Linear Regression:** Fits a 4-parameter logistic (4PL) model to standard curve concentrations (`STD`) to calculate accurate interpolation limits.
* **Metadata Integration:** Merges plate reads with patient profiles, disease groups, and culture condition metrics.
* **Export Automation:** Automatically exports analytical data frames and high-resolution 300 DPI TIFF visualizations.

## Repository File Structure
```text
elisa-analysis-repo/
├── data/
│   ├── ELISAMASP3006_Results.csv
│   └── PatientSpheroidProfiles.csv
├── output/
│   ├── RAnalysedResults.csv
│   ├── StandardCurve.tiff
│   ├── PLOT_OD_Values.tiff
│   └── PLOTGroupedLog2.tiff
└── elisa_analysis.R
```

## Prerequisites
Ensure the following R packages are installed before running the script:
```R
install.packages(c("dplyr", "ggplot2", "drc", "forcats", "tidyr", "stringr", "ggrepel"))
```

<p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>