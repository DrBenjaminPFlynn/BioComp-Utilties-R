# BioComp-Utilities-R

A curated suite of modular, production-ready R data pipelines built for biological data analysis, immunodiagnostic quantification, and translational clinical research workflows.

## Executive Context
This repository houses standardized end-to-end analytical pipelines for processing raw microplate reader data, executing mathematical regression modeling, integrating high-dimensional clinical metadata, and generating publication-grade vector and raster visualizations. 

All sub-projects adhere to strict software design standards:
* **Strict Separation of Concerns:** Standardized directory hierarchy (`data/`, `output/`) across all tools.
* **Dual-Format Publication Graphics:** Automated rendering of 300 DPI rasters (`.tiff`) alongside fully scalable vector graphics (`.svg`).
* **Mathematical Guardrails:** Boundary-capping and fallback logic for non-convergent, sub-baseline, or out-of-range experimental values.

---

## Featured Pipelines

### 1. High-Throughput ELISA & 4PL Regression Analysis (`elisa-analysis-pipeline/`)
An automated pipeline for quantifying target biomarker concentrations from ELISA plate reader absorbances ($OD_{450}$).

* **Core Methodology:** Non-Linear 4-Parameter Logistic (`LL.4`) regression modeling using `drc::drm()`.
* **Key Features:**
  * Replicate-level background ($OD_{450}$) subtraction prior to averaging to prevent error propagation.
  * Integration with 3D organoid culture phenotypic metadata (e.g., spheroid morphometrics, disease cohorts).
  * Automated interpolation via `drc::ED()` with floor-bounding at $0\text{ ng/mL}$ for non-detects.
* **Outputs:** Quality control distribution plots, 4PL calibration curves, raw/$\log_2$ concentration profiles, and merged analytical data frames.

### 2. Complement Hemolytic Assay (CH50 / AP50) Analysis (`hemolytic-assay-analysis/`)
An automated pipeline for evaluating functional classical (CH50) and alternative (AP50) complement pathway activity in human plasma.

* **Core Methodology:** Ordinary Least Squares (OLS) log-linear regression ($Haemolysis \sim \ln(PlasmaFraction)$).
* **Key Features:**
  * Dynamic assay detection (CH50 vs. AP50) based on minimum dilution thresholds with automated coordinate scaling.
  * Baseline control normalization (`NegCon` = 0% lysis, `PosCon` = 100% lysis) to relative hemolysis fractions.
  * Interpolation of $HAEM_{50}$ (serum fraction required for 50% erythrocyte lysis) and clinical status classification (*Normal*, *Low/Borderline*, *Deficient*).
  * Cohort visualization with inverted Y-axis scaling (`scale_y_reverse`), accurately reflecting that lower required plasma volumes denote higher biological activity.
* **Outputs:** Inverted cohort stem/lollipop plots, individual sample dose-response curves with target crosshairs, and extracted model parameter tables.

---

## Repository Architecture

```text
BioComp-Utilities-R/
├── elisa-analysis-pipeline/
│   ├── data/
│   │   ├── Sample_ELISA_Results.csv           # Raw ELISA OD values & dilutions
│   │   └── Sample_3D_Organoid_Phenotypes.csv # Organoid phenotype metadata
│   ├── output/
│   │   ├── RAnalysedResults.csv               # Merged analytical dataset
│   │   ├── StandardCurve.tiff                # 4PL calibration plot
│   │   ├── PLOT_OD_Values.tiff / .svg        # QC OD boxplots
│   │   ├── PLOT_EstConc.tiff / .svg          # Interpolated concentrations
│   │   └── PLOTConcLog2.tiff / .svg          # Log2 concentration profiles
│   ├── elisa_4pl_standard_curve_analysis.R    # Primary execution script
│   └── elisa_4pl_standard_curve_analysis.Rproj# RStudio project file
│
├── hemolytic-assay-analysis/
│   ├── data/
│   │   └── CH50_assay_data.csv             # Raw hemolytic plate reader data
│   ├── output/
│   │   ├── Rresults.csv                    # Extracted HAEM_50 values & model stats
│   │   ├── HAEM_50Results.tiff / .svg      # Cohort summary (inverted Y-axis)
│   │   ├── CH50Curve_[ID]_[Sample].tiff    # Sample-level dose-response curves (TIFF)
│   │   └── CH50Curve_[ID]_[Sample].svg     # Sample-level dose-response curves (SVG)
│   ├── haemolysis_analysis.R               # Primary execution script
│   └── hemolysis_analysis.Rproj            # RStudio project file
│
└── README.md
```

## Prerequisites & Dependencies
Both pipelines require R (≥ 4.0.0). You can install all required packages across both pipelines with the following command:

```R
install.packages(c(
  "tidyverse", # Core suite: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats
  "drc",       # Non-linear dose-response curve fitting (4PL)
  "ggrepel",   # Non-overlapping plot text labels
  "svglite"    # High-performance SVG graphics device
))
```
Note: Detailed execution steps, mathematical formulas, and input data schemas are documented within the respective README.md files located in each pipeline directory.

## License
<p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>