# Kinetic Protease Peptide Probe Assay Pipeline

An automated R pipeline for processing kinetic spectrophotometric microplate reader data from chromogenic protease cleavage assays. This workflow executes dual-stage baseline corrections, fits ordinary least squares (OLS) linear progress curves to quantify reaction velocities ($\Delta OD_{405} / \Delta t$), calculates molar specific enzymatic activity ($\text{pmol}\cdot\text{min}^{-1}\cdot\mu\text{g}^{-1}$), and generates publication-grade cohort visualizations.

## Executive Context
Designed for translational enzymology and clinical biomarker discovery, this repository provides analytical rigor for quantifying target protease cleavage rates (e.g., enzymatic activity using thioester substrate #Z-L-Lys-SBzl_HCl and DTNB/Ellman's reagent). By linking raw kinetic absorbance profiles with plate layout metadata and physical reaction constants, the pipeline automates kinetic rate calculations and exports dual-format vector and raster figures.

## Key Features & Analytical Capabilities
* **Dual-Stage Baseline Adjustments:** Corrects raw kinetic absorbance ($OD_{405}$) by subtracting time-matched blank control wells (`BLANK`) and optional time-zero ($T_0$) background baselines across technical replicates.
* **OLS Kinetic Rate Modeling:** Fits ordinary least squares (OLS) linear progress curves ($OD_{405} \sim Time$) per sample to extract steady-state initial reaction velocities ($\text{Slope} = \Delta OD / \Delta t$).
* **Molar Specific Activity Quantitation:** Converts absorbance slopes ($\text{OD/min}$) directly into molar specific activity ($\text{pmol}\cdot\text{min}^{-1}\cdot\mu\text{g}^{-1}$) using Beer-Lambert optical constants (extinction coefficient $\epsilon$, path length $l$, volume $V$, and enzyme mass $m$).
* **Non-Linear & Linear Curve Fitting:** Renders LOESS progress curves with integrated text labels (`geomtextpath`) alongside linear velocity slopes for kinetic diagnostic inspection.
* **Automated Dual-Format Visual Export:** Automatically generates and saves cohort-level boxplots and kinetic progress curves as high-resolution 300 DPI rasters (`.tiff`) and scalable vectors (`.svg`).

## Input Data Schema Definition
The pipeline processes two CSV files located in the `data/` directory: raw kinetic plate reader measurements (`Data.csv`) and sample annotation maps (`PlateMap.csv`).

### 1. Plate Mapping Schema (`data/PlateMap.csv`)

| Column Name | Data Type | Requirement / Description |
| :--- | :--- | :--- |
| `Wells` | Text | Microplate well coordinate (e.g., `A1`, `B2`). Primary key linking mapping to raw OD data. |
| `Tx` | Text | Sample/treatment identifier. **Must include `BLANK`** for background baseline wells. |
| `Group` | Text | Clinical cohort or experimental group designation (e.g., `Control`, `Cohort_A`). |

### 2. Kinetic OD Read Schema (`data/Data.csv`)

| Column Name | Data Type | Requirement / Description |
| :--- | :--- | :--- |
| `Wells` | Text | Microplate well coordinate matching `PlateMap.csv`. |
| `Time` | Numeric | Reaction elapsed time in minutes (e.g., `0`, `0.5`, `1.0`, ...). |
| `N1`, `N2`, `N3` (`N[k]`) | Numeric | Raw optical density ($OD_{405}$) measurements across technical replicates. |

## Repository Architecture
```text
protease-peptide-assay
├── data/
│   ├── Data.csv                                  # Raw kinetic microplate reader OD405 readings
│   └── PlateMap.csv                              # Well-level treatment and cohort metadata
├── output/
│   ├── Calculated_Protease_Activity.csv          # Extracted slopes and specific activity metrics
│   ├── B.Adjusted.OD_CurveOfBestFit.tiff         # LOESS progress curves (Blank-adjusted, TIFF)
│   ├── B.Adjusted.OD_CurveOfBestFit.svg          # LOESS progress curves (Blank-adjusted, SVG)
│   ├── B.Adjusted.OD_LineOfBestFit.tiff          # OLS linear velocity fits (TIFF)
│   ├── B.Adjusted.OD_LineOfBestFit.svg           # OLS linear velocity fits (SVG)
│   ├── Protease_Specific_Activity_BoxPlot.tiff   # Cohort specific activity boxplot (TIFF)
│   └── Protease_Specific_Activity_BoxPlot.svg    # Cohort specific activity boxplot (SVG)
├── protease_peptide_assay_analysis.R             # Primary processing script
└── protease_peptide_assay.Rproj                  # RStudio environment project file
```

## Prerequisites & Installation
Ensure R (≥ 4.0.0) is installed along with the required CRAN packages:
```R
install.packages(c(
  "dplyr",         # Data manipulation
  "tidyr",         # Data reshaping
  "ggplot2",       # Publication graphics
  "ggrepel",       # Repulsive text annotations
  "geomtextpath",  # Direct text labeling along ggplot paths
  "stringr",       # String utilities
  "svglite"        # Vector graphics device
))
```

## Pipeline Execution Workflow
```text
[1. Config & Inputs] ────────► [2. Dual Baseline Subtraction] ──► [3. Long-Format Reshaping]
  • Reaction constants          • Time-matched BLANK correction  • Replicate expansion
  • Data & PlateMap CSVs        • Time-Zero (T0) adjustment      • Tidy dataframe generation

                                                                         │
[6. Dual Graphics Export] ◄─── [5. Molar Specific Activity] ◄─── [4. OLS Kinetic Modeling]
  • 300 DPI TIFF & SVG          • Beer-Lambert transformation    • Slope (ΔOD/Δt) estimation
  • Progress & Cohort plots     • Output CSV generation          • Linear rate extraction
```

## Detailed Execution Steps
### 1. Configuration & Reaction Parameter Initialization
  * Verifies workspace paths and sets physical constants: well volume ($V = 0.1\text{ mL}$), TNB extinction coefficient ($\epsilon_{405} = 13,260\text{ M}^{-1}\text{cm}^{-1}$), path length ($l = 0.32\text{ cm}$), and loaded enzyme mass ($m = 0.05\text{ }\mu\text{g}$).

### 2. Data Ingestion & Baseline Corrections
  * Merges raw time-series absorbance ($OD_{405}$) with plate annotations.
  * Performs replicate-wise subtraction of time-matched BLANK absorbance values.
  * Subtracts time-zero ($T_0$) background absorbances to eliminate initial substrate background noise.

### 3. Data Tidying & Reshaping
  * Transforms wide-format replicate columns ($N_1, N_2, N_3$) into a tidy long-format dataset (df_ggplot) optimized for linear modeling and ggplot2 rendering routines.
  
### 4. OLS Kinetic Rate Modeling
  * Fits linear regression models ($OD_{405} \sim Time$) per treatment sample to extract reaction rates ($\text{Slope} = \Delta OD_{405} / \Delta t$).
  
### 5. Specific Activity Quantification
  * Converts reaction slopes to specific enzymatic activity ($\text{pmol}\cdot\text{min}^{-1}\cdot\mu\text{g}^{-1}$) using the Beer-Lambert relationship:
  
  \
  $$\text{Specific Activity} = \frac{\text{Slope (OD/min)} \cdot V \cdot 10^{12}}{\epsilon \cdot l \cdot m}$$
  \
  
  * Exports extracted rates and calculated activity metrics to Calculated_Protease_Activity.csv.
  
### 6. Publication Plot Generation & Export
  * Renders non-linear LOESS diagnostic progress curves (geomtextpath) and OLS linear velocity plots.
  * Visualizes cohort-level specific activity using repelled sample labels (ggrepel) and boxplots.
  * Saves graphics simultaneously in high-resolution raster (.tiff, 300 DPI) and scalable vector (.svg) formats.
  
## License
<p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>