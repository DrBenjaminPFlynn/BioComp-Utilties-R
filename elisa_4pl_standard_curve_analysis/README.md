# High-Throughput ELISA & 4PL Regression Analysis Pipeline

An automated R pipeline designed for high-throughput enzyme-linked immunosorbent assay (ELISA) processing. This workflow executes individual replicate-level background corrections, fits non-linear 4-parameter logistic (4PL) standard curves, interpolates target protein concentrations, and integrates multi-variable clinical and organoid metadata.

## Executive Context
Built for multi-cohort biological and clinical research, this repository provides end-to-end analytical rigor for quantifying protein biomarkers (e.g., plasma enzyme levels). By joining raw plate reader measurements with high-throughput 3D organoid phenotypes and cohort metadata, the pipeline automates quality control (QC) and generates publication-grade vector and raster visualizations.

## Key Features & Analytical Capabilities
* **Replicate-Level Background Subtraction:** Applies blank correction ($OD_{450}$) to individual replicates prior to averaging, preventing background noise propagation across downstream metrics.
* **Non-Linear 4PL Modeling:** Leverages `drc` (Dose-Response Curves) to fit a 4-Parameter Logistic (`LL.4`) regression model to standard controls, ensuring accurate calibration across non-linear detection limits.
* **Metadata & Phenotype Fusion:** Integrates optical density metrics with high-dimensional sample metadata (e.g., patient cohorts, disease status, spheroid morphometrics).
* **Automated Dual-Format Graphics Export:** Generates publication-ready figures saved simultaneously as high-resolution 300 DPI rasters (`.tiff`) and fully scalable vectors (`.svg`).
* **Out-of-Range Guardrails:** Handles sub-baseline optical densities and non-convergent estimations by safely flooring non-detectable concentrations at zero.

## Input Data Schema Definition
The pipeline integrates two input CSV files located in the `data/` directory: raw plate reader results (`Sample_ELISA_Results.csv`) and phenotypic metadata (`Sample_3D_Organoid_Phenotypes.csv`).

### 1. Raw ELISA Results Schema (`data/Sample_ELISA_Results.csv`)

| Column Name | Data Type | Requirement / Description |
| :--- | :--- | :--- |
| `Sample` | Text | Primary key / sample designation. **Must include `STD`** for calibration controls. |
| `Conc` | Numeric | Nominal protein concentration (ng/mL) for standard controls (`Sample == "STD"`). Set to `0` for blank control wells. |
| `Dilution.Factor` | Numeric | Fold-dilution multiplier applied during interpolation to restore initial sample concentration. |
| `N1`, `N2` (`N[k]`) | Numeric | Raw absorbance ($OD_{450}$) measurements across technical replicates (e.g., `N1`, `N2`). |

### 2. Phenotypic Metadata Schema (`data/Sample_3D_Organoid_Phenotypes.csv`)

| Column Name | Data Type | Requirement / Description |
| :--- | :--- | :--- |
| `Sample` | Text | Primary key linking metadata records to raw ELISA results (`Sample`). |
| `Spheroid.Profile` | Text | Morphometric or culture profile tag (e.g., growth characteristics). Missing values are automatically assigned `"Unknown"`. |
| `Disease.Group` | Text | Clinical cohort or disease classification tag for downstream stratification. |

## Repository Architecture
```text
elisa-analysis-repo/
├── data/
│   ├── Sample_ELISA_Results.csv              # Raw plate reader OD values & dilutions
│   └── Sample_3D_Organoid_Phenotypes.csv     # Patient cohort & spheroid metadata
├── output/
│   ├── RAnalysedResults.csv                  # Quantified concentrations & merged dataset
│   ├── StandardCurve.tiff                    # 4PL non-linear calibration curve plot
│   ├── PLOT_OD_Values.tiff / .svg            # Quality control OD distributions
│   ├── PLOT_EstConc.tiff / .svg              # Interpolated protein concentrations
│   └── PLOTConcLog2.tiff / .svg              # Log2-transformed concentration profiles
├── elisa_4pl_standard_curve_analysis.R       # Primary execution pipeline
└── elisa_4pl_standard_curve_analysis.Rproj   # RStudio environment project file
```

## Prerequisites & Installation
Ensure R (≥ 4.0.0) is installed along with the required CRAN dependencies:

```R
install.packages(c(
  "dplyr",     # Data manipulation
  "ggplot2",   # Publication graphics
  "drc",       # Dose-response curve fitting (4PL)
  "forcats",   # Factor manipulation
  "tidyr",     # Data tidying
  "stringr",   # String manipulation
  "ggrepel",   # Non-overlapping plot labels
  "svglite"    # High-quality SVG graphics device
))
```

## Pipeline Execution Workflow
```text
[1. Config & Setup] ────────► [2. Data Ingestion] ────────► [3. Blank Correction]
  • Output directory            • Raw ELISA CSV               • Mean blank calculated
  • Reps & target variables     • Metadata CSV                • Subtracted per-replicate

                                                                       │
[6. Export Results] ◄──────── [5. Interpolation] ◄───────── [4. 4PL Regression]
  • CSV dataset                 • ED() limit mapping          • DRM curve fitting
  • 300 DPI TIFF + SVG          • Dilution multiplier         • Calibration plot
```

## Detailed Execution Steps
### 1. Environment Initialization & Configuration
  * Configures analysis parameters (NO.OF.REPS, ELISA.Target) and initializes destination directory structures (output/).

### 2. Data Ingestion & Blank Adjustment
  * Parses raw optical density ($OD_{450}$) values.
  * Calculates the mean background absorbance from blank controls (STD at 0 ng/mL).
  * Subtracts blank values at the individual replicate level ($N_1, N_2, \dots, N_k$) prior to taking mean values to guarantee statistical consistency.

### 3. Metadata Merging & Data Reshaping
  * Merges 3D organoid phenotypic metadata (Sample_3D_Organoid_Phenotypes.csv) using primary key Sample.
  * Reshapes wide-format replicate columns into tidy, long-format data frames (df.ggplot) optimized for ggplot2 plotting routines.

### 4. Quality Control Assessment
  * Generates sample-wise $OD_{450}$ boxplots alongside standard calibration lines (PLOT_OD_Values).
  * Applies viridis color mapping to stratify samples by spheroid morphometrics and detect potential contamination or outliers.

### 5. Non-Linear 4PL Standard Curve Calibration
  * Extracts standard curve control data and fits a 4-Parameter Logistic model using drc::drm():
  
\
$$
y = c + \frac{d - c}{1 + \exp\left(b(\log(x) - \log(e))\right)}$$
<div style="text-align: center;">
  (where $b$ = Hill slope, $c$ = lower asymptote, $d$ = upper asymptote, $e$ = $EC_{50}$).
</div>\


  * Renders and exports the fitted calibration curve (StandardCurve.tiff).

### 6. Concentration Interpolation & Floor Bounding
  * Evaluates unknown sample responses against the fitted 4PL model using drc::ED().
  * Adjusts estimated values by sample-specific dilution factors.Replaces negative or non-estimable values with 0 ng/mL to maintain downstream data integrity.

### 7. Downstream Visualization & Automated Export
  * Generates raw (PLOT_EstConc) and $\log_2(x + 1)$ transformed (PLOTConcLog2) concentration boxplots grouped by clinical cohort.
  * Exports analytical tables (RAnalysedResults.csv) and paired vector (.svg) and high-resolution raster (.tiff, 300 DPI) visual artifacts.License

## License
<p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>