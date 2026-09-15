# Complement Hemolytic Assay (CH50 / AP50) Automated Analysis Pipeline

An automated R pipeline for processing raw spectrophotometric plate reader data from complement-mediated hemolytic assays. This workflow normalizes erythrocyte lysis against baseline controls, fits sample-level log-linear regression models, quantifies the target plasma fraction required to achieve 50% hemolysis ($HAEM_{50}$), and generates publication-ready quality control and cohort visualizations.

## Executive Context
Built for functional immunology and clinical complement diagnostics (classical pathway CH50 & alternative pathway AP50), this repository provides end-to-end analytical rigor for evaluating complement system activity. By coupling raw optical density ($OD$) measurements with automated control normalization and clinical reference stratification, the pipeline streamlines quality control and generates publication-grade vector and raster visualizations with inverted activity scaling.

## Key Features & Analytical Capabilities
* **Dynamic Assay Recognition:** Automatically detects whether input data represents a Classical (CH50) or Alternative (AP50) pathway assay based on minimum dilution fraction thresholds, dynamically adjusting plot coordinates and axis breaks.
* **Control Normalization & Scaling:** Normalizes raw absorbance readings ($OD$) against experimental baseline controls (`NegCon` = 0% lysis, `PosCon` = 100% lysis) to yield standardized relative hemolysis fractions ($0.0 - 1.0$).
* **Log-Linear Regression Modeling:** Fits ordinary least squares (OLS) linear models ($Haemolysis \sim \ln(PlasmaFraction)$) per sample to interpolate exact $HAEM_{50}$ values.
* **Clinical Activity Stratification:** Evaluates calculated $HAEM_{50}$ percentages against patient reference controls (`Normal` vs. `Low`) to automatically classify samples into *Normal*, *Low/Borderline*, or *Deficient* functional complement states.
* **Biological Activity Axis Inversion:** Visualizes cohort comparisons on an inverted Y-axis (`scale_y_reverse`), accurately reflecting that lower required serum volumes correspond to higher biological complement activity.
* **Dual-Format Graphics Export Pipeline:** Automatically renders and exports both cohort-level summary plots and individual sample dose-response curves as high-resolution 300 DPI rasters (`.tiff`) and scalable vectors (`.svg`).

## Input Data Schema Definition
The pipeline processes raw CSV files located in the `data/` directory (e.g., `data/CH50_assay_data.csv`). The input CSV must contain three metadata columns followed by numerical serial dilution headers.

```text
| Column Name | Data Type | Requirement / Description |
| :--- | :--- | :--- |
| `ID` | Numeric / Text | Run or row identifier. |
| `Patient` | Text | Patient or control tag. **Must include `NegCon` (0% baseline) and `PosCon` (100% baseline)**. Optional references: `Normal Ref`, `Low Ref`. |
| `Sample` | Text | Sample subgroup, designation, or replicate tag (e.g., `Control`, `Ref`, `A`, `B`). |
| `<Dilution Headers>` | Numeric | Absorbance ($OD$) values for each serial dilution. Column headers must be numeric fold-dilutions (e.g., `128`, `64`, `32`). |
```

## Repository Architecture

```text
hemolytic-assay-analysis/
├── data/
│   └── CH50_assay_data.csv             # Raw spectrophotometric plate reader dataset
├── output/
│   ├── Rresults.csv                    # Extracted model parameters, HAEM_50 values & clinical status
│   ├── HAEM_50Results.tiff / .svg      # Cohort summary plot (inverted activity axis)
│   ├── CH50Curve_[Patient]_[Sample].tiff # Individual fitted curves with target crosshairs (TIFF)
│   └── CH50Curve_[Patient]_[Sample].svg  # Individual fitted curves (SVG)
├── haemolysis_analysis.R               # Primary processing and visualization script
└── hemolysis_analysis.Rproj            # RStudio environment project file
```

## Prerequisites & Installation
Ensure R (≥ 4.0.0) is installed along with the core tidyverse library:

```R
install.packages("tidyverse") # Includes ggplot2, dplyr, tidyr, readr, and purrr
```

## Pipeline Execution Workflow
```text
[1. Config & Setup] ────────► [2. Data Ingestion & Prep] ────► [3. Control Normalization]
  • Directory verification      • Pivot numeric dilutions       • Baseline mean (NegCon/PosCon)
  • Input CSV validation        • Plasma fraction calculation   • Scaled hemolysis fraction

                                                                       │
[6. Plotting & Export] ◄────── [5. Cohort Stratification] ◄─── [4. Log-Linear Modeling]
  • Inverted activity summary   • Dynamic reference threshold   • OLS regression fitting
  • Individual 50% target curves• Normal/Low/Deficient tag      • HAEM_50 percentage solution
```

## Detailed Execution Steps
### 1. Enviroment & File Validation
  * Verifies required directory structures (data/, output/) and confirms the presence of the input dataset (CH50_assay_data.csv).
### 2. Data Reshaping & Assay Auto-Detection
  * Pivots raw dilution header columns into a tidy long format.Converts fold-dilutions to serum concentration fractions: $PlasmaFraction = \frac{1}{Dilution}$.
  * Evaluates the minimum plasma fraction across the dataset to dynamically set the assay mode to CH50 ($< 0.02$) or AP50 ($\ge 0.02$), automatically adjusting downstream plotting limits and axis breaks.
### 3. Control Normalisation & Hemolysis Scaling
  * Isolates control wells (NegCon and PosCon) to compute baseline mean optical densities.
  * Normalizes raw sample optical density ($OD$) into a relative hemolysis fraction ($0.0 - 1.0$):
  
  $$Haemolysis = \frac{OD - OD_{NegCon}}{OD_{PosCon} - OD_{NegCon}}$$
  
### 4. Log-Linear Regression & $HAEM_{50}$ Calculation
  * Fits an ordinary least squares regression model per sample: $Haemolysis = a \cdot \ln(PlasmaFraction) + b$.
  * Solves for the exact plasma percentage ($HAEM_{50}$) required to induce 50% cell lysis ($y = 0.5$):
  
  $$HAEM_{50} = \exp\left(\frac{0.5 - b}{a}\right) \times 100$$
  
  * Implements guardrails for mathematical edge cases (e.g., negative slopes, infinite or non-convergent calculations) by floor/ceiling capping non-calculable values at $100\%$.
### 5. Cohort Reference Stratification
  * Dynamically extracts reference thresholds (normal_ref_val, low_ref_val) from control samples containing Normal or Low designations (with safe fallbacks of $25\%$ and $55\%$).
  * Stratifies test samples into Normal, Low/Borderline, or Deficient functional complement categories. 
### 6. Automated Visualization & Multi-Format Rendering
  * Cohort Summary Plot: Generates stem/lollipop plots displaying $HAEM_{50}$ values per patient. Uses an inverted Y-axis (scale_y_reverse) so that higher biological activity (lower required plasma fraction) is plotted vertically higher.
  * Individual Dose-Response Curves: Loops through individual samples to render linear regression curves over log plasma fractions, complete with dashed red target crosshairs at $y = 0.5$ and annotated $HAEM_{50}$ callout badges.
  * Exports analytical summary tables (Rresults.csv) alongside paired high-resolution raster (.tiff, 300 DPI) and vector (.svg) image files.
  
## License
<p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>