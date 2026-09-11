# BioComp-Utilities-R

A collection of clean, modular, and production-ready R scripts and data pipelines developed for biological data analysis, translational clinical research workflows, and portfolio presentation. 

## Featured Projects

### ELISA 4PL Analysis Pipeline (`elisa-analysis-pipeline/`)
An automated R pipeline designed for multi-variable clinical studies to quantify immunological enzyme levels in patient plasma. 
* **Key Capabilities:** Performs background adjustments, fits non-linear 4-parameter logistic (4PL) standard curves using the `drc` package, and integrates high-throughput 3D organoid cell culture metadata.
* **Visualization:** Utilizes custom aesthetic mapping, shape-based contamination tracking, and advanced labeling (`ggrepel`) to produce richly stratified figures, exported automatically in both high-resolution raster (`.tiff`) and vector (`.svg`) formats.

## Repository File Structure
```text
BioComp-Utilities-R/
├── elisa-analysis-pipeline/
│   ├── data/
│   │   ├── ELISAMASP3006_Results.csv
│   │   └── PatientSpheroidProfiles.csv
│   ├── output/
│   │   ├── RAnalysedResults.csv
│   │   ├── StandardCurve.tiff
│   │   ├── StandardCurve.svg
│   │   ├── PLOT_OD_Values.tiff
│   │   ├── PLOT_OD_Values.svg
│   │   ├── PLOTGroupedLog2.tiff
│   │   └── PLOTGroupedLog2.svg
│   ├── elisa_analysis_pipeline.Rproj
│   └── elisa_analysis.R
└── README.md
```

## Prerequisites
Please see README.md within subdirectories for details.

<p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="https://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""></a></p>