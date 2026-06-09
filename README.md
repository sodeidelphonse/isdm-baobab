

# Reproducibility Package: Integrated Species Distribution Modelling of the African Baobab in Benin.

This repository contains the data and code required to reproduce the analyses presented in the paper: 

Sode, A.I., Fandohan, A.B., Krainski, E.T., Assogbadjo, A.E., & Glèlè Kakaï, R. (2026). Integrating presence-only and abundance data to predict baobab (*Adansonia digitata L.*) distribution: a Bayesian data fusion framework. *Environmental and Ecological Statistics*. https://doi.org/10.1007/s10651-026-00737-2 

## Table of contents
* [Overview](#overview)
* [Quick Start](#quick-start)
* [Requirements](#requirements)
* [Repository Structure](#repository-structure)
* [Environment Configuration](#environment-configuration)
* [Pipeline Logic](#pipeline-logic)
* [Spatial Data Notes](#spatial-data-notes)
* [Exploratory Analysis Outputs](#exploratory-analysis-outputs)
* [License](#license)
* [Citation](#citation)
* [Contact & Support](#contact-support)

## 📝 Overview 
The study utilises a Bayesian spatial fusion framework with `inlabru` and `isdmtools` to integrate presence-only (GBIF and field records) and structured abundance data of the African baobab. 
The aim is to comprehend and map the spatial variation of this multipurpose agroforestry tree species across the three climatic zones of Benin (West Africa).

## 🚀 Quick Start
1. Open the `.Rproj` file.
2. Run `renv::restore()`.
3. Run `scripts/01_preprocessing.R`.
4. Follow scripts 02–07 in order.

## 🛠 Requirements ️
* **Language:** R (v4.4.1)
* **Key Packages:** `isdmtools`, `inlabru`, `PointedSDMs`, `INLA`, `sf`, `pROC`, and `blockCV`. See the scripts for other required packages. 
* **Environment:** `renv.lock` file is provided to restore the exact library versions used in this analysis.
* **`isdmtools`:** which is not yet available on CRAN can be installed as follows:

```R
install.packages("remotes") 
remotes::install_github("sodeidelphonse/isdmtools@v0.4.0")
```
* **`INLA`:** This workflow requires the package [INLA](https://www.r-inla.org/download/index.html) `v24.06.27` in order to reproduce the outputs presented in the paper.

## 📂 Repository Structure
* **/data**: 
    * `Adansonia_occurrence.csv`: Presence-only records (point pattern data) with 1-km jittered coordinates.
    * `Adansonia_abundance.csv`: Site-level counts (point-referenced data) with 1-km jittered coordinates.
    * `covariates/`: Five clipped `GeoTIFF` files representing the final predictors retained after selection.
    * `shapefile/`: The polygon map used as the study region and the vector lines delineating the three climatic zones.
    * `covariates_pc.rds`: The five final environmental variables stored in a serialised format.
    * `pred_points.rds`: The regular grid points locations stored in a serialised format for the model prediction.

* **/scripts**:
    * `01_preprocessing.R`: Import and clean datasets and generate mesh configurations.
    * `02_EDA.R`: Perform exploratory analysis, particularly spatial dependence assessment, hypothesis testing, and data visualisation.
    * `03_evaluation_pipeline.R`: Set up the ISDM pipeline from spatial blocking, fitting, prediction to evaluation.
    * `04_run_blockCV.R`: Run the block cross-validation strategy for the integrated modelling workflow.
    * `05_model_fitting.R`: Run the INLA-SPDE integrated modelling workflow for the selected model using the full datasets.
    * `06_prediction.R`: Predict the intensity surface from the fitted model and compute relevant target quantities.
    * `07_predictive_check.R`: Perform posterior predictive check (PPC) from the integrated models.
    * `08_utils.R`: Utility functions used in other scripts during the analysis.
    
* **/figures**: The figures generated from the analysis using the raw datasets. 

* **/results**: Contains the pre-computed covariance parameters for `LGCP` and `variofit` models (see `02_EDA.R`)

* **/software**: The archive of the `isdmtools` package used for preparing the manuscript. Both binary and source codes are provided.

* `README.md`: The project documentation (this page).

## 📦 Environment Configuration
To reproduce the environment used for the analysis:

* Make sure you have the `renv` package installed by running `install.packages("renv")`. 

* With the project root as your working directory, run `renv::restore()` to automatically download and install all required library versions.

## ⚙️ Pipeline Logic
This analysis is designed as a sequential pipeline (01–07). Each stage relies on the outputs of the preceding stages and the shared utility functions:

* **Sequential Flow**: To reproduce the full results, it is recommended to run the scripts in numerical order.

* **Cross-Stage Dependencies**:
    * Exploratory Data Analysis (Stage 02): Utilises processed data from stage 01 to perform the initial variogram and spatial point pattern analyses as well as data visualisation. 
    
    * Integrated Model Evaluation (Stage 04): Utilises processed data from 01 and the evaluation pipeline sourced from 03.
    
    * Modelling (Stage 05): Requires the environmental covariates' stack, datasets and mesh objects prepared in 01.

    * Prediction (Stage 06): Requires the fitted model objects from 05 and the prediction locations to be imported.
    
    * Validation(Stage 07): Utilises processed data from 01 to rerun models with `inlabru` and perform residual diagnostics.
    
## 💾 Spatial Data Notes

* Due to serialisation constraints common with high-resolution rasters and complex spatial model outputs (e.g., `terra` and `INLA` objects), we recommend the following approach for full reproducibility:

  * Re-running Pre-processing: It is strongly recommended that you re-run `01_preprocessing.R` step, particularly after restarting your R session. 
  This ensures that the environmental covariates are correctly loaded into your local R memory for subsequent stages.

  * Raw Data Access: The original environmental layers are provided in `data/covariates/` to demonstrate the transition from raw geospatial data to the analysis-ready stacks.

  * Serialised Alternative: Pre-processed covariates are available as `covariates_pc.rds`. If using it, ensure that object names are mapped according to the naming conventions established in all scripts.

* To facilitate public distribution and ensure the privacy of surveyed settlements is protected, the datasets provided in the `data/` folder have been de-identified using a 1-km jitter. 
While this does not significantly alter the scale of the predictions, users may notice minor variations if they re-run the distance and plotting scripts.

## 📊 Exploratory Analysis Outputs

The Exploratory Data Analysis (Stage 02) produces two critical outputs, among others.

- **Bootstrap Replicates**: To assess the uncertainty of covariance parameters (estimated via `kppm` and `variofit`), 
1000 parametric bootstrap replicates are generated.

- **Fast-Track**: To save time, pre-calculated replicates for these parameters are provided in:
  - `results/simulations_lgcp.csv`
  - `results/simulations_variofit.csv`
  
- **Purpose**: These estimates establish the baseline spatial range used to inform the SPDE priors in the integrated models.

## ⚖️ License
This repository is released under the [MIT License](LICENSE).

## 📑 Citation 
If you use this repository, the underlying data, or the specific analysis pipeline in your research, please cite them as follows:

- **The Core Methodology & Application Paper**: 
Our foundational study on the African baobab datasets in Benin is accepted for publication in *Environmental and Ecological Statistics*:

  Sode, A.I., Fandohan, A.B., Krainski, E.T., Assogbadjo, A.E., & Glèlè Kakaï, R. (2026). Integrating presence-only and abundance data to predict baobab (*Adansonia digitata* L.) distribution: a Bayesian data fusion framework. *Environmental and Ecological Statistics*. https://doi.org/10.1007/s10651-026-00737-2

- **The Project Archive**: 
The complete reproducible archive for this study—including the 1-km jittered datasets, repository export, `isdmtools v0.4.0` source code, pre-computed outputs, and figures—is permanently archived on Zenodo:

  Sode, A.I., Fandohan, A.B., Krainski, E.T., Assogbadjo, A.E., & Glèlè Kakaï, R. (2026). Research compendium for integrating and predicting baobab (*Adansonia digitata* L.) distribution in Benin using a Bayesian data fusion framework: Data, pipelines, and isdmtools v0.4.0 source code (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.19227943

- **The `isdmtools` R Package**: 
The underlying computational infrastructure used for multisource spatial data resampling and ISDM evaluation is part of the `isdmtools` R package. A dedicated software manuscript is currently in preparation. To cite the package software itself with its latest version, please use:

```R
citation("isdmtools")
```
## 📧 Contact & Support
The best way to get help with `isdmtools` or the data analysis pipeline is to open an issue on the package [issue tracker](https://github.com/sodeidelphonse/isdmtools/issues). This allows the community to benefit from the discussion.

For private inquiries or data sharing questions that cannot be posted publicly, you may contact the corresponding author at [sdidelphonse@gmail.com](mailto:sdidelphonse@gmail.com).

