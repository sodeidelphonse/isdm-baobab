

# Reproducibility Package: Integrated Species Distribution Modeling of the African Baobab in Benin.

This repository contains the data and code required to reproduce the analyses presented in the manuscript: *"Integrating Presence-only and Abundance Data to Predict Baobab (*Adansonia digitata* L.) Distribution: A Bayesian Data Fusion Framework"*. 

## 1. Overview
This study utilises a spatial fusion framework via `inlabru` (via `PointedSDMs`) and `isdmtools` to integrate presence-only (e.g., GBIF) and structured abundance data. 
The analysis is conducted across the three climatic zones of Benin (West Africa).

## 2. Requirements
* **Language:** R (v4.4.1)
* **Key Packages:** `isdmtools`, `inlabru`, `INLA`, `PointedSDMs`, `sf`, `pROC`, and `blockCV`. See the scripts. 
* **Environment:** An `renv.lock` file is provided to restore the exact library versions used in this analysis.
* The package `isdmtools` which is not yet available on CRAN can be installed as follows:

```R
install.packages("remotes") 
remotes::install_github("sodeidelphonse/isdmtools@v0.4.0")
```
* The proposed reproducible workflow require the package [`INLA`](https://www.r-inla.org/download-install) version `v24.06.27` 
before expecting reproduce the outputs presented in the paper.
* Though some scripts can work independently, they may require at list the pre-processing step.

## 3. Repository Structure
* **/data**: 
    * `Adansonia_occurrence.csv`: Presence-only records (point pattern data).
    * `Adansonia_abundance.csv`: Village-level count data (point-referenced data).
    * `covariates/`: 5 clipped GeoTIFFs (Final predictors retained after selection).
    * `shapefile/`: The polygon map used as the study region.
* **/scripts**:
    * `01_preprocessing.R`: Importing and cleaning datasets and generating mesh settings.
    * `02_EDA.R`: Perform exploratory analysis, hypothesis testing and data visualisation.
    * `03_evaluation_pipeline.R`: Generate the ISDM pipeline from spatial blocking, fitting, prediction to evaluation.
    * `04_run_blockCV.R`: Run the block cross-validation process for the integrated modeling workflow.
    * `05_model_fitting.R`: INLA-SPDE integrated modeling workflow for the selected model.
    * `06_prediction.R`: Make prediction from the fitted model and compute relevant target quantities
    * `07_predictive_check.R`: Posterior predictive check from the integrated modelling workflow.
* **/figures**: The figures generated from the analysis.
* **/results**: Contains the pre-computed covariance parameters for `LGCP` and `variofit` models (voir `02_EDA.R`)
* **/software**: The archived version of `isdmtools` package used for preparing the manuscript.
* `README.md`: Project documentation.
* `renv.lock`: An environment file provided to restore the exact library versions used in this analysis.

## 4. Data Notes
* **Covariates (covariates_pc.rds):** We provide the five final environmental variables used in the model (clipped to the study area) as a serialised format.
* **Prediction grids (pred_points.rds):** They are regular grid points locations generated for the model prediction.

## 5. Exploratory Data Analysis (EDA) 📊 
The script **`02-EDA.R`** performs the initial variogram and spatial point pattern analyses and requires the script **`01-preprocessing.R`**. 

- **Parametric Bootstrapping**: To assess the uncertainty of covariance parameters (estimated via `kppm` and `variofit`), 
1,000 parametric bootstrap replicates are generated.

- **Fast-Track**: To save time, pre-calculated replicates are provided in:
  - `outputs/simulations_lgcp.csv`
  - `outputs/simulations_variofit.csv`
- **Purpose**: These replicates establish the baseline spatial range used to inform the SPDE priors in the integrated models.

## 6. Contact
For questions regarding the `isdmtools` implementation or data processing, please contact the corresponding author or open an issue on the package repository.

## 7. License
This repository is released under the [MIT License](LICENSE).
