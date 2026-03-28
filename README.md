

# Reproducibility Package: Integrated Species Distribution Modeling of the African Baobab in Benin.

This repository contains the data and code required to reproduce the analyses presented in the manuscript: *"Integrating Presence-only and Abundance Data to Predict Baobab (*Adansonia digitata* L.) Distribution: A Bayesian Data Fusion Framework"*. 

## 📝 Overview 
The study utilises a Bayesian spatial fusion framework with `inlabru` (via `PointedSDMs`) and `isdmtools` to integrate presence-only (GBIF and field records) and structured abundance data of the African baobab. 
The aim is to comprehend and map the spatial variation of this multipurpose agroforestry tree species across the three climatic zones of Benin (West Africa).

## 🛠 Requirements ️
* **Language:** R (v4.4.1)
* **Key Packages:** `isdmtools`, `inlabru`, `INLA`, `PointedSDMs`, `sf`, `pROC`, and `blockCV`. See the scripts for other required packages. 
* **Environment:** `renv.lock` file is provided to restore the exact library versions used in this analysis.
* The package `isdmtools` which is not yet available on CRAN can be installed as follows:

```R
install.packages("remotes") 
remotes::install_github("sodeidelphonse/isdmtools@v0.4.0")
```
* The proposed workflow requires the package [`INLA`](https://www.r-inla.org/download-install) `v24.06.27` 
in order to reproduce the outputs presented in the paper.
* Although some scripts can work independently, they may require at least the pre-processing step.

## 📂 Repository Structure
* **/data**: 
    * `Adansonia_occurrence.csv`: Presence-only records (point pattern data).
    * `Adansonia_abundance.csv`: Site-level counts (point-referenced data).
    * `covariates/`: Five clipped GeoTIFF files representing the final predictors retained after selection.
    * `shapefile/`: The polygon map used as the study region and the vector lines delineating the three climatic zones.
    * `covariates_pc.rds`: The five final environmental variables stored in a serialized format.
    * `pred_points.rds`: The regular grid points locations stored in a serialized format for the model prediction.

* **/scripts**:
    * `01_preprocessing.R`: Import and clean datasets and generate the mesh configurations.
    * `02_EDA.R`: Perform exploratory analysis, particularly spatial dependence assessment, hypothesis testing, and data visualisation.
    * `03_evaluation_pipeline.R`: Set up the ISDM pipeline from spatial blocking, fitting, prediction to evaluation.
    * `04_run_blockCV.R`: Run the block cross-validation strategy for the integrated modeling workflow.
    * `05_model_fitting.R`: INLA-SPDE integrated modeling workflow for the selected model using the original datasets.
    * `06_prediction.R`: Make prediction from the fitted model and compute relevant target quantities.
    * `07_predictive_check.R`: Perform posterior predictive check (PPC) from the integrated models.
    * `08_utils.R`: Utility functions used in other scripts.
    
* **/figures**: The figures generated from the analysis.

* **/results**: Contains the pre-computed covariance parameters for `LGCP` and `variofit` models (see `02_EDA.R`)

* **/software**: The archive of `isdmtools` package used for preparing the manuscript.

* `README.md`: The project documentation.

* `renv.lock`: A virtual environment file provided to restore the exact library versions used in this analysis.

## 📦 Reproducing the virtual environment

- Make sure you have the `renv` package installed:

```R
install.packages("renv")
```

- With the project directory as your working directory, use `renv` to install all packages listed in the `renv.lock` file:

```R
renv::restore()
```

## 📊 Exploratory Data Analysis (EDA) 
The script **`02-EDA.R`** performs the initial variogram and spatial point pattern analyses and requires the script **`01-preprocessing.R`**. 

- **Parametric Bootstrapping**: To assess the uncertainty of covariance parameters (estimated via `kppm` and `variofit`), 
1,000 parametric bootstrap replicates are generated.

- **Fast-Track**: To save time, pre-calculated replicates for these parameters are provided in:
  - `results/simulations_lgcp.csv`
  - `results/simulations_variofit.csv`
- **Purpose**: These replicates establish the baseline spatial range used to inform the SPDE priors in the integrated models.

## 📧 Contact
For questions regarding the [`isdmtools`](https://sodeidelphonse.github.io/isdmtools/) implementation or data processing, 
please contact the corresponding author or open an issue on the package [repository](https://github.com/sodeidelphonse/isdmtools/issues).

## ⚖️ License
This repository is released under the [MIT License](LICENSE).
