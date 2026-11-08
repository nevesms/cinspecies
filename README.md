# cinspecies: Automated Ecological Niche Modeling and Realized Habitat Assessment

## Overview

The `cinspecies` package provides an end-to-end, reproducible framework for Ecological Niche Modeling (ENM) and Species Distribution Modeling (SDM) in R. Traditional modeling pipelines frequently terminate at the projection of potential abiotically suitable areas, overlooking anthropogenic landscape conversion and biotic constraints that govern actual species occupancy in the field. This software bridges the gap between potential Grinnellian niches and realized distributions by coupling multi-algorithm machine learning ensembles with automated Land Use and Land Cover (LULC) post-processing, recursive biotic interaction modeling, and geodetic area accounting across baseline and future climate scenarios.

The package standardizes complex spatial and statistical operations into a single master orchestrator function (`cins_pipeline`) while preserving complete user control over ecological assumptions, data filtering, and hyperparameter tuning.

## Databases and Coupled R Architecture

The internal architecture of `cinspecies` integrates specialized spatial, ecological, and statistical R packages with global environmental repositories to execute five sequential modeling stages:

### 1. Biodiversity Records and Spatial Filtering
Occurrence data are retrieved dynamically from the Global Biodiversity Information Facility (GBIF) via the `rgbif` package. Raw coordinates undergo automated quality control using `CoordinateCleaner` to eliminate fossil records, institutional coordinates (museums, herbaria, universities), country/province centroids, zero-coordinate artifacts, and sea-land mismatches. To reduce spatial autocorrelation and sampling bias, records are filtered using geographic distance thinning via `spThin`. When altitudinal limits are specified, the pipeline downloads Digital Elevation Models (DEM) derived from the Shuttle Radar Topography Mission (SRTM) through `geodata` and excludes records falling outside the focal species elevation envelope.

### 2. Baseline and Future Environmental Predictors
For terrestrial organisms, baseline bioclimatic variables (BIO1 through BIO19, representing annual trends, seasonality, and extreme temperature and precipitation conditions) are retrieved from WorldClim v2.1 via `geodata` at spatial resolutions of 10, 5, 2.5, or 0.5 arc-minutes (approximately 18 km, 9 km, 4.5 km, and 1 km at the equator). For marine organisms, benthic and surface oceanographic predictors are accessed from Bio-ORACLE and MARSPEC repositories via `sdmpredictors`.

Future climate projections incorporate Coupled Model Intercomparison Project Phase 6 (CMIP6) downscaled datasets across four Shared Socioeconomic Pathways (`ssp126`, `ssp245`, `ssp370`, and `ssp585`) and four time horizons (`2030` [2021-2040], `2050` [2041-2060], `2070` [2061-2080], and `2090` [2081-2100]). Rather than relying on a single climate trajectory, `cinspecies` downloads multiple user-defined Global Climate Models (such as `MIROC6`, `IPSL-CM6A-LR`, `ACCESS-CM2`, `HadGEM3-GC31-LL`, `MRI-ESM2-0`, `CNRM-CM6-1`, `EC-Earth3-Veg`, or `BCC-CSM2-MR`) and computes a multi-GCM cell-by-cell ensemble mean prior to niche projection, minimizing inter-model uncertainty.

### 3. Collinearity Control and Spatial Processing
All raster manipulation, coordinate reference system (CRS) alignment, bilinear resampling, polygon masking, and ellipsoidal geodetic area calculations (`expanse`) are executed through the `terra` engine. Political and administrative vector boundaries at continental, national, and state/provincial scales are sourced from the Natural Earth repository via `rnaturalearth` and `rnaturalearthhires`. Prior to model fitting, environmental values are extracted at occurrence points and evaluated for multicollinearity using a stepwise Variance Inflation Factor (`vifstep`) algorithm from the `usdm` package, discarding redundant covariates that exceed the user-defined threshold.

### 4. Land Use and Land Cover (LULC) Integration
To transition from potential climatic suitability to realized habitat availability, the post-processing module intersects binary climate projections with empirical land cover data. Users can trigger automated global downloads from the European Space Agency (ESA WorldCover 10m database, aggregated to 30 arc-seconds via `geodata`) for fractional tree cover (`esa_trees`), shrubland (`esa_shrubs`), grassland (`esa_grassland`), or cropland (`esa_cropland`). Alternatively, users may supply a local file path to high-resolution regional categorical rasters, such as MapBiomas collections, specifying the exact integer class codes corresponding to native habitat formations.

## Supported Modeling Algorithms (AI and Statistical Engines)

Model calibration, cross-validation, and ensemble forecasting are powered by the `sdm` package coupled with parallel processing (`parallel`). Users can combine any of the following machine learning, artificial intelligence, and statistical algorithms within the `algorithms` parameter:

| Algorithm Code | Full Method Name | Category | Analytical Profile |
| :--- | :--- | :--- | :--- |
| `rf` | Random Forest | Machine Learning | High-performance ensemble of decision trees; robust against non-linear interactions and overfitting. |
| `maxnet` | Maxent (via glmnet) | Machine Learning | Presence-background regularised maximum entropy modeling; standard in ENM studies. |
| `maxent` | Maxent (Java engine) | Machine Learning | Classic Java implementation of Maximum Entropy (requires `dismo` and Java installed). |
| `brt` / `gbm` | Boosted Regression Trees | Machine Learning | Gradient boosting combining sequential trees; highly effective for complex ecological responses. |
| `svm` | Support Vector Machines | Machine Learning | Kernel-based classification that constructs optimal hyperplanes in high-dimensional niche space. |
| `mlp` / `nnet` | Multi-Layer Perceptron | Artificial Neural Network | Feed-forward neural networks capable of modeling highly complex, non-linear relationships. |
| `mars` | Multivariate Adaptive Regression Splines | Machine Learning | Flexible non-parametric regression using piecewise linear splines. |
| `cart` | Classification and Regression Trees | Machine Learning | Single recursive partitioning tree structure. |
| `glm` | Generalized Linear Models | Parametric Statistical | Classic regression framework fitting linear and quadratic environmental responses; high transferability. |
| `gam` | Generalized Additive Models | Semi-Parametric Statistical | Extends GLMs by fitting smooth non-parametric spline functions to environmental predictors. |
| `mahal` | Mahalanobis Distance | Environmental Envelope | Presence-only distance metric accounting for covariance among environmental variables. |
| `bioclim` | BIOCLIM | Environmental Envelope | Classic rectilinear envelope model based on percentile distributions of presence records. |

Individual model replicates are evaluated via cross-validation and synthesized into a final consensus map weighted by the True Skill Statistic (TSS = Sensitivity + Specificity - 1), ensuring that algorithms or replicates with poor predictive performance do not degrade the final ensemble.

## Comprehensive Parameter Reference

The master function `cins_pipeline()` accepts 23 arguments that control every phase of the modeling workflow:

| Parameter | Default | Accepted Values / Type | Description and Usage |
| :--- | :--- | :--- | :--- |
| `species` | Required | Character string | Scientific binomial name of the focal taxon (e.g., `"Ilex paraguariensis"`). |
| `region` | `"South America"` | Character vector | Calibration area (background/M area). Accepts continents, countries, or `"Country:State"` (e.g., `"Brazil:Parana"`). |
| `proj_region` | `NULL` | Character vector | Geographic target for spatial transferability. If `NULL`, projects onto `region`. |
| `dependent_species` | `NULL` | Character vector | Names of host, prey, or canopy species. Recursively models these taxa first and includes their suitability as biotic predictors. |
| `future_scenario` | `NULL` | Character string | CMIP6 emission pathway: `"ssp126"`, `"ssp245"`, `"ssp370"`, or `"ssp585"`. |
| `future_year` | `2050` | Numeric | CMIP6 time horizon: `2030`, `2050`, `2070`, or `2090`. |
| `gcms` | `c("MIROC6", "IPSL-CM6A-LR")` | Character vector | Global Climate Models to download and average for future projections. |
| `ecosystem` | `"terrestrial"` | Character string | `"terrestrial"` (WorldClim + land mask) or `"marine"` (Bio-ORACLE/MARSPEC + ocean mask). |
| `resolution` | `5` | Numeric | Baseline climate spatial resolution in arc-minutes: `10`, `5`, `2.5`, or `0.5`. |
| `fine_resolution` | `NULL` | Numeric | Optional target resolution in decimal degrees for bilinear downscaling of output rasters (e.g., `0.01`). |
| `binary_threshold` | `"maxSSS"` | `"maxSSS"` or Numeric | Cutoff for presence/absence conversion. Accepts `"maxSSS"` or a strict probability value between `0` and `1` (e.g., `0.5`). |
| `lulc_source` | `NULL` | Character string | Land cover filter: `"esa_trees"`, `"esa_shrubs"`, `"esa_grassland"`, `"esa_cropland"`, or path to a local `.tif` file. |
| `lulc_classes` | `NULL` | Numeric vector | Integer class IDs to retain when using a categorical local raster (e.g., `c(3, 4)` for MapBiomas forest/savanna). |
| `area_report` | `NULL` | Character string | Geodetic area stratification level: `"continent"`, `"country"`, or `"state"`. Exports CSV tables. |
| `vif_cutoff` | `10` | Numeric | Maximum Variance Inflation Factor allowed during stepwise collinearity filtering (recommended: `5` to `10`). |
| `thin_km` | `10` | Numeric | Minimum geographic distance (in kilometers) between occurrence records during spatial thinning. |
| `n_pseudo` | `1000` | Numeric | Number of random background/pseudo-absence points generated for model calibration. |
| `alt_min` | `NULL` | Numeric | Minimum elevation threshold (in meters above sea level) used to filter occurrence points via SRTM DEM. |
| `alt_max` | `NULL` | Numeric | Maximum elevation threshold (in meters above sea level) used to filter occurrence points via SRTM DEM. |
| `algorithms` | `c("rf", "maxnet", "glm")` | Character vector | Suite of statistical and machine learning algorithms to train and ensemble. |
| `replicates` | `1` | Numeric | Number of cross-validation subsampling replicates per algorithm (`1` for single split; `>= 5` recommended for publication). |
| `ncores` | `1` | Numeric | Number of CPU processor cores allocated for parallel model training and raster projection. |
| `plot_results` | `TRUE` | Logical | If `TRUE`, automatically renders a multi-panel comparative visualization of continuous, binary, LULC, and future maps. |
| `output_dir` | `getwd()` | Character string | Directory path where `.tif` rasters and `.csv` area reports will be written. |

## Installation

Install the package directly from GitHub using `remotes`:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("nevesms/cinspecies")
```

## Reproducible Research Workflows

### 1. Standard Calibration with Automated ESA Land Cover and CMIP6 Projections

The following workflow calibrates an ensemble model across South America, transfers the projection to Brazil, applies a 50% suitability threshold, filters the potential niche by existing tree cover from the European Space Agency, projects into 2050 under SSP5-8.5 across two GCMs, and computes state-level habitat loss:

```r
library(cinspecies)

model_esa <- cins_pipeline(
  species = "Ilex paraguariensis",
  region = "South America",
  proj_region = "Brazil",
  future_scenario = "ssp585",
  future_year = 2050,
  gcms = c("MIROC6", "IPSL-CM6A-LR"),
  resolution = 5,
  binary_threshold = 0.5,
  lulc_source = "esa_trees",
  area_report = "state",
  vif_cutoff = 5,
  thin_km = 15,
  n_pseudo = 1500,
  alt_min = 305,
  algorithms = c("rf", "maxnet", "glm"),
  replicates = 5,
  ncores = 4,
  plot_results = TRUE
)

# Inspect geodetic area (km2) and anthropogenic habitat loss per state
print(model_esa$area_table_projected)
print(model_esa$area_table_future)
```

### 2. Regional Calibration with Local MapBiomas Raster and Biotic Interactions

For regional studies requiring high-resolution categorical land use data (such as MapBiomas) and biotic dependencies, users can constrain the training region to specific states or countries, supply the local `.tif` path, and declare the associated species:

```r
model_regional <- cins_pipeline(
  species = "Oncidium baueri",
  dependent_species = c("Araucaria angustifolia"),
  region = c("Brazil:Parana", "Brazil:Santa Catarina", "Brazil:Rio Grande do Sul"),
  proj_region = "Brazil",
  binary_threshold = 0.6,
  lulc_source = "C:/data/mapbiomas_brazil_collection.tif",
  lulc_classes = c(3, 4),
  area_report = "state",
  vif_cutoff = 5,
  algorithms = c("rf", "brt", "maxnet", "gam"),
  replicates = 5,
  ncores = 4
)
```

### Returned Object Structure

The `cins_pipeline()` function writes all raster layers (`.tif`) and summary tables (`.csv`) directly to `output_dir` and returns an R `list` containing:
* `occurrences`: Cleaned and thinned data frame of species presence coordinates (`lon`, `lat`).
* `model`: Full `sdmModels` object containing evaluation metrics (AUC, TSS, COR, Deviance) and variable importance scores.
* `cutoff_used`: Exact numeric probability value applied during binary conversion.
* `map_native` & `map_native_binary` & `map_native_lulc`: `SpatRaster` objects of continuous, binary, and LULC-filtered suitability in the calibration region.
* `map_projected` & `map_projected_binary` & `map_projected_lulc`: `SpatRaster` objects for the target transfer region.
* `map_future` & `map_future_binary` & `map_future_lulc`: `SpatRaster` objects for the CMIP6 future scenario.
* `area_table_native`, `area_table_projected`, & `area_table_future`: Data frames reporting `Potential_Climatic_km2`, `Realized_LULC_km2`, and `Habitat_Loss_pct` per administrative unit.

## Methodological Considerations and Known Limitations

Statistical and machine learning algorithms identify environmental correlations based strictly on the input coordinates and background extent provided. Researchers utilizing `cinspecies` for ecological inference or conservation planning must account for three primary methodological considerations during model calibration and interpretation.

The first consideration involves calibration extent and background sampling bias. Selecting an overly broad geographic background for model training (for instance, setting `region = "South America"` for a species restricted to the Atlantic Forest biome) causes random pseudo-absence sampling to draw points from extreme environmental zones such as the Atacama Desert, the Amazon Basin, or high Andean peaks. When algorithms contrast presence records against such extreme backgrounds, evaluation metrics become artificially inflated, and the fitted response curves become overly permissive. Consequently, transitional regions with moderate seasonal stress (such as the central Brazilian Cerrado in Mato Grosso or Minas Gerais) may be misclassified as highly suitable. Researchers should constrain the `region` parameter to countries or state combinations that reflect the historically accessible area (the M area in the BAM framework) of the focal taxon, alongside enforcing stricter multicollinearity control (e.g., `vif_cutoff = 5`).

The second consideration pertains to binary thresholding and overprediction. By default, `cinspecies` employs the `maxSSS` threshold to convert continuous probabilities into binary maps. Although widely validated in the literature, `maxSSS` assigns equal weight to sensitivity and specificity. In datasets with broad environmental gradients or sparse occurrences, this mathematical balance often yields low probability cutoffs (frequently between 0.10 and 0.20). Classifying pixels with 15% suitability as presence inevitably overestimates the potential distribution across marginal habitats. When modeling specialist species or when strict commission error control is required, users should override the automated metric by supplying a direct numeric threshold (such as `binary_threshold = 0.5` or `0.6`) and inspect the continuous suitability maps alongside the binary outputs.

The third consideration involves raster-vector boundary effects in geodetic area reports. When generating stratified area tables (`area_report = "state"` or `"country"`), users may observe minor area values assigned to administrative units immediately adjacent to the target projection region (for example, small area fragments attributed to departments in Uruguay, Paraguay, or Bolivia when projecting onto Brazil). This phenomenon is not a projection error, but a geometric artifact of spatial resolution: environmental rasters consist of rectangular grid cells (e.g., 5 arc-minutes, approximately 9 km at the equator), whereas administrative boundaries from `rnaturalearth` are high-precision vector polygons. Grid cells intersecting international borders retain their climatic values during regional masking and are subsequently intersected by neighboring state polygons during zonal extraction. These boundary fragments are typically orders of magnitude smaller than core state areas and should be filtered out or disregarded during final tabular reporting.

## Authorship and License

Developed and maintained by Mateus dos Santos Neves. Source code and issue tracking are available at https://github.com/nevesms/cinspecies. This software is released under the MIT License.
