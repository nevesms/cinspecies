# cinspecies: niche and realized-habitat modelling of forest and grassland plants

`cinspecies` (Climate Intelligence for Niche Simulation of Species) is an R package that automates ecological niche / species distribution modelling for **terrestrial plants**: trees, shrubs, natural-grassland (campo) species, lianas and epiphytes. It couples multi-algorithm ensemble models with **land use and land cover (LULC)** so that the output is not only the potential climatic niche but the **realized habitat at 30 m**, in the present and under CMIP6 futures, with geodetic area reports of habitat loss.

A single call to `cins_model()` runs every step. All defaults are literature-based and every step is referenced (`cins_references()`).

## Installation

```r
install.packages("remotes")
remotes::install_github("nevesms/cinspecies")

# algorithm packages used by sdm (install those you need)
install.packages(c("randomForest", "maxnet", "glmnet", "mgcv", "gbm", "e1071", "earth", "rpart"))
```

## Quick start

```r
library(cinspecies)

m <- cins_model(
  species = "Araucaria angustifolia",
  life_form = "tree",            # tree | shrub | grassland | liana | epiphyte
  savanna_habitat = FALSE,       # must be set: does Savanna Formation count as habitat?
  region = "auto",               # calibration area (M) from the records
  future_scenario = c("ssp245", "ssp585"),
  future_year = c(2050, 2090),
  area_report = "state",
  replicates = 5,
  seed = 42
)
m                      # summary
m$native$area_table    # potential vs realized km2 and habitat loss per state
cins_references()      # what to cite
```

Epiphytes and lianas can use their hosts (phorophytes) as biotic predictors:

```r
m_epi <- cins_model(
  species = "Vriesea incurvata", life_form = "epiphyte", savanna_habitat = FALSE,
  host_species = c("Araucaria angustifolia"), host_life_form = "tree",
  region = c("Brazil:Parana", "Brazil:Santa Catarina", "Brazil:Rio Grande do Sul"),
  algorithms = c("rf", "maxnet", "gam"), replicates = 5, seed = 1
)
```

Records from a previous study can be supplied for full reproducibility (GBIF content changes over time):

```r
m <- cins_model("Ocotea porosa", life_form = "tree", savanna_habitat = FALSE,
                occurrences = "my_records.csv", region = "Brazil")
```

## Workflow

| Step | What is done | Main tools |
| :--- | :--- | :--- |
| 1. Records | GBIF query by accepted name (synonyms resolved, kingdom Plantae checked); removal of fossils, living collections, cultivated/managed records, records before `year_min` (1970) and with coordinate uncertainty > 10 km; CoordinateCleaner tests (capitals, centroids, equal, GBIF HQ, institutions, seas, zeros); optional elevation filter (30″ DEM); spatial thinning; removal of records inside forest plantations (MapBiomas class 9). | `rgbif`, `CoordinateCleaner`, `spThin`, `geodata` |
| 2. Calibration area (M) | `region = "auto"`: convex hull of the records + `m_buffer_km` (300 km), clipped to land. Or continents, countries, `"Country:State"`, polygons. | `terra`, GADM |
| 3. Predictors | WorldClim 2.1 BIO1–BIO19; optional SoilGrids 2.0 (`soil_vars`) and elevation/slope (`topography = TRUE`); stepwise VIF. | `geodata`, `usdm` |
| 4. Hosts | Each host is modelled first; host suitability × host habitat fraction enters as a predictor. | — |
| 5. Ensemble | Random background, sub-sampling cross-validation, TSS-weighted ensemble; evaluation (AUC, TSS, Kappa, COR) and variable importance exported. | `sdm` |
| 6. Threshold | **maxTSS** (max sensitivity + specificity) computed on the ensemble, by default; any numeric value given by the user overrides it. | `sdm` |
| 7. LULC | Suitability downscaled to **30 m** and combined with native-habitat classes chosen from `life_form` (MapBiomas in Brazil, ESA WorldCover elsewhere). | `terra` |
| 8. Transfer | Other regions (`proj_region`) and CMIP6 futures (every SSP × year, mean of several GCMs). The calibration threshold is reused so that maps are comparable. | `geodata` |
| 9. Areas | Geodetic area of potential and realized habitat and % habitat loss per continent, country or state; each cell is counted in a single unit. | `terra`, GADM |

### Habitat classes by life form

| `life_form` | MapBiomas Brazil (Collection 10) | ESA WorldCover |
| :--- | :--- | :--- |
| `tree`, `liana`, `epiphyte` | 3 Forest Formation, 5 Mangrove, 6 Floodable Forest, 49 Wooded Sandbank Vegetation | trees + mangroves |
| `shrub` | 12 Grassland, 29 Rocky Outcrop, 49 Wooded Sandbank Vegetation, 50 Herbaceous Sandbank Vegetation | shrubs |
| `grassland` | 11 Wetland, 12 Grassland, 29 Rocky Outcrop, 50 Herbaceous Sandbank Vegetation | grassland + wetland |
| + `savanna_habitat = TRUE` | adds 4 Savanna Formation | adds trees to shrubs |

`savanna_habitat` has **no default on purpose**: whether the Savanna Formation is habitat is an ecological decision (it usually is for shrubs and for Cerrado species of any life form, and usually is not for forest species), so the function stops with an error if it is not set while MapBiomas classes are derived from `life_form`. Classes can always be given explicitly with `lulc_classes`.

### Resolution of realized-habitat maps

Realized maps are produced at **30 m by default** (`lulc_resolution = 30`, the MapBiomas native resolution). The continuous suitability is bilinearly resampled to that grid and multiplied by habitat. Larger values (e.g. `lulc_resolution = 250`) aggregate habitat into the fraction of habitat per cell. That runs faster, but it blurs habitat edges and fragmentation and **reduces the spatial accuracy of realized maps and areas**, and the function warns about this. ESA WorldCover (as distributed by `geodata`) is limited to ~1 km. At 30 m, large regions produce very large rasters: the function reports the grid size before processing, and all outputs are written as compressed GeoTIFFs.

## Main arguments of `cins_model()`

| Argument | Default | Description |
| :--- | :--- | :--- |
| `species` | required | Scientific name. |
| `life_form` | required | `"tree"`, `"shrub"`, `"grassland"`, `"liana"`, `"epiphyte"`. |
| `savanna_habitat` | required with MapBiomas | Savanna Formation (class 4) as habitat (`TRUE`/`FALSE`). |
| `occurrences` | `NULL` | User records (data frame / CSV); otherwise GBIF. |
| `region` | `"auto"` | Calibration area: `"auto"`, continent, country, `"Country:State"`, polygon. |
| `m_buffer_km` | `300` | Buffer for the automatic M. |
| `proj_region` | `NULL` | Region for spatial transfer. |
| `host_species` / `host_life_form` | `NULL` / `"tree"` | Hosts (phorophytes) as biotic predictors. |
| `future_scenario` / `future_year` | `NULL` / `2050` | One or more SSPs and years (2021–2100). |
| `gcms` | `"default"` | ACCESS-CM2, IPSL-CM6A-LR, MIROC6, MPI-ESM1-2-HR, MRI-ESM2-0; or `"all"` or a vector. |
| `resolution` | `5` | Climate resolution (arc-min): 10, 5, 2.5, 0.5. |
| `soil_vars`, `soil_depth`, `topography` | `NULL`, `5`, `FALSE` | Optional static predictors. |
| `lulc_source` | `"auto"` | `"auto"`, `"mapbiomas"`, `"esa"`, raster path/`SpatRaster`, or `NULL`. |
| `lulc_year` | `2024` | MapBiomas year. |
| `lulc_classes` | `NULL` | Override the classes derived from `life_form`. |
| `lulc_resolution` | `30` | Cell size (m) of realized maps; coarser = less accurate. |
| `habitat_fraction_min` | `0.5` | Minimum habitat fraction for coarse cells to be mapped as realized. |
| `exclude_cultivated` / `native_only` | `TRUE` / `FALSE` | Record filters. |
| `year_min` / `max_uncertainty_m` | `1970` / `10000` | Record filters. |
| `thin_km` | `10` | Thinning distance. |
| `alt_min` / `alt_max` | `NULL` | Elevation limits (m). |
| `binary_threshold` | `"maxTSS"` | Or `"sp=se"`, `"maxKappa"`, `"P10"`, … or a number in (0, 1). |
| `vif_cutoff` | `10` | VIF threshold (`NULL` to skip). |
| `n_pseudo` | `1000` | Background points. |
| `algorithms` | `c("rf", "maxnet", "glm")` | Any of `rf`, `maxnet`, `maxent`, `brt`, `svm`, `mlp`, `mars`, `cart`, `glm`, `gam`, `mahal`, `bioclim`. |
| `replicates` / `test_percent` | `1` / `30` | Sub-sampling replicates (≥ 5 recommended) and test share. |
| `seed` | `NULL` | Set for reproducible runs. |
| `area_report` | `"none"` | `"total"`, `"continent"`, `"country"`, `"state"`. |
| `output_dir` / `cache_dir` | `getwd()` / `cins_cache_dir()` | Outputs and persistent download cache. |

## Outputs

Files in `output_dir` (prefix = species name):
`*_occurrences.csv`, `*_evaluation.csv`, `*_variable_importance.csv`, `*_native_suitability.tif`, `*_native_binary.tif`, `*_native_habitat.tif`, `*_native_lulc_stack.tif` (realized suitability, potential binary, realized binary and realized weight at 30 m), equivalent `*_projected_*` and `*_future_<ssp>_<year>_*` files, `*_area_<level>_<period>.csv`, and `*_settings.txt` (all settings and package versions).

The returned `cins_model` object holds `occurrences`, `region`, `model`, `threshold`, `selected_vars`, `evaluation`, `var_importance`, `native`, `projected`, `future` (named `ssp_year`) and `settings`. Each of `native`, `projected` and `future[[i]]` contains `continuous`, `binary`, `realized_continuous`, `realized_binary`, `habitat` and `area_table`.

## Methodological notes

* **Calibration area.** Overly broad backgrounds inflate evaluation metrics and make response curves too permissive (Barve et al. 2011). Prefer `region = "auto"` or a region reflecting the accessible area of the species.
* **Threshold.** maxTSS is the default because it performs well with presence-background data (Liu et al. 2013, 2016). For specialists, or when commission errors must be minimised, give a numeric threshold and inspect the continuous maps.
* **Static land cover in futures.** Future realized maps combine future climate with the *current* land cover (`lulc_year`); they show climate-driven change within today's remaining habitat.
* **Downscaling.** Resampling suitability to 30 m refines where habitat is delineated; it does not add climatic information beyond that of the climate grid.
* **Boundaries.** Areas are computed on ~1 km blocks of the 30 m grid and each block is assigned to a single administrative unit (GADM 4.1) by its centre, which avoids double counting along borders.

## Migrating from 0.1.0

`cins_pipeline()` still works but is deprecated: use `cins_model()`. `dependent_species` became `host_species`, `"esa_trees"`-type sources became `lulc_source = "esa"` (classes now come from `life_form`), `fine_resolution` was replaced by `lulc_resolution` (metres), and the marine module was removed. See `NEWS.md`.

## Citation

```r
citation("cinspecies")
cins_references()                       # all references, by step
cins_references(algorithms = c("rf", "maxnet", "glm"))
```

## References

**R and R packages**
* R Core Team. *R: A language and environment for statistical computing*. R Foundation for Statistical Computing, Vienna. https://www.R-project.org/
* Naimi, B. & Araújo, M.B. (2016). sdm: a reproducible and extensible R platform for species distribution modelling. *Ecography*, 39, 368–375. https://doi.org/10.1111/ecog.01881
* Hijmans, R.J. et al. terra: Spatial Data Analysis. R package. https://doi.org/10.32614/CRAN.package.terra
* Hijmans, R.J. et al. geodata: Download Geographic Data. R package. https://doi.org/10.32614/CRAN.package.geodata
* Chamberlain, S., Oldoni, D. & Waller, J. rgbif: Interface to the Global Biodiversity Information Facility API. R package. https://doi.org/10.32614/CRAN.package.rgbif
* Zizka, A. et al. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. *Methods in Ecology and Evolution*, 10, 744–751. https://doi.org/10.1111/2041-210X.13152
* Aiello-Lammens, M.E. et al. (2015). spThin: an R package for spatial thinning of species occurrence records for use in ecological niche models. *Ecography*, 38, 541–545. https://doi.org/10.1111/ecog.01132
* Naimi, B. et al. (2014). Where is positional uncertainty a problem for species distribution modelling? *Ecography*, 37, 191–203 (usdm). https://doi.org/10.1111/j.1600-0587.2013.00205.x

**Data**
* GBIF.org. GBIF Occurrence Data. https://www.gbif.org
* Fick, S.E. & Hijmans, R.J. (2017). WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37, 4302–4315. https://doi.org/10.1002/joc.5086
* Eyring, V. et al. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. *Geoscientific Model Development*, 9, 1937–1958. https://doi.org/10.5194/gmd-9-1937-2016
* Poggio, L. et al. (2021). SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty. *SOIL*, 7, 217–240. https://doi.org/10.5194/soil-7-217-2021
* Farr, T.G. et al. (2007). The Shuttle Radar Topography Mission. *Reviews of Geophysics*, 45, RG2004. https://doi.org/10.1029/2005RG000183
* Souza, C.M. Jr. et al. (2020). Reconstructing three decades of land use and land cover changes in Brazilian biomes with Landsat archive and Earth Engine. *Remote Sensing*, 12, 2735 (MapBiomas). https://doi.org/10.3390/rs12172735
* Zanaga, D. et al. (2022). ESA WorldCover 10 m 2021 v200. Zenodo. https://doi.org/10.5281/zenodo.7254221
* GADM (2022). Database of Global Administrative Areas, version 4.1. https://gadm.org

**Methods**
* Soberón, J. & Peterson, A.T. (2005). Interpretation of models of fundamental ecological niches and species' distributional areas. *Biodiversity Informatics*, 2, 1–10. https://doi.org/10.17161/bi.v2i0.4
* Barve, N. et al. (2011). The crucial role of the accessible area in ecological niche modeling and species distribution modeling. *Ecological Modelling*, 222, 1810–1819. https://doi.org/10.1016/j.ecolmodel.2011.02.011
* Dormann, C.F. et al. (2013). Collinearity: a review of methods to deal with it and a simulation study evaluating their performance. *Ecography*, 36, 27–46. https://doi.org/10.1111/j.1600-0587.2012.07348.x
* Barbet-Massin, M. et al. (2012). Selecting pseudo-absences for species distribution models: how, where and how many? *Methods in Ecology and Evolution*, 3, 327–338. https://doi.org/10.1111/j.2041-210X.2011.00172.x
* Araújo, M.B. & New, M. (2007). Ensemble forecasting of species distributions. *Trends in Ecology & Evolution*, 22, 42–47. https://doi.org/10.1016/j.tree.2006.09.010
* Allouche, O., Tsoar, A. & Kadmon, R. (2006). Assessing the accuracy of species distribution models: prevalence, kappa and the true skill statistic (TSS). *Journal of Applied Ecology*, 43, 1223–1232. https://doi.org/10.1111/j.1365-2664.2006.01214.x
* Liu, C., White, M. & Newell, G. (2013). Selecting thresholds for the prediction of species occurrence with presence-only data. *Journal of Biogeography*, 40, 778–789. https://doi.org/10.1111/jbi.12058
* Liu, C., Newell, G. & White, M. (2016). On the selection of thresholds for predicting species occurrence with presence-only data. *Ecology and Evolution*, 6, 337–348. https://doi.org/10.1002/ece3.1878

**Algorithms**
* Breiman, L. (2001). Random Forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324
* Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190, 231–259. https://doi.org/10.1016/j.ecolmodel.2005.03.026
* Phillips, S.J. et al. (2017). Opening the black box: an open-source release of Maxent. *Ecography*, 40, 887–893. https://doi.org/10.1111/ecog.03049
* Friedman, J., Hastie, T. & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software*, 33(1), 1–22. https://doi.org/10.18637/jss.v033.i01
* Elith, J., Leathwick, J.R. & Hastie, T. (2008). A working guide to boosted regression trees. *Journal of Animal Ecology*, 77, 802–813. https://doi.org/10.1111/j.1365-2656.2008.01390.x
* McCullagh, P. & Nelder, J.A. (1989). *Generalized Linear Models*, 2nd ed. Chapman & Hall.
* Wood, S.N. (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. *Journal of the Royal Statistical Society B*, 73, 3–36. https://doi.org/10.1111/j.1467-9868.2010.00749.x
* Cortes, C. & Vapnik, V. (1995). Support-vector networks. *Machine Learning*, 20, 273–297. https://doi.org/10.1007/BF00994018
* Friedman, J.H. (1991). Multivariate adaptive regression splines. *The Annals of Statistics*, 19, 1–67. https://doi.org/10.1214/aos/1176347963
* Farber, O. & Kadmon, R. (2003). Assessment of alternative approaches for bioclimatic modeling with special emphasis on the Mahalanobis distance. *Ecological Modelling*, 160, 115–130. https://doi.org/10.1016/S0304-3800(02)00327-7

## License

GPL (>= 3). Developed and maintained by Mateus dos Santos Neves (UFPR). Issues: https://github.com/nevesms/cinspecies/issues
