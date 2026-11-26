#' cinspecies: niche and realized-habitat modelling of terrestrial plants
#'
#' Automated workflow for ecological niche / species distribution modelling of
#' forest and natural-grassland plants (trees, shrubs, grassland species,
#' lianas and epiphytes) that couples ensemble models with land use and land
#' cover to map realized habitat at 30 m. The main function is
#' [cins_model()]. Each step relies on published methods and R packages; see
#' [cins_references()] for what to cite.
#'
#' @keywords internal
#' @importFrom sdm sdm sdmData ensemble
"_PACKAGE"

#' @noRd
.cins_refs <- data.frame(
  step = c(
    "Software", "Software", "Software", "Software", "Software", "Software", "Software", "Software",
    "Records", "Records", "Records",
    "Calibration area", "Calibration area",
    "Predictors", "Predictors", "Predictors", "Predictors", "Predictors",
    "Boundaries",
    "Modelling", "Modelling", "Modelling",
    "Algorithms", "Algorithms", "Algorithms", "Algorithms", "Algorithms", "Algorithms",
    "Algorithms", "Algorithms", "Algorithms", "Algorithms",
    "Algorithms", "Algorithms", "Algorithms", "Algorithms", "Algorithms", "Algorithms", "Algorithms",
    "Evaluation and threshold", "Evaluation and threshold", "Evaluation and threshold",
    "Land use and land cover", "Land use and land cover"
  ),
  reference = c(
    "R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/",
    "Naimi, B. & Araujo, M.B. (2016). sdm: a reproducible and extensible R platform for species distribution modelling. Ecography, 39, 368-375.",
    "Hijmans, R.J. terra: Spatial Data Analysis. R package.",
    "Hijmans, R.J., Barbosa, M., Ghosh, A. & Mandel, A. geodata: Download Geographic Data. R package.",
    "Chamberlain, S., Oldoni, D. & Waller, J. rgbif: Interface to the Global Biodiversity Information Facility API. R package.",
    "Zizka, A. et al. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10, 744-751.",
    "Aiello-Lammens, M.E., Boria, R.A., Radosavljevic, A., Vilela, B. & Anderson, R.P. (2015). spThin: an R package for spatial thinning of species occurrence records for use in ecological niche models. Ecography, 38, 541-545.",
    "Naimi, B., Hamm, N.A.S., Groen, T.A., Skidmore, A.K. & Toxopeus, A.G. (2014). Where is positional uncertainty a problem for species distribution modelling? Ecography, 37, 191-203. (usdm package)",
    "GBIF.org. GBIF Occurrence Data. Global Biodiversity Information Facility. https://www.gbif.org (cite the download DOI or query date).",
    "Zizka, A. et al. (2019). Methods in Ecology and Evolution, 10, 744-751 (coordinate cleaning).",
    "Farr, T.G. et al. (2007). The Shuttle Radar Topography Mission. Reviews of Geophysics, 45, RG2004 (elevation filter / topography).",
    "Soberon, J. & Peterson, A.T. (2005). Interpretation of models of fundamental ecological niches and species' distributional areas. Biodiversity Informatics, 2, 1-10.",
    "Barve, N. et al. (2011). The crucial role of the accessible area in ecological niche modeling and species distribution modeling. Ecological Modelling, 222, 1810-1819.",
    "Fick, S.E. & Hijmans, R.J. (2017). WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology, 37, 4302-4315.",
    "Eyring, V. et al. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. Geoscientific Model Development, 9, 1937-1958.",
    "Poggio, L. et al. (2021). SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty. SOIL, 7, 217-240.",
    "Dormann, C.F. et al. (2013). Collinearity: a review of methods to deal with it and a simulation study evaluating their performance. Ecography, 36, 27-46.",
    "Naimi, B. et al. (2014). Ecography, 37, 191-203 (stepwise VIF, usdm).",
    "GADM (2022). Database of Global Administrative Areas, version 4.1. https://gadm.org",
    "Araujo, M.B. & New, M. (2007). Ensemble forecasting of species distributions. Trends in Ecology & Evolution, 22, 42-47.",
    "Barbet-Massin, M., Jiguet, F., Albert, C.H. & Thuiller, W. (2012). Selecting pseudo-absences for species distribution models: how, where and how many? Methods in Ecology and Evolution, 3, 327-338.",
    "Naimi, B. & Araujo, M.B. (2016). Ecography, 39, 368-375 (sdm: fitting, replication and ensemble).",
    "rf - Breiman, L. (2001). Random Forests. Machine Learning, 45, 5-32. (randomForest: Liaw, A. & Wiener, M. (2002). R News, 2(3), 18-22.)",
    "maxnet - Phillips, S.J., Anderson, R.P., Dudik, M., Schapire, R.E. & Blair, M.E. (2017). Opening the black box: an open-source release of Maxent. Ecography, 40, 887-893.",
    "maxent - Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006). Maximum entropy modeling of species geographic distributions. Ecological Modelling, 190, 231-259.",
    "brt/gbm - Elith, J., Leathwick, J.R. & Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology, 77, 802-813.",
    "glm/glmpoly - McCullagh, P. & Nelder, J.A. (1989). Generalized Linear Models, 2nd ed. Chapman & Hall, London.",
    "maxnet/glmnet - Friedman, J., Hastie, T. & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33(1), 1-22.",
    "gam - Wood, S.N. (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society B, 73, 3-36.",
    "svm - Cortes, C. & Vapnik, V. (1995). Support-vector networks. Machine Learning, 20, 273-297.",
    "mars - Friedman, J.H. (1991). Multivariate adaptive regression splines. The Annals of Statistics, 19, 1-67.",
    "mahal.dismo - Farber, O. & Kadmon, R. (2003). Assessment of alternative approaches for bioclimatic modeling with special emphasis on the Mahalanobis distance. Ecological Modelling, 160, 115-130.",
    "mlp/rbf - Bergmeir, C. & Benitez, J.M. (2012). Neural networks in R using the Stuttgart Neural Network Simulator: RSNNS. Journal of Statistical Software, 46(7), 1-26.",
    "svm (kernlab engine) - Karatzoglou, A., Smola, A., Hornik, K. & Zeileis, A. (2004). kernlab - an S4 package for kernel methods in R. Journal of Statistical Software, 11(9), 1-20.",
    "fda/mda - Hastie, T., Tibshirani, R. & Buja, A. (1994). Flexible discriminant analysis by optimal scoring. Journal of the American Statistical Association, 89, 1255-1270.",
    "maxlike - Royle, J.A., Chandler, R.B., Yackulic, C. & Nichols, J.D. (2012). Likelihood analysis of species occurrence probability from presence-only data for modelling species distributions. Methods in Ecology and Evolution, 3, 545-554.",
    "cart/rpart - Breiman, L., Friedman, J.H., Olshen, R.A. & Stone, C.J. (1984). Classification and Regression Trees. Wadsworth, Belmont.",
    "bioclim/bioclim.dismo - Busby, J.R. (1991). BIOCLIM - a bioclimate analysis and prediction system. In: Margules, C.R. & Austin, M.P. (eds) Nature Conservation: Cost Effective Biological Surveys and Data Analysis. CSIRO, Melbourne, 64-68.",
    "domain.dismo - Carpenter, G., Gillison, A.N. & Winter, J. (1993). DOMAIN: a flexible modelling procedure for mapping potential distributions of plants and animals. Biodiversity and Conservation, 2, 667-680.",
    "Allouche, O., Tsoar, A. & Kadmon, R. (2006). Assessing the accuracy of species distribution models: prevalence, kappa and the true skill statistic (TSS). Journal of Applied Ecology, 43, 1223-1232.",
    "Liu, C., White, M. & Newell, G. (2013). Selecting thresholds for the prediction of species occurrence with presence-only data. Journal of Biogeography, 40, 778-789.",
    "Liu, C., Newell, G. & White, M. (2016). On the selection of thresholds for predicting species occurrence with presence-only data. Ecology and Evolution, 6, 337-348.",
    "Souza, C.M. Jr. et al. (2020). Reconstructing three decades of land use and land cover changes in Brazilian biomes with Landsat archive and Earth Engine. Remote Sensing, 12, 2735. (MapBiomas Collection 10)",
    "Zanaga, D. et al. (2022). ESA WorldCover 10 m 2021 v200. Zenodo."
  ),
  doi = c(
    NA, "10.1111/ecog.01881", "10.32614/CRAN.package.terra", "10.32614/CRAN.package.geodata",
    "10.32614/CRAN.package.rgbif", "10.1111/2041-210X.13152", "10.1111/ecog.01132",
    "10.1111/j.1600-0587.2013.00205.x",
    NA, "10.1111/2041-210X.13152", "10.1029/2005RG000183",
    "10.17161/bi.v2i0.4", "10.1016/j.ecolmodel.2011.02.011",
    "10.1002/joc.5086", "10.5194/gmd-9-1937-2016", "10.5194/soil-7-217-2021",
    "10.1111/j.1600-0587.2012.07348.x", "10.1111/j.1600-0587.2013.00205.x",
    NA,
    "10.1016/j.tree.2006.09.010", "10.1111/j.2041-210X.2011.00172.x", "10.1111/ecog.01881",
    "10.1023/A:1010933404324", "10.1111/ecog.03049", "10.1016/j.ecolmodel.2005.03.026",
    "10.1111/j.1365-2656.2008.01390.x", NA, "10.18637/jss.v033.i01", "10.1111/j.1467-9868.2010.00749.x",
    "10.1007/BF00994018", "10.1214/aos/1176347963", "10.1016/S0304-3800(02)00327-7",
    "10.18637/jss.v046.i07", "10.18637/jss.v011.i09", "10.1080/01621459.1994.10476866",
    "10.1111/j.2041-210X.2011.00182.x", NA, NA, "10.1007/BF00051966",
    "10.1111/j.1365-2664.2006.01214.x", "10.1111/jbi.12058", "10.1002/ece3.1878",
    "10.3390/rs12172735", "10.5281/zenodo.7254221"
  ),
  stringsAsFactors = FALSE
)

#' References used by cinspecies
#'
#' Returns (and prints) the references for every step of the workflow: R
#' packages, data sources and methods. Because cinspecies chains many tools,
#' studies using it should cite the original works behind the steps they used,
#' plus the package itself (`citation("cinspecies")`).
#'
#' @param step Optional character to filter by step (e.g. `"Algorithms"`).
#' @param algorithms Optional character vector of algorithm codes; only those
#'   algorithms are listed under "Algorithms".
#' @param print Logical. Print a formatted list.
#' @return Invisibly, a data frame with `step`, `reference` and `doi`.
#' @examples
#' cins_references("Land use and land cover")
#' cins_references(algorithms = c("rf", "maxnet", "glm"))
#' @export
cins_references <- function(step = NULL, algorithms = NULL, print = TRUE) {
  r <- .cins_refs
  if (!is.null(step)) r <- r[tolower(r$step) %in% tolower(step), ]
  if (!is.null(algorithms)) {
    alg <- r$step == "Algorithms"
    code <- sub(" .*", "", r$reference)
    keep_alg <- vapply(strsplit(code, "/"), function(x) any(x %in% algorithms), logical(1))
    r <- r[!alg | keep_alg, ]
  }
  if (print) {
    for (s in unique(r$step)) {
      cat("\n== ", s, " ==\n", sep = "")
      x <- r[r$step == s, ]
      for (i in seq_len(nrow(x))) {
        cat("- ", x$reference[i], if (!is.na(x$doi[i])) paste0(" https://doi.org/", x$doi[i]), "\n", sep = "")
      }
    }
  }
  invisible(r)
}
