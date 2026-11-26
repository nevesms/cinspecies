# Model fitting, ensemble and thresholds --------------------------------------

#' @noRd
.threshold_codes <- c("sp=se" = 1, "maxtss" = 2, "maxsss" = 2, "max(se+sp)" = 2, "mincost" = 3,
                      "minrocdist" = 4, "maxkappa" = 5, "maxppv+npv" = 6, "ppv=npv" = 7,
                      "maxnmi" = 8, "maxccr" = 9, "prevalence" = 10, "maxmcc" = 11,
                      "p10" = 12, "p5" = 13, "p1" = 14, "p0" = 15)

#' Compute the binary threshold of the ensemble
#'
#' By default the threshold that maximises the True Skill Statistic
#' (`"maxTSS"`, equivalent to max(sensitivity + specificity), also called
#' maxSSS) is computed on the **ensemble** predictions, as recommended for
#' presence-background data (Liu et al. 2013, 2016). A fixed numeric value in
#' (0, 1) supplied by the user always overrides the automatic criterion.
#'
#' @param model `sdmModels` object.
#' @param binary_threshold `"maxTSS"` (default), another criterion
#'   (`"sp=se"`, `"maxKappa"`, `"minROCdist"`, `"P10"`, `"P5"`, `"maxMCC"`,
#'   `"prevalence"`), or a number between 0 and 1.
#' @return Numeric threshold, with attribute `"method"`.
#' @references
#' Liu, C., White, M. & Newell, G. (2013). Selecting thresholds for the
#' prediction of species occurrence with presence-only data. Journal of
#' Biogeography, 40, 778-789. \doi{10.1111/jbi.12058}
#'
#' Liu, C., Newell, G. & White, M. (2016). On the selection of thresholds for
#' predicting species occurrence with presence-only data. Ecology and
#' Evolution, 6, 337-348. \doi{10.1002/ece3.1878}
#'
#' Allouche, O., Tsoar, A. & Kadmon, R. (2006). Assessing the accuracy of
#' species distribution models: prevalence, kappa and the true skill statistic
#' (TSS). Journal of Applied Ecology, 43, 1223-1232.
#' \doi{10.1111/j.1365-2664.2006.01214.x}
#' @export
model_threshold <- function(model, binary_threshold = "maxTSS") {
  if (is.numeric(binary_threshold)) {
    if (binary_threshold <= 0 || binary_threshold >= 1) stop("Numeric `binary_threshold` must be in (0, 1).", call. = FALSE)
    return(structure(binary_threshold, method = "user-defined"))
  }
  key <- gsub("[ _]", "", tolower(binary_threshold))
  if (!key %in% names(.threshold_codes)) {
    stop("Unknown `binary_threshold`: ", binary_threshold, ". Use 'maxTSS', 'sp=se', 'maxKappa', ",
         "'minROCdist', 'P10', 'P5', 'maxMCC', 'prevalence' or a number in (0, 1).", call. = FALSE)
  }
  opt <- unname(.threshold_codes[key])
  th <- tryCatch(
    sdm::threshold(model, id = "ensemble", opt = opt,
                   setting = list(method = "weighted", stat = "TSS", opt = 2)),
    error = function(e) NULL)
  if (is.null(th) || !is.finite(th)) {
    ev <- sdm::getEvaluation(model, stat = "threshold", opt = opt)
    th <- mean(ev$threshold, na.rm = TRUE)
    .msg("Threshold", "Ensemble threshold unavailable; using the mean of individual-model thresholds.")
  }
  structure(as.numeric(th), method = binary_threshold)
}

#' Convert a continuous suitability map to binary
#'
#' @param raster_map `SpatRaster` of continuous suitability.
#' @param threshold Numeric threshold (see [model_threshold()]).
#' @return `SpatRaster` with 1 (suitable) and 0.
#' @export
post_process_raster <- function(raster_map, threshold) {
  if (!is.numeric(threshold)) stop("`threshold` must be numeric; see model_threshold().", call. = FALSE)
  b <- terra::ifel(raster_map >= threshold, 1L, 0L)
  names(b) <- "binary_suitability"
  b
}

#' Fit a multi-algorithm ensemble model
#'
#' Builds presence-background data (random background within the calibration
#' area; Barbet-Massin et al. 2012), fits the selected algorithms with `sdm`
#' using random sub-sampling cross-validation, and predicts a TSS-weighted
#' ensemble (Araujo & New 2007). Model evaluation (AUC, TSS, Kappa, COR,
#' threshold) and variable importance are returned as tables.
#'
#' @param occurrences Data frame with `decimalLongitude`, `decimalLatitude`.
#' @param predictors `SpatRaster` of predictors for the calibration area.
#' @param algorithms Character. sdm method codes (see README).
#' @param n_pseudo Numeric. Number of background points.
#' @param replicates Numeric. Sub-sampling replicates per algorithm (>= 5 recommended).
#' @param test_percent Numeric. Percentage of records held out for testing.
#' @param ncores Numeric. CPU cores.
#' @param binary_threshold See [model_threshold()].
#' @param seed Integer or `NULL`.
#' @param output_dir Directory for outputs.
#' @param species_name File prefix.
#' @return List with `model`, `ensemble_map`, `binary_map`, `threshold`,
#'   `selected_vars`, `evaluation`, `var_importance`.
#' @references
#' Naimi, B. & Araujo, M.B. (2016). sdm: a reproducible and extensible R
#' platform for species distribution modelling. Ecography, 39, 368-375.
#' \doi{10.1111/ecog.01881}
#'
#' Araujo, M.B. & New, M. (2007). Ensemble forecasting of species
#' distributions. Trends in Ecology & Evolution, 22, 42-47.
#' \doi{10.1016/j.tree.2006.09.010}
#'
#' Barbet-Massin, M. et al. (2012). Selecting pseudo-absences for species
#' distribution models: how, where and how many? Methods in Ecology and
#' Evolution, 3, 327-338. \doi{10.1111/j.2041-210X.2011.00172.x}
#' @export
fit_models <- function(occurrences, predictors, algorithms = c("rf", "maxnet", "glm"),
                       n_pseudo = 1000, replicates = 1, test_percent = 30, ncores = 1,
                       binary_threshold = "maxTSS", seed = NULL, output_dir = getwd(),
                       species_name = "Species") {
  .check_algorithms(algorithms)
  if (!is.null(seed)) set.seed(seed)
  pts <- terra::vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  pts$occurrence <- 1
  pts <- pts[, "occurrence"]

  .msg("Modeling", "Building sdmData with ", n_pseudo, " background points...")
  d <- sdm::sdmData(occurrence ~ ., train = pts, predictors = predictors,
                    bg = list(n = n_pseudo, method = "gRandom", remove = TRUE))

  par_cfg <- if (ncores > 1) list(ncore = ncores, method = "parallel") else NULL
  .msg("Modeling", "Fitting ", paste(algorithms, collapse = ", "), " (", replicates,
       " replicate[s], ", test_percent, "% test, ", ncores, " core[s])...")
  m <- sdm::sdm(occurrence ~ ., data = d, methods = algorithms, replication = "sub",
                test.percent = test_percent, n = replicates, parallelSettings = par_cfg)
  ok <- sdm::getModelInfo(m)
  if (!any(ok$success)) stop("All models failed to fit. Check the algorithm packages and the data.", call. = FALSE)

  .msg("Modeling", "Predicting TSS-weighted ensemble...")
  ens <- sdm::ensemble(m, newdata = predictors, setting = list(method = "weighted", stat = "TSS", opt = 2))
  names(ens) <- "suitability"
  ens <- .write_rast(ens, file.path(output_dir, paste0(species_name, "_native_suitability.tif")))

  th <- model_threshold(m, binary_threshold)
  .msg("Threshold", attr(th, "method"), " = ", round(th, 4))
  bin <- .write_rast(post_process_raster(ens, th),
                     file.path(output_dir, paste0(species_name, "_native_binary.tif")), "INT1U")

  ev <- tryCatch({
    e <- sdm::getEvaluation(m, stat = c("AUC", "TSS", "Kappa", "COR", "sensitivity", "specificity", "threshold"), opt = 2)
    info <- ok[, c("modelID", "method", "replication")]
    merge(info, e, by = "modelID")
  }, error = function(e) NULL)
  if (!is.null(ev)) utils::write.csv(ev, file.path(output_dir, paste0(species_name, "_evaluation.csv")), row.names = FALSE)

  vi <- tryCatch({
    v <- sdm::getVarImp(m)
    v@varImportanceMean$corTest
  }, error = function(e) NULL)
  if (!is.null(vi)) utils::write.csv(vi, file.path(output_dir, paste0(species_name, "_variable_importance.csv")), row.names = FALSE)

  list(model = m, ensemble_map = ens, binary_map = bin, threshold = th,
       selected_vars = names(predictors), evaluation = ev, var_importance = vi)
}
