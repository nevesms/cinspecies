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
#' @param ensemble_setting List passed to `sdm::ensemble()` (see [fit_models()]).
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
model_threshold <- function(model, binary_threshold = "maxTSS",
                            ensemble_setting = list(method = "weighted", stat = "TSS", opt = 2)) {
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
    sdm::threshold(model, id = "ensemble", opt = opt, setting = ensemble_setting),
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

#' Build the sdm ensemble setting
#'
#' @param ensemble_method sdm ensemble method: `"weighted"` (default),
#'   `"unweighted"`, `"median"`, `"mean-weighted"`, `"median-weighted"`, ...
#' @param ensemble_stat Statistic used for weights and selection (`"TSS"`, `"AUC"`, ...).
#' @param ensemble_min Numeric or `NULL`. Only models with
#'   `ensemble_stat >= ensemble_min` enter the ensemble (e.g. 0.5 for TSS).
#' @return A list for the `setting` argument of `sdm::ensemble()`.
#' @export
ensemble_setting <- function(ensemble_method = "weighted", ensemble_stat = "TSS", ensemble_min = NULL) {
  st <- list(method = ensemble_method, stat = ensemble_stat, opt = 2)
  if (!is.null(ensemble_min)) st$expr <- paste(ensemble_stat, ">=", ensemble_min)
  st
}

#' Fit a multi-algorithm ensemble model
#'
#' Builds presence-background data (one or more independent sets of random
#' background points within the calibration area; Barbet-Massin et al. 2012),
#' fits the selected algorithms with `sdm` using sub-sampling, k-fold
#' cross-validation or bootstrapping, and predicts an ensemble (TSS-weighted by
#' default, optionally restricted to models with TSS >= `ensemble_min`;
#' Araujo & New 2007). Evaluation (AUC, TSS, Kappa, COR, sensitivity,
#' specificity, threshold) per model, variable importance and an uncertainty
#' map (coefficient of variation among models) are produced.
#'
#' @param occurrences Data frame with `decimalLongitude`, `decimalLatitude`.
#' @param predictors `SpatRaster` of predictors for the calibration area.
#' @param algorithms Character. Any code from [list_algorithms()].
#' @param n_pseudo Numeric. Background points per set.
#' @param n_pseudo_sets Numeric. Number of independent background sets.
#' @param bg_method `"gRandom"` (geographic random, default) or `"eRandom"`
#'   (environmentally stratified random).
#' @param cv_method `"sub"` (random sub-sampling), `"cv"` (k-fold) or `"boot"`.
#' @param replicates Numeric. Repetitions of the partitioning per algorithm and set.
#' @param test_percent Numeric. Percentage held out (sub-sampling/bootstrap).
#' @param cv_folds Numeric. Folds for `cv_method = "cv"`.
#' @param ensemble_method,ensemble_stat,ensemble_min See [ensemble_setting()].
#' @param uncertainty Logical. Produce the coefficient-of-variation map.
#' @param ncores Numeric. CPU cores.
#' @param binary_threshold See [model_threshold()].
#' @param seed Integer or `NULL`.
#' @param output_dir Directory for outputs.
#' @param species_name File prefix.
#' @return List with `model`, `ensemble_map`, `binary_map`, `uncertainty_map`,
#'   `threshold`, `ensemble_setting`, `selected_vars`, `evaluation`,
#'   `var_importance`.
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
                       n_pseudo = 1000, n_pseudo_sets = 1, bg_method = "gRandom",
                       cv_method = "sub", replicates = 1, test_percent = 30, cv_folds = 5,
                       ensemble_method = "weighted", ensemble_stat = "TSS", ensemble_min = NULL,
                       uncertainty = TRUE, ncores = 1, binary_threshold = "maxTSS", seed = NULL,
                       output_dir = getwd(), species_name = "Species") {
  .check_algorithms(algorithms)
  if (!cv_method %in% c("sub", "cv", "boot")) stop("`cv_method` must be 'sub', 'cv' or 'boot'.", call. = FALSE)
  if (!bg_method %in% c("gRandom", "eRandom")) stop("`bg_method` must be 'gRandom' or 'eRandom'.", call. = FALSE)
  pts <- terra::vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  pts$occurrence <- 1
  pts <- pts[, "occurrence"]
  par_cfg <- if (ncores > 1) list(ncore = ncores, method = "parallel") else NULL

  m <- NULL
  set_id <- integer(0)
  for (k in seq_len(n_pseudo_sets)) {
    if (!is.null(seed)) set.seed(seed + k - 1)
    .msg("Modeling", "Background set ", k, "/", n_pseudo_sets, ": ", n_pseudo, " points (", bg_method, ")...")
    d <- sdm::sdmData(occurrence ~ ., train = pts, predictors = predictors,
                      bg = list(n = n_pseudo, method = bg_method, remove = TRUE))
    .msg("Modeling", "Fitting ", paste(algorithms, collapse = ", "), " (", cv_method, ", ", replicates,
         " replicate[s]", if (cv_method == "cv") paste0(", ", cv_folds, " folds") else paste0(", ", test_percent, "% test"), ")...")
    args <- list(occurrence ~ ., data = d, methods = algorithms, replication = cv_method,
                 n = replicates, parallelSettings = par_cfg)
    if (cv_method == "cv") args$cv.folds <- cv_folds else args$test.percent <- test_percent
    mk <- do.call(sdm::sdm, args)
    n_new <- nrow(sdm::getModelInfo(mk))
    set_id <- c(set_id, rep(k, n_new))
    m <- if (is.null(m)) mk else m + mk
  }
  info <- sdm::getModelInfo(m)
  info$pa_set <- set_id[seq_len(nrow(info))]
  if (!any(info$success)) stop("All models failed to fit. Check the algorithm packages and the data.", call. = FALSE)

  ev <- tryCatch({
    e <- sdm::getEvaluation(m, stat = c("AUC", "TSS", "Kappa", "COR", "sensitivity", "specificity", "threshold"), opt = 2)
    merge(info[, c("modelID", "method", "replication", "pa_set")], e, by = "modelID")
  }, error = function(e) NULL)
  setting <- ensemble_setting(ensemble_method, ensemble_stat, ensemble_min)
  if (!is.null(ensemble_min) && !is.null(ev)) {
    ok <- sum(ev[[ensemble_stat]] >= ensemble_min, na.rm = TRUE)
    .msg("Modeling", ok, " of ", nrow(ev), " models have ", ensemble_stat, " >= ", ensemble_min, " and enter the ensemble.")
    if (ok == 0) {
      warning("No model reached ", ensemble_stat, " >= ", ensemble_min, "; all models were used.", call. = FALSE)
      setting$expr <- NULL
    }
    ev$in_ensemble <- ev[[ensemble_stat]] >= ensemble_min
  }
  if (!is.null(ev)) utils::write.csv(ev, file.path(output_dir, paste0(species_name, "_evaluation.csv")), row.names = FALSE)

  .msg("Modeling", "Predicting ensemble (", setting$method, ", ", setting$stat,
       if (!is.null(setting$expr)) paste0(", ", setting$expr), ")...")
  ens <- .predict_ensemble(m, predictors, setting)
  names(ens) <- "suitability"
  ens <- .write_rast(ens, file.path(output_dir, paste0(species_name, "_native_suitability.tif")))

  unc <- NULL
  if (uncertainty) unc <- .uncertainty_map(m, predictors, setting,
                                           file.path(output_dir, paste0(species_name, "_native_uncertainty_cv.tif")))

  th <- model_threshold(m, binary_threshold, setting)
  .msg("Threshold", attr(th, "method"), " = ", round(th, 4))
  bin <- .write_rast(post_process_raster(ens, th),
                     file.path(output_dir, paste0(species_name, "_native_binary.tif")), "INT1U")

  vi <- tryCatch(sdm::getVarImp(m)@varImportanceMean$corTest, error = function(e) NULL)
  if (!is.null(vi)) utils::write.csv(vi, file.path(output_dir, paste0(species_name, "_variable_importance.csv")), row.names = FALSE)

  list(model = m, ensemble_map = ens, binary_map = bin, uncertainty_map = unc, threshold = th,
       ensemble_setting = setting, selected_vars = names(predictors), evaluation = ev, var_importance = vi)
}

#' Model IDs entering the ensemble
#' @noRd
.ensemble_ids <- function(model, setting) {
  info <- sdm::getModelInfo(model)
  ids <- info$modelID[info$success]
  if (!is.null(setting$expr)) {
    thr <- as.numeric(sub(".*>=\\s*", "", setting$expr))
    ev <- sdm::getEvaluation(model, stat = setting$stat, opt = 2)
    keep <- ev$modelID[ev[[setting$stat]] >= thr]
    if (length(keep) > 0) ids <- intersect(ids, keep)
  }
  ids
}

#' Ensemble prediction (falls back to a single model when only one is available)
#' @noRd
.predict_ensemble <- function(model, newdata, setting) {
  ids <- .ensemble_ids(model, setting)
  if (length(ids) == 1) {
    .msg("Modeling", "Only one model available; its prediction is used instead of an ensemble.")
    pf <- methods::selectMethod("predict", "sdmModels")
    return(pf(model, newdata = newdata, id = ids)[[1]])
  }
  sdm::ensemble(model, newdata = newdata, setting = setting)
}

#' Coefficient of variation among the models of the ensemble (in %)
#' @noRd
.uncertainty_map <- function(model, newdata, setting, filename) {
  if (length(.ensemble_ids(model, setting)) < 2) return(NULL)
  st <- setting
  st$method <- "cv"
  u <- tryCatch(sdm::ensemble(model, newdata = newdata, setting = st), error = function(e) NULL)
  if (is.null(u)) return(NULL)
  u <- u * 100
  names(u) <- "cv_percent"
  u <- .write_rast(u, filename)
  mv <- terra::global(u, "mean", na.rm = TRUE)[1, 1]
  .msg("Uncertainty", "Mean coefficient of variation among models: ", round(mv, 2))
  u
}
