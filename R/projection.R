# Spatial and temporal transfer -----------------------------------------------

#' Project a calibrated ensemble to another region and/or a future scenario
#'
#' Current climate (WorldClim 2.1) or a multi-GCM CMIP6 mean is loaded for the
#' target region, static predictors (soil, topography) and host layers are
#' added, and the ensemble is predicted. The binary map uses the threshold of
#' the calibrated model, so that native, projected and future maps are
#' directly comparable.
#'
#' @param cins_result Output of [fit_models()] (or [cins_model()]).
#' @param target_region Region definition or `SpatVector`.
#' @param resolution Climate resolution (arc-min).
#' @param future_scenario `NULL` (current climate) or an SSP code.
#' @param future_year Numeric.
#' @param gcms GCMs (see [get_future_climate()]).
#' @param static_layers `SpatRaster` of static predictors for the target region, or `NULL`.
#' @param biotic_layers `SpatRaster` of host layers for the target region/scenario, or `NULL`.
#' @param new_climate Optional user climate for the target (a `SpatRaster`, a
#'   file path, or a list of them - e.g. several GCMs from any source - which
#'   are averaged cell by cell). Layer names must match the calibration
#'   variables (`bio1` ... `bio19`; 19-layer stacks are renamed automatically).
#'   When given, nothing is downloaded.
#' @param uncertainty Logical. Also write the coefficient-of-variation map.
#' @param output_dir,species_name,tag Output location and file naming.
#' @param cache_dir Directory where data are cached.
#' @return List with `continuous`, `binary`, `threshold`.
#' @export
project_niche <- function(cins_result, target_region, resolution = 5, future_scenario = NULL,
                          future_year = 2050, gcms = "default", static_layers = NULL,
                          biotic_layers = NULL, new_climate = NULL, uncertainty = FALSE,
                          output_dir = getwd(), species_name = "Species",
                          tag = "projected", cache_dir = cins_cache_dir()) {
  poly <- if (inherits(target_region, "SpatVector")) target_region else get_region_polygon(target_region, cache_dir)
  vars <- cins_result$selected_vars
  clim_vars <- grep("^bio[0-9]+$", vars, value = TRUE)

  if (!is.null(new_climate)) {
    .msg("Projection", "Using user-supplied climate (", tag, ")...")
    new <- .load_user_climate(new_climate, poly, clim_vars)
  } else if (is.null(future_scenario)) {
    .msg("Projection", "Current climate for the target region...")
    new <- get_climate(region = poly, resolution = resolution, vif_cutoff = NULL, cache_dir = cache_dir)[[clim_vars]]
  } else {
    new <- get_future_climate(region = poly, resolution = resolution, future_scenario = future_scenario,
                              future_year = future_year, gcms = gcms, selected_vars = clim_vars,
                              cache_dir = cache_dir)
  }
  if (!is.null(static_layers)) {
    st <- static_layers[[intersect(names(static_layers), vars)]]
    if (terra::nlyr(st) > 0) new <- c(new, terra::resample(st, new, method = "bilinear"))
  }
  if (!is.null(biotic_layers)) {
    bl <- biotic_layers[[intersect(names(biotic_layers), vars)]]
    if (terra::nlyr(bl) > 0) new <- c(new, terra::resample(bl, new, method = "bilinear"))
  }
  missing <- setdiff(vars, names(new))
  if (length(missing) > 0) stop("Predictors missing for projection: ", paste(missing, collapse = ", "), call. = FALSE)
  new <- new[[vars]]

  .msg("Projection", "Predicting ensemble (", tag, ")...")
  setting <- if (!is.null(cins_result$ensemble_setting)) cins_result$ensemble_setting else ensemble_setting()
  p <- .predict_ensemble(cins_result$model, new, setting)
  names(p) <- "suitability"
  p <- .write_rast(p, file.path(output_dir, paste0(species_name, "_", tag, "_suitability.tif")))
  b <- .write_rast(post_process_raster(p, cins_result$threshold),
                   file.path(output_dir, paste0(species_name, "_", tag, "_binary.tif")), "INT1U")
  u <- if (uncertainty) .uncertainty_map(cins_result$model, new, setting,
                                         file.path(output_dir, paste0(species_name, "_", tag, "_uncertainty_cv.tif"))) else NULL
  list(continuous = p, binary = b, uncertainty = u, threshold = cins_result$threshold)
}

#' Load user climate (one raster or several, averaged)
#' @noRd
.load_user_climate <- function(x, poly, vars) {
  if (!is.list(x) || inherits(x, "SpatRaster")) x <- list(x)
  lays <- lapply(x, function(r) {
    r <- if (inherits(r, "SpatRaster")) r else terra::rast(r)
    if (!terra::same.crs(r, "EPSG:4326")) r <- terra::project(r, "EPSG:4326")
    if (terra::nlyr(r) == 19 && !all(vars %in% names(r))) r <- .std_bio_names(r)
    miss <- setdiff(vars, names(r))
    if (length(miss) > 0) stop("User climate lacks layers: ", paste(miss, collapse = ", "), call. = FALSE)
    .crop_region(r[[vars]], poly)
  })
  if (length(lays) == 1) return(lays[[1]])
  out <- terra::rast(lapply(vars, function(v) terra::app(terra::rast(lapply(lays, `[[`, v)), mean, na.rm = TRUE)))
  names(out) <- vars
  out
}
