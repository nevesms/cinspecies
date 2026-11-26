# Main workflow ---------------------------------------------------------------

#' @noRd
.validate_lulc_args <- function(lulc_source, lulc_classes, savanna_habitat, lulc_resolution) {
  if (is.null(lulc_source)) return(invisible(TRUE))
  needs_savanna <- !identical(lulc_source, "esa") && is.null(lulc_classes)
  if (needs_savanna && (is.null(savanna_habitat) || !is.logical(savanna_habitat))) {
    stop("Please set `savanna_habitat = TRUE` or `FALSE`: should MapBiomas Savanna Formation (class 4) ",
         "count as habitat for this species? Usually TRUE for shrubs and Cerrado species, FALSE for ",
         "forest species. (Alternatively, give `lulc_classes` explicitly.)", call. = FALSE)
  }
  if (!is.numeric(lulc_resolution) || lulc_resolution < 10) stop("`lulc_resolution` must be a number of metres >= 10.", call. = FALSE)
  if (lulc_resolution > 30) {
    warning("lulc_resolution = ", lulc_resolution, " m is coarser than the 30 m default. Habitat will be ",
            "aggregated into fractions per cell, which blurs habitat boundaries and fragmentation and lowers ",
            "the spatial accuracy of realized-habitat maps and areas. Use 30 m whenever computing resources allow.",
            call. = FALSE)
  }
  invisible(TRUE)
}

#' Resolve LULC source + habitat classes for a polygon
#' @noRd
.lulc_setup <- function(lulc_source, poly, life_form, lulc_classes, savanna_habitat, lulc_year, cache_dir) {
  if (is.null(lulc_source)) return(NULL)
  src <- lulc_source
  if (identical(src, "auto")) {
    share <- .share_in_brazil(poly, cache_dir)
    src <- if (share >= 0.95) "mapbiomas" else "esa"
    .msg("Land use", "lulc_source = 'auto' -> ", src, " (", round(100 * share), "% of the region in Brazil).")
  }
  if (identical(src, "esa")) {
    vars <- habitat_classes(life_form, "esa", savanna_habitat)
    lulc <- get_lulc("esa", poly, esa_vars = vars, cache_dir = cache_dir)
    return(list(lulc = lulc, classes = NULL, mapbiomas_legend = FALSE))
  }
  classes <- if (!is.null(lulc_classes)) lulc_classes else habitat_classes(life_form, "mapbiomas", savanna_habitat)
  lulc <- get_lulc(src, poly, year = lulc_year, cache_dir = cache_dir)
  .msg("Land use", "Habitat classes: ", paste(classes, collapse = ", "))
  list(lulc = lulc, classes = classes,
       mapbiomas_legend = identical(src, "mapbiomas") || is.null(lulc_classes))
}

#' Realized maps + area report for one map (native / projected / future)
#' @noRd
.realize <- function(res, setup, lulc_resolution, habitat_fraction_min, area_report, output_dir,
                     sp, tag, cache_dir, habitat_cache = NULL) {
  out <- list(continuous = res$continuous, binary = res$binary)
  if (is.null(setup)) {
    out$area_table <- calculate_area_report(res$binary, NULL, area_report, output_dir, sp, tag, cache_dir)
    return(out)
  }
  hab <- if (!is.null(habitat_cache)) habitat_cache else
    habitat_layer(setup$lulc, setup$classes, lulc_resolution,
                  filename = file.path(output_dir, paste0(sp, "_", tag, "_habitat.tif")))
  .lulc_size_message(hab, max(lulc_resolution, setup$lulc$native_res_m))
  rl <- apply_lulc_mask(res$continuous, res$threshold, hab, habitat_fraction_min,
                        file_prefix = file.path(output_dir, paste0(sp, "_", tag)))
  out$realized_continuous <- rl$continuous
  out$realized_binary <- rl$binary
  out$habitat <- hab
  out$area_table <- calculate_area_report(rl$potential_binary, rl$weight, area_report, output_dir, sp, tag, cache_dir)
  out
}

#' Niche and realized-habitat model of a terrestrial plant species
#'
#' `cins_model()` runs the whole cinspecies workflow for trees, shrubs,
#' natural-grassland species, lianas and epiphytes:
#'
#' 1. **Records** - GBIF (or user) records; kingdom Plantae check; removal of
#'    fossil, cultivated/living-collection and managed records; temporal and
#'    positional-uncertainty filters; CoordinateCleaner; optional elevation
#'    filter; spatial thinning ([get_occurrences()]). Records in forest
#'    plantations (MapBiomas class 9) are also removed.
#' 2. **Calibration area (M)** - automatic (convex hull + buffer,
#'    [define_calibration_area()]) or user-defined ([get_region_polygon()]).
#' 3. **Predictors** - WorldClim 2.1 bioclimatic variables, optional SoilGrids
#'    soil properties and topography, stepwise VIF ([select_predictors()]).
#' 4. **Host species** (lianas, epiphytes, hemiparasites ...) - each host is
#'    modelled first and its suitability x habitat on the climate grid enters
#'    as a biotic predictor.
#' 5. **Ensemble model** - multiple algorithms with `sdm`, TSS-weighted
#'    ensemble, binary map with the maxTSS threshold by default ([fit_models()]).
#' 6. **Land use and land cover** - the suitability map is downscaled to
#'    30 m (default) and combined with native habitat from MapBiomas (Brazil)
#'    or ESA WorldCover, chosen from `life_form` ([habitat_classes()],
#'    [apply_lulc_mask()]).
#' 7. **Transfer** to another region and to CMIP6 futures (several SSPs,
#'    years and GCMs) ([project_niche()]).
#' 8. **Geodetic area reports** of potential and realized habitat and habitat
#'    loss per continent, country or state ([calculate_area_report()]).
#'
#' All rasters (GeoTIFF), tables (CSV) and a settings file are written to
#' `output_dir`. Run [cins_references()] for the full list of references to
#' cite for each step.
#'
#' @param species Character. Scientific name.
#' @param life_form Character (required). `"tree"`, `"shrub"`, `"grassland"`,
#'   `"liana"` or `"epiphyte"`. Defines the land-cover classes used as habitat.
#' @param occurrences Optional user records (data frame, `SpatVector` or CSV
#'   path). If `NULL`, GBIF is queried. Supplying the records used in a previous
#'   study makes results reproducible (GBIF content changes over time).
#' @param region Calibration area: `"auto"` (default, from the records),
#'   a continent, country, `"Country:State"` vector, `SpatVector`, `sf` or file path.
#' @param m_buffer_km Numeric. Buffer (km) for `region = "auto"`.
#' @param proj_region Optional region to transfer the model to (same formats as `region`).
#' @param host_species Optional character vector of host / phorophyte species.
#' @param host_life_form Life form of the hosts (default `"tree"`).
#' @param future_scenario `NULL` or one or more of `"ssp126"`, `"ssp245"`,
#'   `"ssp370"`, `"ssp585"`.
#' @param future_year One or more years (2021-2100); every scenario x year
#'   combination is projected.
#' @param gcms GCMs to average: `"default"` (5 models), `"all"` or a vector.
#' @param resolution Climate resolution in arc-minutes: 10, 5, 2.5 or 0.5.
#' @param soil_vars `NULL` or SoilGrids variables (e.g. `c("clay", "sand", "phh2o", "soc")`).
#' @param soil_depth SoilGrids depth interval (see [get_static_predictors()]).
#' @param topography Logical. Add elevation and slope as predictors.
#' @param lulc_source `"auto"` (MapBiomas in Brazil, ESA elsewhere),
#'   `"mapbiomas"`, `"esa"`, a raster path, a `SpatRaster`, or `NULL` (no LULC step).
#' @param lulc_year Integer. MapBiomas year (Collection 10: 1985-2024).
#' @param lulc_classes Optional integer classes overriding those derived from `life_form`.
#' @param savanna_habitat Logical. Should MapBiomas Savanna Formation (4) count
#'   as habitat? **Required** whenever MapBiomas classes are derived from
#'   `life_form` (no default, by design).
#' @param lulc_resolution Numeric. Cell size (m) of realized-habitat maps.
#'   Default 30 m (MapBiomas native). Larger values aggregate habitat into
#'   fractions and reduce spatial accuracy (a warning is issued). ESA
#'   WorldCover via `geodata` is limited to ~1 km.
#' @param habitat_fraction_min Numeric. Minimum habitat fraction of a coarse
#'   cell to be mapped as realized (only when `lulc_resolution` > native).
#' @param exclude_cultivated,native_only,year_min,max_uncertainty_m,gbif_limit,thin_km,alt_min,alt_max
#'   Record filters, see [get_occurrences()].
#' @param binary_threshold `"maxTSS"` (default) or another criterion or a
#'   number in (0, 1) (see [model_threshold()]).
#' @param vif_cutoff Numeric or `NULL`. VIF threshold.
#' @param n_pseudo Numeric. Number of background points.
#' @param algorithms Character. sdm algorithms (see README).
#' @param replicates,test_percent Cross-validation settings.
#' @param ncores Numeric. CPU cores.
#' @param seed Integer or `NULL`. Seed for thinning, VIF sampling, background
#'   and data partitioning (set it for reproducible runs).
#' @param area_report `"none"`, `"total"`, `"continent"`, `"country"` or `"state"`.
#' @param plot_results Logical. Plot the main maps.
#' @param output_dir Directory for outputs.
#' @param cache_dir Directory for downloaded data ([cins_cache_dir()]).
#' @return An object of class `cins_model` (a list) with `occurrences`,
#'   `region`, `model`, `threshold`, `selected_vars`, `evaluation`,
#'   `var_importance`, `native`, `projected`, `future` (named by
#'   `ssp_year`) and `settings`. Each of `native`, `projected` and the elements
#'   of `future` contains `continuous`, `binary`, `realized_continuous`,
#'   `realized_binary`, `habitat` and `area_table`.
#' @references See [cins_references()].
#' @examples
#' \dontrun{
#' m <- cins_model(
#'   species = "Araucaria angustifolia", life_form = "tree",
#'   savanna_habitat = FALSE, region = "auto",
#'   future_scenario = c("ssp245", "ssp585"), future_year = 2050,
#'   area_report = "state", replicates = 5, seed = 42
#' )
#' m
#' }
#' @export
cins_model <- function(species, life_form, occurrences = NULL,
                       region = "auto", m_buffer_km = 300, proj_region = NULL,
                       host_species = NULL, host_life_form = "tree",
                       future_scenario = NULL, future_year = 2050, gcms = "default",
                       resolution = 5, soil_vars = NULL, soil_depth = 5, topography = FALSE,
                       lulc_source = "auto", lulc_year = 2024, lulc_classes = NULL,
                       savanna_habitat = NULL, lulc_resolution = 30, habitat_fraction_min = 0.5,
                       exclude_cultivated = TRUE, native_only = FALSE, year_min = 1970,
                       max_uncertainty_m = 10000, gbif_limit = 20000, thin_km = 10,
                       alt_min = NULL, alt_max = NULL,
                       binary_threshold = "maxTSS", vif_cutoff = 10, n_pseudo = 1000,
                       algorithms = c("rf", "maxnet", "glm"), replicates = 1, test_percent = 30,
                       ncores = 1, seed = NULL, area_report = "none", plot_results = TRUE,
                       output_dir = getwd(), cache_dir = cins_cache_dir()) {

  # 0. Checks (fail fast, before any download) --------------------------------
  if (missing(life_form) || is.null(life_form) || !life_form %in% .life_forms) {
    stop("`life_form` is required: one of ", paste(.life_forms, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(host_species) && !host_life_form %in% .life_forms) stop("Invalid `host_life_form`.", call. = FALSE)
  .check_resolution(resolution)
  .check_algorithms(algorithms)
  .validate_lulc_args(lulc_source, lulc_classes, savanna_habitat, lulc_resolution)
  if (!is.null(future_scenario)) { invisible(lapply(future_year, .future_period)); .resolve_gcms(gcms) }
  if (!area_report %in% c("none", "total", "continent", "country", "state")) stop("Invalid `area_report`.", call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  sp <- .sp_slug(species)
  settings <- as.list(environment())[setdiff(names(formals()), c("occurrences"))]

  # 1. Records ----------------------------------------------------------------
  pts <- get_occurrences(species, occurrences = occurrences, gbif_limit = gbif_limit, year_min = year_min,
                         max_uncertainty_m = max_uncertainty_m, exclude_cultivated = exclude_cultivated,
                         native_only = native_only, thin_km = thin_km, alt_min = alt_min,
                         alt_max = alt_max, seed = seed, cache_dir = cache_dir)

  # 2. Calibration area -------------------------------------------------------
  poly <- if (identical(region, "auto")) define_calibration_area(pts, m_buffer_km, cache_dir)
  else get_region_polygon(region, cache_dir)
  if (!is.null(poly)) {
    inside <- terra::is.related(terra::vect(pts, geom = c("decimalLongitude", "decimalLatitude"),
                                            crs = "EPSG:4326"), poly, "intersects")
    if (any(!inside)) .msg("Occurrences", sum(!inside), " records outside the calibration area removed.")
    pts <- pts[inside, ]
  }

  # 3. LULC of the calibration area (also removes records in plantations) -----
  setup_native <- .lulc_setup(lulc_source, poly, life_form, lulc_classes, savanna_habitat, lulc_year, cache_dir)
  if (exclude_cultivated && isTRUE(setup_native$mapbiomas_legend)) pts <- .drop_plantation_points(pts, setup_native$lulc)
  if (nrow(pts) < 5) stop("Fewer than 5 records remain in the calibration area.", call. = FALSE)
  utils::write.csv(pts, file.path(output_dir, paste0(sp, "_occurrences.csv")), row.names = FALSE)

  # 4. Predictors -------------------------------------------------------------
  clim <- get_climate(poly, resolution, vif_cutoff = NULL, cache_dir = cache_dir)
  static <- get_static_predictors(clim, soil_vars, soil_depth, topography, cache_dir)
  preds <- if (is.null(static)) clim else c(clim, static)
  if (!is.null(vif_cutoff)) preds <- select_predictors(preds, vif_cutoff, seed = seed)

  proj_poly <- if (!is.null(proj_region)) get_region_polygon(proj_region, cache_dir) else NULL
  fut_poly <- if (!is.null(proj_poly)) proj_poly else poly
  setup_proj <- if (!is.null(proj_poly)) .lulc_setup(lulc_source, proj_poly, life_form, lulc_classes,
                                                    savanna_habitat, lulc_year, cache_dir) else NULL
  setup_fut <- if (!is.null(proj_poly)) setup_proj else setup_native

  # 5. Host species (biotic predictors) ---------------------------------------
  host_native <- host_proj <- NULL
  host_future <- list()
  if (!is.null(host_species)) {
    .msg("Hosts", "Modelling host species first: ", paste(host_species, collapse = ", "))
    # Habitat fraction of the host life form on the climate grid (current land cover)
    host_frac <- function(setup, poly_t) {
      if (is.null(setup)) return(NULL)
      template <- get_climate(poly_t, resolution, vif_cutoff = NULL, cache_dir = cache_dir)[[1]]
      if (setup$lulc$type == "categorical") {
        h <- habitat_layer(setup$lulc, habitat_classes(host_life_form, "mapbiomas", isTRUE(savanna_habitat)),
                           resolution_m = setup$lulc$native_res_m)
      } else {
        lu <- get_lulc("esa", poly_t, esa_vars = habitat_classes(host_life_form, "esa", isTRUE(savanna_habitat)),
                       cache_dir = cache_dir)
        h <- habitat_layer(lu, NULL, resolution_m = lu$native_res_m)
      }
      .habitat_on_grid(h, template)
    }
    fr_native <- host_frac(setup_native, poly)
    fr_proj <- if (!is.null(proj_poly)) host_frac(setup_proj, proj_poly) else NULL
    fr_fut <- if (!is.null(proj_poly)) fr_proj else fr_native
    hn <- list(); hp <- list(); hf <- list()
    for (h in host_species) {
      nm <- paste0("host_", .sp_slug(h))
      hm <- cins_model(h, life_form = host_life_form, region = poly, proj_region = proj_poly,
                       future_scenario = future_scenario, future_year = future_year, gcms = gcms,
                       resolution = resolution, soil_vars = soil_vars, soil_depth = soil_depth,
                       topography = topography, lulc_source = NULL,
                       exclude_cultivated = exclude_cultivated, native_only = native_only,
                       year_min = year_min, max_uncertainty_m = max_uncertainty_m, gbif_limit = gbif_limit,
                       thin_km = thin_km, binary_threshold = binary_threshold, vif_cutoff = vif_cutoff,
                       n_pseudo = n_pseudo, algorithms = algorithms, replicates = replicates,
                       test_percent = test_percent, ncores = ncores, seed = seed, area_report = "none",
                       plot_results = FALSE, output_dir = file.path(output_dir, "hosts"), cache_dir = cache_dir)
      weight <- function(r, fr) {
        if (!is.null(fr)) r <- terra::resample(r, fr, method = "bilinear") * fr
        names(r) <- nm
        r
      }
      hn[[nm]] <- weight(hm$native$continuous, fr_native)
      if (!is.null(hm$projected)) hp[[nm]] <- weight(hm$projected$continuous, fr_proj)
      for (k in names(hm$future)) hf[[k]][[nm]] <- weight(hm$future[[k]]$continuous, fr_fut)
    }
    host_native <- terra::rast(unname(hn))
    preds <- c(preds, terra::mask(terra::resample(host_native, preds[[1]], method = "bilinear"), preds[[1]]))
    if (length(hp) > 0) host_proj <- terra::rast(unname(hp))
    host_future <- lapply(hf, function(x) terra::rast(unname(x)))
  }

  # 6. Ensemble model ---------------------------------------------------------
  fit <- fit_models(pts, preds, algorithms = algorithms, n_pseudo = n_pseudo, replicates = replicates,
                    test_percent = test_percent, ncores = ncores, binary_threshold = binary_threshold,
                    seed = seed, output_dir = output_dir, species_name = sp)

  # 7. Native realized habitat ------------------------------------------------
  native <- .realize(list(continuous = fit$ensemble_map, binary = fit$binary_map, threshold = fit$threshold),
                     setup_native, lulc_resolution, habitat_fraction_min, area_report, output_dir, sp,
                     "native", cache_dir)

  # 8. Spatial transfer -------------------------------------------------------
  projected <- NULL
  if (!is.null(proj_poly)) {
    st_p <- if (!is.null(static)) get_static_predictors(get_climate(proj_poly, resolution, NULL, cache_dir = cache_dir),
                                                       soil_vars, soil_depth, topography, cache_dir) else NULL
    pr <- project_niche(fit, proj_poly, resolution, NULL, static_layers = st_p, biotic_layers = host_proj,
                        output_dir = output_dir, species_name = sp, tag = "projected", cache_dir = cache_dir)
    projected <- .realize(pr, setup_proj, lulc_resolution, habitat_fraction_min, area_report, output_dir,
                          sp, "projected", cache_dir)
  }

  # 9. Future scenarios -------------------------------------------------------
  future <- list()
  if (!is.null(future_scenario)) {
    st_f <- if (!is.null(static)) get_static_predictors(get_climate(fut_poly, resolution, NULL, cache_dir = cache_dir),
                                                       soil_vars, soil_depth, topography, cache_dir) else NULL
    hab_cache <- if (!is.null(setup_fut)) (if (!is.null(proj_poly)) projected$habitat else native$habitat) else NULL
    for (ssp in future_scenario) for (yr in future_year) {
      key <- paste0(tolower(ssp), "_", yr)
      fr <- project_niche(fit, fut_poly, resolution, ssp, yr, gcms, static_layers = st_f,
                          biotic_layers = host_future[[key]], output_dir = output_dir,
                          species_name = sp, tag = paste0("future_", key), cache_dir = cache_dir)
      future[[key]] <- .realize(fr, setup_fut, lulc_resolution, habitat_fraction_min, area_report,
                                output_dir, sp, paste0("future_", key), cache_dir, habitat_cache = hab_cache)
    }
  }

  # 10. Settings file (reproducibility) --------------------------------------
  pk <- unique(c("cinspecies", "sdm", "terra", "geodata", "rgbif", "CoordinateCleaner", "spThin", "usdm",
                 unname(.algorithm_pkgs[algorithms])))
  ver <- vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), "")
  writeLines(c(paste("cinspecies run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               paste("Species:", species), paste("Records used:", nrow(pts)),
               paste("Threshold:", attr(fit$threshold, "method"), "=", round(fit$threshold, 4)),
               paste("Selected predictors:", paste(fit$selected_vars, collapse = ", ")),
               "", "Settings:", utils::capture.output(utils::str(settings[!vapply(settings, is.function, TRUE)])),
               "", "Package versions:", paste0("  ", names(ver), " ", ver),
               "", "Cite: run cinspecies::cins_references()"),
             file.path(output_dir, paste0(sp, "_settings.txt")))

  res <- structure(list(species = species, life_form = life_form, occurrences = pts, region = poly,
                        model = fit$model, threshold = fit$threshold, selected_vars = fit$selected_vars,
                        evaluation = fit$evaluation, var_importance = fit$var_importance,
                        native = native, projected = projected, future = future,
                        output_dir = normalizePath(output_dir), settings = settings),
                   class = c("cins_model", "list"))
  if (plot_results) plot(res)
  .msg("Done", "cinspecies finished. Outputs in: ", normalizePath(output_dir))
  res
}

#' @export
print.cins_model <- function(x, ...) {
  cat("<cins_model> ", x$species, " (", x$life_form, ")\n", sep = "")
  cat("  Records used      : ", nrow(x$occurrences), "\n", sep = "")
  cat("  Predictors        : ", paste(x$selected_vars, collapse = ", "), "\n", sep = "")
  cat("  Threshold         : ", attr(x$threshold, "method"), " = ", round(x$threshold, 4), "\n", sep = "")
  if (!is.null(x$evaluation)) {
    e <- x$evaluation
    cat("  Mean AUC / TSS    : ", round(mean(e$AUC, na.rm = TRUE), 3), " / ", round(mean(e$TSS, na.rm = TRUE), 3), "\n", sep = "")
  }
  cat("  Realized (LULC)   : ", !is.null(x$native$realized_binary), "\n", sep = "")
  cat("  Projected region  : ", !is.null(x$projected), "\n", sep = "")
  cat("  Future scenarios  : ", if (length(x$future)) paste(names(x$future), collapse = ", ") else "none", "\n", sep = "")
  cat("  Outputs           : ", x$output_dir, "\n", sep = "")
  invisible(x)
}

#' Plot the main maps of a cins_model
#' @param x A `cins_model` object.
#' @param ... Not used.
#' @export
plot.cins_model <- function(x, ...) {
  maps <- list("Native suitability" = x$native$continuous,
               "Native realized (LULC)" = x$native$realized_continuous)
  if (!is.null(x$projected)) {
    maps[["Projected suitability"]] <- x$projected$continuous
    maps[["Projected realized (LULC)"]] <- x$projected$realized_continuous
  }
  for (k in names(x$future)) {
    maps[[paste("Future", k)]] <- x$future[[k]]$continuous
    maps[[paste("Future", k, "realized")]] <- x$future[[k]]$realized_continuous
  }
  maps <- maps[!vapply(maps, is.null, logical(1))]
  n <- length(maps)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(ceiling(n / 2), min(n, 2)))
  for (nm in names(maps)) terra::plot(maps[[nm]], main = paste0(x$species, " - ", nm), range = c(0, 1))
  invisible(x)
}

#' Deprecated: use cins_model()
#'
#' `cins_pipeline()` is kept for backward compatibility and forwards to
#' [cins_model()]. `dependent_species` is passed as `host_species`, the old
#' ESA options (`"esa_trees"`, ...) become `lulc_source = "esa"`, and
#' `fine_resolution` is replaced by `lulc_resolution`. Marine modelling was
#' removed in version 0.2.0.
#'
#' @param species Scientific name.
#' @param ... Arguments passed to [cins_model()].
#' @param dependent_species Old name of `host_species`.
#' @param ecosystem Defunct.
#' @param fine_resolution Defunct.
#' @return See [cins_model()].
#' @export
cins_pipeline <- function(species, ..., dependent_species = NULL, ecosystem = NULL, fine_resolution = NULL) {
  .Deprecated("cins_model", package = "cinspecies")
  args <- list(...)
  if (identical(ecosystem, "marine")) stop("Marine modelling was removed in cinspecies 0.2.0 (terrestrial plants only).", call. = FALSE)
  if (!is.null(fine_resolution)) warning("`fine_resolution` is defunct; use `lulc_resolution` (metres).", call. = FALSE)
  if (!is.null(dependent_species)) args$host_species <- dependent_species
  if (is.character(args$lulc_source) && grepl("^esa_", args$lulc_source[1])) args$lulc_source <- "esa"
  if (identical(args$gcms, "all")) args$gcms <- "default"
  if (identical(args$binary_threshold, "maxSSS")) args$binary_threshold <- "maxTSS"
  if (is.null(args$life_form)) stop("cins_model() needs `life_form` ('tree', 'shrub', 'grassland', 'liana' or 'epiphyte').", call. = FALSE)
  do.call(cins_model, c(list(species = species), args))
}
