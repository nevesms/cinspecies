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
.lulc_setup <- function(lulc_source, poly, life_form, lulc_classes, savanna_habitat, lulc_year, cache_dir,
                        collection = 11) {
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
  lulc <- get_lulc(src, poly, year = lulc_year, collection = collection, cache_dir = cache_dir)
  .msg("Land use", "Habitat classes: ", paste(classes, collapse = ", "))
  list(lulc = lulc, classes = classes,
       mapbiomas_legend = identical(src, "mapbiomas") || is.null(lulc_classes))
}

#' Realized maps, conflict table and area report for one period
#' @noRd
.realize <- function(res, setup, opts, sp, tag, habitat_cache = NULL) {
  out <- list(continuous = res$continuous, binary = res$binary, uncertainty = res$uncertainty)
  if (is.null(setup)) {
    out$area_table <- calculate_area_report(res$binary, NULL, opts$area_report, opts$output_dir, sp, tag, opts$cache_dir)
    out$change_weight <- res$binary
    return(out)
  }
  na_cl <- if (opts$mask_water && isTRUE(setup$mapbiomas_legend)) .water_classes else NULL
  hab <- if (!is.null(habitat_cache)) habitat_cache else
    habitat_layer(setup$lulc, setup$classes, opts$lulc_resolution, na_classes = na_cl,
                  filename = file.path(opts$output_dir, paste0(sp, "_", tag, "_habitat.tif")))
  .lulc_size_message(hab, max(opts$lulc_resolution, setup$lulc$native_res_m))
  rl <- apply_lulc_mask(res$continuous, res$threshold, hab, opts$habitat_fraction_min,
                        file_prefix = file.path(opts$output_dir, paste0(sp, "_", tag)))
  out$realized_continuous <- rl$continuous
  out$realized_binary <- rl$binary
  out$potential_binary_fine <- rl$potential_binary
  out$change_weight <- rl$weight
  out$habitat <- hab
  out$area_table <- calculate_area_report(rl$potential_binary, rl$weight, opts$area_report, opts$output_dir,
                                          sp, tag, opts$cache_dir)
  if (opts$conflict_table && setup$lulc$type == "categorical" && isTRUE(setup$mapbiomas_legend)) {
    .msg("Land use", "Land-use conflict table (", tag, ")...")
    groups <- if (is.null(opts$lulc_groups)) mapbiomas_groups() else opts$lulc_groups
    out$conflict_table <- lulc_conflict_table(rl$potential_binary, setup$lulc, setup$classes, groups,
                                              opts$mask_water,
                                              file.path(opts$output_dir, paste0(sp, "_", tag, "_lulc_conflict.csv")))
  }
  out
}

#' Niche and realized-habitat model of terrestrial plant species
#'
#' `cins_model()` runs the whole cinspecies workflow for one or several
#' species of trees, shrubs, natural-grassland species, lianas and epiphytes.
#' Defaults reproduce a published workflow (VIF < 5, three background sets,
#' three 80/20 sub-sampling repetitions, TSS-weighted ensemble of models with
#' TSS >= 0.5, maxTSS threshold, climate downscaled to the 30 m MapBiomas grid,
#' water masked, land-use conflict table and altitudinal shift), and every
#' step can be changed through the arguments below.
#'
#' 1. **Records** - GBIF (or user) records; kingdom Plantae check; removal of
#'    fossil, cultivated/living-collection and managed records; temporal and
#'    positional-uncertainty filters; CoordinateCleaner; optional elevation
#'    filter; spatial thinning ([get_occurrences()]); records in forest
#'    plantations (MapBiomas class 9) removed.
#' 2. **Calibration area (M)** - automatic ([define_calibration_area()]) or
#'    user-defined ([get_region_polygon()]).
#' 3. **Predictors** - WorldClim 2.1, optional SoilGrids and topography,
#'    stepwise VIF ([select_predictors()]).
#' 4. **Host species** - modelled first; suitability x habitat enters as a
#'    biotic predictor (epiphytes, lianas, dependent species).
#' 5. **Ensemble** - any algorithms of [list_algorithms()] (machine learning,
#'    neural networks, maximum entropy, regression, envelopes), several
#'    background sets and partitions, ensemble with optional performance
#'    filter and uncertainty map ([fit_models()]).
#' 6. **Land use and land cover** - suitability downscaled to 30 m and combined
#'    with native habitat chosen from `life_form` ([habitat_classes()],
#'    [apply_lulc_mask()]); land-use conflict table ([lulc_conflict_table()]).
#' 7. **Transfer** to another region and to CMIP6 futures (any WorldClim GCM,
#'    or user climate through `future_climate`) ([project_niche()]).
#' 8. **Reports** - areas per unit ([calculate_area_report()]), net change,
#'    gains and losses ([area_change_table()]) and altitudinal shift
#'    ([elevation_shift()]).
#'
#' **Several species.** If `species` has more than one name, each species is
#' modelled in its own sub-folder and a `cins_multi` object is returned with a
#' summary table and richness maps (sum of binary maps). `life_form` and
#' `savanna_habitat` may be vectors with one value per species; `occurrences`,
#' `host_species` and `lulc_classes` may be lists with one element per species
#' (named by species or in the same order). A single value applies to all.
#'
#' @param species Character. One or more scientific names.
#' @param life_form Character (required). `"tree"`, `"shrub"`, `"grassland"`,
#'   `"liana"` or `"epiphyte"` (one per species, or one for all).
#' @param occurrences Optional user records (data frame, `SpatVector`, CSV path;
#'   a list for several species).
#' @param region Calibration area: `"auto"`, continent, country,
#'   `"Country:State"` vector, `SpatVector`, `sf` or file path.
#' @param m_buffer_km Numeric. Buffer (km) for `region = "auto"`.
#' @param proj_region Optional region to transfer the model to.
#' @param host_species Optional host / phorophyte species (character; a list for several species).
#' @param host_life_form Life form of the hosts.
#' @param future_scenario `NULL` or SSP codes (`"ssp126"`, `"ssp245"`, `"ssp370"`, `"ssp585"`).
#' @param future_year One or more years (2021-2100).
#' @param gcms `"default"` (5 GCMs), `"all"` (26 GCMs) or a vector from [list_gcms()].
#' @param future_climate Optional named list of user climate for any GCM or
#'   source not distributed by WorldClim: each element (named e.g.
#'   `"ssp245_2070"`) is a `SpatRaster`, file path, or list of them (averaged).
#'   Used instead of `future_scenario`/`gcms`.
#' @param resolution Climate resolution (arc-min): 10, 5, 2.5 or 0.5.
#' @param soil_vars,soil_depth,topography Static predictors (see [get_static_predictors()]).
#' @param lulc_source `"auto"`, `"mapbiomas"`, `"esa"`, raster path/`SpatRaster`, or `NULL`.
#' @param lulc_year Integer. MapBiomas year.
#' @param mapbiomas_collection Integer. 11 (default) or 10.
#' @param lulc_classes Optional habitat classes overriding `life_form`.
#' @param savanna_habitat Logical, required with MapBiomas: Savanna Formation (4) as habitat?
#' @param lulc_resolution Numeric. Cell size (m) of realized maps (default 30;
#'   coarser values reduce accuracy and trigger a warning).
#' @param habitat_fraction_min Numeric. Minimum habitat fraction of coarse cells.
#' @param mask_water Logical. Exclude water bodies from areas and tables.
#' @param conflict_table Logical. Land-use conflict table per period.
#' @param lulc_groups Optional named list grouping LULC classes (default [mapbiomas_groups()]).
#' @param exclude_cultivated,native_only,year_min,max_uncertainty_m,gbif_limit,cc_tests,thin_km,alt_min,alt_max
#'   Record filters (see [get_occurrences()]).
#' @param binary_threshold `"maxTSS"` (default), another criterion or a number in (0, 1).
#' @param vif_cutoff Numeric or `NULL`.
#' @param algorithms Character. Codes from [list_algorithms()].
#' @param n_pseudo,n_pseudo_sets,bg_method Background settings.
#' @param cv_method,replicates,test_percent,cv_folds Partitioning settings.
#' @param ensemble_method,ensemble_stat,ensemble_min Ensemble settings (see [ensemble_setting()]).
#' @param uncertainty Logical. Coefficient-of-variation maps.
#' @param elevation_shift Logical. Altitudinal-shift table (with futures).
#' @param dem DEM for `elevation_shift` (`"worldclim"`, path or `SpatRaster`).
#' @param elevation_basis `"potential"` or `"realized"`.
#' @param ncores Numeric. CPU cores.
#' @param seed Integer or `NULL`.
#' @param area_report `"none"`, `"total"`, `"continent"`, `"country"` or `"state"`.
#' @param plot_results Logical.
#' @param output_dir Directory for outputs.
#' @param cache_dir Directory for downloaded data.
#' @return A `cins_model` object (one species) or a `cins_multi` object
#'   (several species). A `cins_model` contains `occurrences`, `region`,
#'   `model`, `threshold`, `selected_vars`, `evaluation`, `var_importance`,
#'   `native`, `projected`, `future`, `area_change`, `elevation` and
#'   `settings`; `native`, `projected` and each `future[[i]]` contain
#'   `continuous`, `binary`, `uncertainty`, `realized_continuous`,
#'   `realized_binary`, `habitat`, `area_table` and `conflict_table`.
#' @references See [cins_references()].
#' @examples
#' \dontrun{
#' m <- cins_model(
#'   species = "Ilex paraguariensis", life_form = "tree", savanna_habitat = TRUE,
#'   region = c("Brazil", "Argentina", "Paraguay", "Uruguay"),
#'   proj_region = c("Brazil:Parana", "Brazil:Santa Catarina", "Brazil:Rio Grande do Sul"),
#'   future_scenario = c("ssp245", "ssp585"), future_year = 2070,
#'   gcms = c("MIROC6", "IPSL-CM6A-LR", "MRI-ESM2-0"),
#'   algorithms = c("rf", "maxnet", "glm", "mlp"), area_report = "state", seed = 42
#' )
#' m$native$conflict_table
#' m$area_change
#' m$elevation
#'
#' # several species at once
#' mm <- cins_model(c("Araucaria angustifolia", "Ocotea porosa"), life_form = "tree",
#'                  savanna_habitat = FALSE, region = "auto")
#' mm$summary
#' }
#' @export
cins_model <- function(species, life_form, occurrences = NULL,
                       region = "auto", m_buffer_km = 300, proj_region = NULL,
                       host_species = NULL, host_life_form = "tree",
                       future_scenario = NULL, future_year = 2070, gcms = "default",
                       future_climate = NULL,
                       resolution = 5, soil_vars = NULL, soil_depth = 5, topography = FALSE,
                       lulc_source = "auto", lulc_year = 2024, mapbiomas_collection = 11,
                       lulc_classes = NULL, savanna_habitat = NULL, lulc_resolution = 30,
                       habitat_fraction_min = 0.5, mask_water = TRUE, conflict_table = TRUE,
                       lulc_groups = NULL,
                       exclude_cultivated = TRUE, native_only = FALSE, year_min = 1970,
                       max_uncertainty_m = 10000, gbif_limit = 20000,
                       cc_tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"),
                       thin_km = 5, alt_min = NULL, alt_max = NULL,
                       binary_threshold = "maxTSS", vif_cutoff = 5,
                       algorithms = c("rf", "maxnet", "glm"),
                       n_pseudo = 1000, n_pseudo_sets = 3, bg_method = "gRandom",
                       cv_method = "sub", replicates = 3, test_percent = 20, cv_folds = 5,
                       ensemble_method = "weighted", ensemble_stat = "TSS", ensemble_min = 0.5,
                       uncertainty = TRUE, elevation_shift = TRUE, dem = "worldclim",
                       elevation_basis = "potential",
                       ncores = 1, seed = NULL, area_report = "none", plot_results = TRUE,
                       output_dir = getwd(), cache_dir = cins_cache_dir()) {

  # 0. Checks (fail fast, before any download) --------------------------------
  if (missing(life_form) || is.null(life_form) || !all(life_form %in% .life_forms)) {
    stop("`life_form` is required: one of ", paste(.life_forms, collapse = ", "), ".", call. = FALSE)
  }
  if (length(species) > 1) {
    args <- mget(names(formals()), envir = environment())
    return(.cins_multi(args))
  }
  if (length(life_form) > 1) stop("Give one `life_form` per species.", call. = FALSE)
  if (!is.null(host_species) && !host_life_form %in% .life_forms) stop("Invalid `host_life_form`.", call. = FALSE)
  .check_resolution(resolution)
  .check_algorithms(algorithms)
  .validate_lulc_args(lulc_source, lulc_classes, savanna_habitat, lulc_resolution)
  if (!is.null(future_scenario) && is.null(future_climate)) {
    invisible(lapply(future_year, .future_period)); .resolve_gcms(gcms)
  }
  if (!is.null(future_climate) && (!is.list(future_climate) || is.null(names(future_climate)) ||
                                   any(!nzchar(names(future_climate))))) {
    stop("`future_climate` must be a named list (e.g. list(ssp245_2070 = c('gcm1.tif', 'gcm2.tif'))).", call. = FALSE)
  }
  if (!area_report %in% c("none", "total", "continent", "country", "state")) stop("Invalid `area_report`.", call. = FALSE)
  if (!elevation_basis %in% c("potential", "realized")) stop("Invalid `elevation_basis`.", call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  sp <- .sp_slug(species)
  settings <- mget(setdiff(names(formals()), c("occurrences", "future_climate")), envir = environment())
  opts <- list(area_report = area_report, output_dir = output_dir, cache_dir = cache_dir,
               mask_water = mask_water, lulc_resolution = lulc_resolution,
               habitat_fraction_min = habitat_fraction_min, conflict_table = conflict_table,
               lulc_groups = lulc_groups)

  # 1. Records ----------------------------------------------------------------
  pts <- get_occurrences(species, occurrences = occurrences, gbif_limit = gbif_limit, year_min = year_min,
                         max_uncertainty_m = max_uncertainty_m, exclude_cultivated = exclude_cultivated,
                         native_only = native_only, cc_tests = cc_tests, thin_km = thin_km,
                         alt_min = alt_min, alt_max = alt_max, seed = seed, cache_dir = cache_dir)

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
  lsetup <- function(p) .lulc_setup(lulc_source, p, life_form, lulc_classes, savanna_habitat, lulc_year,
                                    cache_dir, mapbiomas_collection)
  setup_native <- lsetup(poly)
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
  setup_proj <- if (!is.null(proj_poly)) lsetup(proj_poly) else NULL
  setup_fut <- if (!is.null(proj_poly)) setup_proj else setup_native
  fut_keys <- if (!is.null(future_climate)) names(future_climate)
              else if (!is.null(future_scenario)) as.vector(outer(tolower(future_scenario), future_year, paste, sep = "_"))
              else character(0)

  # 5. Host species (biotic predictors) ---------------------------------------
  host_proj <- NULL
  host_future <- list()
  if (!is.null(host_species)) {
    .msg("Hosts", "Modelling host species first: ", paste(host_species, collapse = ", "))
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
                       future_climate = future_climate, resolution = resolution, soil_vars = soil_vars,
                       soil_depth = soil_depth, topography = topography, lulc_source = NULL,
                       exclude_cultivated = exclude_cultivated, native_only = native_only,
                       year_min = year_min, max_uncertainty_m = max_uncertainty_m, gbif_limit = gbif_limit,
                       cc_tests = cc_tests, thin_km = thin_km, binary_threshold = binary_threshold,
                       vif_cutoff = vif_cutoff, algorithms = algorithms, n_pseudo = n_pseudo,
                       n_pseudo_sets = n_pseudo_sets, bg_method = bg_method, cv_method = cv_method,
                       replicates = replicates, test_percent = test_percent, cv_folds = cv_folds,
                       ensemble_method = ensemble_method, ensemble_stat = ensemble_stat,
                       ensemble_min = ensemble_min, uncertainty = FALSE, elevation_shift = FALSE,
                       ncores = ncores, seed = seed, area_report = "none", plot_results = FALSE,
                       output_dir = file.path(output_dir, "hosts"), cache_dir = cache_dir)
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
  fit <- fit_models(pts, preds, algorithms = algorithms, n_pseudo = n_pseudo, n_pseudo_sets = n_pseudo_sets,
                    bg_method = bg_method, cv_method = cv_method, replicates = replicates,
                    test_percent = test_percent, cv_folds = cv_folds, ensemble_method = ensemble_method,
                    ensemble_stat = ensemble_stat, ensemble_min = ensemble_min, uncertainty = uncertainty,
                    ncores = ncores, binary_threshold = binary_threshold, seed = seed,
                    output_dir = output_dir, species_name = sp)

  # 7. Native realized habitat ------------------------------------------------
  native <- .realize(list(continuous = fit$ensemble_map, binary = fit$binary_map,
                          uncertainty = fit$uncertainty_map, threshold = fit$threshold),
                     setup_native, opts, sp, "native")

  # 8. Spatial transfer -------------------------------------------------------
  projected <- NULL
  if (!is.null(proj_poly)) {
    st_p <- if (!is.null(static)) get_static_predictors(get_climate(proj_poly, resolution, NULL, cache_dir = cache_dir),
                                                       soil_vars, soil_depth, topography, cache_dir) else NULL
    pr <- project_niche(fit, proj_poly, resolution, NULL, static_layers = st_p, biotic_layers = host_proj,
                        uncertainty = uncertainty, output_dir = output_dir, species_name = sp,
                        tag = "projected", cache_dir = cache_dir)
    projected <- .realize(pr, setup_proj, opts, sp, "projected")
  }

  # 9. Future scenarios -------------------------------------------------------
  future <- list()
  if (length(fut_keys) > 0) {
    st_f <- if (!is.null(static)) get_static_predictors(get_climate(fut_poly, resolution, NULL, cache_dir = cache_dir),
                                                       soil_vars, soil_depth, topography, cache_dir) else NULL
    ref <- if (!is.null(proj_poly)) projected else native
    for (key in fut_keys) {
      parts <- strsplit(key, "_")[[1]]
      fr <- project_niche(fit, fut_poly, resolution,
                          future_scenario = if (is.null(future_climate)) parts[1] else NULL,
                          future_year = if (is.null(future_climate)) as.numeric(parts[2]) else 2070,
                          gcms = gcms, static_layers = st_f, biotic_layers = host_future[[key]],
                          new_climate = if (!is.null(future_climate)) future_climate[[key]] else NULL,
                          uncertainty = uncertainty, output_dir = output_dir, species_name = sp,
                          tag = paste0("future_", key), cache_dir = cache_dir)
      future[[key]] <- .realize(fr, setup_fut, opts, sp, paste0("future_", key), habitat_cache = ref$habitat)
    }
  }

  # 10. Change and elevation reports ------------------------------------------
  area_change <- NULL
  elevation <- NULL
  if (length(future) > 0) {
    ref <- if (!is.null(proj_poly)) projected else native
    area_change <- area_change_table(ref$change_weight, lapply(future, `[[`, "change_weight"),
                                     area_report, file.path(output_dir, paste0(sp, "_area_change.csv")), cache_dir)
    if (elevation_shift) {
      suit <- c(list(current = ref$continuous), lapply(future, `[[`, "continuous"))
      real <- if (!is.null(ref$realized_binary)) c(list(current = ref$realized_binary), lapply(future, `[[`, "realized_binary")) else NULL
      basis <- if (elevation_basis == "realized" && is.null(real)) "potential" else elevation_basis
      elevation <- elevation_shift(suit, fit$threshold, real, dem, basis,
                                   file.path(output_dir, paste0(sp, "_elevation_shift.csv")), cache_dir)
    }
  }

  # 11. Settings file (reproducibility) --------------------------------------
  pk <- unique(c("cinspecies", "sdm", "terra", "geodata", "rgbif", "CoordinateCleaner", "spThin", "usdm",
                 .algorithm_package(algorithms)))
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
                        area_change = area_change, elevation = elevation,
                        output_dir = normalizePath(output_dir), settings = settings),
                   class = c("cins_model", "list"))
  if (plot_results) plot(res)
  .msg("Done", "cinspecies finished for ", species, ". Outputs in: ", normalizePath(output_dir))
  res
}

#' Run several species and combine them
#' @noRd
.cins_multi <- function(args) {
  spp <- args$species
  n <- length(spp)
  pick <- function(x, i, atomic) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && !is.data.frame(x) && !inherits(x, c("SpatVector", "SpatRaster"))) {
      if (!is.null(names(x)) && spp[i] %in% names(x)) return(x[[spp[i]]])
      if (length(x) != n) stop("Per-species lists must have one element per species.", call. = FALSE)
      return(x[[i]])
    }
    if (atomic && length(x) == n) return(x[i])
    x
  }
  per_atomic <- c("life_form", "savanna_habitat")
  per_list <- c("occurrences", "host_species", "lulc_classes")
  if (length(args$life_form) > 1 && length(args$life_form) != n) stop("`life_form` must have length 1 or one per species.", call. = FALSE)
  out <- list()
  for (i in seq_len(n)) {
    a <- args
    a$species <- spp[i]
    for (p in per_atomic) a[p] <- list(pick(args[[p]], i, TRUE))
    for (p in per_list) a[p] <- list(pick(args[[p]], i, FALSE))
    a$output_dir <- file.path(args$output_dir, .sp_slug(spp[i]))
    a$plot_results <- FALSE
    .msg("Species", i, "/", n, ": ", spp[i])
    out[[spp[i]]] <- do.call(cins_model, a)
  }
  combine_species(out, output_dir = args$output_dir, cache_dir = args$cache_dir)
}

#' Combine several cins_model results
#'
#' Builds a summary table (records, predictors, mean AUC/TSS, threshold and
#' areas per period) and richness maps - the number of species predicted
#' suitable per cell - for the potential (climate grid) and realized (LULC grid)
#' binary maps of every period shared by the models.
#'
#' @param models Named list of `cins_model` objects.
#' @param output_dir Directory for the richness rasters and summary CSV.
#' @param cache_dir Not used directly; kept for consistency.
#' @return A `cins_multi` object: `models`, `summary`, `richness`.
#' @export
combine_species <- function(models, output_dir = getwd(), cache_dir = cins_cache_dir()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tot <- function(r) if (is.null(r)) NA else round(.area_by_zone(r, NULL)$Area_km2, 2)
  rows <- lapply(names(models), function(k) {
    m <- models[[k]]
    e <- m$evaluation
    ens <- if (!is.null(e) && "in_ensemble" %in% names(e)) e[e$in_ensemble %in% TRUE, ] else e
    r <- data.frame(species = k, life_form = m$life_form, records = nrow(m$occurrences),
                    predictors = paste(m$selected_vars, collapse = ";"),
                    AUC = if (is.null(ens)) NA else round(mean(ens$AUC, na.rm = TRUE), 3),
                    TSS = if (is.null(ens)) NA else round(mean(ens$TSS, na.rm = TRUE), 3),
                    threshold = round(as.numeric(m$threshold), 3),
                    potential_current_km2 = tot(if (!is.null(m$native$potential_binary_fine)) m$native$potential_binary_fine else m$native$binary),
                    realized_current_km2 = tot(m$native$change_weight))
    for (f in names(m$future)) r[[paste0("realized_", f, "_km2")]] <- tot(m$future[[f]]$change_weight)
    r
  })
  cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] })
  summary <- do.call(rbind, rows)
  utils::write.csv(summary, file.path(output_dir, "species_summary.csv"), row.names = FALSE)

  rich <- function(layers, file) {
    layers <- layers[!vapply(layers, is.null, logical(1))]
    if (length(layers) < 2) return(NULL)
    r <- tryCatch(terra::mosaic(terra::sprc(lapply(layers, function(x) terra::classify(x, cbind(NA, 0)))), fun = "sum"),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    names(r) <- "richness"
    .write_rast(r, file, "INT2U")
  }
  .msg("Species", "Richness maps...")
  periods <- c("native", Reduce(intersect, lapply(models, function(m) names(m$future))))
  richness <- list()
  for (p in periods) {
    get <- function(m, what) if (p == "native") m$native[[what]] else m$future[[p]][[what]]
    richness[[p]] <- list(
      potential = rich(lapply(models, get, "binary"), file.path(output_dir, paste0("richness_", p, "_potential.tif"))),
      realized = rich(lapply(models, get, "realized_binary"), file.path(output_dir, paste0("richness_", p, "_realized.tif"))))
  }
  structure(list(models = models, summary = summary, richness = richness, output_dir = normalizePath(output_dir)),
            class = c("cins_multi", "list"))
}

#' @export
print.cins_multi <- function(x, ...) {
  cat("<cins_multi> ", length(x$models), " species\n\n", sep = "")
  print(x$summary[, intersect(c("species", "life_form", "records", "AUC", "TSS", "threshold",
                                "potential_current_km2", "realized_current_km2"), names(x$summary))],
        row.names = FALSE)
  cat("\nOutputs: ", x$output_dir, "\n", sep = "")
  invisible(x)
}

#' Plot richness maps of a cins_multi
#' @param x A `cins_multi` object.
#' @param ... Not used.
#' @export
plot.cins_multi <- function(x, ...) {
  maps <- unlist(lapply(names(x$richness), function(p) {
    r <- x$richness[[p]]; r <- r[!vapply(r, is.null, logical(1))]
    stats::setNames(r, paste(p, names(r)))
  }), recursive = FALSE)
  if (length(maps) == 0) return(invisible(x))
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(ceiling(length(maps) / 2), min(length(maps), 2)))
  for (nm in names(maps)) terra::plot(maps[[nm]], main = paste("Richness -", nm))
  invisible(x)
}

#' @export
print.cins_model <- function(x, ...) {
  cat("<cins_model> ", x$species, " (", x$life_form, ")\n", sep = "")
  cat("  Records used      : ", nrow(x$occurrences), "\n", sep = "")
  cat("  Predictors        : ", paste(x$selected_vars, collapse = ", "), "\n", sep = "")
  cat("  Threshold         : ", attr(x$threshold, "method"), " = ", round(x$threshold, 4), "\n", sep = "")
  if (!is.null(x$evaluation)) {
    e <- x$evaluation
    ens <- if ("in_ensemble" %in% names(e)) e[e$in_ensemble %in% TRUE, ] else e
    cat("  Models (ensemble) : ", nrow(e), " fitted, ", nrow(ens), " in the ensemble\n", sep = "")
    cat("  Mean AUC / TSS    : ", round(mean(ens$AUC, na.rm = TRUE), 3), " / ", round(mean(ens$TSS, na.rm = TRUE), 3),
        " (ensemble members)\n", sep = "")
  }
  cat("  Realized (LULC)   : ", !is.null(x$native$realized_binary), "\n", sep = "")
  cat("  Projected region  : ", !is.null(x$projected), "\n", sep = "")
  cat("  Future scenarios  : ", if (length(x$future)) paste(names(x$future), collapse = ", ") else "none", "\n", sep = "")
  if (!is.null(x$native$conflict_table)) {
    ct <- x$native$conflict_table
    av <- ct$Area_pct[ct$Category == "Available habitat"]
    cat("  Available habitat : ", if (length(av)) av else 0, "% of the current suitable area\n", sep = "")
  }
  if (!is.null(x$elevation)) {
    cat("  Elevation shift   : ", paste0(x$elevation$Period[-1], " ", sprintf("%+.1f", x$elevation$Shift_m[-1]), " m",
                                         collapse = "; "), "\n", sep = "")
  }
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
