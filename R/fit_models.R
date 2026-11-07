#' Post-Process Raster: Resampling and Binary Thresholding
#'
#' @export
post_process_raster <- function(raster_map, sdm_model, fine_resolution = NULL, binary_threshold = "maxSSS") {

  if (!is.null(fine_resolution)) {
    message(">>> [Post-Process] Resampling map to fine resolution (", fine_resolution, " arc-min)...")
    target_deg <- fine_resolution / 60
    template <- terra::rast(terra::ext(raster_map), resolution = target_deg, crs = terra::crs(raster_map))
    raster_map <- terra::resample(raster_map, template, method = "bilinear")
  }

  bin_map <- NULL
  th_val <- 0.5

  if (!is.null(binary_threshold)) {
    evals <- tryCatch(sdm::getEvaluation(sdm_model, stat = c("TSS", "threshold"), opt = 2), error = function(e) NULL)
    if (!is.null(evals) && "threshold" %in% colnames(evals)) {
      th_val <- mean(evals$threshold, na.rm = TRUE)
    } else if (is.numeric(binary_threshold)) {
      th_val <- binary_threshold
    }
    message(">>> [Post-Process] Applying binary threshold (cutoff = ", round(th_val, 3), ")...")
    bin_map <- terra::ifel(raster_map >= th_val, 1, 0)
    names(bin_map) <- "binary_suitability"
  }

  return(list(continuous = raster_map, binary = bin_map, threshold = th_val))
}

#' Fit Multi-Algorithm SDM and Generate Weighted Ensemble
#'
#' @export
fit_models <- function(occurrences, climate, algorithms = c("rf", "maxnet", "glm"), n_pseudo = 1000, replicates = 1, ncores = 1, fine_resolution = NULL, binary_threshold = "maxSSS", output_dir = getwd(), species_name = "Species") {

  require(sdm)
  pts_vect <- terra::vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  pts_vect$Presenca <- 1
  pts_vect <- pts_vect[, "Presenca"]

  message(">>> [Modeling] Building sdmData with ", n_pseudo, " pseudo-absences...")
  d <- sdm::sdmData(formula = Presenca ~ ., train = pts_vect, predictors = climate, bg = list(n = n_pseudo, method = "gRandom"))

  parallel_cfg <- if (ncores > 1) list(ncore = ncores, method = "parallel") else NULL

  message(">>> [Modeling] Fitting ", paste(algorithms, collapse = ", "), " (", replicates, " replicate[s], ", ncores, " core[s])...")
  if (replicates > 1) {
    m <- sdm::sdm(Presenca ~ ., data = d, methods = algorithms, replication = "sub", test.percent = 30, n = replicates, parallelSettings = parallel_cfg)
  } else {
    m <- sdm::sdm(Presenca ~ ., data = d, methods = algorithms, replication = "sub", test.percent = 30, n = 1, parallelSettings = parallel_cfg)
  }

  out_file <- file.path(output_dir, paste0(species_name, "_ensemble_native.tif"))
  message(">>> [Modeling] Calculating TSS-weighted ensemble...")
  ens_map <- sdm::ensemble(m, newdata = climate, filename = out_file, setting = list(method = "weighted", stat = "TSS", opt = 2), overwrite = TRUE)

  processed <- post_process_raster(ens_map, m, fine_resolution = fine_resolution, binary_threshold = binary_threshold)

  if (!is.null(processed$binary)) {
    bin_file <- file.path(output_dir, paste0(species_name, "_binary_native.tif"))
    terra::writeRaster(processed$binary, filename = bin_file, overwrite = TRUE)
  }

  return(list(model = m, ensemble_map = processed$continuous, binary_map = processed$binary, threshold = processed$threshold, selected_vars = names(climate)))
}
#' Apply Land Use and Land Cover (LULC) Filter to Suitability Maps
#'
#' @param cont_raster SpatRaster. Continuous suitability map.
#' @param bin_raster SpatRaster. Binary suitability map.
#' @param lulc_source String or SpatRaster. Path to local raster (e.g., MapBiomas .tif), SpatRaster object, or "esa_trees", "esa_shrubs", "esa_grassland".
#' @param lulc_classes Numeric vector. Target class codes to keep if using a categorical raster like MapBiomas (e.g., c(3, 4) for forest).
#'
#' @export
apply_lulc_mask <- function(cont_raster, bin_raster = NULL, lulc_source = NULL, lulc_classes = NULL) {
  if (is.null(lulc_source)) return(list(continuous = cont_raster, binary = bin_raster))

  message(">>> [Land Use] Integrating Land Use / Land Cover (LULC) filter...")

  # Opção 1: Download automático global (ESA WorldCover via geodata)
  if (is.character(lulc_source) && lulc_source %in% c("esa_trees", "esa_shrubs", "esa_grassland", "esa_cropland")) {
    var_type <- gsub("esa_", "", lulc_source)
    lulc_raw <- geodata::landcover(var = var_type, path = tempdir())
    lulc_crop <- terra::crop(lulc_raw, terra::ext(cont_raster))
    lulc_aligned <- terra::resample(lulc_crop, cont_raster, method = "bilinear")

    # Mantém áreas com pelo menos 10% de cobertura vegetacional no pixel
    lulc_mask <- terra::ifel(lulc_aligned >= 0.10, 1, 0)

  } else {
    # Opção 2: Raster customizado do usuário (ex: MapBiomas 30m)
    lulc_raw <- if (is.character(lulc_source)) terra::rast(lulc_source) else lulc_source
    lulc_crop <- terra::crop(lulc_raw, terra::ext(cont_raster))

    # Se o MapBiomas for mais fino que o clima, refina o mapa climático para a resolução do MapBiomas
    if (terra::res(lulc_crop)[1] < terra::res(cont_raster)[1]) {
      message(">>> [Land Use] Aligning climate model to high-resolution LULC grid...")
      cont_raster <- terra::resample(cont_raster, lulc_crop, method = "bilinear")
      if (!is.null(bin_raster)) bin_raster <- terra::resample(bin_raster, lulc_crop, method = "near")
      lulc_aligned <- lulc_crop
    } else {
      lulc_aligned <- terra::resample(lulc_crop, cont_raster, method = "near")
    }

    # Filtra os códigos específicos do MapBiomas (ex: 3 = Formação Florestal)
    if (!is.null(lulc_classes)) {
      lulc_mask <- terra::ifel(lulc_aligned %in% lulc_classes, 1, 0)
    } else {
      lulc_mask <- terra::ifel(lulc_aligned > 0, 1, 0)
    }
  }

  # Multiplica o clima adequado pela presença de vegetação/uso válido
  cont_masked <- cont_raster * lulc_mask
  bin_masked <- if (!is.null(bin_raster)) bin_raster * lulc_mask else NULL

  return(list(continuous = cont_masked, binary = bin_masked))
}
