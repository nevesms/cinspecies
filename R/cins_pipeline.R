#' Complete CINS Pipeline with Full Control, Biotic Interactions, LULC, and Future Scenarios
#'
#' @param species String. Scientific name of the target species.
#' @param region String/Vector. Native region to train the model (e.g., "South America").
#' @param proj_region String/Vector. Optional. New region to project the model.
#' @param dependent_species String/Vector. Optional. Species that the target species depends on.
#' @param future_scenario String. Optional. Shared Socioeconomic Pathway (e.g., "ssp245", "ssp585").
#' @param future_year Numeric. Optional. Projection year (e.g., 2050, 2070).
#' @param gcms Vector/String. Optional. Global Climate Models (e.g., "all", c("MIROC6", "IPSL-CM6A-LR")).
#' @param ecosystem String. "auto", "terrestrial", or "marine".
#' @param resolution Numeric. Spatial resolution in arc-minutes (e.g., 5).
#' @param fine_resolution Numeric. Optional. Final resolution in arc-minutes for resample.
#' @param binary_threshold String/Numeric. Metric for suitability cutoff ("maxSSS" or numeric 0-1).
#' @param lulc_source String or SpatRaster. Path to local LULC raster (e.g., MapBiomas .tif), SpatRaster, or "esa_trees", "esa_shrubs", "esa_grassland", "esa_cropland".
#' @param lulc_classes Numeric vector. Target class codes to keep if using a categorical raster like MapBiomas (e.g., c(3, 4) for forest).
#' @param area_report String. "none", "total", "continent", "country", or "state".
#' @param vif_cutoff Numeric. Threshold for Variance Inflation Factor.
#' @param thin_km Numeric. Minimum distance in km to thin occurrences.
#' @param n_pseudo Numeric. Number of pseudo-absences/background points.
#' @param alt_min Numeric. Minimum altitude of occurrence (m). NULL to ignore.
#' @param alt_max Numeric. Maximum altitude of occurrence (m). NULL to ignore.
#' @param algorithms Vector. Machine Learning algorithms (e.g., c("rf", "maxnet", "glm")).
#' @param replicates Numeric. Number of cross-validation replicates.
#' @param ncores Numeric. Number of CPU cores for parallel processing.
#' @param plot_results Logical. If TRUE, automatically plots the generated maps in a multi-panel layout.
#' @param output_dir String. Path to save the final maps and tables.
#'
#' @export
cins_pipeline <- function(species,
                          region = "South America",
                          proj_region = NULL,
                          dependent_species = NULL,
                          future_scenario = NULL,
                          future_year = 2050,
                          gcms = "all",
                          ecosystem = "auto",
                          resolution = 5,
                          fine_resolution = NULL,
                          binary_threshold = "maxSSS",
                          lulc_source = NULL,
                          lulc_classes = NULL,
                          area_report = "none",
                          vif_cutoff = 10,
                          thin_km = 10,
                          n_pseudo = 1000,
                          alt_min = NULL,
                          alt_max = NULL,
                          algorithms = c("rf", "maxnet", "glm"),
                          replicates = 1,
                          ncores = 1,
                          plot_results = TRUE,
                          output_dir = getwd()) {

  sp_clean <- gsub(" ", "_", species)

  # 1. Ocorrências e Clima Base
  occ_res <- get_occurrences(species = species, ecosystem = ecosystem, thin_km = thin_km, alt_min = alt_min, alt_max = alt_max)
  pts <- occ_res$data
  eco_detected <- occ_res$ecosystem

  clim_native <- get_climate(ecosystem = eco_detected, region = region, resolution = resolution, vif_cutoff = vif_cutoff)

  # 2. Modelos Dependentes (Interações Bióticas)
  biotic_native <- NULL
  biotic_proj <- NULL
  biotic_future <- NULL

  if (!is.null(dependent_species)) {
    message(">>> [Biotic Interaction] Modeling dependent species first: ", paste(dependent_species, collapse = ", "))
    dep_native_list <- list()
    dep_proj_list <- list()
    dep_fut_list <- list()

    for (dep_sp in dependent_species) {
      dep_clean <- paste0("biotic_", gsub(" ", "_", dep_sp))
      dep_mod <- cins_pipeline(
        species = dep_sp, region = region, proj_region = proj_region,
        dependent_species = NULL, future_scenario = future_scenario,
        future_year = future_year, gcms = gcms, ecosystem = eco_detected,
        resolution = resolution, fine_resolution = NULL, binary_threshold = NULL,
        lulc_source = lulc_source, lulc_classes = lulc_classes,
        area_report = "none", vif_cutoff = vif_cutoff, thin_km = thin_km,
        n_pseudo = n_pseudo, algorithms = algorithms, replicates = replicates,
        ncores = ncores, plot_results = FALSE, output_dir = output_dir
      )

      r_nat <- terra::resample(dep_mod$ensemble_map, clim_native, method = "bilinear")
      names(r_nat) <- dep_clean
      dep_native_list[[dep_clean]] <- r_nat

      if (!is.null(proj_region) && !is.null(dep_mod$projected_map)) {
        r_prj <- dep_mod$projected_map
        names(r_prj) <- dep_clean
        dep_proj_list[[dep_clean]] <- r_prj
      }
      if (!is.null(future_scenario) && !is.null(dep_mod$future_map)) {
        r_fut <- dep_mod$future_map
        names(r_fut) <- dep_clean
        dep_fut_list[[dep_clean]] <- r_fut
      }
    }
    biotic_native <- terra::rast(dep_native_list)
    clim_native <- c(clim_native, biotic_native)
    if (length(dep_proj_list) > 0) biotic_proj <- terra::rast(dep_proj_list)
    if (length(dep_fut_list) > 0) biotic_future <- terra::rast(dep_fut_list)
  }

  # 3. Ajuste do Modelo Nativo + Filtro de Uso do Solo (MapBiomas / ESA)
  res_native <- fit_models(
    occurrences = pts, climate = clim_native, algorithms = algorithms,
    n_pseudo = n_pseudo, replicates = replicates, ncores = ncores,
    fine_resolution = fine_resolution, binary_threshold = binary_threshold,
    output_dir = output_dir, species_name = sp_clean
  )

  lulc_native <- apply_lulc_mask(
    cont_raster = res_native$ensemble_map,
    bin_raster = res_native$binary_map,
    lulc_source = lulc_source,
    lulc_classes = lulc_classes
  )
  res_native$realized_ensemble_map <- lulc_native$continuous
  res_native$realized_binary_map <- lulc_native$binary

  res_native$area_table_native <- calculate_area_report(
    binary_raster = res_native$binary_map,
    lulc_binary_raster = if (!is.null(lulc_source)) res_native$realized_binary_map else NULL,
    area_report = area_report,
    output_dir = output_dir,
    species_name = sp_clean,
    suffix = "native"
  )

  # 4. Transferabilidade Espacial (Clima Atual em Nova Região) + Filtro LULC
  if (!is.null(proj_region)) {
    proj_res <- project_niche(
      cins_result = res_native, target_region = proj_region, ecosystem = eco_detected,
      resolution = resolution, fine_resolution = fine_resolution, binary_threshold = binary_threshold,
      future_scenario = NULL, biotic_layers = biotic_proj, output_dir = output_dir,
      species_name = sp_clean, tag = "projected_current"
    )
    res_native$projected_map <- proj_res$continuous
    res_native$projected_binary <- proj_res$binary

    lulc_proj <- apply_lulc_mask(
      cont_raster = proj_res$continuous,
      bin_raster = proj_res$binary,
      lulc_source = lulc_source,
      lulc_classes = lulc_classes
    )
    res_native$projected_realized_map <- lulc_proj$continuous
    res_native$projected_realized_binary <- lulc_proj$binary

    res_native$area_table_projected <- calculate_area_report(
      binary_raster = proj_res$binary,
      lulc_binary_raster = if (!is.null(lulc_source)) res_native$projected_realized_binary else NULL,
      area_report = area_report,
      output_dir = output_dir,
      species_name = sp_clean,
      suffix = "projected"
    )
  }

  # 5. Projeção de Cenários Futuros (CMIP6 / GCMs) + Filtro LULC
  if (!is.null(future_scenario)) {
    target_fut_reg <- if (is.null(proj_region)) region else proj_region
    fut_res <- project_niche(
      cins_result = res_native, target_region = target_fut_reg, ecosystem = eco_detected,
      resolution = resolution, fine_resolution = fine_resolution, binary_threshold = binary_threshold,
      future_scenario = future_scenario, future_year = future_year, gcms = gcms,
      biotic_layers = if (!is.null(proj_region)) biotic_proj else biotic_future,
      output_dir = output_dir, species_name = sp_clean, tag = paste0("future_", future_scenario, "_", future_year)
    )
    res_native$future_map <- fut_res$continuous
    res_native$future_binary <- fut_res$binary

    lulc_fut <- apply_lulc_mask(
      cont_raster = fut_res$continuous,
      bin_raster = fut_res$binary,
      lulc_source = lulc_source,
      lulc_classes = lulc_classes
    )
    res_native$future_realized_map <- lulc_fut$continuous
    res_native$future_realized_binary <- lulc_fut$binary

    res_native$area_table_future <- calculate_area_report(
      binary_raster = fut_res$binary,
      lulc_binary_raster = if (!is.null(lulc_source)) res_native$future_realized_binary else NULL,
      area_report = area_report,
      output_dir = output_dir,
      species_name = sp_clean,
      suffix = paste0("future_", future_scenario, "_", future_year)
    )
  }

  # 6. Painel Visual Automático (Suporta extensões geográficas diferentes)
  if (plot_results) {
    plots_list <- list()
    plots_list[["1. Native Climate"]] <- res_native$ensemble_map
    if (!is.null(lulc_source) && !is.null(res_native$realized_ensemble_map)) {
      plots_list[["2. Native Realized (LULC)"]] <- res_native$realized_ensemble_map
    }
    if (!is.null(res_native$projected_map)) {
      plots_list[["3. Projected Climate"]] <- res_native$projected_map
      if (!is.null(lulc_source) && !is.null(res_native$projected_realized_map)) {
        plots_list[["4. Projected Realized (LULC)"]] <- res_native$projected_realized_map
      }
    }
    if (!is.null(res_native$future_map)) {
      plots_list[[paste0("5. Future (", future_scenario, " - ", future_year, ")")]] <- res_native$future_map
      if (!is.null(lulc_source) && !is.null(res_native$future_realized_map)) {
        plots_list[[paste0("6. Future Realized (LULC)")]] <- res_native$future_realized_map
      }
    }

    n_p <- length(plots_list)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(ceiling(n_p / 2), ifelse(n_p > 1, 2, 1)))
    for (nm in names(plots_list)) {
      terra::plot(plots_list[[nm]], main = paste0(species, " - ", nm))
    }
  }

  terra::tmpFiles(current = TRUE, remove = TRUE)
  message(">>> [SUCCESS] CINS Pipeline completed! Outputs saved in: ", output_dir)
  return(res_native)
}
