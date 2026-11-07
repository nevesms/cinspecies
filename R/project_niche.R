#' Project Trained SDM to New Geographic Region or Future Scenario
#'
#' @export
project_niche <- function(cins_result, target_region, ecosystem = "terrestrial", resolution = 5, fine_resolution = NULL, binary_threshold = "maxSSS", future_scenario = NULL, future_year = 2050, gcms = "all", biotic_layers = NULL, output_dir = getwd(), species_name = "Species", tag = "projected") {

  require(sdm)
  base_vars <- setdiff(cins_result$selected_vars, if (!is.null(biotic_layers)) names(biotic_layers) else character(0))

  if (is.null(future_scenario)) {
    message(">>> [Projection] Loading current climate for target region: ", paste(target_region, collapse = ", "))
    if (ecosystem == "terrestrial") {
      new_clim <- geodata::worldclim_global(var = "bio", res = resolution, path = tempdir())
    } else {
      marine_layers <- c("BO_sstmean", "BO_salinity", "BO_bathymean")
      new_clim <- terra::rast(sdmpredictors::load_layers(marine_layers, datadir = tempdir()))
    }
    poly <- get_region_polygon(target_region)
    if (!is.null(poly)) new_clim <- terra::crop(new_clim, poly, mask = TRUE)
    new_clim <- new_clim[[base_vars]]
  } else {
    new_clim <- get_future_climate(region = target_region, resolution = resolution, future_scenario = future_scenario, future_year = future_year, gcms = gcms, selected_vars = base_vars)
  }

  if (!is.null(biotic_layers)) {
    biotic_res <- terra::resample(biotic_layers, new_clim, method = "bilinear")
    new_clim <- c(new_clim, biotic_res)
  }

  out_file <- file.path(output_dir, paste0(species_name, "_", tag, "_", paste(gsub(" ", "_", target_region), collapse = "_"), ".tif"))

  proj_map <- sdm::ensemble(x = cins_result$model, newdata = new_clim, filename = out_file, setting = list(method = "weighted", stat = "TSS", opt = 2), overwrite = TRUE)

  processed <- post_process_raster(proj_map, cins_result$model, fine_resolution = fine_resolution, binary_threshold = binary_threshold)

  if (!is.null(processed$binary)) {
    bin_file <- file.path(output_dir, paste0(species_name, "_", tag, "_binary_", paste(gsub(" ", "_", target_region), collapse = "_"), ".tif"))
    terra::writeRaster(processed$binary, filename = bin_file, overwrite = TRUE)
  }

  return(processed)
}
