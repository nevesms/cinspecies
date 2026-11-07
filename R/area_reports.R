#' Calculate Suitable and Realized Area (km2) with Habitat Loss Metrics
#'
#' @export
calculate_area_report <- function(binary_raster, lulc_binary_raster = NULL, area_report = "none", output_dir = getwd(), species_name = "Species", suffix = "native") {

  if (area_report == "none" || is.null(binary_raster)) return(NULL)

  message(">>> [Area Report] Computing geodetic area (km2) by level: ", area_report)

  # Função interna auxiliar para somar área
  calc_zonal <- function(r, bounds, col_nm) {
    suite <- terra::ifel(r == 1, 1, NA)
    if (is.null(bounds)) {
      val <- terra::expanse(suite, unit = "km")
      return(data.frame(Unit = "Total Study Area", Area_km2 = round(sum(val$area, na.rm = TRUE), 2)))
    }
    px_km2 <- terra::cellSize(suite, unit = "km") * suite
    z_sum <- terra::extract(px_km2, bounds, fun = sum, na.rm = TRUE)
    df_r <- data.frame(Unit = bounds[[col_nm]][, 1], Area_km2 = z_sum[, 2])
    df_agg <- stats::aggregate(Area_km2 ~ Unit, data = df_r, FUN = sum, na.rm = TRUE)
    df_agg$Area_km2 <- round(df_agg$Area_km2, 2)
    return(df_agg)
  }

  bounds <- NULL
  col_name <- NULL
  if (area_report %in% c("continent", "country")) {
    bounds <- terra::crop(rnaturalearth::ne_countries(scale = "medium", returnclass = "sv"), terra::ext(binary_raster))
    col_name <- ifelse(area_report == "continent", "continent", "name")
  } else if (area_report == "state") {
    bounds <- terra::crop(rnaturalearth::ne_states(returnclass = "sv"), terra::ext(binary_raster))
    col_name <- "name"
  }

  df_out <- calc_zonal(binary_raster, bounds, col_name)
  colnames(df_out)[2] <- "Potential_Climatic_km2"

  # Se usou filtro de uso do solo (MapBiomas/ESA), calcula a perda de habitat
  if (!is.null(lulc_binary_raster)) {
    df_lulc <- calc_zonal(lulc_binary_raster, bounds, col_name)
    colnames(df_lulc)[2] <- "Realized_LULC_km2"
    df_out <- merge(df_out, df_lulc, by = "Unit", all.x = TRUE)
    df_out$Realized_LULC_km2[is.na(df_out$Realized_LULC_km2)] <- 0
    df_out$Habitat_Loss_pct <- round((1 - (df_out$Realized_LULC_km2 / df_out$Potential_Climatic_km2)) * 100, 2)
  }

  df_out <- df_out[df_out$Potential_Climatic_km2 > 0, ]
  df_out <- df_out[order(-df_out$Potential_Climatic_km2), ]
  rownames(df_out) <- NULL

  csv_file <- file.path(output_dir, paste0(species_name, "_area_report_", area_report, "_", suffix, ".csv"))
  utils::write.csv(df_out, file = csv_file, row.names = FALSE)
  message(">>> [Area Report] Saved CSV table to: ", csv_file)

  return(df_out)
}
