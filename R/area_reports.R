# Geodetic area reports -------------------------------------------------------

#' @noRd
.report_bounds <- function(level, ext, cache_dir) {
  if (level %in% c("total", "none")) return(NULL)
  world <- terra::crop(.world_polygons(cache_dir), ext)
  if (level == "continent") {
    world$Unit <- world$continent
    return(terra::aggregate(world, by = "Unit"))
  }
  if (level == "country") { world$Unit <- world$NAME_0; return(world[, "Unit"]) }
  if (level == "state") {
    st <- terra::crop(.states(unique(world$GID_0), cache_dir), ext)
    st$Unit <- paste0(st$NAME_1, " (", st$GID_0, ")")
    return(st[, "Unit"])
  }
  stop("`area_report` must be 'none', 'total', 'continent', 'country' or 'state'.", call. = FALSE)
}

#' Area (km2) of a 0-1 weight raster per zone, memory-safe for 30 m grids
#'
#' The weight raster is summed into blocks of ~1 km, multiplied by the exact
#' geodetic area of the fine cells, and summed per zone. Zones are rasterised
#' by cell centre, so every cell is assigned to a single unit (no double
#' counting along borders).
#' @noRd
.area_by_zone <- function(w, bounds) {
  res_m <- terra::res(w)[1] * 111320
  f <- max(1, floor(1000 / res_m))
  agg <- if (f > 1) terra::aggregate(w, fact = f, fun = "sum", na.rm = TRUE) else w
  km2 <- agg * terra::cellSize(agg, unit = "km") / (f^2)
  if (is.null(bounds)) {
    return(data.frame(Unit = "Total", Area_km2 = sum(terra::global(km2, "sum", na.rm = TRUE)[, 1])))
  }
  bounds$zone_id <- seq_len(nrow(bounds))
  z <- terra::rasterize(bounds, km2, field = "zone_id")
  zs <- terra::zonal(km2, z, fun = "sum", na.rm = TRUE)
  data.frame(Unit = bounds$Unit[zs[, 1]], Area_km2 = zs[, 2])
}

#' Potential (climatic) and realized (LULC) suitable area per administrative unit
#'
#' @param binary_raster `SpatRaster`. Potential binary map (1 = suitable). When
#'   LULC is used, pass the potential binary at the habitat grid so that both
#'   areas share the same grid.
#' @param realized_weight `SpatRaster` or `NULL`. Realized habitat weight
#'   (binary x habitat fraction) from [apply_lulc_mask()].
#' @param area_report `"none"`, `"total"`, `"continent"`, `"country"` or `"state"`.
#' @param output_dir,species_name,suffix Output file naming.
#' @param cache_dir Directory where boundaries are cached.
#' @return Data frame with `Unit`, `Potential_Climatic_km2`, and, with LULC,
#'   `Realized_LULC_km2` and `Habitat_Loss_pct`; also written to CSV.
#' @export
calculate_area_report <- function(binary_raster, realized_weight = NULL, area_report = "none",
                                  output_dir = getwd(), species_name = "Species", suffix = "native",
                                  cache_dir = cins_cache_dir()) {
  if (identical(area_report, "none") || is.null(binary_raster)) return(NULL)
  .msg("Area report", "Geodetic area (km2) by ", area_report, " (", suffix, ")...")
  bounds <- .report_bounds(area_report, terra::ext(binary_raster), cache_dir)
  out <- .area_by_zone(binary_raster, bounds)
  names(out)[2] <- "Potential_Climatic_km2"
  if (!is.null(realized_weight)) {
    rl <- .area_by_zone(realized_weight, bounds)
    names(rl)[2] <- "Realized_LULC_km2"
    out <- merge(out, rl, by = "Unit", all.x = TRUE)
    out$Realized_LULC_km2[is.na(out$Realized_LULC_km2)] <- 0
    out$Habitat_Loss_pct <- ifelse(out$Potential_Climatic_km2 > 0,
                                   100 * (1 - out$Realized_LULC_km2 / out$Potential_Climatic_km2), NA)
  }
  out <- out[out$Potential_Climatic_km2 > 0, , drop = FALSE]
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], round, 2)
  out <- out[order(-out$Potential_Climatic_km2), , drop = FALSE]
  rownames(out) <- NULL
  f <- file.path(output_dir, paste0(species_name, "_area_", area_report, "_", suffix, ".csv"))
  utils::write.csv(out, f, row.names = FALSE)
  .msg("Area report", "Saved: ", f)
  out
}
