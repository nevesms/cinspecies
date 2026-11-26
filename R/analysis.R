# Post-modelling analyses: altitudinal shift and area change ------------------

#' @noRd
.load_dem <- function(dem, ext, cache_dir) {
  d <- if (identical(dem, "worldclim")) geodata::elevation_global(res = 0.5, path = cache_dir)
       else if (inherits(dem, "SpatRaster")) dem else terra::rast(dem)
  if (is.null(d)) stop("DEM could not be loaded.", call. = FALSE)
  if (!terra::same.crs(d, "EPSG:4326")) d <- terra::project(d, "EPSG:4326")
  terra::crop(d, ext, snap = "out")
}

#' Elevation statistics of suitable areas and altitudinal shift between periods
#'
#' The mean elevation of all cells classified as suitable is computed for each
#' period and compared with the current (first) map, quantifying vertical
#' niche shifts without imposing a fixed altitudinal filter. Continuous
#' suitability is resampled (bilinear) onto the DEM grid and binarised with
#' the model threshold (`basis = "potential"`), or the realized binary map is
#' summarised as a habitat fraction on the DEM grid (`basis = "realized"`,
#' fraction-weighted mean).
#'
#' @param suitability Named list of continuous suitability `SpatRaster`s; the
#'   first element is the reference (current) period.
#' @param threshold Numeric model threshold.
#' @param realized Optional named list of realized binary `SpatRaster`s (same names).
#' @param dem `"worldclim"` (SRTM-derived WorldClim 2.1 DEM, 30 arc-s), a file
#'   path or a `SpatRaster` (e.g. a finer DEM from AWS Terrain Tiles).
#' @param basis `"potential"` or `"realized"`.
#' @param file Optional CSV path.
#' @param cache_dir Directory where data are cached.
#' @return Data frame with `Period`, `Mean_elevation_m`, `SD_m`, `Min_m`,
#'   `Max_m`, `Shift_m` (difference to the first period).
#' @references
#' Farr, T.G. et al. (2007). The Shuttle Radar Topography Mission. Reviews of
#' Geophysics, 45, RG2004. \doi{10.1029/2005RG000183}
#' @export
elevation_shift <- function(suitability, threshold, realized = NULL, dem = "worldclim",
                            basis = c("potential", "realized"), file = NULL,
                            cache_dir = cins_cache_dir()) {
  basis <- match.arg(basis)
  if (basis == "realized" && is.null(realized)) stop("`realized` maps are needed for basis = 'realized'.", call. = FALSE)
  .msg("Elevation", "Elevation of suitable areas (", basis, ")...")
  rows <- lapply(names(suitability), function(k) {
    s <- suitability[[k]]
    d <- .load_dem(dem, terra::ext(s), cache_dir)
    if (basis == "potential") {
      sd_ <- terra::resample(s, d, method = "bilinear")
      w <- terra::ifel(sd_ >= threshold, 1, NA)
    } else {
      r <- realized[[k]]
      f <- max(1, round(terra::res(d)[1] / terra::res(r)[1]))
      fr <- if (f > 1) terra::aggregate(r, fact = f, fun = "mean", na.rm = TRUE) else r
      w <- terra::resample(fr, d, method = "bilinear")
      w <- terra::ifel(w > 0, w, NA)
    }
    z <- terra::mask(d, w)
    sw <- terra::global(w, "sum", na.rm = TRUE)[1, 1]
    if (!is.finite(sw) || sw == 0) {
      return(data.frame(Period = k, Mean_elevation_m = NA, SD_m = NA, Min_m = NA, Max_m = NA))
    }
    mu <- terra::global(z * w, "sum", na.rm = TRUE)[1, 1] / sw
    v <- terra::global((z - mu)^2 * w, "sum", na.rm = TRUE)[1, 1] / sw
    mm <- terra::global(z, c("min", "max"), na.rm = TRUE)
    data.frame(Period = k, Mean_elevation_m = round(mu, 1), SD_m = round(sqrt(v), 1),
               Min_m = round(mm[1, 1], 1), Max_m = round(mm[1, 2], 1))
  })
  out <- do.call(rbind, rows)
  out$Shift_m <- round(out$Mean_elevation_m - out$Mean_elevation_m[1], 1)
  if (!is.null(file)) utils::write.csv(out, file, row.names = FALSE)
  out
}

#' Net change, gains and losses of suitable area between periods
#'
#' Compares the current map with each future map cell by cell and reports,
#' per administrative unit, the current and future areas, the net change (%)
#' and the gained and lost areas (km2). Works with realized maps (habitat
#' weights at 30 m or fractions) or with potential binary maps.
#'
#' @param current `SpatRaster` (0-1 weights or binary) for the current period.
#' @param futures Named list of `SpatRaster`s on the same grid as `current`.
#' @param area_report `"total"`, `"continent"`, `"country"` or `"state"`.
#' @param file Optional CSV path.
#' @param cache_dir Directory where boundaries are cached.
#' @return Data frame, one row per unit.
#' @export
area_change_table <- function(current, futures, area_report = "total", file = NULL,
                              cache_dir = cins_cache_dir()) {
  if (identical(area_report, "none")) area_report <- "total"
  .msg("Area change", "Net change, gains and losses by ", area_report, "...")
  bounds <- .report_bounds(area_report, terra::ext(current), cache_dir)
  out <- .area_by_zone(current, bounds)
  names(out)[2] <- "Area_current_km2"
  for (k in names(futures)) {
    f <- futures[[k]]
    if (!terra::compareGeom(f, current, stopOnError = FALSE)) f <- terra::resample(f, current, method = "near")
    gl <- terra::lapp(c(current, f), fun = function(n, x) {
      n[is.na(n)] <- 0; x[is.na(x)] <- 0
      cbind(x, pmax(x - n, 0), pmax(n - x, 0))
    })
    a <- lapply(1:3, function(i) .area_by_zone(gl[[i]], bounds))
    tab <- data.frame(Unit = a[[1]]$Unit, A = a[[1]]$Area_km2, G = a[[2]]$Area_km2, L = a[[3]]$Area_km2)
    names(tab)[2:4] <- paste0(c("Area_", "Gain_", "Loss_"), k, "_km2")
    out <- merge(out, tab, by = "Unit", all = TRUE)
    out[[paste0("Net_change_", k, "_pct")]] <- ifelse(out$Area_current_km2 > 0,
      100 * (out[[paste0("Area_", k, "_km2")]] - out$Area_current_km2) / out$Area_current_km2, NA)
  }
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(x) round(x, 2))
  out <- out[rowSums(out[num], na.rm = TRUE) > 0, , drop = FALSE]
  out <- out[order(-out$Area_current_km2), , drop = FALSE]
  if (nrow(out) > 1) {
    tot <- out[1, , drop = FALSE]; tot$Unit <- "TOTAL"
    for (nm in names(out)[num]) tot[[nm]] <- sum(out[[nm]], na.rm = TRUE)
    for (k in names(futures)) {
      tot[[paste0("Net_change_", k, "_pct")]] <- round(100 * (tot[[paste0("Area_", k, "_km2")]] -
                                                         tot$Area_current_km2) / tot$Area_current_km2, 2)
    }
    out <- rbind(out, tot)
  }
  rownames(out) <- NULL
  if (!is.null(file)) utils::write.csv(out, file, row.names = FALSE)
  out
}
