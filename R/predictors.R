# Environmental predictors ----------------------------------------------------

#' @noRd
.crop_region <- function(r, poly) {
  if (is.null(poly)) return(r)
  terra::crop(r, poly, mask = TRUE)
}

#' Download baseline bioclimatic variables (WorldClim 2.1)
#'
#' @param region Region definition (see [get_region_polygon()]) or a `SpatVector`.
#' @param resolution Numeric. 10, 5, 2.5 or 0.5 arc-minutes.
#' @param vif_cutoff Numeric or `NULL`. If not `NULL`, collinear variables are
#'   removed with a stepwise VIF procedure ([select_predictors()]).
#' @param seed Integer or `NULL`. Seed for the VIF background sample.
#' @param cache_dir Directory where data are cached.
#' @return `SpatRaster` with layers `bio1` ... `bio19` (or the retained subset).
#' @references
#' Fick, S.E. & Hijmans, R.J. (2017). WorldClim 2: new 1-km spatial resolution
#' climate surfaces for global land areas. International Journal of
#' Climatology, 37, 4302-4315. \doi{10.1002/joc.5086}
#' @export
get_climate <- function(region = "global", resolution = 2.5, vif_cutoff = 5, seed = NULL,
                        cache_dir = cins_cache_dir()) {
  .check_resolution(resolution)
  .msg("Climate", "Loading WorldClim 2.1 bioclimatic variables (", resolution, " arc-min)...")
  clim <- geodata::worldclim_global(var = "bio", res = resolution, path = cache_dir)
  if (is.null(clim)) stop("WorldClim download failed.", call. = FALSE)
  clim <- .std_bio_names(clim)
  poly <- if (inherits(region, "SpatVector")) region else get_region_polygon(region, cache_dir)
  clim <- .crop_region(clim, poly)
  if (!is.null(vif_cutoff)) clim <- select_predictors(clim, vif_cutoff = vif_cutoff, seed = seed)
  clim
}

#' Remove collinear predictors with stepwise VIF
#'
#' Variance inflation factors are computed on a random sample of up to 10000
#' background cells of the calibration area (`usdm::vifstep`), and the variable
#' with the highest VIF is removed iteratively until all VIF < `vif_cutoff`.
#'
#' @param predictors `SpatRaster`.
#' @param vif_cutoff Numeric. Maximum VIF (5-10 recommended; Dormann et al. 2013).
#' @param keep Character. Layers never removed (e.g. host-suitability layers).
#' @param seed Integer or `NULL`.
#' @return The `SpatRaster` with the retained layers.
#' @references
#' Naimi, B. et al. (2014). Where is positional uncertainty a problem for
#' species distribution modelling? Ecography, 37, 191-203.
#' \doi{10.1111/j.1600-0587.2013.00205.x}
#'
#' Dormann, C.F. et al. (2013). Collinearity: a review of methods to deal with
#' it and a simulation study evaluating their performance. Ecography, 36, 27-46.
#' \doi{10.1111/j.1600-0587.2012.07348.x}
#' @export
select_predictors <- function(predictors, vif_cutoff = 5, keep = NULL, seed = NULL) {
  test <- setdiff(names(predictors), keep)
  if (length(test) < 2) return(predictors)
  .msg("Predictors", "Stepwise VIF (cutoff = ", vif_cutoff, ")...")
  if (!is.null(seed)) set.seed(seed)
  df <- terra::spatSample(predictors[[test]], size = 10000, method = "random", na.rm = TRUE, as.df = TRUE)
  v <- usdm::vifstep(df, th = vif_cutoff)
  sel <- c(as.character(v@results$Variables), keep)
  .msg("Predictors", "Retained: ", paste(sel, collapse = ", "))
  predictors[[sel]]
}

#' Download CMIP6 bioclimatic projections and average several GCMs
#'
#' @param region Region definition or `SpatVector`.
#' @param resolution Numeric. 10, 5, 2.5 or 0.5 arc-minutes.
#' @param future_scenario Character. `"ssp126"`, `"ssp245"`, `"ssp370"` or `"ssp585"`.
#' @param future_year Numeric (2021-2100); mapped to the WorldClim periods
#'   2021-2040, 2041-2060, 2061-2080, 2081-2100.
#' @param gcms Character vector of GCMs, `"default"` (ACCESS-CM2, IPSL-CM6A-LR,
#'   MIROC6, MPI-ESM1-2-HR, MRI-ESM2-0) or `"all"`.
#' @param selected_vars Character. Layers to keep (e.g. from the calibrated model).
#' @param cache_dir Directory where data are cached.
#' @return `SpatRaster` with the cell-wise mean across GCMs.
#' @references
#' Eyring, V. et al. (2016). Overview of the Coupled Model Intercomparison
#' Project Phase 6 (CMIP6) experimental design and organization. Geoscientific
#' Model Development, 9, 1937-1958. \doi{10.5194/gmd-9-1937-2016}
#' @export
get_future_climate <- function(region = "global", resolution = 2.5, future_scenario = "ssp585",
                               future_year = 2050, gcms = "default", selected_vars = NULL,
                               cache_dir = cins_cache_dir()) {
  .check_resolution(resolution)
  gcms <- .resolve_gcms(gcms)
  period <- .future_period(future_year)
  ssp <- gsub("ssp", "", tolower(future_scenario))
  if (!ssp %in% c("126", "245", "370", "585")) stop("Invalid `future_scenario`: ", future_scenario, call. = FALSE)
  poly <- if (inherits(region, "SpatVector")) region else get_region_polygon(region, cache_dir)

  stack <- list()
  for (g in gcms) {
    .msg("Future climate", "GCM ", g, " | SSP", ssp, " | ", period)
    r <- tryCatch(geodata::cmip6_world(model = g, ssp = ssp, time = period, var = "bioc",
                                       res = resolution, path = cache_dir),
                  error = function(e) { .msg("Future climate", "skipped ", g, ": ", conditionMessage(e)); NULL })
    if (is.null(r)) next
    names(r) <- paste0("bio", seq_len(terra::nlyr(r)))
    r <- .crop_region(r, poly)
    if (!is.null(selected_vars)) r <- r[[intersect(selected_vars, names(r))]]
    stack[[g]] <- r
  }
  if (length(stack) == 0) stop("No GCM could be loaded for SSP", ssp, " ", period, ".", call. = FALSE)
  if (length(stack) == 1) return(stack[[1]])
  .msg("Future climate", "Averaging ", length(stack), " GCMs cell by cell...")
  vars <- names(stack[[1]])
  out <- terra::rast(lapply(vars, function(v) terra::app(terra::rast(lapply(stack, `[[`, v)), mean, na.rm = TRUE)))
  names(out) <- vars
  out
}

#' Static edaphic and topographic predictors
#'
#' Soil properties from SoilGrids 2.0 (30 arc-s) and elevation/slope from the
#' SRTM-derived WorldClim DEM, aggregated (mean) to the climate grid. These
#' layers are held constant in future projections.
#'
#' @param template `SpatRaster` defining the target grid (climate layers).
#' @param soil_vars Character or `NULL`. SoilGrids variables, e.g.
#'   `c("clay", "sand", "phh2o", "soc", "cec", "nitrogen", "bdod")`.
#' @param soil_depth Numeric. Lower limit (cm) of the SoilGrids depth interval:
#'   5 (0-5 cm), 15 (5-15 cm), 30 (15-30 cm), 60, 100 or 200.
#' @param topography Logical. Add `elevation` and `slope` (degrees).
#' @param cache_dir Directory where data are cached.
#' @return `SpatRaster` aligned to `template`, or `NULL`.
#' @references
#' Poggio, L. et al. (2021). SoilGrids 2.0: producing soil information for the
#' globe with quantified spatial uncertainty. SOIL, 7, 217-240.
#' \doi{10.5194/soil-7-217-2021}
#'
#' Farr, T.G. et al. (2007). The Shuttle Radar Topography Mission. Reviews of
#' Geophysics, 45, RG2004. \doi{10.1029/2005RG000183}
#' @export
get_static_predictors <- function(template, soil_vars = NULL, soil_depth = 5, topography = FALSE,
                                  cache_dir = cins_cache_dir()) {
  out <- list()
  to_grid <- function(r) {
    r <- terra::crop(r, terra::ext(template), snap = "out")
    f <- max(1, round(terra::res(template)[1] / terra::res(r)[1]))
    if (f > 1) r <- terra::aggregate(r, fact = f, fun = "mean", na.rm = TRUE)
    terra::mask(terra::resample(r, template, method = "bilinear"), template[[1]])
  }
  if (!is.null(soil_vars)) {
    .msg("Predictors", "Loading SoilGrids 2.0: ", paste(soil_vars, collapse = ", "), " (", soil_depth, " cm)...")
    for (v in soil_vars) {
      s <- geodata::soil_world(var = v, depth = soil_depth, stat = "mean", path = cache_dir)
      if (is.null(s)) stop("SoilGrids download failed for ", v, call. = FALSE)
      s <- to_grid(s)
      names(s) <- paste0("soil_", v)
      out[[v]] <- s
    }
  }
  if (topography) {
    .msg("Predictors", "Loading elevation and slope (30 arc-s DEM)...")
    dem <- geodata::elevation_global(res = 0.5, path = cache_dir)
    dem <- terra::crop(dem, terra::ext(template), snap = "out")
    slope <- terra::terrain(dem, v = "slope", unit = "degrees")
    e <- to_grid(dem); names(e) <- "elevation"
    s <- to_grid(slope); names(s) <- "slope"
    out[["elevation"]] <- e
    out[["slope"]] <- s
  }
  if (length(out) == 0) return(NULL)
  terra::rast(unname(out))
}
