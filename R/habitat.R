# Land use and land cover (LULC): habitat of terrestrial plants ---------------

#' @noRd
.life_forms <- c("tree", "shrub", "grassland", "liana", "epiphyte")

#' @noRd
.mapbiomas_default_url <- paste0("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/",
                                 "collection_10/lulc/coverage/brazil_coverage_%d.tif")

#' Native-habitat land-cover classes for a plant life form
#'
#' Translates a life form into the land-cover classes treated as available
#' habitat.
#'
#' **MapBiomas Brazil (Collection 10 legend)**
#' * `tree`, `liana`, `epiphyte`: Forest Formation (3), Mangrove (5),
#'   Floodable Forest (6), Wooded Sandbank Vegetation / restinga arborea (49).
#' * `shrub`: Grassland (12), Rocky Outcrop (29), Wooded Sandbank Vegetation (49),
#'   Herbaceous Sandbank Vegetation (50).
#' * `grassland` (natural grasslands, campos): Wetland (11), Grassland (12),
#'   Rocky Outcrop (29), Herbaceous Sandbank Vegetation (50).
#' * Savanna Formation (4) is added to any life form when
#'   `savanna_habitat = TRUE`. Because the decision is ecological (e.g. Cerrado
#'   trees vs. Atlantic Forest trees), it must be set explicitly by the user:
#'   an error is raised if it is left `NULL`. Typical choice: `TRUE` for shrubs
#'   and for Cerrado species of any life form, `FALSE` for forest species.
#'
#' **ESA WorldCover (fractional cover, 30 arc-s, via `geodata::landcover`)**
#' * `tree`, `liana`, `epiphyte`: trees + mangroves.
#' * `shrub`: shrubs (+ trees if `savanna_habitat = TRUE`).
#' * `grassland`: grassland + wetland.
#'
#' Forest Plantation (MapBiomas class 9) is never habitat; records falling in
#' it are removed when `exclude_cultivated = TRUE` in [cins_model()].
#'
#' @param life_form One of `"tree"`, `"shrub"`, `"grassland"`, `"liana"`, `"epiphyte"`.
#' @param lulc_source `"mapbiomas"` (or any categorical raster with the
#'   MapBiomas legend) or `"esa"`.
#' @param savanna_habitat Logical. Treat Savanna Formation as habitat. Required
#'   for MapBiomas.
#' @return Integer class codes (MapBiomas) or ESA layer names.
#' @references
#' Souza, C.M. Jr. et al. (2020). Reconstructing three decades of land use and
#' land cover changes in Brazilian biomes with Landsat archive and Earth Engine.
#' Remote Sensing, 12, 2735. \doi{10.3390/rs12172735}
#'
#' Zanaga, D. et al. (2022). ESA WorldCover 10 m 2021 v200. Zenodo.
#' \doi{10.5281/zenodo.7254221}
#' @examples
#' habitat_classes("tree", "mapbiomas", savanna_habitat = FALSE)
#' habitat_classes("shrub", "mapbiomas", savanna_habitat = TRUE)
#' habitat_classes("grassland", "esa")
#' @export
habitat_classes <- function(life_form, lulc_source = "mapbiomas", savanna_habitat = NULL) {
  if (is.null(life_form) || !life_form %in% .life_forms) {
    stop("`life_form` must be one of: ", paste(.life_forms, collapse = ", "), call. = FALSE)
  }
  forest <- c("tree", "liana", "epiphyte")
  if (identical(lulc_source, "esa")) {
    v <- if (life_form %in% forest) c("trees", "mangroves")
    else if (life_form == "shrub") c("shrubs", if (isTRUE(savanna_habitat)) "trees")
    else c("grassland", "wetland")
    return(v)
  }
  if (is.null(savanna_habitat) || !is.logical(savanna_habitat) || is.na(savanna_habitat)) {
    stop("`savanna_habitat` must be set to TRUE or FALSE when MapBiomas classes are derived from ",
         "`life_form`: should Savanna Formation (class 4) count as habitat? Usually TRUE for shrubs ",
         "and Cerrado species, FALSE for forest species. Alternatively give `lulc_classes` directly.",
         call. = FALSE)
  }
  cl <- if (life_form %in% forest) c(3L, 5L, 6L, 49L)
  else if (life_form == "shrub") c(12L, 29L, 49L, 50L)
  else c(11L, 12L, 29L, 50L)
  sort(unique(c(cl, if (savanna_habitat) 4L)))
}

#' Load land use / land cover for a region
#'
#' * `"mapbiomas"`: MapBiomas Brazil Collection 10 (30 m, 1985-2024). The
#'   region window is read directly from the cloud-optimised GeoTIFF through
#'   GDAL `/vsicurl/` (only the needed tiles are transferred) and cached; if
#'   that fails the whole annual map is downloaded. Another collection can be
#'   used with `options(cinspecies.mapbiomas_url = "<url with %d for the year>")`.
#' * `"esa"`: ESA WorldCover 2021 fractional cover aggregated to 30 arc-s by
#'   `geodata::landcover()` (global).
#' * `"auto"`: MapBiomas if at least 95% of the region is in Brazil, ESA otherwise.
#' * a file path or `SpatRaster`: user raster. Integer rasters are treated as
#'   categorical (MapBiomas legend unless `lulc_classes` is given); float rasters
#'   with values in 0-1 as habitat fractions.
#'
#' @param lulc_source See above.
#' @param region `SpatVector` polygon (or region definition).
#' @param year Integer. MapBiomas year (1985-2024 for Collection 10).
#' @param esa_vars Character. ESA layers to sum (see [habitat_classes()]).
#' @param cache_dir Directory where data are cached.
#' @return A list with `raster`, `type` (`"categorical"` or `"fraction"`),
#'   `source` and `native_res_m`.
#' @references
#' Souza, C.M. Jr. et al. (2020). Remote Sensing, 12, 2735. \doi{10.3390/rs12172735}
#'
#' Zanaga, D. et al. (2022). ESA WorldCover 10 m 2021 v200. \doi{10.5281/zenodo.7254221}
#' @export
get_lulc <- function(lulc_source = "auto", region, year = 2024, esa_vars = NULL,
                     cache_dir = cins_cache_dir()) {
  poly <- if (inherits(region, "SpatVector")) region else get_region_polygon(region, cache_dir)
  if (identical(lulc_source, "auto")) {
    share <- .share_in_brazil(poly, cache_dir)
    lulc_source <- if (share >= 0.95) "mapbiomas" else "esa"
    .msg("Land use", "lulc_source = 'auto' -> ", lulc_source, " (", round(100 * share), "% of region in Brazil).")
  }

  if (identical(lulc_source, "mapbiomas")) {
    url <- sprintf(getOption("cinspecies.mapbiomas_url", .mapbiomas_default_url), as.integer(year))
    e <- round(as.vector(terra::ext(poly)), 3)
    key <- gsub("[^0-9a-zA-Z]+", "_", paste(c(year, e), collapse = "_"))
    f <- file.path(cache_dir, "mapbiomas", paste0("mapbiomas_", key, ".tif"))
    if (!file.exists(f)) {
      dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
      .msg("Land use", "Reading MapBiomas ", year, " window (30 m) from ", url)
      src <- tryCatch(terra::rast(paste0("/vsicurl/", url)), error = function(e) NULL)
      if (is.null(src)) {
        full <- file.path(cache_dir, "mapbiomas", basename(url))
        if (!file.exists(full)) {
          .msg("Land use", "Remote read failed; downloading the full annual map (several GB)...")
          old <- options(timeout = max(7200, getOption("timeout"))); on.exit(options(old), add = TRUE)
          utils::download.file(url, full, mode = "wb")
        }
        src <- terra::rast(full)
      }
      terra::crop(src, poly, mask = TRUE, filename = f, datatype = "INT1U", overwrite = TRUE,
                  gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"))
    }
    r <- terra::rast(f)
    return(list(raster = r, type = "categorical", source = "mapbiomas",
                native_res_m = terra::res(r)[1] * 111320))
  }

  if (identical(lulc_source, "esa")) {
    if (is.null(esa_vars)) stop("`esa_vars` is required for ESA WorldCover.", call. = FALSE)
    .msg("Land use", "Loading ESA WorldCover fractional cover: ", paste(esa_vars, collapse = " + "))
    lays <- lapply(esa_vars, function(v) {
      x <- geodata::landcover(var = v, path = cache_dir)
      if (is.null(x)) stop("ESA WorldCover download failed for ", v, call. = FALSE)
      terra::crop(x, poly, mask = TRUE)
    })
    r <- if (length(lays) == 1) lays[[1]] else terra::app(terra::rast(lays), sum, na.rm = TRUE)
    r <- terra::clamp(r, 0, 1)
    names(r) <- "habitat_fraction"
    return(list(raster = r, type = "fraction", source = "esa", native_res_m = terra::res(r)[1] * 111320))
  }

  r <- if (inherits(lulc_source, "SpatRaster")) lulc_source else terra::rast(lulc_source)
  if (!terra::same.crs(r, "EPSG:4326")) {
    .msg("Land use", "Projecting user LULC raster to EPSG:4326 (nearest neighbour)...")
    r <- terra::project(r, "EPSG:4326", method = "near")
  }
  r <- terra::crop(r, poly, mask = TRUE)
  type <- if (grepl("^INT", terra::datatype(r)[1]) || isTRUE(terra::is.int(r)[1])) "categorical" else "fraction"
  list(raster = r, type = type, source = "user", native_res_m = terra::res(r)[1] * 111320)
}

#' Habitat layer at the requested resolution
#'
#' Converts land cover into habitat: 1/0 at the native resolution of a
#' categorical map, or the fraction of habitat per cell when aggregated to a
#' coarser `resolution_m`.
#'
#' @param lulc Output of [get_lulc()].
#' @param classes Class codes (categorical maps).
#' @param resolution_m Numeric. Target cell size in metres (default 30, the
#'   MapBiomas native resolution).
#' @param filename Optional output file.
#' @return `SpatRaster` with values 0-1.
#' @export
habitat_layer <- function(lulc, classes = NULL, resolution_m = 30, filename = "") {
  r <- lulc$raster
  if (lulc$type == "categorical") {
    if (is.null(classes)) stop("Habitat classes are required for categorical LULC.", call. = FALSE)
    h <- terra::mask(terra::classify(r, cbind(classes, 1L), others = 0L), r)
  } else {
    h <- r
  }
  native <- lulc$native_res_m
  if (resolution_m < native * 0.9) {
    warning("The finest resolution available for this LULC source is ~", round(native), " m; ",
            "`lulc_resolution = ", resolution_m, "` cannot be honoured and the native grid is used.",
            call. = FALSE)
  }
  f <- round(resolution_m / native)
  if (f > 1) h <- terra::aggregate(h, fact = f, fun = "mean", na.rm = TRUE)
  names(h) <- "habitat"
  if (nzchar(filename)) h <- .write_rast(h, filename, if (f > 1 || lulc$type == "fraction") "FLT4S" else "INT1U")
  h
}

#' @noRd
.lulc_size_message <- function(habitat, resolution_m) {
  n <- terra::ncell(habitat)
  gb <- n * 4 / 1024^3
  .msg("Land use", "Realized-habitat grid: ", format(round(n), big.mark = ","), " cells at ",
       round(resolution_m), " m (~", round(gb, 1), " GB per float layer before compression).")
  if (n > 2e9) {
    warning("Very large 30 m grid (", format(round(n), big.mark = ","), " cells): expect long runtimes and ",
            "large disk use. Consider a smaller region, or a coarser `lulc_resolution` (at the cost of ",
            "spatial accuracy).", call. = FALSE)
  }
}

#' Apply land-use / land-cover habitat to suitability maps
#'
#' The continuous climatic suitability is resampled (bilinear) onto the habitat
#' grid (30 m by default) and combined with habitat:
#' * realized suitability = suitability x habitat;
#' * potential binary = suitability >= threshold (at the habitat grid);
#' * realized binary = potential binary AND habitat >= `habitat_fraction_min`;
#' * realized weight = potential binary x habitat (used for fraction-weighted
#'   areas when habitat is a fraction; equal to realized binary at 30 m).
#'
#' @param cont_raster `SpatRaster`. Continuous suitability (0-1).
#' @param threshold Numeric. Binary threshold of the model.
#' @param habitat `SpatRaster` from [habitat_layer()].
#' @param habitat_fraction_min Numeric. Minimum habitat fraction for a coarse
#'   cell to be mapped as realized habitat (ignored at native 30 m).
#' @param file_prefix Character. Path prefix for the output GeoTIFFs.
#' @return List of `SpatRaster`: `continuous`, `potential_binary`, `binary`, `weight`.
#' @export
apply_lulc_mask <- function(cont_raster, threshold, habitat, habitat_fraction_min = 0.5,
                            file_prefix = tempfile("cins_")) {
  .msg("Land use", "Downscaling suitability to the habitat grid and applying LULC...")
  s <- terra::resample(cont_raster, habitat, method = "bilinear")
  th <- threshold
  hmin <- habitat_fraction_min
  out <- terra::lapp(c(s, habitat), fun = function(s, h) {
    pb <- as.numeric(s >= th)
    cbind(s * h, pb, as.numeric(pb == 1 & h >= hmin & h > 0), pb * h)
  }, filename = paste0(file_prefix, "_lulc_stack.tif"), overwrite = TRUE,
  wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER")))
  names(out) <- c("realized_suitability", "potential_binary", "realized_binary", "realized_weight")
  list(continuous = out[[1]], potential_binary = out[[2]], binary = out[[3]], weight = out[[4]])
}

#' Habitat fraction on a coarse grid (e.g. climate grid, for host layers)
#' @noRd
.habitat_on_grid <- function(habitat, template) {
  f <- max(1, round(terra::res(template)[1] / terra::res(habitat)[1]))
  h <- if (f > 1) terra::aggregate(habitat, fact = f, fun = "mean", na.rm = TRUE) else habitat
  terra::mask(terra::resample(h, template, method = "bilinear"), template[[1]])
}

#' Remove records located in forest plantations (MapBiomas class 9)
#' @noRd
.drop_plantation_points <- function(pts, lulc) {
  if (is.null(lulc) || lulc$type != "categorical") return(pts)
  cls <- terra::extract(lulc$raster, as.matrix(pts[, c("decimalLongitude", "decimalLatitude")]))[, 1]
  drop <- !is.na(cls) & cls == 9
  if (any(drop)) .msg("Occurrences", sum(drop), " records in Forest Plantation (MapBiomas class 9) removed.")
  pts[!drop, ]
}
