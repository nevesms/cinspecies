# Geographic regions, administrative boundaries and calibration area (M) ------

#' @noRd
.world_polygons <- function(cache_dir = cins_cache_dir()) {
  w <- geodata::world(resolution = 2, level = 0, path = cache_dir)
  if (is.null(w)) stop("Could not download GADM world boundaries (check your internet connection).", call. = FALSE)
  cc <- geodata::country_codes()
  w$continent <- cc$continent[match(w$GID_0, cc$ISO3)]
  w
}

#' @noRd
.country_iso3 <- function(name) {
  cc <- geodata::country_codes()
  key <- .norm_txt(name)
  hit <- which(.norm_txt(cc$NAME) == key | .norm_txt(cc$ISO3) == key |
                 .norm_txt(cc$NAME_LOCAL) == key | .norm_txt(cc$NAME_ISO) == key |
                 .norm_txt(cc$NAME_FAO) == key)
  if (length(hit) == 0) return(NA_character_)
  cc$ISO3[hit[1]]
}

#' @noRd
.states <- function(iso3, cache_dir = cins_cache_dir()) {
  s <- geodata::gadm(iso3, level = 1, path = cache_dir, resolution = 2)
  if (is.null(s)) stop("Could not download GADM level-1 boundaries for ", paste(iso3, collapse = ", "), call. = FALSE)
  if (is.list(s) && !inherits(s, "SpatVector")) s <- do.call(rbind, s)
  s
}

#' Build a polygon for a named region
#'
#' Accepts continents (e.g. `"South America"`), countries (`"Brazil"`, `"BRA"`),
#' states written as `"Country:State"` (accents optional, e.g. `"Brazil:Parana"`),
#' `"global"`, a `SpatVector`/`sf` polygon, or a path to a vector file
#' (`.shp`, `.gpkg`, `.geojson`). Boundaries come from GADM via `geodata`.
#'
#' @param region Character vector, `SpatVector`, `sf` object or file path.
#' @param cache_dir Directory where boundaries are cached.
#' @return A single dissolved `SpatVector` polygon (or `NULL` for `"global"`).
#' @references
#' GADM (2022). Database of Global Administrative Areas, version 4.1. <https://gadm.org>
#'
#' Hijmans, R.J., Barbosa, M., Ghosh, A. & Mandel, A. (2024). geodata: Download
#' Geographic Data. R package. \doi{10.32614/CRAN.package.geodata}
#' @export
get_region_polygon <- function(region, cache_dir = cins_cache_dir()) {
  if (is.null(region) || identical(region, "global")) return(NULL)
  if (inherits(region, "sf")) region <- terra::vect(region)
  if (inherits(region, "SpatVector")) return(terra::aggregate(terra::project(region, "EPSG:4326")))
  if (is.character(region) && length(region) == 1 && file.exists(region)) {
    return(terra::aggregate(terra::project(terra::vect(region), "EPSG:4326")))
  }
  if (identical(region, "auto")) stop("region = 'auto' needs occurrences; use define_calibration_area().", call. = FALSE)

  pieces <- list()
  world <- NULL
  for (item in region) {
    if (grepl(":", item, fixed = TRUE)) {
      parts <- strsplit(item, ":", fixed = TRUE)[[1]]
      iso <- .country_iso3(parts[1])
      if (is.na(iso)) stop("Country not recognised in '", item, "'.", call. = FALSE)
      st <- .states(iso, cache_dir)
      hit <- st[.norm_txt(st$NAME_1) == .norm_txt(parts[2]), ]
      if (nrow(hit) == 0) {
        stop("State '", parts[2], "' not found in ", parts[1], ". Available: ",
             paste(sort(st$NAME_1), collapse = ", "), call. = FALSE)
      }
      pieces[[item]] <- hit
      next
    }
    if (is.null(world)) world <- .world_polygons(cache_dir)
    cont <- world[.norm_txt(world$continent) == .norm_txt(item), ]
    if (nrow(cont) > 0) { pieces[[item]] <- cont; next }
    iso <- .country_iso3(item)
    if (!is.na(iso)) { pieces[[item]] <- world[world$GID_0 == iso, ]; next }
    stop("Region not recognised: '", item, "'. Use a continent, a country, or 'Country:State' ",
         "(e.g. 'Brazil:Parana').", call. = FALSE)
  }
  pieces <- lapply(pieces, function(p) p[, character(0)])
  terra::aggregate(do.call(rbind, unname(pieces)))
}

#' Define the calibration area (accessible area, M) automatically
#'
#' Builds the accessible area M of the BAM framework (Soberon & Peterson 2005;
#' Barve et al. 2011) from the occurrence records: a minimum convex polygon
#' around the records buffered by `buffer_km`, clipped to land. This avoids
#' calibrating models against environmentally extreme backgrounds that the
#' species never had access to.
#'
#' @param occurrences Data frame with `decimalLongitude` and `decimalLatitude`.
#' @param buffer_km Numeric. Buffer distance (km) around the convex hull.
#' @param cache_dir Directory where boundaries are cached.
#' @return A dissolved `SpatVector` polygon.
#' @references
#' Soberon, J. & Peterson, A.T. (2005). Interpretation of models of fundamental
#' ecological niches and species' distributional areas. Biodiversity
#' Informatics, 2, 1-10. \doi{10.17161/bi.v2i0.4}
#'
#' Barve, N. et al. (2011). The crucial role of the accessible area in
#' ecological niche modeling and species distribution modeling. Ecological
#' Modelling, 222, 1810-1819. \doi{10.1016/j.ecolmodel.2011.02.011}
#' @export
define_calibration_area <- function(occurrences, buffer_km = 300, cache_dir = cins_cache_dir()) {
  pts <- terra::vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  hull <- if (nrow(pts) >= 3) terra::convHull(pts) else terra::aggregate(terra::buffer(pts, 1))
  m <- terra::buffer(hull, width = buffer_km * 1000)
  land <- terra::aggregate(.world_polygons(cache_dir)[, character(0)])
  m <- terra::intersect(terra::aggregate(m), terra::crop(land, terra::ext(m)))
  .msg("Calibration area", "M defined as convex hull + ", buffer_km, " km buffer (",
       round(sum(terra::expanse(m, unit = "km"))), " km2).")
  terra::aggregate(m)
}

#' Share of a polygon that lies inside Brazil (used to choose MapBiomas)
#' @noRd
.share_in_brazil <- function(poly, cache_dir = cins_cache_dir()) {
  if (is.null(poly)) return(0)
  w <- .world_polygons(cache_dir)
  bra <- w[w$GID_0 == "BRA", ]
  inter <- tryCatch(terra::intersect(poly, bra), error = function(e) NULL)
  if (is.null(inter) || nrow(inter) == 0) return(0)
  sum(terra::expanse(inter, unit = "km")) / sum(terra::expanse(poly, unit = "km"))
}
