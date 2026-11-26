# Occurrence records for terrestrial plants ----------------------------------

#' @noRd
.find_col <- function(df, candidates) {
  hit <- names(df)[tolower(names(df)) %in% tolower(candidates)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

#' Read and standardise user-supplied occurrences
#' @noRd
.read_user_occurrences <- function(occurrences, species) {
  if (is.character(occurrences)) occurrences <- utils::read.csv(occurrences, stringsAsFactors = FALSE)
  if (inherits(occurrences, "SpatVector")) {
    xy <- terra::crds(occurrences)
    occurrences <- data.frame(decimalLongitude = xy[, 1], decimalLatitude = xy[, 2])
  }
  occurrences <- as.data.frame(occurrences)
  lon <- .find_col(occurrences, c("decimalLongitude", "longitude", "lon", "long", "x"))
  lat <- .find_col(occurrences, c("decimalLatitude", "latitude", "lat", "y"))
  if (is.na(lon) || is.na(lat)) stop("Could not find longitude/latitude columns in `occurrences`.", call. = FALSE)
  out <- data.frame(species = species,
                    decimalLongitude = as.numeric(occurrences[[lon]]),
                    decimalLatitude = as.numeric(occurrences[[lat]]))
  out[stats::complete.cases(out), ]
}

#' Download GBIF records of a terrestrial plant and apply plant-specific filters
#' @noRd
.download_gbif_plants <- function(species, limit, year_min, max_uncertainty_m,
                                  exclude_cultivated, native_only) {
  bb <- rgbif::name_backbone(name = species)
  if (is.null(bb$usageKey) || is.na(bb$usageKey[1])) stop("Taxon not found in the GBIF backbone: ", species, call. = FALSE)
  if (!identical(as.character(bb$kingdom[1]), "Plantae")) {
    stop("cinspecies models terrestrial plants only; GBIF classifies '", species, "' in kingdom ",
         bb$kingdom[1], ".", call. = FALSE)
  }
  key <- if (!is.null(bb$acceptedUsageKey) && !is.na(bb$acceptedUsageKey[1])) bb$acceptedUsageKey[1] else bb$usageKey[1]
  if (!identical(toupper(as.character(bb$status[1])), "ACCEPTED")) {
    .msg("Occurrences", "'", species, "' is a ", tolower(bb$status[1]),
         " in GBIF; records of the accepted name (taxonKey ", key, ") are used.")
  }

  .msg("Occurrences", "Downloading GBIF records (taxonKey ", key, ")...")
  res <- rgbif::occ_data(taxonKey = key, hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                         occurrenceStatus = "PRESENT", limit = limit)
  d <- res$data
  if (is.null(d) || nrow(d) == 0) stop("No georeferenced GBIF records for: ", species, call. = FALSE)
  if (!is.null(res$meta$count) && res$meta$count > limit) {
    .msg("Occurrences", "GBIF holds ", res$meta$count, " records; only the first ", limit,
         " were retrieved (increase `gbif_limit`).")
  }
  n0 <- nrow(d)
  col <- function(nm) if (nm %in% names(d)) as.character(d[[nm]]) else rep(NA_character_, nrow(d))

  keep <- rep(TRUE, nrow(d))
  bor <- toupper(col("basisOfRecord"))
  keep <- keep & !(bor %in% "FOSSIL_SPECIMEN")
  if (exclude_cultivated) {
    est <- toupper(col("establishmentMeans"))
    deg <- toupper(col("degreeOfEstablishment"))
    keep <- keep & !(bor %in% "LIVING_SPECIMEN") & !(est %in% "MANAGED") &
      !grepl("CULTIVATED|CAPTIVE|MANAGED", deg)
  }
  if (native_only) {
    est <- toupper(col("establishmentMeans"))
    keep <- keep & !(est %in% c("INTRODUCED", "INVASIVE", "NATURALISED", "VAGRANT"))
  }
  if (!is.null(year_min)) {
    yr <- suppressWarnings(as.integer(col("year")))
    keep <- keep & (is.na(yr) | yr >= year_min)
  }
  if (!is.null(max_uncertainty_m)) {
    unc <- suppressWarnings(as.numeric(col("coordinateUncertaintyInMeters")))
    keep <- keep & (is.na(unc) | unc <= max_uncertainty_m)
  }
  d <- d[keep, ]
  .msg("Occurrences", n0 - nrow(d), " records removed by plant filters (fossil",
       if (exclude_cultivated) ", cultivated/living collections" else "",
       if (native_only) ", non-native" else "",
       if (!is.null(year_min)) paste0(", before ", year_min) else "",
       if (!is.null(max_uncertainty_m)) paste0(", uncertainty > ", max_uncertainty_m, " m") else "", ").")
  out <- data.frame(species = species,
                    decimalLongitude = as.numeric(d$decimalLongitude),
                    decimalLatitude = as.numeric(d$decimalLatitude))
  unique(out[stats::complete.cases(out), ])
}

#' Download and clean occurrence records of a terrestrial plant
#'
#' Retrieves GBIF records (via `rgbif`) for vascular plants and other
#' terrestrial plants (trees, shrubs, grassland species, lianas, epiphytes), or
#' uses records supplied by the user, and applies:
#' (1) taxonomic check (kingdom Plantae; synonyms resolved to the accepted name),
#' (2) removal of fossil records and, by default, cultivated / living-collection
#' and managed records (important for planted forest species),
#' (3) temporal and positional-uncertainty filters,
#' (4) automated coordinate cleaning with `CoordinateCleaner`,
#' (5) optional altitudinal filter on an SRTM-derived DEM (WorldClim 2.1, 30 arc-s),
#' (6) spatial thinning with `spThin`.
#'
#' @param species Character. Scientific name.
#' @param occurrences Optional data frame, `SpatVector` or CSV path with the
#'   user's own records (longitude/latitude columns are detected automatically).
#'   When supplied, GBIF is not queried.
#' @param clean Logical. Apply coordinate cleaning and thinning to user records.
#' @param gbif_limit Numeric. Maximum number of GBIF records (max 100000).
#' @param year_min Numeric or `NULL`. Oldest record year kept (default 1970,
#'   matching the WorldClim 2.1 baseline 1970-2000).
#' @param max_uncertainty_m Numeric or `NULL`. Maximum coordinate uncertainty (m).
#' @param exclude_cultivated Logical. Remove living specimens, managed and
#'   cultivated records.
#' @param native_only Logical. Also remove records flagged as introduced,
#'   invasive or naturalised.
#' @param cc_tests Character. CoordinateCleaner tests (see
#'   `CoordinateCleaner::clean_coordinates()`), e.g. add `"outliers"` or `"urban"`.
#' @param thin_km Numeric. Minimum distance (km) between records.
#' @param alt_min,alt_max Numeric or `NULL`. Elevation limits (m a.s.l.).
#' @param seed Integer or `NULL`. Random seed for the thinning step.
#' @param cache_dir Directory where the DEM is cached.
#' @return Data frame with `species`, `decimalLongitude`, `decimalLatitude`.
#' @references
#' Chamberlain, S. et al. rgbif: Interface to the Global Biodiversity
#' Information Facility API. R package. \doi{10.32614/CRAN.package.rgbif}
#'
#' Zizka, A. et al. (2019). CoordinateCleaner: Standardized cleaning of
#' occurrence records from biological collection databases. Methods in Ecology
#' and Evolution, 10, 744-751. \doi{10.1111/2041-210X.13152}
#'
#' Aiello-Lammens, M.E. et al. (2015). spThin: an R package for spatial thinning
#' of species occurrence records for use in ecological niche models. Ecography,
#' 38, 541-545. \doi{10.1111/ecog.01132}
#'
#' Farr, T.G. et al. (2007). The Shuttle Radar Topography Mission. Reviews of
#' Geophysics, 45, RG2004. \doi{10.1029/2005RG000183}
#' @export
get_occurrences <- function(species, occurrences = NULL, clean = TRUE, gbif_limit = 20000,
                            year_min = 1970, max_uncertainty_m = 10000,
                            exclude_cultivated = TRUE, native_only = FALSE,
                            cc_tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"),
                            thin_km = 5, alt_min = NULL, alt_max = NULL, seed = NULL,
                            cache_dir = cins_cache_dir()) {

  if (!is.null(occurrences)) {
    .msg("Occurrences", "Using user-supplied records for ", species)
    pts <- .read_user_occurrences(occurrences, species)
  } else {
    pts <- .download_gbif_plants(species, gbif_limit, year_min, max_uncertainty_m,
                                 exclude_cultivated, native_only)
    clean <- TRUE
  }
  pts <- unique(pts)

  if (clean) {
    .msg("Occurrences", "Cleaning coordinates with CoordinateCleaner (", nrow(pts), " records)...")
    pts <- CoordinateCleaner::clean_coordinates(
      x = pts, lon = "decimalLongitude", lat = "decimalLatitude", species = "species",
      tests = cc_tests,
      value = "clean", verbose = FALSE
    )
  }
  if (nrow(pts) < 5) stop("Fewer than 5 records remain after cleaning for ", species, ".", call. = FALSE)

  if (!is.null(alt_min) || !is.null(alt_max)) {
    .msg("Occurrences", "Applying altitudinal filter (DEM 30 arc-s)...")
    elev <- geodata::elevation_global(res = 0.5, path = cache_dir)
    z <- terra::extract(elev, as.matrix(pts[, c("decimalLongitude", "decimalLatitude")]))[, 1]
    keep <- !is.na(z)
    if (!is.null(alt_min)) keep <- keep & z >= alt_min
    if (!is.null(alt_max)) keep <- keep & z <= alt_max
    pts <- pts[keep, ]
  }

  if (clean && !is.null(thin_km) && thin_km > 0 && nrow(pts) > 1) {
    .msg("Occurrences", "Spatial thinning (", thin_km, " km)...")
    if (!is.null(seed)) set.seed(seed)
    th <- spThin::thin(loc.data = pts, lat.col = "decimalLatitude", long.col = "decimalLongitude",
                       spec.col = "species", thin.par = thin_km, reps = 1,
                       locs.thinned.list.return = TRUE, write.files = FALSE,
                       write.log.file = FALSE, verbose = FALSE)[[1]]
    pts <- data.frame(species = species, decimalLongitude = th$Longitude, decimalLatitude = th$Latitude)
  }
  if (nrow(pts) < 5) stop("Fewer than 5 records remain after filtering for ", species, ".", call. = FALSE)
  rownames(pts) <- NULL
  .msg("Occurrences", "Final records for ", species, ": ", nrow(pts))
  pts
}
