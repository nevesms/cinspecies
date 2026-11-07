#' Download and Clean GBIF Occurrences with Altitudinal Filter
#'
#' @param species String. Scientific name of the species.
#' @param ecosystem String. "auto", "terrestrial", or "marine".
#' @param limit Numeric. Max records to download from GBIF.
#' @param thin_km Numeric. Minimum distance (km) for spatial thinning.
#' @param alt_min Numeric. Minimum altitude (m). NULL to ignore.
#' @param alt_max Numeric. Maximum altitude (m). NULL to ignore.
#'
#' @export
get_occurrences <- function(species, ecosystem = "auto", limit = 5000, thin_km = 10, alt_min = NULL, alt_max = NULL) {

  message(">>> [Occurrences] Downloading GBIF records for: ", species)
  search_res <- rgbif::occ_data(scientificName = species, hasCoordinate = TRUE, limit = limit)
  if (is.null(search_res$data)) stop("No coordinates found in GBIF for: ", species)

  pts <- unique(stats::na.omit(search_res$data[, c("species", "decimalLongitude", "decimalLatitude")]))

  world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sv")
  pts_vect <- terra::vect(pts, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

  if (ecosystem == "auto") {
    overlap <- terra::is.related(pts_vect, world_map, "intersects")
    ecosystem <- ifelse((sum(overlap) / nrow(pts)) >= 0.5, "terrestrial", "marine")
  }

  message(">>> [Occurrences] Cleaning coordinate errors (", ecosystem, ")...")
  if (ecosystem == "terrestrial") {
    pts_clean <- CoordinateCleaner::clean_coordinates(
      x = pts, lon = "decimalLongitude", lat = "decimalLatitude", species = "species",
      tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"), value = "clean", verbose = FALSE
    )
  } else {
    pts_clean <- pts[!terra::is.related(pts_vect, world_map, "intersects"), ]
  }

  if (nrow(pts_clean) < 5) stop("Insufficient points remaining after coordinate cleaning.")

  message(">>> [Occurrences] Applying spatial thinning (", thin_km, " km)...")
  thin_res <- spThin::thin(
    loc.data = pts_clean, lat.col = "decimalLatitude", long.col = "decimalLongitude",
    spec.col = "species", thin.par = thin_km, reps = 1, locs.thinned.list.return = TRUE,
    write.files = FALSE, write.log.file = FALSE
  )
  pts_final <- thin_res[[1]]
  colnames(pts_final) <- c("decimalLongitude", "decimalLatitude")
  pts_final$species <- species

  if (!is.null(alt_min) || !is.null(alt_max)) {
    message(">>> [Occurrences] Applying DEM altitudinal filter...")
    elev <- geodata::elevation_global(res = 5, path = tempdir())
    pts_vect_final <- terra::vect(pts_final, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    pts_alt <- terra::extract(elev, pts_vect_final)[, 2]

    keep <- !is.na(pts_alt)
    if (!is.null(alt_min)) keep <- keep & (pts_alt >= alt_min)
    if (!is.null(alt_max)) keep <- keep & (pts_alt <= alt_max)

    pts_final <- pts_final[which(keep), ]
  }

  message(">>> [Occurrences] Final ready points for ", species, ": ", nrow(pts_final))
  return(list(data = pts_final, ecosystem = ecosystem))
}
