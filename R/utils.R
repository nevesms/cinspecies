# Internal helpers -----------------------------------------------------------

#' Persistent cache directory for downloaded datasets
#'
#' All downloads (WorldClim, CMIP6, SoilGrids, ESA WorldCover, MapBiomas, GADM
#' boundaries and elevation) are stored in a persistent cache so that they are
#' only downloaded once. The location can be changed with
#' `options(cinspecies.cache = "path")` or with the `CINSPECIES_CACHE`
#' environment variable.
#'
#' @param create Logical. Create the directory if it does not exist.
#' @return Character path of the cache directory.
#' @examples
#' cins_cache_dir(create = FALSE)
#' @export
cins_cache_dir <- function(create = TRUE) {
  path <- getOption("cinspecies.cache", Sys.getenv("CINSPECIES_CACHE", ""))
  if (!nzchar(path)) path <- tools::R_user_dir("cinspecies", which = "cache")
  if (create && !dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

#' @noRd
.msg <- function(stage, ...) {
  message(">>> [", stage, "] ", ...)
}

#' @noRd
.norm_txt <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  tolower(trimws(gsub("[^A-Za-z0-9 ]", "", x)))
}

#' @noRd
.sp_slug <- function(species) gsub("[^A-Za-z0-9]+", "_", trimws(species))

#' Write a raster with compression (keeps returned objects file-backed)
#' @noRd
.write_rast <- function(r, file, datatype = "FLT4S") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(r, filename = file, overwrite = TRUE, datatype = datatype,
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"))
}

#' Standardise WorldClim/CMIP6 bioclimatic layer names to bio1 ... bio19
#' @noRd
.std_bio_names <- function(r) {
  nm <- names(r)
  num <- suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", nm)))
  if (terra::nlyr(r) == 19 && (any(is.na(num)) || anyDuplicated(num))) num <- 1:19
  names(r) <- paste0("bio", num)
  r
}

#' Approximate conversion of metres to decimal degrees (equator)
#' @noRd
.m_to_deg <- function(m) m / 111320

#' Valid CMIP6 GCMs distributed by WorldClim 2.1 / geodata
#' @noRd
.valid_gcms <- c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CanESM5",
                 "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1",
                 "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G",
                 "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR",
                 "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")

#' @noRd
.default_gcms <- c("ACCESS-CM2", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")

#' @noRd
.resolve_gcms <- function(gcms) {
  if (identical(gcms, "default")) return(.default_gcms)
  if (identical(gcms, "all")) return(.valid_gcms)
  bad <- setdiff(gcms, .valid_gcms)
  if (length(bad) > 0) {
    stop("Invalid GCM(s): ", paste(bad, collapse = ", "), ". Valid options: ",
         paste(.valid_gcms, collapse = ", "), call. = FALSE)
  }
  gcms
}

#' @noRd
.check_resolution <- function(resolution) {
  if (!as.character(resolution) %in% c("10", "5", "2.5", "0.5")) {
    stop("`resolution` must be one of 10, 5, 2.5 or 0.5 arc-minutes.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.future_period <- function(future_year) {
  if (!is.numeric(future_year) || future_year < 2021 || future_year > 2100) {
    stop("`future_year` must be between 2021 and 2100.", call. = FALSE)
  }
  if (future_year <= 2040) "2021-2040" else if (future_year <= 2060) "2041-2060"
  else if (future_year <= 2080) "2061-2080" else "2081-2100"
}

#' @noRd
.check_suggested <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required ", why, ". Install it with install.packages('", pkg, "').",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Map an algorithm code to the package it needs (for early, readable errors)
#' @noRd
.algorithm_pkgs <- c(rf = "randomForest", maxnet = "maxnet", maxent = "dismo", brt = "gbm",
                     gbm = "gbm", svm = "e1071", mlp = "RSNNS", nnet = "nnet", mars = "earth",
                     cart = "rpart", glm = "stats", gam = "mgcv", mahal = "stats", bioclim = "stats")

#' @noRd
.check_algorithms <- function(algorithms) {
  unknown <- setdiff(algorithms, names(.algorithm_pkgs))
  if (length(unknown) > 0) stop("Unknown algorithm(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  pk <- unique(.algorithm_pkgs[algorithms])
  miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop("Algorithm package(s) not installed: ", paste(miss, collapse = ", "),
         ". Install them (or run sdm::installAll()).", call. = FALSE)
  }
  invisible(TRUE)
}
