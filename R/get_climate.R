#' Helper to Crop Spatial Polygon by Country, Continent, or State
#' @noRd
get_region_polygon <- function(region_names) {
  if (identical(region_names, "global")) return(NULL)

  world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sv")
  poly <- world_map[world_map$name %in% region_names | world_map$continent %in% region_names, ]

  if (nrow(poly) == 0) {
    states_map <- rnaturalearth::ne_states(returnclass = "sv")
    poly <- states_map[states_map$name %in% region_names, ]
  }
  if (nrow(poly) == 0) stop("Region or State not found: ", paste(region_names, collapse = ", "))
  return(poly)
}

#' Download Current Climate and Apply VIF Filter
#'
#' @export
get_climate <- function(ecosystem = "terrestrial", region = "global", resolution = 5, vif_cutoff = 10) {
  message(">>> [Climate] Loading baseline environmental layers...")

  if (ecosystem == "terrestrial") {
    clim <- geodata::worldclim_global(var = "bio", res = resolution, path = tempdir())
  } else {
    marine_layers <- c("BO_sstmean", "BO_salinity", "BO_bathymean")
    clim <- terra::rast(sdmpredictors::load_layers(marine_layers, datadir = tempdir()))
  }

  poly <- get_region_polygon(region)
  if (!is.null(poly)) {
    message(">>> [Climate] Cropping to training region: ", paste(region, collapse = ", "))
    clim <- terra::crop(clim, poly, mask = TRUE)
  }

  if (!is.null(vif_cutoff)) {
    message(">>> [Climate] Calculating VIF (cutoff = ", vif_cutoff, ")...")
    df_sample <- terra::spatSample(clim, size = 5000, method = "random", na.rm = TRUE, as.df = TRUE)
    vif_res <- usdm::vifstep(df_sample, th = vif_cutoff)
    selected_vars <- vif_res@results$Variables
    clim <- clim[[selected_vars]]
    message(">>> [Climate] Retained non-collinear variables: ", paste(selected_vars, collapse = ", "))
  }

  return(clim)
}

#' Download CMIP6 Future Climate and Build Multi-GCM Consensus
#'
#' @export
get_future_climate <- function(region = "global", resolution = 5, future_scenario = "ssp585", future_year = 2050, gcms = "all", selected_vars = NULL) {

  valid_gcms <- c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2", "CNRM-CM6-1", "GISS-E2-1-G", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
  if (identical(gcms, "all")) {
    gcms_to_use <- c("MIROC6", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
  } else {
    gcms_to_use <- intersect(gcms, valid_gcms)
    if (length(gcms_to_use) == 0) stop("Invalid GCM specified. Choose from CMIP6 standard models.")
  }

  period <- ifelse(future_year <= 2040, "2021-2040",
                   ifelse(future_year <= 2060, "2041-2060",
                          ifelse(future_year <= 2080, "2061-2080", "2081-2100")))

  ssp_code <- gsub("ssp", "", tolower(future_scenario))
  poly <- get_region_polygon(region)

  gcm_list <- list()
  for (g in gcms_to_use) {
    message(">>> [Future Climate] Downloading GCM: ", g, " (SSP", ssp_code, " | ", period, ")...")
    r <- geodata::cmip6_world(model = g, ssp = ssp_code, time = period, var = "bioc", res = resolution, path = tempdir())
    # Padroniza nomes das camadas com o WorldClim atual
    names(r) <- paste0("wc2.1_", resolution, "m_bio_", 1:19)

    if (!is.null(poly)) r <- terra::crop(r, poly, mask = TRUE)
    if (!is.null(selected_vars)) r <- r[[intersect(selected_vars, names(r))]]
    gcm_list[[g]] <- r
  }

  if (length(gcm_list) == 1) return(gcm_list[[1]])

  message(">>> [Future Climate] Averaging multi-GCM ensemble stack across ", length(gcm_list), " models...")
  vars <- names(gcm_list[[1]])
  avg_layers <- lapply(vars, function(v) {
    s <- terra::rast(lapply(gcm_list, function(x) x[[v]]))
    terra::app(s, fun = mean, na.rm = TRUE)
  })
  avg_stack <- terra::rast(avg_layers)
  names(avg_stack) <- vars
  return(avg_stack)
}
