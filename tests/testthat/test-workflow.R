skip_if_not_installed("randomForest")

ext_poly <- function(xmin, xmax, ymin, ymax) terra::as.polygons(terra::ext(xmin, xmax, ymin, ymax), crs = "EPSG:4326")

mock_downloads <- function(env = parent.frame()) {
  p <- make_predictors()
  lulc_fake <- make_lulc(p)
  local_mocked_bindings(
    get_occurrences = function(species, ...) {
      d <- make_occurrences(); d$species <- species
      if (grepl("other", species)) d$decimalLatitude <- d$decimalLatitude + 0.5
      d
    },
    define_calibration_area = function(...) ext_poly(-52, -48, -27, -23),
    get_region_polygon = function(region, ...) if (inherits(region, "SpatVector")) region else ext_poly(-51, -48, -26, -23),
    get_climate = function(region, resolution = 5, vif_cutoff = 10, ...) terra::crop(p, region, mask = TRUE),
    get_future_climate = function(region, selected_vars = NULL, ...) {
      f <- terra::crop(p, region, mask = TRUE); f$bio1 <- f$bio1 + 2
      if (!is.null(selected_vars)) f <- f[[intersect(selected_vars, names(f))]]
      f
    },
    .lulc_setup = function(lulc_source, poly, ...) {
      if (is.null(lulc_source)) return(NULL)
      list(lulc = list(raster = terra::crop(lulc_fake$raster, poly), type = "categorical",
                       source = "mapbiomas", native_res_m = lulc_fake$native_res_m),
           classes = c(3L, 4L), mapbiomas_legend = TRUE)
    },
    .report_bounds = function(level, ext, cache_dir) {
      if (level == "total") return(NULL)
      z <- rbind(ext_poly(-52, -50, -27, -23), ext_poly(-50, -48, -27, -23)); z$Unit <- c("West", "East"); z
    },
    .package = "cinspecies", .env = env
  )
}

make_dem <- function() {
  r <- terra::rast(nrows = 80, ncols = 80, xmin = -52, xmax = -48, ymin = -27, ymax = -23, crs = "EPSG:4326")
  terra::setValues(r, (terra::crds(r, na.rm = FALSE)[, 2] + 27) * 200)
}

test_that("cins_model runs the full SBSR-style workflow (downloads mocked)", {
  mock_downloads()
  out <- tempfile("cins_wf")
  m <- suppressMessages(suppressWarnings(cins_model(
    "Testus epiphyticus", life_form = "epiphyte", savanna_habitat = TRUE,
    host_species = "Testus arboreus", proj_region = "Brazil:Parana",
    future_scenario = c("ssp245", "ssp585"), future_year = 2070,
    algorithms = c("glm", "rf"), n_pseudo = 300, n_pseudo_sets = 2, replicates = 2,
    vif_cutoff = NULL, ensemble_min = 0.3, dem = make_dem(),
    area_report = "state", seed = 3, plot_results = FALSE, output_dir = out)))
  expect_s3_class(m, "cins_model")
  expect_true("host_Testus_arboreus" %in% m$selected_vars)
  expect_named(m$future, c("ssp245_2070", "ssp585_2070"))
  # background sets and ensemble filter
  expect_setequal(unique(m$evaluation$pa_set), 1:2)
  expect_true("in_ensemble" %in% names(m$evaluation))
  # uncertainty, conflict, change and elevation outputs
  expect_s4_class(m$native$uncertainty, "SpatRaster")
  ct <- m$projected$conflict_table
  expect_true(all(c("Available habitat", "Pasture", "Forest plantation") %in% ct$Category))
  expect_equal(sum(ct$Area_pct), 100, tolerance = 0.1)
  expect_true(all(c("Area_current_km2", "Net_change_ssp585_2070_pct", "Gain_ssp245_2070_km2") %in% names(m$area_change)))
  expect_equal(m$area_change$Unit[nrow(m$area_change)], "TOTAL")
  expect_equal(m$elevation$Period, c("current", "ssp245_2070", "ssp585_2070"))
  expect_equal(m$elevation$Shift_m[1], 0)
  expect_true(all(m$native$area_table$Realized_LULC_km2 <= m$native$area_table$Potential_Climatic_km2))
  expect_true(file.exists(file.path(out, "Testus_epiphyticus_elevation_shift.csv")))
  expect_true(file.exists(file.path(out, "Testus_epiphyticus_native_lulc_conflict.csv")))
  expect_output(print(m), "Elevation shift")
  pdf(NULL); expect_silent(plot(m)); dev.off()
})

test_that("user-supplied future climate (custom GCMs) is used without downloads", {
  mock_downloads()
  p <- make_predictors()
  g1 <- p; g1$bio1 <- g1$bio1 + 1
  g2 <- p; g2$bio1 <- g2$bio1 + 3
  m <- suppressMessages(suppressWarnings(cins_model(
    "Testus arboreus", life_form = "tree", savanna_habitat = FALSE,
    future_climate = list(myssp_2080 = list(g1, g2)),
    algorithms = c("glm"), n_pseudo = 300, n_pseudo_sets = 1, replicates = 1, vif_cutoff = NULL,
    ensemble_min = NULL, uncertainty = FALSE, elevation_shift = FALSE, conflict_table = FALSE,
    seed = 1, plot_results = FALSE, output_dir = tempfile("cins_fc"))))
  expect_named(m$future, "myssp_2080")
  expect_false(is.null(m$area_change))
})

test_that("several species are modelled and combined", {
  mock_downloads()
  out <- tempfile("cins_multi")
  mm <- suppressMessages(suppressWarnings(cins_model(
    c("Testus arboreus", "Testus other"), life_form = c("tree", "shrub"), savanna_habitat = c(FALSE, TRUE),
    future_scenario = "ssp585", future_year = 2050,
    algorithms = "glm", n_pseudo = 300, n_pseudo_sets = 1, replicates = 1, vif_cutoff = NULL,
    ensemble_min = NULL, elevation_shift = FALSE, seed = 2, plot_results = FALSE, output_dir = out)))
  expect_s3_class(mm, "cins_multi")
  expect_equal(nrow(mm$summary), 2)
  expect_equal(mm$summary$life_form, c("tree", "shrub"))
  expect_true("realized_ssp585_2050_km2" %in% names(mm$summary))
  rv <- terra::values(mm$richness$native$realized)[, 1]
  expect_true(max(rv, na.rm = TRUE) <= 2 && min(rv, na.rm = TRUE) >= 0)
  expect_true(file.exists(file.path(out, "species_summary.csv")))
  expect_true(dir.exists(file.path(out, "Testus_other")))
  expect_output(print(mm), "2 species")
})
