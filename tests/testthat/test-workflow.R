skip_if_not_installed("randomForest")

test_that("cins_model runs the full workflow (downloads mocked)", {
  p <- make_predictors()
  ext_poly <- function(xmin, xmax, ymin, ymax) terra::as.polygons(terra::ext(xmin, xmax, ymin, ymax), crs = "EPSG:4326")
  lulc_fake <- make_lulc(p)
  local_mocked_bindings(
    get_occurrences = function(species, ...) { d <- make_occurrences(); d$species <- species; d },
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
      z <- rbind(ext_poly(-52, -50, -27, -23), ext_poly(-50, -48, -27, -23)); z$Unit <- c("West", "East"); z
    }
  )
  out <- tempfile("cins_wf")
  m <- suppressMessages(suppressWarnings(cins_model(
    "Testus epiphyticus", life_form = "epiphyte", savanna_habitat = FALSE,
    host_species = "Testus arboreus", proj_region = "Brazil:Parana",
    future_scenario = c("ssp245", "ssp585"), future_year = 2050,
    algorithms = c("glm", "rf"), n_pseudo = 300, replicates = 2, vif_cutoff = NULL,
    area_report = "state", seed = 3, plot_results = FALSE, output_dir = out)))
  expect_s3_class(m, "cins_model")
  expect_true("host_Testus_arboreus" %in% m$selected_vars)
  expect_named(m$future, c("ssp245_2050", "ssp585_2050"))
  expect_true(all(c("Potential_Climatic_km2", "Realized_LULC_km2", "Habitat_Loss_pct") %in% names(m$native$area_table)))
  expect_true(all(m$native$area_table$Realized_LULC_km2 <= m$native$area_table$Potential_Climatic_km2))
  expect_false(is.null(m$projected$realized_binary))
  expect_false(is.null(m$future$ssp585_2050$realized_binary))
  expect_true(file.exists(file.path(out, "Testus_epiphyticus_settings.txt")))
  expect_true(file.exists(file.path(out, "Testus_epiphyticus_occurrences.csv")))
  expect_output(print(m), "cins_model")
  pdf(NULL); expect_silent(plot(m)); dev.off()
})
