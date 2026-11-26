skip_if_not_installed("randomForest")

test_that("ensemble, maxTSS threshold, LULC and areas work end-to-end on synthetic data", {
  p <- make_predictors()
  pts <- make_occurrences()
  out <- withr_tmp <- tempfile("cins_test"); dir.create(out)
  fit <- suppressMessages(suppressWarnings(
    fit_models(pts, p, algorithms = c("glm", "rf"), n_pseudo = 300, replicates = 2, seed = 7,
               output_dir = out, species_name = "Testus")))
  expect_s4_class(fit$model, "sdmModels")
  expect_true(is.numeric(fit$threshold) && fit$threshold > 0 && fit$threshold < 1)
  expect_equal(attr(fit$threshold, "method"), "maxTSS")
  expect_true(file.exists(file.path(out, "Testus_native_suitability.tif")))
  expect_true(file.exists(file.path(out, "Testus_evaluation.csv")))
  expect_true(all(c("AUC", "TSS") %in% names(fit$evaluation)))

  # user-defined numeric threshold overrides maxTSS (bug in 0.1.0)
  th <- model_threshold(fit$model, 0.5)
  expect_equal(as.numeric(th), 0.5)
  b <- post_process_raster(fit$ensemble_map, th)
  expect_setequal(stats::na.omit(unique(terra::values(b)[, 1])), c(0, 1))

  lu <- make_lulc(p)
  hab <- habitat_layer(lu, c(3L), resolution_m = lu$native_res_m)
  rl <- suppressMessages(apply_lulc_mask(fit$ensemble_map, fit$threshold, hab, file_prefix = file.path(out, "t")))
  expect_equal(terra::res(rl$continuous), terra::res(hab))
  vr <- terra::values(rl$binary)[, 1]; vp <- terra::values(rl$potential_binary)[, 1]
  expect_true(all(vr <= vp, na.rm = TRUE))

  a_pot <- cinspecies:::.area_by_zone(rl$potential_binary, NULL)$Area_km2
  a_real <- cinspecies:::.area_by_zone(rl$weight, NULL)$Area_km2
  expect_true(a_real <= a_pot)
  exact <- sum(terra::global(terra::cellSize(rl$potential_binary, unit = "km") * rl$potential_binary,
                             "sum", na.rm = TRUE)[, 1])
  expect_lt(abs(a_pot - exact) / exact, 0.01)

  zones <- terra::as.polygons(terra::ext(-52, -50, -27, -23), crs = "EPSG:4326")
  zones <- rbind(zones, terra::as.polygons(terra::ext(-50, -48, -27, -23), crs = "EPSG:4326"))
  zones$Unit <- c("West", "East")
  z <- cinspecies:::.area_by_zone(rl$potential_binary, zones)
  expect_lt(abs(sum(z$Area_km2) - a_pot) / a_pot, 0.01)
})
