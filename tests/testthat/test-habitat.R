test_that("life-form classes follow the MapBiomas legend and require savanna choice", {
  expect_equal(habitat_classes("tree", "mapbiomas", FALSE), c(3L, 5L, 6L, 49L))
  expect_true(4L %in% habitat_classes("tree", "mapbiomas", TRUE))
  expect_equal(habitat_classes("grassland", "mapbiomas", FALSE), c(11L, 12L, 29L, 50L))
  expect_true(4L %in% habitat_classes("shrub", "mapbiomas", TRUE))
  expect_equal(habitat_classes("epiphyte", "mapbiomas", FALSE), habitat_classes("tree", "mapbiomas", FALSE))
  expect_error(habitat_classes("tree", "mapbiomas", NULL), "savanna_habitat")
  expect_error(habitat_classes("alga", "mapbiomas", FALSE), "life_form")
  expect_equal(habitat_classes("tree", "esa"), c("trees", "mangroves"))
})

test_that("habitat layer is binary at native grid and a fraction when aggregated", {
  p <- make_predictors()
  lu <- make_lulc(p)
  h30 <- habitat_layer(lu, c(3L), resolution_m = lu$native_res_m)
  expect_setequal(stats::na.omit(unique(terra::values(h30)[, 1])), c(0, 1))
  hc <- habitat_layer(lu, c(3L), resolution_m = lu$native_res_m * 5)
  v <- terra::values(hc)[, 1]
  expect_true(all(v >= 0 & v <= 1))
  expect_true(any(v > 0 & v < 1))
  expect_warning(habitat_layer(lu, 3L, resolution_m = lu$native_res_m / 3), "finest resolution")
})

test_that("cins_model validates life form and savanna before downloading", {
  expect_error(cins_model("Araucaria angustifolia"), "life_form")
  expect_error(cins_model("Araucaria angustifolia", life_form = "tree"), "savanna_habitat")
  expect_error(cins_model("x", life_form = "tree", savanna_habitat = FALSE, future_scenario = "ssp585",
                          future_year = 2200), "future_year")
})

test_that("coarser lulc_resolution warns about accuracy", {
  expect_warning(cinspecies:::.validate_lulc_args("mapbiomas", NULL, FALSE, 1000), "coarser")
  expect_silent(cinspecies:::.validate_lulc_args("mapbiomas", NULL, FALSE, 30))
  expect_silent(cinspecies:::.validate_lulc_args("esa", NULL, NULL, 30))
})

test_that("records in forest plantations are removed", {
  p <- make_predictors()
  lu <- make_lulc(p)
  pts <- make_occurrences()
  cls <- terra::extract(lu$raster, as.matrix(pts[, 2:3]))[, 1]
  out <- suppressMessages(cinspecies:::.drop_plantation_points(pts, lu))
  expect_equal(nrow(out), sum(cls != 9))
})

test_that("cells outside the region stay NA (not 0) in the habitat layer", {
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(3L, NA, 9L, 4L), crs = "EPSG:4326")
  lu <- list(raster = r, type = "categorical", source = "mapbiomas", native_res_m = terra::res(r)[1] * 111320)
  expect_equal(terra::values(habitat_layer(lu, c(3L, 4L), lu$native_res_m))[, 1], c(1, NA, 0, 1))
})

test_that("water classes are excluded and conflict table sums to 100%", {
  p <- make_predictors()
  lu <- make_lulc(p)
  v <- terra::values(lu$raster)[, 1]; v[seq(1, length(v), 7)] <- 33L
  lu$raster <- terra::setValues(lu$raster, v)
  h <- habitat_layer(lu, c(3L, 4L), lu$native_res_m, na_classes = c(26L, 31L, 33L))
  expect_true(all(is.na(terra::values(h)[v == 33, 1])))
  pb <- terra::ifel(terra::init(lu$raster, "y") > -25, 1, 0)
  ct <- lulc_conflict_table(pb, lu, c(3L, 4L))
  expect_false("Water" %in% ct$Category)
  expect_equal(sum(ct$Area_pct), 100, tolerance = 0.05)
  expect_equal(ct$Category[1], "Available habitat")
})

test_that("elevation shift and area change compute expected directions", {
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:4326")
  y <- terra::crds(r, na.rm = FALSE)[, 2]
  dem <- terra::setValues(r, y * 1000)
  cur <- terra::setValues(r, as.numeric(y < 1))
  fut <- terra::setValues(r, as.numeric(y > 1))
  e <- suppressMessages(elevation_shift(list(current = cur, f = fut), 0.5, dem = dem))
  expect_gt(e$Shift_m[2], 900)
  ch <- suppressMessages(area_change_table(cur, list(f = fut), "total"))
  expect_equal(ch$Net_change_f_pct, 0, tolerance = 1)
  expect_gt(ch$Gain_f_km2, 0)
  expect_equal(ch$Gain_f_km2, ch$Loss_f_km2, tolerance = 1)
})
