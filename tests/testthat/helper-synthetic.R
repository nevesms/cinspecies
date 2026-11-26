# Synthetic landscape: two climatic gradients over a small area in southern Brazil
make_predictors <- function() {
  r <- terra::rast(nrows = 40, ncols = 40, xmin = -52, xmax = -48, ymin = -27, ymax = -23, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  b1 <- terra::setValues(r, 25 - (xy[, 2] + 27) * 2)           # temperature-like (lat)
  b12 <- terra::setValues(r, 1200 + (xy[, 1] + 52) * 150)      # precipitation-like (lon)
  b4 <- terra::setValues(r, 300 + sin(xy[, 1]) * 50 + stats::runif(terra::ncell(r), 0, 5))
  s <- c(b1, b12, b4); names(s) <- c("bio1", "bio12", "bio4"); s
}
make_occurrences <- function(n = 60) {
  set.seed(1)
  data.frame(species = "Testus arboreus",
             decimalLongitude = stats::runif(n, -50.5, -48.5),
             decimalLatitude = stats::runif(n, -26.5, -24.5))
}
make_lulc <- function(template) {
  f <- terra::disagg(template[[1]], fact = 10)   # 10x finer "MapBiomas"
  set.seed(2)
  v <- sample(c(3L, 4L, 9L, 12L, 15L), terra::ncell(f), replace = TRUE)
  f <- terra::setValues(f, v)
  list(raster = f, type = "categorical", source = "mapbiomas", native_res_m = terra::res(f)[1] * 111320)
}
