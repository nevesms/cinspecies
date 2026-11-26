# cinspecies 0.2.0

## Focus: terrestrial plants (forest and natural-grassland species)

* New main function `cins_model()`; `cins_pipeline()` is deprecated and forwards to it.
* New required argument `life_form` (`"tree"`, `"shrub"`, `"grassland"`, `"liana"`, `"epiphyte"`) that selects the native-habitat land-cover classes automatically (`habitat_classes()`).
* New argument `savanna_habitat` (no default by design, error if unset with MapBiomas): whether Savanna Formation counts as habitat.
* Marine modelling (Bio-ORACLE/MARSPEC, `ecosystem`) removed.
* GBIF records: taxonomic check (Plantae, accepted name), removal of fossil, living-collection, cultivated and managed records, `year_min`, `max_uncertainty_m`, `native_only`, CoordinateCleaner `seas` test, removal of records in forest plantations (MapBiomas class 9). User records via `occurrences`.
* `dependent_species` renamed `host_species` (phorophytes/hosts); host suitability is weighted by host habitat.

## Land use and land cover

* Realized maps at 30 m by default (`lulc_resolution`); coarser values aggregate habitat into fractions and trigger an accuracy warning.
* Automatic MapBiomas Collection 10 (window read from the cloud GeoTIFF and cached) or ESA WorldCover (`lulc_source = "auto"`).
* Area reports computed from fraction-weighted habitat, each cell assigned to one administrative unit (GADM 4.1); memory-safe for 30 m grids.

## Automation

* `region = "auto"`: calibration area (M) from the records (`define_calibration_area()`).
* Optional SoilGrids and topography predictors (`get_static_predictors()`).
* Several SSPs and years in one call; GCM list validated against WorldClim 2.1.
* Persistent download cache (`cins_cache_dir()`), `seed` for reproducibility, settings and package versions saved to `*_settings.txt`, evaluation and variable-importance tables exported.

## Bug fixes

* Numeric `binary_threshold` was ignored whenever sdm returned a threshold; it now always overrides. The default is maxTSS computed on the ensemble (previously the mean of individual-model thresholds).
* `"Country:State"` regions (e.g. `"Brazil:Parana"`) now work, accents optional.
* `area_report = "state"` no longer requires the undeclared `rnaturalearthhires`.
* Future layers at 0.5 arc-min had wrong names; layer names are now standardised (`bio1`...`bio19`).
* `terra::tmpFiles(remove = TRUE)` could delete files backing the returned rasters; outputs are now written to `output_dir`.
* `fine_resolution` units were inconsistent between code and docs (replaced by `lulc_resolution` in metres).
* `sdm` must be attached to work, so it moved to `Depends`; `require()` calls inside functions removed.
* README and DESCRIPTION brought in line with the code (license GPL-3, defaults, GCM list).

## References

* References added to every function (`@references`), to the README, to `inst/CITATION` and through `cins_references()`.
