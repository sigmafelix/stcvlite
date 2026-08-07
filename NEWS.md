# 0.3.0 (unreleased)
- **Breaking:** dropped the `stdt` S3 wrapper class (`list(stdt = <table>,
  crs = <chr>)`). All functions (`generate_cv_index()`,
  `generate_spt_index()`, `generate_block_sp_index()`, `plot_cv_folds()`,
  `convert_cv_index_rset()`) now take a plain `data.frame`-like object
  (e.g. the output of `prep_input()`) with `lon`, `lat`, and `time`
  columns directly, and always return a fold vector matching that
  object's row count. `spdat` is now a `data.table` with a `"crs"`
  attribute instead of a wrapped list.
- `generate_block_sp_index()` now clusters/joins on *unique* `lon`/`lat`
  locations and broadcasts the result back onto all rows, instead of
  operating on every row. This keeps spatial blocking correct for
  irregular/unbalanced spacetime panels (e.g. stations observed a
  different number of times) and fixes the sf/SpatVector polygon-block
  path, which previously received the whole wrapped `stdt` object instead
  of a table `sf::st_as_sf()` could consume.

# 0.2
- Added multiple backend engines for spatial blocking
- Added plot function
- Added preparation function

# 0.1
- Subset the indexing functions from amadeus