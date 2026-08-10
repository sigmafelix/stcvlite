# Spatiotemporal CV on gstat's Air Quality Data with xgboost and BART

## Introduction

This vignette works through a complete `stcvlite` workflow on a real
spatiotemporal dataset: `DE_RB_2005`, the daily rural-background PM10
readings from 69 German monitoring stations in 2005, shipped with
**gstat** (it is the same series used throughout gstat’s spatio-temporal
kriging tutorials, originally distributed as `air`/`dates`/`stations` in
**spacetime**). We will

1.  reshape the station-day records into the flat `lon`/`lat`/`time`
    table `stcvlite` expects, using the station coordinates and the
    daily time stamps,
2.  build spatiotemporal indices with
    [`generate_spt_index()`](https://sigmafelix.github.io/stcvlite/reference/generate_spt_index.md)
    and inspect the spatial blocks that
    [`generate_block_sp_index()`](https://sigmafelix.github.io/stcvlite/reference/generate_block_sp_index.md)
    derives from the station coordinates,
3.  generate every cross-validation scheme `stcvlite` offers and see how
    fold structure differs across them, and
4.  fit two very different prediction models – **xgboost** (gradient
    boosted trees) and **BART** (Bayesian Additive Regression Trees, via
    **dbarts**) – on the `rset` objects `stcvlite` produces, and compare
    how each cross-validation scheme rates their accuracy.

The point of the last step is not to declare a “winner” model, but to
show *why* the scheme matters: naive random CV lets a model interpolate
between neighboring stations and nearby days, which inflates apparent
accuracy. Leaving out whole stations or whole blocks of time is a
stricter, more realistic test of how these models would perform at a new
location or in a future period.

``` r

library(stcvlite)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(rsample)
library(xgboost)
library(dbarts)

set.seed(2026)
```

## 1. The data: gstat’s rural-background PM10 stations

`DE_RB_2005` is a
[`spacetime::STSDF`](https://rdrr.io/pkg/spacetime/man/STSDF-class.html)
– a sparse space-time object holding one `PM10` measurement per
station-day, plus per-station metadata such as elevation.

``` r

data("DE_RB_2005", package = "gstat")
```

We flatten it to a long table with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), drop
missing readings, and keep the station identifier, its coordinates
(originally UTM zone 32N), the date, station altitude (a genuine
predictor of background PM10), and the response.

``` r

air_raw <- as.data.frame(DE_RB_2005)
setDT(air_raw)
air_raw <- air_raw[!is.na(PM10)]
air_raw <- air_raw[, .(
  station_id = sp.ID,
  coords.x1, coords.x2,
  time = as.Date(time),
  station_altitude,
  PM10
)]

cat(sprintf(
  "%d station-day observations across %d stations\n",
  nrow(air_raw), uniqueN(air_raw$station_id)
))
#> 23230 station-day observations across 69 stations
```

### Station coordinates

The station coordinates live once per station in `DE_RB_2005@sp`; after
flattening they are repeated on every row. We reproject them from UTM32N
to longitude/latitude (EPSG:4326), matching the convention `stcvlite`
uses for its own sample data (`spdat`), and join them back onto the long
table.

``` r

stations_sf <- sf::st_as_sf(
  unique(air_raw[, .(station_id, coords.x1, coords.x2)]),
  coords = c("coords.x1", "coords.x2"), crs = 32632
)
stations_sf <- sf::st_transform(stations_sf, 4326)
station_coords <- data.table(
  station_id = stations_sf$station_id,
  lon = sf::st_coordinates(stations_sf)[, 1],
  lat = sf::st_coordinates(stations_sf)[, 2]
)

air_dt <- merge(air_raw, station_coords, by = "station_id")
```

We add a pair of day-of-year harmonics so the prediction models have a
smooth seasonal signal to work with, alongside the station coordinates
and altitude.

``` r

air_dt[, doy := as.integer(format(time, "%j"))]
air_dt[, doy_sin := sin(2 * pi * doy / 365)]
air_dt[, doy_cos := cos(2 * pi * doy / 365)]
air_dt[, c("coords.x1", "coords.x2", "doy") := NULL]
setcolorder(air_dt, c("station_id", "lon", "lat", "time"))

str(air_dt)
#> Classes 'data.table' and 'data.frame':   23230 obs. of  8 variables:
#>  $ station_id      : chr  "DEBB053" "DEBB053" "DEBB053" "DEBB053" ...
#>  $ lon             : num  14 14 14 14 14 ...
#>  $ lat             : num  52.6 52.6 52.6 52.6 52.6 ...
#>  $ time            : Date, format: "2005-01-01" "2005-01-02" ...
#>  $ station_altitude: int  88 88 88 88 88 88 88 88 88 88 ...
#>  $ PM10            : num  27.17 10.21 6.75 12.56 13.67 ...
#>  $ doy_sin         : num  0.0172 0.0344 0.0516 0.0688 0.086 ...
#>  $ doy_cos         : num  1 0.999 0.999 0.998 0.996 ...
#>  - attr(*, ".internal.selfref")=<pointer: 0x55c01b5ccf20> 
#>  - attr(*, "sorted")= chr "station_id"
```

## 2. Preparing the covariate table

`stcvlite`’s index generators expect a plain `data.frame`-like object
with `lon`, `lat`, and `time` columns (plus whatever covariates you
like) – no special wrapper class.
[`prep_input()`](https://sigmafelix.github.io/stcvlite/reference/prep_input.md)
does the coordinate/time extraction for `sf`, `sftime`, `SpatVector`,
and `SpatRasterDataset` inputs and returns exactly that shape, so we
hand it an `sf` version of `air_dt` and use its output as-is.

``` r

air_sf <- sf::st_as_sf(air_dt, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
air_prepped <- suppressWarnings(stcvlite::prep_input(air_sf))

class(air_prepped)
#> [1] "data.table" "data.frame"
air_prepped
#>             lon      lat       time station_id station_altitude   PM10
#>           <num>    <num>     <Date>     <char>            <int>  <num>
#>     1: 14.01525 52.56383 2005-01-01    DEBB053               88 27.167
#>     2: 14.01525 52.56383 2005-01-02    DEBB053               88 10.208
#>     3: 14.01525 52.56383 2005-01-03    DEBB053               88  6.750
#>     4: 14.01525 52.56383 2005-01-04    DEBB053               88 12.565
#>     5: 14.01525 52.56383 2005-01-05    DEBB053               88 13.667
#>    ---                                                                
#> 23226: 13.64492 52.97184 2005-10-20    DEUB040               70 25.083
#> 23227: 13.64492 52.97184 2005-10-21    DEUB040               70 23.417
#> 23228: 13.64492 52.97184 2005-10-22    DEUB040               70 20.955
#> 23229: 13.64492 52.97184 2005-10-23    DEUB040               70  5.708
#> 23230: 13.64492 52.97184 2005-10-24    DEUB040               70  8.375
#>            doy_sin   doy_cos
#>              <num>     <num>
#>     1:  0.01721336 0.9998518
#>     2:  0.03442161 0.9994074
#>     3:  0.05161967 0.9986668
#>     4:  0.06880243 0.9976303
#>     5:  0.08596480 0.9962982
#>    ---                      
#> 23226: -0.94559639 0.3253421
#> 23227: -0.93985606 0.3415708
#> 23228: -0.93383723 0.3576982
#> 23229: -0.92754168 0.3737197
#> 23230: -0.92097129 0.3896304
```

([`prep_input()`](https://sigmafelix.github.io/stcvlite/reference/prep_input.md)
calls
[`sf::st_point_on_surface()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
internally, which warns for point geometries even though the coordinates
it returns are exact – safe to ignore here.)

## 3. Spatiotemporal indices and station-coordinate blocks

[`generate_spt_index()`](https://sigmafelix.github.io/stcvlite/reference/generate_spt_index.md)
turns the station coordinates (and, optionally, time) into a single
grouping key. In `"spatial"` mode it collapses to one key per station;
in `"spatiotemporal"` mode every station-day gets its own key.

``` r

air_spatial <- generate_spt_index(air_prepped, mode = "spatial")
air_spt <- generate_spt_index(air_prepped, mode = "spatiotemporal")

cat(sprintf(
  "unique spatial keys: %d (= number of stations)\n",
  uniqueN(air_spatial$sp_index)
))
#> unique spatial keys: 69 (= number of stations)
cat(sprintf(
  "unique spatiotemporal keys: %d (= number of rows)\n",
  uniqueN(air_spt$sp_index)
))
#> unique spatiotemporal keys: 23230 (= number of rows)
```

[`generate_block_sp_index()`](https://sigmafelix.github.io/stcvlite/reference/generate_block_sp_index.md)
uses the same `lon`/`lat` columns to group stations into contiguous
spatial blocks (DBSCAN by default), clustering each *unique* station
location once regardless of how many days it was observed. This is the
mechanism behind the `"lblo"` and `"lblto"` schemes below – six blocks
of stations, chosen purely from their coordinates:

``` r

air_blocks <- generate_block_sp_index(air_prepped, cv_fold = 6)
table(air_blocks$sp_index)
#> 
#>    1    2    3    4    5    6 
#> 3380 2873 5443 3615 5194 2725
```

``` r

plot_cv_folds(air_blocks, air_blocks$sp_index, max_points = 5000)
```

## 4. Every cross-validation scheme, side by side

[`generate_cv_index()`](https://sigmafelix.github.io/stcvlite/reference/generate_cv_index.md)
implements seven schemes (see the package README for the full
description of each acronym). We build all of them here to compare the
fold structure they induce on the same 23,230-row dataset:
leave-one-station-out (`lolo`) and leave-one-day-out (`loto`) each hold
out one coordinate at a time; `lolto` is plain leave-one-row-out;
`lblo`/`lbto`/`lblto` group stations and/or days into contiguous blocks;
and `random` is the naive baseline.

``` r

cv_indices <- list(
  lolo   = generate_cv_index(air_prepped, cv_mode = "lolo"),
  loto   = generate_cv_index(air_prepped, cv_mode = "loto"),
  lolto  = generate_cv_index(air_prepped, cv_mode = "lolto"),
  lblo   = generate_cv_index(air_prepped, cv_mode = "lblo", cv_fold = 6),
  lbto   = generate_cv_index(air_prepped, cv_mode = "lbto", cv_fold = 6),
  lblto  = generate_cv_index(air_prepped, cv_mode = "lblto", sp_fold = 4, t_fold = 3),
  random = generate_cv_index(air_prepped, cv_mode = "random", cv_fold = 5)
)

fold_summary <- rbindlist(lapply(names(cv_indices), function(m) {
  tb <- table(cv_indices[[m]])
  data.table(
    cv_mode = m,
    n_folds = length(tb),
    min_fold_size = min(tb),
    max_fold_size = max(tb)
  )
}))
knitr::kable(fold_summary)
```

| cv_mode | n_folds | min_fold_size | max_fold_size |
|:--------|--------:|--------------:|--------------:|
| lolo    |      69 |            79 |           365 |
| loto    |     365 |            54 |            68 |
| lolto   |   23230 |             1 |             1 |
| lblo    |       6 |          2725 |          5443 |
| lbto    |       6 |          3569 |          3973 |
| lblto   |      12 |           832 |          3957 |
| random  |       5 |          4543 |          4734 |

`lolo` and `lblo`/`lbto`/`lblto` hold out *whole stations or blocks*, so
every fold tests extrapolation in space and/or time. `loto` and `lolto`
produce far more, much smaller folds (365 and 23,230 respectively) –
faithful to the “leave-one-out” definition, but too many to refit two
models on for this vignette. We generate their indices to show the fold
structure, and evaluate models on the other five schemes below, using
the same helper code that would apply identically to `loto`/`lolto`
given enough time.

## 5. Prediction models: xgboost and BART

Both models use the same feature set – station coordinates, altitude,
and the seasonal harmonics – and are refit independently on the analysis
set of every fold, then scored on that fold’s assessment set.

``` r

feature_cols <- c("lon", "lat", "station_altitude", "doy_sin", "doy_cos")
target_col <- "PM10"

fit_predict_xgb <- function(train, test) {
  dtrain <- xgboost::xgb.DMatrix(
    data = as.matrix(train[, ..feature_cols]),
    label = train[[target_col]]
  )
  fit <- xgboost::xgb.train(
    data = dtrain,
    nrounds = 150,
    params = list(
      objective = "reg:squarederror",
      max_depth = 4, eta = 0.1, subsample = 0.8
    ),
    verbose = 0
  )
  predict(fit, as.matrix(test[, ..feature_cols]))
}

fit_predict_bart <- function(train, test) {
  fit <- dbarts::bart(
    x.train = as.matrix(train[, ..feature_cols]),
    y.train = train[[target_col]],
    x.test = as.matrix(test[, ..feature_cols]),
    ntree = 50, ndpost = 300, nskip = 100,
    verbose = FALSE, keeptrees = FALSE
  )
  fit$yhat.test.mean
}
```

[`convert_cv_index_rset()`](https://sigmafelix.github.io/stcvlite/reference/convert_cv_index_rset.md)
turns a fold-assignment vector into a `tidymodels`-compatible `rset`;
[`rsample::analysis()`](https://rsample.tidymodels.org/reference/as.data.frame.rsplit.html)
and
[`rsample::assessment()`](https://rsample.tidymodels.org/reference/as.data.frame.rsplit.html)
pull the train/test tables out of each split. `evaluate_rset()` loops
over every fold of a given `rset`, fits both models, and records
RMSE/MAE:

``` r

rmse <- function(actual, pred) sqrt(mean((actual - pred)^2))
mae <- function(actual, pred) mean(abs(actual - pred))

evaluate_rset <- function(rset, cv_mode) {
  per_fold <- lapply(seq_len(nrow(rset)), function(i) {
    split <- rset$splits[[i]]
    train <- as.data.table(rsample::analysis(split))
    test <- as.data.table(rsample::assessment(split))

    pred_xgb <- fit_predict_xgb(train, test)
    pred_bart <- fit_predict_bart(train, test)

    data.table(
      cv_mode = cv_mode,
      fold = rset$id[i],
      n_test = nrow(test),
      engine = c("xgboost", "bart"),
      rmse = c(rmse(test[[target_col]], pred_xgb), rmse(test[[target_col]], pred_bart)),
      mae = c(mae(test[[target_col]], pred_xgb), mae(test[[target_col]], pred_bart))
    )
  })
  rbindlist(per_fold)
}
```

We run this on five schemes: `random` as the (optimistic) baseline,
`lolo` for pure spatial extrapolation to unseen stations, `lblo` for
block-wise spatial extrapolation, `lbto` for pure temporal extrapolation
to unseen time blocks, and `lblto` for combined spatiotemporal block
extrapolation. `lolo` alone means refitting both models 69 times; on a
laptop this section takes roughly two to three minutes.

``` r

rset_random <- convert_cv_index_rset(cv_indices$random, air_prepped, cv_mode = "random")
rset_lolo   <- convert_cv_index_rset(cv_indices$lolo, air_prepped, cv_mode = "lolo")
rset_lblo   <- convert_cv_index_rset(cv_indices$lblo, air_prepped, cv_mode = "lblo")
rset_lbto   <- convert_cv_index_rset(cv_indices$lbto, air_prepped, cv_mode = "lbto")
rset_lblto  <- convert_cv_index_rset(cv_indices$lblto, air_prepped, cv_mode = "lblto")

results <- rbindlist(list(
  evaluate_rset(rset_random, "random"),
  evaluate_rset(rset_lolo, "lolo"),
  evaluate_rset(rset_lblo, "lblo"),
  evaluate_rset(rset_lbto, "lbto"),
  evaluate_rset(rset_lblto, "lblto")
))
```

## 6. Comparing the schemes

Averaging RMSE within each scheme/engine combination shows how much the
CV scheme itself changes the apparent accuracy:

``` r

summary_tbl <- results[, .(
  n_folds = uniqueN(fold),
  mean_rmse = mean(rmse),
  mean_mae = mean(mae)
), by = .(cv_mode, engine)]

knitr::kable(summary_tbl[order(cv_mode, engine)], digits = 2)
```

| cv_mode | engine  | n_folds | mean_rmse | mean_mae |
|:--------|:--------|--------:|----------:|---------:|
| lblo    | bart    |       6 |      9.04 |     6.65 |
| lblo    | xgboost |       6 |      8.17 |     5.73 |
| lblto   | bart    |      12 |      8.28 |     6.03 |
| lblto   | xgboost |      12 |      7.70 |     5.57 |
| lbto    | bart    |       6 |     17.19 |    13.29 |
| lbto    | xgboost |       6 |     11.36 |     8.44 |
| lolo    | bart    |      69 |      7.66 |     5.68 |
| lolo    | xgboost |      69 |      7.08 |     5.16 |
| random  | bart    |       5 |      6.98 |     4.83 |
| random  | xgboost |       5 |      6.77 |     4.66 |

``` r

plotly::plot_ly(
  data = results,
  x = ~cv_mode, y = ~rmse, color = ~engine,
  type = "box"
) |>
  plotly::layout(
    yaxis = list(title = "Fold RMSE (PM10, ug/m3)"),
    xaxis = list(title = "CV scheme"),
    boxmode = "group"
  )
```

## Takeaways

- `random` CV lets both models lean on neighboring stations and adjacent
  days, so it reads as the easiest scheme – and it is: mean RMSE for
  both engines sits lowest here.
- `lolo` and `lblo` remove entire stations from training. RMSE rises
  only modestly over `random`, which makes sense given the feature set:
  longitude, latitude, and altitude are smooth surfaces, so a model can
  extrapolate them reasonably well to a station it has never seen.
- `lbto` removes contiguous *blocks of days* instead, and RMSE roughly
  doubles relative to `random` for both engines. The seasonal harmonics
  (`doy_sin`/`doy_cos`) repeat every year, so once a block of
  consecutive days is entirely held out, neither model has anything in
  the feature set that encodes the actual meteorology of that stretch of
  time – only its position in the yearly cycle. `lblto` combines both
  stresses and lands in between.
- xgboost and BART track each other closely within every scheme; the gap
  between schemes is driven far more by *what* is held out than by
  *which* model is used. That gap – not the model choice – is what
  `stcvlite`’s CV schemes are built to expose, and it is invisible if
  you only ever run random k-fold.

The `rset` objects used above (`rset_random`, `rset_lolo`, …) are
ordinary `rsample` objects, so they also drop directly into
`parsnip`/`workflows`/`tune` if you prefer that interface – see the
README’s “Integration with tidymodels” section for the `fit_resamples()`
form of the same workflow.
