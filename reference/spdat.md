# Sample full spatiotemporal point data

This dataset contains 12,000 spatiotemporal points with `lon`, `lat`,
`time`, and `value` columns. The points can be used to demonstrate the
functions in `stcvlite`.

## Usage

``` r
spdat
```

## Format

A `data.table` with 12,000 rows (100 locations \* 120 times) and six
variables, plus a `"crs"` attribute (`attr(spdat, "crs")`):

- time:

  Unique point identifier. Arbitrarily generated.

- lon:

  Longitude

- lat:

  Latitude

- x1:

  covariate 1

- x2:

  covariate 2

- y:

  response variable

## Note

Coordinates are in EPSG:4326

## Examples

``` r
data("spdat", package = "stcvlite")
```
