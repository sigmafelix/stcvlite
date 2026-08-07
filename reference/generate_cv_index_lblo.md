# Generate spatio-temporal cross-validation index (leave-block-location-out)

Generate spatio-temporal cross-validation index
(leave-block-location-out)

## Usage

``` r
generate_cv_index_lblo(covars, cv_fold = NULL, blocks = NULL, block_id = NULL)
```

## Arguments

- covars:

  data.frame-like object with `lon`, `lat`, and `time` columns. See
  [`prep_input()`](https://sigmafelix.github.io/stcvlite/reference/prep_input.md)
  for details.

- cv_fold:

  integer(1). Number of folds for cross-validation.

- blocks:

  integer(2)/sf/SpatVector object.

- block_id:

  character(1). The unique identifier of each block.

## Value

An integer vector.

## Author

Insang Song
