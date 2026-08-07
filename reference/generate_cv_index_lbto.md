# Generate spatio-temporal cross-validation index (leave-block-time-out)

Generate spatio-temporal cross-validation index (leave-block-time-out)

## Usage

``` r
generate_cv_index_lbto(covars, cv_fold = NULL)
```

## Arguments

- covars:

  data.frame-like object with `lon`, `lat`, and `time` columns. See
  [`prep_input()`](https://sigmafelix.github.io/stcvlite/reference/prep_input.md)
  for details.

- cv_fold:

  integer(1). Number of folds for cross-validation.

## Value

An integer vector.

## Author

Insang Song
