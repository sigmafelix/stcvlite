# Generate spatio-temporal cross-validation index (random)

Generate spatio-temporal cross-validation index (random)

## Usage

``` r
generate_cv_index_random(covars, cv_fold = NULL)
```

## Arguments

- covars:

  stdt. See
  [`prep_input()`](https://sigmafelix.github.io/stcvlite/reference/prep_input.md)
  for details.

- cv_fold:

  integer(1). Number of folds for cross-validation.

## Value

An integer vector with unique values of `seq(1, cv_fold)`

## Author

Insang Song
