# Generate manual rset object from spatiotemporal cross-validation indices

Generate manual rset object from spatiotemporal cross-validation indices

## Usage

``` r
convert_cv_index_rset(cvindex, data, cv_mode)
```

## Arguments

- cvindex:

  integer. Output of
  [`generate_cv_index`](https://sigmafelix.github.io/stcvlite/reference/generate_cv_index.md).

- data:

  data.frame from `stdt`. Should be the same object as what was used for
  `covars` argument of
  [`generate_cv_index`](https://sigmafelix.github.io/stcvlite/reference/generate_cv_index.md)

- cv_mode:

  character(1). Spatiotemporal cross-validation indexing method. See
  `cv_mode` description in
  [`generate_cv_index`](https://sigmafelix.github.io/stcvlite/reference/generate_cv_index.md).

## Value

rset object of `rsample` package. A tibble with a list column of
training-test data.frames and a column of labels.

## Author

Insang Song
