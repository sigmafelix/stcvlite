# Generate spatio-temporal cross-validation index (leave-block-location-time-out)

Generate spatio-temporal cross-validation index
(leave-block-location-time-out)

## Usage

``` r
generate_cv_index_lblto(covars, sp_fold, t_fold, blocks, block_id = NULL)
```

## Arguments

- covars:

  stdt. See
  [`prep_input`](https://sigmafelix.github.io/stcvlite/reference/prep_input.md)
  for details.

- sp_fold:

  integer(1). Number of subfolds for spatial blocks.

- t_fold:

  integer(1). Number of subfolds for temporal blocks.

- blocks:

  integer(2)/sf/SpatVector object.

- block_id:

  character(1). The unique identifier of each block.

## Value

An integer vector.

## Details

The maximum of results is `sp_fold * t_fold`.

## Author

Insang Song
