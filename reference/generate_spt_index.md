# Generate unique spatiotemporal identifier from stdt object

It generates unique spatiotemporal identifier in the input stdt object.
Regardless of mode values (should be one of `"spatial"` or
`"spatiotemporal"`), the index column will be named "sp_index".

## Usage

``` r
generate_spt_index(covars, mode = c("spatial", "spatiotemporal"))
```

## Arguments

- covars:

  stdt. See
  [`prep_input()`](https://sigmafelix.github.io/stcvlite/reference/prep_input.md)

- mode:

  One of `"spatial"` or `"spatiotemporal"`

## Value

stdt with a new index named "sp_index"

## Author

Insang Song
