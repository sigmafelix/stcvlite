#' Sample full spatiotemporal point data
#'
#' This dataset contains 12,000 spatiotemporal points with
#' `lon`, `lat`, `time`, and `value` columns. The points can be used
#' to demonstrate the functions in `stcvlite`.
#'
#' @family Dataset
#' @format A `data.table` with 12,000 rows (100 locations * 120 times) and
#'   six variables, plus a `"crs"` attribute (`attr(spdat, "crs")`):
#' \describe{
#' \item{time}{Unique point identifier. Arbitrarily generated.}
#' \item{lon}{Longitude}
#' \item{lat}{Latitude}
#' \item{x1}{covariate 1}
#' \item{x2}{covariate 2}
#' \item{y}{response variable}
#' }
#' @note Coordinates are in EPSG:4326
#' @importFrom data.table data.table
#' @examples
#' data("spdat", package = "stcvlite")
"spdat"
