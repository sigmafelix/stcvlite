make_covars_fixture <- function() {
  dt <- data.table::data.table(
    lon = c(0, 0, 1, 1, 2, 2),
    lat = c(0, 0, 1, 1, 2, 2),
    time = as.Date("2020-01-01") + c(0, 1, 0, 1, 0, 1),
    value = 1:6
  )

  structure(
    list(stdt = dt, crs = "EPSG:4326"),
    class = "stdt"
  )
}
