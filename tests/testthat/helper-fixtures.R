make_covars_fixture <- function() {
  dt <- data.table::data.table(
    lon = c(0, 0, 1, 1, 2, 2),
    lat = c(0, 0, 1, 1, 2, 2),
    time = as.Date("2020-01-01") + c(0, 1, 0, 1, 0, 1),
    value = 1:6
  )

  attr(dt, "crs") <- "EPSG:4326"
  dt
}


make_covars_fixture_scaled <- function(
    n_locations = 100L,
    n_times = 120L,
    seed = 2026L) {
  set.seed(seed)

  loc_dt <- data.table::data.table(
    loc_id = seq_len(n_locations),
    lon = runif(n_locations, min = -120, max = -70),
    lat = runif(n_locations, min = 30, max = 50),
    x1 = rnorm(n_locations),
    x2 = runif(n_locations, min = -1, max = 1)
  )

  time_vec <- as.Date("2021-01-01") + seq_len(n_times) - 1L
  time_dt <- data.table::data.table(time = time_vec)

  dt <- data.table::CJ(loc_id = loc_dt$loc_id, time = time_dt$time)
  dt <- merge(dt, loc_dt, by = "loc_id", all.x = TRUE, sort = FALSE)

  time_index <- as.integer(dt$time - min(dt$time))
  noise <- rnorm(nrow(dt), sd = 0.2)
  dt[, y := 0.03 * lon - 0.02 * lat + 0.6 * x1 - 0.4 * x2 + 0.15 * time_index + noise]
  dt[, loc_id := NULL]

  data.table::setcolorder(dt, c("lon", "lat", "time", "x1", "x2", "y"))

  attr(dt, "crs") <- "EPSG:4326"
  dt
}

make_blocks_fixture <- function(input) {
  lon_range <- range(input$lon)
  lat_range <- range(input$lat)

  lon_breaks <- seq(lon_range[1], lon_range[2], length.out = 4)
  lat_breaks <- seq(lat_range[1], lat_range[2], length.out = 4)

  lonlat_grid <- expand.grid(
    lon_block = lon_breaks,
    lat_block = lat_breaks
  )
  lonlat_grid_sf <-
    lonlat_grid |>
    sf::st_as_sf(coords = c("lon_block", "lat_block"), crs = sf::st_crs(4326)) |>
    sf::st_make_grid(n = c(4, 4)) |>
    sf::st_as_sf()
  lonlat_grid_sf[["block_id"]] <- seq_len(nrow(lonlat_grid_sf))
  lonlat_grid_sf
}

# blockgrid <- make_blocks_fixture(make_covars_fixture_scaled()$stdt)
