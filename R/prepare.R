# data preparation for generate_index

#' Convert sf and sftime objects to stdt format
#' @param input An sf, sftime, SpatVector or SpatRasterDataset object
#' @return A data.table object with `lon`, `lat`, `time` and a `crs` attribute
#' @export
prep_input <- function(input) {
  detect_time_col <- function(nms) {
    nms_low <- tolower(nms)
    preferred <- c("time", "datetime", "date", "timestamp", "ts")
    matched <- nms[nms_low %in% preferred]

    if (length(matched) == 0L) {
      stop("input does not contain a temporal column")
    }

    matched_low <- tolower(matched)
    for (cand in preferred) {
      idx <- which(matched_low == cand)
      if (length(idx) > 0L) {
        return(matched[[idx[1]]])
      }
    }

    return(matched[[1]])
  }

  melt_raster_var <- function(rast, value_name) {
    df <- as.data.frame(rast, xy = TRUE)
    dt <- data.table::as.data.table(df)
    data.table::setnames(dt, old = names(dt)[1:2], new = c("lon", "lat"))
    measure_cols <- names(dt)[3:ncol(dt)]

    if (length(measure_cols) == 0L) {
      stop("input SpatRasterDataset does not contain temporal layers")
    }

    dt_long <- data.table::melt(
      dt,
      id.vars = c("lon", "lat"),
      measure.vars = measure_cols,
      variable.name = ".time_src",
      value.name = value_name,
      variable.factor = FALSE
    )

    ts_vec <- terra::time(rast)
    if (!is.null(ts_vec) && length(ts_vec) == length(measure_cols)) {
      time_map <- data.table::data.table(
        .time_src = measure_cols,
        time = as.character(ts_vec)
      )
    } else {
      time_map <- data.table::data.table(
        .time_src = measure_cols,
        time = measure_cols
      )
    }

    dt_long <- data.table::merge.data.table(
      dt_long,
      time_map,
      by = ".time_src",
      all.x = TRUE,
      sort = FALSE
    )

    dt_long[[".time_src"]] <- NULL
    return(dt_long)
  }

  if (inherits(input, "sf") || inherits(input, "sftime")) {
    crs_dt <- as.character(sf::st_crs(input))[1]
    stdt <- data.table::as.data.table(sf::st_drop_geometry(input))
    time_col <- detect_time_col(names(stdt))

    rep_points <- sf::st_point_on_surface(sf::st_geometry(input))
    coords <- sf::st_coordinates(rep_points)
    if (nrow(coords) != nrow(stdt)) {
      coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(input)))
    }

    stdt[["lon"]] <- coords[, 1]
    stdt[["lat"]] <- coords[, 2]
    data.table::setnames(stdt, old = time_col, new = "time")
  } else if (inherits(input, "SpatVector")) {
    crs_dt <- terra::crs(input)
    stdt <- data.table::as.data.table(as.data.frame(input, geom = FALSE))
    time_col <- detect_time_col(names(stdt))

    rep_points <- terra::centroids(input)
    coords <- as.data.frame(rep_points, geom = "XY")
    stdt[["lon"]] <- coords$x
    stdt[["lat"]] <- coords$y
    data.table::setnames(stdt, old = time_col, new = "time")
  } else if (inherits(input, "SpatRasterDataset")) {
    crs_dt <- terra::crs(input)
    vars <- names(input)

    if (length(vars) == 0L) {
      stop("input SpatRasterDataset has no variables")
    }

    stdt <- melt_raster_var(input[1], vars[1])
    if (length(vars) > 1L) {
      for (var in vars[2:length(vars)]) {
        dt_var <- melt_raster_var(input[var], var)
        stdt <- data.table::merge.data.table(
          stdt,
          dt_var,
          by = c("lon", "lat", "time"),
          all = TRUE,
          sort = FALSE
        )
      }
    }
  } else {
    stop("stobj class not accepted")
  }

  data.table::setcolorder(
    stdt,
    c("lon", "lat", "time", setdiff(names(stdt), c("lon", "lat", "time")))
  )
  attr(stdt, "crs") <- crs_dt
  return(stdt)
}


from_stdt <- function(input, to = "sf") {
  match.arg(to, choices = c("sf", "sftime", "SpatVector"))

  if (to == "sf") {
    sf::st_as_sf(
      input,
      coords = c("lon", "lat"),
      crs = attr(input, "crs"),
      remove = FALSE
    )
  } else if (to == "sftime") {
    sf::st_as_sf(
      input,
      coords = c("lon", "lat"),
      crs = attr(input, "crs"),
      remove = FALSE
    )
  } else if (to == "SpatVector") {
    terra::vect(
      input,
      geom = c("lon", "lat"),
      crs = attr(input, "crs")
    )
  }
}