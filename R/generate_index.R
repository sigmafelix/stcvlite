#' Generate spatio-temporal cross-validation index
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()] for details.
#' @param cv_mode character(1). One of `"lolo"` (leave-one-location-out),
#' `"loto"` (leave-one-time-out),
#' `"lolto"` (leave-one-location-time-out),
#' `"lblo"` (leave-block-location-out),
#' `"lbto"` (leave-block-time-out),
#' `"loo"`, `"lblto"` (leave-block-location-time-out)
#' `"random"` (full random selection)
#' @param cv_fold integer(1). Number of folds for cross-validation.
#' @param sp_fold integer(1). Number of subfolds for spatial blocks.
#' @param t_fold integer(1). Number of subfolds for temporal blocks.
#' @param blocks integer(2)/sf/SpatVector object.
#' @param block_id character(1). The unique identifier of each block.
#' @details \code{blocks} is NULL as default
#' (i.e., no block structure is considered); then "lb*" cv_mode
#' is not working.
#' \code{blocks} should be one of integer vector (length 2),
#' sf, sftime, or SpatVector object.
#' Please be advised that we cannot provide error messages
#' related to the \code{blocks} type invalidity
#' due to the cyclometric complexity.
#' If any valid object is assigned here,
#' cv_fold will not make sense and be ignored.
#' Definition of unit grid (horizontal, vertical) size,
#' arbitrary shape of blocks (sf/SpatVector case)
#' where users can consider common spatial hierarchies such as
#' states, counties, hydrological unit code areas, etc.
#' \code{lblto} mode will accept \code{sp_fold} and \code{t_fold}.
#' The maximum of results in this case is \code{sp_fold * t_fold}.
#' @returns A numeric vector with length of the number of input rows
#' @author Insang Song
#' @examples
#' data(spdat)
#'
#' lolo <- generate_cv_index(spdat, cv_mode = "lolo")
#' loto <- generate_cv_index(spdat, cv_mode = "loto")
#' lolto <- generate_cv_index(spdat, cv_mode = "lolto")
#' lblo <- generate_cv_index(spdat, cv_mode = "lblo", cv_fold = 5)
#' lbto <- generate_cv_index(spdat, cv_mode = "lbto", cv_fold = 5)
#' lblto <- generate_cv_index(spdat, cv_mode = "lblto", sp_fold = 5, t_fold = 4)
#' random <- generate_cv_index(spdat, cv_mode = "random", cv_fold = 5)
#' plot_cv_folds(spdat, lblto)
#' @export
generate_cv_index <- function(
    covars,
    cv_mode = c("lolo", "loto", "random", "lblo", "lbto", "lblto", "lolto", "loo"),
    cv_fold = 5L,
    sp_fold = NULL,
    t_fold = NULL,
    blocks = NULL,
    block_id = NULL) {
  # type check
  if (!is.data.frame(covars) || !all(c("lon", "lat", "time") %in% names(covars))) {
    stop("`covars` must be a data.frame-like object with 'lon', 'lat',
    and 'time' columns. Please consider using prep_input().")
  }
  cv_mode <- match.arg(cv_mode)
  # no block check
  # sensible cv_mode and cv_fold
  if (startsWith(cv_mode, "lb") && is.null(cv_fold)) {
    stop("Inputs for blocks argument are invalid.
    Please revisit the help page of this function
    for the details of proper setting of blocks.\n")
  }

  if (cv_mode != "lblto") {
    if ((!is.null(sp_fold) || !is.null(t_fold))) {
      stop("sp_fold and t_fold values are only applicable to
        cv_mode == 'lblto'")
    }
  }

  index_cv <- switch(cv_mode,
    lolo = generate_cv_index_lolo(covars),
    loto = generate_cv_index_loto(covars),
    lolto = generate_cv_index_lolto(covars),
    loo = generate_cv_index_lolto(covars),
    lblo = generate_cv_index_lblo(covars, cv_fold, blocks, block_id),
    lbto = generate_cv_index_lbto(covars, cv_fold),
    lblto = generate_cv_index_lblto(covars,
      sp_fold = sp_fold, t_fold = t_fold,
      blocks = blocks, block_id = block_id
    ),
    random = generate_cv_index_random(covars, cv_fold)
  )

  return(index_cv)
}


#' Generate unique spatiotemporal identifier from stdt object
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()]
#' @param mode One of `"spatial"` or `"spatiotemporal"`
#' @description It generates unique spatiotemporal identifier in the
#' input stdt object. Regardless of mode values
#' (should be one of `"spatial"` or `"spatiotemporal"`),
#' the index column will be named "sp_index".
#' @returns stdt with a new index named "sp_index"
#' @author Insang Song
#' @export
generate_spt_index <- function(
    covars,
    mode = c("spatial", "spatiotemporal")) {
  mode <- match.arg(mode)
  # generate unique sp_index
  if ("sp_index" %in% colnames(covars)) {
    return(covars)
  }
  covars[["sp_index"]] <-
    paste0(covars[["lon"]], "_", covars[["lat"]])

  if (mode == "spatiotemporal") {
    covars[["sp_index"]] <-
      paste0(covars[["sp_index"]], "_", covars[["time"]])
  }

  return(covars)
}


#' Broadcast a per-location `sp_index` onto every row of `covars`
#' @param covars data.frame-like object with `lon`/`lat` columns (one row
#' per observation; locations may repeat).
#' @param locs data.frame-like object with unique `lon`/`lat` pairs and a
#' `sp_index` column, as produced by clustering/joining on deduplicated
#' locations.
#' @details Matching is by exact `lon`/`lat` equality, which is correct
#' whenever coordinates are literal repeated copies (e.g. joined once per
#' site and broadcast across rows), keeping the assignment invariant to how
#' many times each location was observed.
#' @return `covars` with a `sp_index` column added, in the original row order.
#' @noRd
.join_sp_index_by_location <- function(covars, locs) {
  key_covars <- paste(covars[["lon"]], covars[["lat"]], sep = "_")
  key_locs <- paste(locs[["lon"]], locs[["lat"]], sep = "_")
  covars[["sp_index"]] <- locs[["sp_index"]][match(key_covars, key_locs)]
  covars
}


# nolint start
#' Generate blocked spatial index
#' @param covars data.frame-like object with `lon`, `lat`, and `time` columns.
#' @param cv_fold integer(1). Number of folds for cross-validation.
#' @param blocks numeric(2)/sf/SpatVector configuration of blocks.
#' @param block_id character(1). The unique identifier of each block.
#' @details "Block" in this function refers to a group of contiguous spatial/temporal entities.
#' Regardless of `mode` value (should be one of `"spatial"` or `"spatiotemporal"`),
#' the index column will be named "sp_index". \code{block_id} should be set with a proper field name
#' when blocks is a sf or SpatVector object. When \code{cv_fold} is an integer, then the *unique*
#' `lon`/`lat` coordinates in \code{covars} are clustered (not every row), so the result is
#' unaffected by how many times each location was observed — this keeps the block
#' assignment correct for irregular/unbalanced spacetime panels. Location matching when
#' broadcasting the cluster/block id back onto \code{covars} is by exact `lon`/`lat` equality;
#' coordinates that jitter slightly across repeated observations of "the same" site should be
#' rounded/snapped beforehand. When \code{cv_fold} is an integer, the result object will include
#' attributes (accessible with \code{attr} function) about "dbscan_eps" and "dbscan_minPts".
#' @return data.frame-like object (same class as \code{covars}) with a `sp_index` column
#' @author Insang Song
#' @examples
#' data(spdat)
#' covars <- spdat
#' covars_block <- generate_block_sp_index(covars, cv_fold = 5)
#' covars_block
#'
#' plot_cv_folds(covars_block, covars_block$sp_index)
#' @importFrom methods is
#' @importFrom sf st_join
#' @importFrom terra intersect
#' @importFrom data.table .SD
#' @importFrom data.table copy
#' @export
# nolint end
generate_block_sp_index <- function(
    covars,
    cv_fold = NULL,
    blocks = NULL,
    block_id = NULL) {
  detected_class <- class(blocks)[1]

  if (!is.null(cv_fold)) {
    locs <- unique(data.table::data.table(lon = covars[["lon"]], lat = covars[["lat"]]))
    locs <- .assign_spindex(
      locs,
      sp_cols = c("lon", "lat"),
      nclusters = cv_fold,
      engine = "dbscan"
    )
    covars <- .join_sp_index_by_location(covars, locs)
    if (!is.null(attr(locs, "dbscan_eps"))) {
      attr(covars, "dbscan_eps") <- attr(locs, "dbscan_eps")
      attr(covars, "dbscan_minPts") <- attr(locs, "dbscan_minPts")
    }
  }

  if (inherits(blocks, "sf") ||
        inherits(blocks, "sftime") ||
        inherits(blocks, "SpatVector")) {
    if (is.null(block_id)) {
      stop("block_id must be set for this type of argument blocks.")
    }
    if (any(duplicated(unlist(blocks[[block_id]])))) {
      stop("block_id has duplicates.
      Please make sure 'blocks' argument consists of unique identifiers.")
    }

    # spatial join
    fun_stjoin <- switch(detected_class,
      sf = sf::st_join,
      sftime = sf::st_join,
      SpatVector = terra::intersect
    )

    locs <- unique(data.table::data.table(lon = covars[["lon"]], lat = covars[["lat"]]))
    attr(locs, "crs") <- attr(covars, "crs")
    locs_recov <- from_stdt(locs, to = detected_class)
    locs_recov_id <- fun_stjoin(locs_recov, blocks[, block_id])
    locs[["sp_index"]] <- unlist(locs_recov_id[[block_id]])
    covars <- .join_sp_index_by_location(covars, locs)
  }

  if (is.numeric(blocks)) {
    step_lon <- blocks[1]
    step_lat <- blocks[2]

    locs <- unique(data.table::data.table(lon = covars[["lon"]], lat = covars[["lat"]]))

    vlon <- locs[["lon"]]
    vlat <- locs[["lat"]]
    vlon_cuts <- seq(min(vlon), max(vlon), step_lon)
    vlat_cuts <- seq(min(vlat), max(vlat), step_lat)

    x_range <-
      cut(vlon, vlon_cuts, include.lowest = TRUE)
    y_range <-
      cut(vlat, vlat_cuts, include.lowest = TRUE)
    xy_range <- paste(
      as.character(x_range),
      as.character(y_range),
      sep = "|"
    )
    locs[["sp_index"]] <- as.numeric(factor(xy_range))
    covars <- .join_sp_index_by_location(covars, locs)
  }

  return(covars)
}



#' Generate spatio-temporal cross-validation index (leave-one-time-out)
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()] for details.
#' @author Insang Song
#' @return An integer vector.
#' @export
generate_cv_index_loto <-
  function(
      covars) {
    origin_ts <- covars$time
    sorted_ts <- sort(unique(origin_ts))
    cv_index <- as.numeric(factor(origin_ts, levels = sorted_ts))
    return(cv_index)
  }

#' Generate spatio-temporal cross-validation index (leave-one-location-out)
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()] for details.
#' @author Insang Song
#' @return An integer vector.
#' @export
generate_cv_index_lolo <-
  function(
    covars
  ) {
    covars_sp_index <- generate_spt_index(covars, mode = "spatial")
    sp_index_origin <- unlist(covars_sp_index[["sp_index"]])
    sp_index_unique <- sort(unique(sp_index_origin))
    cv_index <- as.numeric(factor(sp_index_origin, levels = sp_index_unique))
    return(cv_index)
  }


#' Generate spatio-temporal cross-validation index leave-one-location-time-out)
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()] for details.
#' @author Insang Song
#' @return An integer vector.
#' @export
generate_cv_index_lolto <-
  function(
      covars
    ) {
    rows <- nrow(covars)
    cv_index <- seq(1, rows)
    return(cv_index)
  }


# nolint start
#' Generate spatio-temporal cross-validation index (leave-block-location-out)
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()] for details.
#' @param cv_fold integer(1). Number of folds for cross-validation.
#' @param blocks integer(2)/sf/SpatVector object.
#' @param block_id character(1). The unique identifier of each block.
#' @author Insang Song
#' @return An integer vector.
#' @export
# nolint end
generate_cv_index_lblo <-
  function(
      covars,
      cv_fold = NULL,
      blocks = NULL,
      block_id = NULL) {
    if (is.null(cv_fold) && is.null(blocks)) {
      stop("Argument cv_fold cannot be NULL unless
      valid argument for blocks is entered.
      Please set a proper number.")
    }
    covars_sp_index <- generate_block_sp_index(
      covars,
      cv_fold = cv_fold, blocks, block_id
    )

    cv_index <- covars_sp_index$sp_index

    # if cv_index is character (when block vector is entered)
    # convert cv_index to factor
    if (is.character(cv_index)) {
      cv_index <- factor(cv_index)
    }
    return(cv_index)
  }

# nolint start
#' Generate spatio-temporal cross-validation index (leave-block-time-out)
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()] for details.
#' @param cv_fold integer(1). Number of folds for cross-validation.
#' @author Insang Song
#' @return An integer vector.
#' @export
# nolint end
generate_cv_index_lbto <- function(
    covars,
    cv_fold = NULL) {
  if (is.null(cv_fold)) {
    stop("Argument cv_fold cannot be NULL. Please set a proper number.\n")
  }

  origin_ts <- covars$time
  origin_ts_min <- min(origin_ts)
  origin_ts_diff <- as.integer(origin_ts - origin_ts_min)

  # We assume that we use daily data...
  # Since POSIX* is based on seconds, we divide
  # 86400 (1 day) from diff
  if (startsWith(class(origin_ts_min)[1], "POSIX")) {
    origin_ts_diff <- origin_ts_diff / 86400
  }
  origin_ts_diff <- origin_ts_diff + 1
  sorted_ts <- sort(unique(origin_ts))
  length_ts <- length(sorted_ts)

  if (length_ts < cv_fold) {
    stop(sprintf("The total length of time series in the input is
    shorter than cv_fold.\n
    Length of the input time series: %d\n
    cv_fold: %d\n", length_ts, cv_fold))
  }
  unit_split <- ceiling(length_ts / cv_fold)
  cv_index <- ceiling(origin_ts_diff / unit_split)
  return(cv_index)
}

# nolint start
#' Generate spatio-temporal cross-validation index (leave-block-location-time-out)
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See \code{\link{prep_input}} for details.
#' @param sp_fold integer(1). Number of subfolds for spatial blocks.
#' @param t_fold integer(1). Number of subfolds for temporal blocks.
#' @param blocks integer(2)/sf/SpatVector object.
#' @param block_id character(1). The unique identifier of each block.
#' @author Insang Song
#' @return An integer vector.
#' @details The maximum of results is \code{sp_fold * t_fold}.
#' @export
# nolint end
generate_cv_index_lblto <- function(
    covars,
    sp_fold,
    t_fold,
    blocks,
    block_id = NULL) {
  covars_sp_index <- generate_block_sp_index(
    covars,
    cv_fold = sp_fold, blocks, block_id
  )
  ts_index <- generate_cv_index_lbto(
    covars,
    cv_fold = t_fold
  )
  spt_index <- sprintf(
    "S%04d-T%04d",
    covars_sp_index$sp_index, ts_index
  )
  cv_index <- as.numeric(factor(spt_index))
  return(cv_index)
}


#' Generate spatio-temporal cross-validation index (random)
#' @param covars data.frame-like object with `lon`, `lat`, and `time`
#' columns. See [prep_input()] for details.
#' @param cv_fold integer(1). Number of folds for cross-validation.
#' @author Insang Song
#' @return An integer vector with unique values of \code{seq(1, cv_fold)}
#' @export
generate_cv_index_random <- function(
    covars,
    cv_fold = NULL) {
  if (is.null(cv_fold)) {
    stop("Argument cv_fold cannot be NULL. Please set a proper number.\n")
  }
  rows <- nrow(covars)
  cv_index <- sample(seq(1, cv_fold), rows, replace = TRUE)
  return(cv_index)
}
