# data preparation for generate_index

#' Convert sf and sftime objects to stdt format
#' @param stobj An sf, sftime, SpatVector or SpatRasterDataset object
#' @return A stdt object
#' @export
convert_stobj_to_stdt <- function(stobj) {
  geometry <- NULL
  format <- class(stobj)[[1]]
  if (format == "sf" || format == "sftime") {
    if (any(!(c("geometry", "time") %in% colnames(stobj)))) {
      stop("stobj does not contain geometry and time columns")
    }
    crs_dt <- as.character(sf::st_crs(stobj))[1]
    stobj$lon <- sf::st_coordinates(stobj)[, 1]
    stobj$lat <- sf::st_coordinates(stobj)[, 2]
    stdt <- data.table::as.data.table(stobj)
    stdt <- stdt[, geometry := NULL]
  } else if (format == "SpatVector") {
    if (!("time") %in% names(stobj)) {
      stop("stobj does not contain time column")
    }
    crs_dt <- terra::crs(stobj)
    stdf <- as.data.frame(stobj, geom = "XY")
    names(stdf)[names(stdf) == "x"] <- "lon"
    names(stdf)[names(stdf) == "y"] <- "lat"
    stdt <- data.table::as.data.table(stdf)
  } else if (format == "SpatRasterDataset") {
    crs_dt <- terra::crs(stobj)
    stdf <- as.data.frame(stobj[1], xy = TRUE)
    colnames(stdf)[1] <- "lon"
    colnames(stdf)[2] <- "lat"
    # -- tranform from wide to long format
    stdf <- stdf %>% tidyr::pivot_longer(
      cols = 3:ncol(stdf),
      names_to = "time",
      values_to = names(stobj)[1]
    )
    for (var in names(stobj)[2:length(names(stobj))]) {
      # test that the ts is identical to the ts of the 1st variable
      if (!(identical(names(stobj[var]), names(stobj[1])))) {
        stop("Error in SpatRastDataset: timeserie is different for at least
             2 variables - or not ordered for one of these.")
      }
      df_var <- as.data.frame(stobj[var], xy = TRUE)
      # -- tranform from wide to long format
      df_var <- df_var %>% tidyr::pivot_longer(
        cols = 3:ncol(df_var),
        names_to = "time",
        values_to = var
      )
      stdf[, var] <- df_var[, var]
    }
    stdt <- data.table::as.data.table(stdf)
  } else {
    stop("stobj class not accepted")
  }
  stdtobj <- create_stdtobj(stdt, crs_dt)
  return(stdtobj)
}
