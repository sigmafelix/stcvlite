testthat::test_that("prep_input converts sf input to data.table with coordinates and time", {
  testthat::skip_if_not_installed("sf")

  df <- data.frame(
    Date = as.Date(c("2020-01-01", "2020-01-02")),
    val = c(10, 20),
    x = c(10, 11),
    y = c(45, 46)
  )
  obj_sf <- sf::st_as_sf(df, coords = c("x", "y"), crs = 4326)

  out <- suppressWarnings(stcvlite::prep_input(obj_sf))

  testthat::expect_s3_class(out, "data.table")
  testthat::expect_true(all(c("lon", "lat", "time") %in% names(out)))
  testthat::expect_equal(nrow(out), nrow(df))
  testthat::expect_equal(as.character(out$time), as.character(df$Date))
  testthat::expect_true(!is.null(attr(out, "crs")))
})

testthat::test_that("prep_input errors when class is unsupported", {
  testthat::expect_error(
    stcvlite::prep_input(list(a = 1)),
    "stobj class not accepted"
  )
})

testthat::test_that("convert_cv_index_rset builds an rset with expected split count", {
  covars <- make_covars_fixture()
  cvindex <- c(1, 2, 1, 2, 1, 2)

  out <- stcvlite::convert_cv_index_rset(
    cvindex = cvindex,
    data = covars,
    cv_mode = "lolo"
  )

  testthat::expect_s3_class(out, "rset")
  testthat::expect_equal(nrow(out), 2)
  testthat::expect_true(all(grepl("^cvfold_lolo_", out$id)))
})
