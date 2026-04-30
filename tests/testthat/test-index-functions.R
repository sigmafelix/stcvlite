testthat::test_that("generate_spt_index creates spatial and spatiotemporal identifiers", {
  covars <- make_covars_fixture()

  out_sp <- stcvlite::generate_spt_index(covars, mode = "spatial")
  out_spt <- stcvlite::generate_spt_index(covars, mode = "spatiotemporal")

  testthat::expect_true("sp_index" %in% names(out_sp$stdt))
  testthat::expect_equal(length(unique(out_sp$stdt$sp_index)), 3)
  testthat::expect_equal(length(unique(out_spt$stdt$sp_index)), nrow(covars$stdt))
})

testthat::test_that("generate_block_sp_index supports kmeans and numeric blocks", {
  covars <- make_covars_fixture()

  out_km <- stcvlite::generate_block_sp_index(covars, cv_fold = 2)
  testthat::expect_true("sp_index" %in% names(out_km$stdt))

  out_grid <- stcvlite::generate_block_sp_index(covars, blocks = c(1, 1))
  testthat::expect_true("sp_index" %in% names(out_grid$stdt))
  testthat::expect_equal(length(out_grid$stdt$sp_index), nrow(covars$stdt))
})

testthat::test_that("index generators return expected fold structures", {
  covars <- make_covars_fixture()

  idx_loto <- stcvlite::generate_cv_index_loto(covars)
  idx_lolo <- stcvlite::generate_cv_index_lolo(covars)
  idx_lolto <- stcvlite::generate_cv_index_lolto(covars)
  idx_lblo <- stcvlite::generate_cv_index_lblo(covars, cv_fold = 2)
  idx_lbto <- stcvlite::generate_cv_index_lbto(covars, cv_fold = 2)
  idx_lblto <- stcvlite::generate_cv_index_lblto(covars, sp_fold = 2, t_fold = 2, blocks = NULL)

  testthat::expect_equal(length(unique(idx_loto)), 2)
  testthat::expect_equal(length(unique(idx_lolo)), 3)
  testthat::expect_equal(idx_lolto, seq_len(nrow(covars$stdt)))
  testthat::expect_equal(length(idx_lblo), nrow(covars$stdt))
  testthat::expect_true(all(idx_lblo >= 1 & idx_lblo <= 2))
  testthat::expect_true(all(idx_lbto >= 1 & idx_lbto <= 2))
  testthat::expect_true(all(idx_lblto >= 1))
})

testthat::test_that("lblo cases are consistent across generate_cv_index* functions", {
  covars <- make_covars_fixture()

  set.seed(123)
  idx_lblo_direct <- stcvlite::generate_cv_index_lblo(covars, cv_fold = 2)
  set.seed(123)
  idx_lblo_dispatch <- stcvlite::generate_cv_index(covars, cv_mode = "lblo", cv_fold = 2)

  testthat::expect_equal(idx_lblo_dispatch, idx_lblo_direct)

  idx_lblo_blocks_direct <- stcvlite::generate_cv_index_lblo(
    covars,
    cv_fold = NULL,
    blocks = c(1, 1)
  )
  idx_lblo_blocks_dispatch <- stcvlite::generate_cv_index(
    covars,
    cv_mode = "lblo",
    cv_fold = 2,
    blocks = c(1, 1)
  )

  testthat::expect_equal(idx_lblo_blocks_dispatch, idx_lblo_blocks_direct)
  testthat::expect_equal(length(idx_lblo_blocks_direct), nrow(covars$stdt))

  testthat::expect_error(
    stcvlite::generate_cv_index_lblo(covars, cv_fold = NULL),
    "cannot be NULL"
  )
  testthat::expect_error(
    stcvlite::generate_cv_index(covars, cv_mode = "lblo", cv_fold = NULL),
    "invalid"
  )
})

testthat::test_that("generate_cv_index dispatches correctly for all available modes", {
  covars <- make_covars_fixture()

  idx_lolo <- stcvlite::generate_cv_index(covars, cv_mode = "lolo")
  idx_loto <- stcvlite::generate_cv_index(covars, cv_mode = "loto")
  idx_lolto <- stcvlite::generate_cv_index(covars, cv_mode = "lolto")
  set.seed(777)
  idx_lblo <- stcvlite::generate_cv_index(covars, cv_mode = "lblo", cv_fold = 2)
  idx_lbto <- stcvlite::generate_cv_index(covars, cv_mode = "lbto", cv_fold = 2)
  set.seed(888)
  idx_lblto <- stcvlite::generate_cv_index(
    covars,
    cv_mode = "lblto",
    sp_fold = 2,
    t_fold = 2,
    blocks = NULL
  )
  set.seed(99)
  idx_random <- stcvlite::generate_cv_index(covars, cv_mode = "random", cv_fold = 3)

  testthat::expect_equal(idx_lolo, stcvlite::generate_cv_index_lolo(covars))
  testthat::expect_equal(idx_loto, stcvlite::generate_cv_index_loto(covars))
  testthat::expect_equal(idx_lolto, stcvlite::generate_cv_index_lolto(covars))
  set.seed(777)
  testthat::expect_equal(idx_lblo, stcvlite::generate_cv_index_lblo(covars, cv_fold = 2))
  testthat::expect_equal(idx_lbto, stcvlite::generate_cv_index_lbto(covars, cv_fold = 2))
  set.seed(888)
  testthat::expect_equal(
    idx_lblto,
    stcvlite::generate_cv_index_lblto(covars, sp_fold = 2, t_fold = 2, blocks = NULL)
  )
  set.seed(99)
  testthat::expect_equal(
    idx_random,
    stcvlite::generate_cv_index_random(covars, cv_fold = 3)
  )
})

testthat::test_that("generate_cv_index validates arguments", {
  covars <- make_covars_fixture()

  testthat::expect_error(
    stcvlite::generate_cv_index(covars$stdt),
    "Only stdt object is acceptable"
  )
  testthat::expect_error(
    stcvlite::generate_cv_index(covars, cv_mode = "loto", sp_fold = 2),
    "only applicable to"
  )
})

testthat::test_that("generate_cv_index_lbto and random validate required cv_fold", {
  covars <- make_covars_fixture()

  testthat::expect_error(
    stcvlite::generate_cv_index_lbto(covars, cv_fold = NULL),
    "cannot be NULL"
  )
  testthat::expect_error(
    stcvlite::generate_cv_index_random(covars, cv_fold = NULL),
    "cannot be NULL"
  )

  set.seed(42)
  idx_random <- stcvlite::generate_cv_index_random(covars, cv_fold = 4)
  testthat::expect_true(all(idx_random %in% 1:4))
})
