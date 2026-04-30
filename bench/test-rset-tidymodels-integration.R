testthat::test_that("scaled rset works with tidymodels workflow across splits", {
  testthat::skip_if_not_installed("workflows")
  testthat::skip_if_not_installed("parsnip")

  covars <- make_covars_fixture_scaled(n_locations = 100, n_times = 120, seed = 2026)
  cvindex <- stcvlite::generate_cv_index(
    covars,
    cv_mode = "lblo",
    cv_fold = 6
  )

  rset_obj <- stcvlite::convert_cv_index_rset(
    cvindex = cvindex,
    data = covars$stdt,
    cv_mode = "lblo"
  )

  spec <- parsnip::linear_reg() |> parsnip::set_engine("lm")
  wf <- workflows::workflow() |>
    workflows::add_model(spec) |>
    workflows::add_formula(y ~ lon + lat + x1 + x2)

  split_metrics <- lapply(rset_obj$splits, function(s) {
    train_df <- rsample::analysis(s)
    test_df <- rsample::assessment(s)
    fit_obj <- workflows::fit(wf, data = train_df)
    preds <- predict(fit_obj, new_data = test_df)

    data.frame(
      n_assessment = nrow(test_df),
      n_pred = nrow(preds),
      rmse = sqrt(mean((preds$.pred - test_df$y)^2)),
      stringsAsFactors = FALSE
    )
  })

  metric_df <- do.call(rbind, split_metrics)

  testthat::expect_s3_class(rset_obj, "rset")
  testthat::expect_equal(nrow(rset_obj), max(cvindex))
  testthat::expect_true(all(metric_df$n_assessment > 0))
  testthat::expect_true(all(metric_df$n_pred == metric_df$n_assessment))
  testthat::expect_true(all(is.finite(metric_df$rmse)))
})
