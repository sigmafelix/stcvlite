#' Plot generated CV folds
#'
#' Creates a plotly scatterplot (x, y, fold index as z) with fold/cluster
#' number shown by color. For large datasets, points are downsampled to
#' at most `max_points` for display.
#'
#' @param covars A `data.frame`/`data.table` (or inheriting object) that
#'   contains `lon` and `lat` columns.
#' @param cv_index Fold assignment vector returned by [generate_cv_index()].
#'   Must have length `nrow(covars)`.
#' @param max_points Maximum number of points to display (default `10000`).
#' @param seed Random seed used when downsampling.
#'
#' @return A `plotly` 3D scatter object.
#' @examples
#' data(spdat)
#' cv_index <-
#'   generate_cv_index(spdat, cv_mode = "lblo", cv_fold = 6)
#'
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   p <- plot_cv_folds(spdat, cv_index)
#'   p
#' }
#' @export
plot_cv_folds <-
  function(
    covars,
    cv_index,
    max_points = 1200L,
    seed = 2026L
  ) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for plot_cv_folds_25d().", call. = FALSE)
    }

    coords <- covars

    if (!is.data.frame(coords) || !all(c("lon", "lat") %in% names(coords))) {
      stop("`covars` must contain `lon` and `lat` columns.", call. = FALSE)
    }

    n <- nrow(coords)
    if (length(cv_index) != n) {
      stop("`cv_index` length must equal number of rows in `covars`.", call. = FALSE)
    }
    fold_id <- as.integer(as.factor(cv_index))

    plot_df <- data.frame(
      x = coords[["lon"]],
      y = coords[["lat"]],
      z = coords[["time"]],
      fold = fold_id,
      stringsAsFactors = FALSE
    )
    plot_df <- plot_df[!is.na(plot_df$fold), , drop = FALSE]

    if (nrow(plot_df) > max_points) {
      set.seed(seed)
      keep <- sample.int(nrow(plot_df), size = max_points)
      plot_df <- plot_df[keep, , drop = FALSE]
    }

    plot_df$fold_f <- as.factor(plot_df$fold)

    plotly::plot_ly(
      data = plot_df,
      x = ~x,
      y = ~y,
      z = ~z,
      type = "scatter3d",
      mode = "markers",
      color = ~fold_f,
      colors = "Set3",
      marker = list(size = 3, opacity = 0.85)
    ) |>
      plotly::layout(
        scene = list(
          xaxis = list(title = "lon"),
          yaxis = list(title = "lat"),
          zaxis = list(title = "Time")
        ),
        legend = list(title = list(text = "Fold / Cluster"))
      )
  }
