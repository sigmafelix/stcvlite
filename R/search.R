.assign_spindex <-
  function(
    data,
    sp_cols = c("lon", "lat"),
    nclusters = 5L,
    engine = c("dbscan", "FNN"),
    ...
  ) {
    stopifnot(inherits(data, "data.frame"))
    # spmat <- as.matrix(sp_subset)
    # if (any(!is.finite(spmat))) {
    #   stop("Spatial columns must contain only finite values.")
    # }

    worker <- switch(engine,
      dbscan = .assign_spindex_dbscan,
      FNN = .assign_spindex_FNN,
      stop("Invalid engine. Choose one of 'dbscan', 'FNN'.")
    )
    worker(data, sp_cols, nclusters, ...)
  }


.assign_spindex_FNN <-
  function(
    data,
    sp_cols,
    nclusters,
    k = 10L
  ) {
    sp_subset <-
      if (data.table::is.data.table(data)) {
        as.data.frame(data)[, sp_cols, drop = TRUE]
      } else {
        data[, sp_cols, drop = FALSE]
      }
    x <- as.matrix(sp_subset)
    storage.mode(x) <- "double"

    knn_index <- FNN::knn.index(x, k = k)

    adj <- matrix(0L, nrow = nrow(x), ncol = nrow(x))
    for (i in seq_len(nrow(knn_index))) {
      adj[i, knn_index[i, ]] <- 1L
    }
    adj <- adj + t(adj)
    adj[adj > 0L] <- 1L

    cl <- integer(nrow(x))
    current_cluster <- 0L
    unvisited <- seq_len(nrow(x))

    while (length(unvisited) > 0L) {
      seed <- unvisited[1]
      current_cluster <- current_cluster + 1L
      queue <- seed
      cl[seed] <- current_cluster
      unvisited <- setdiff(unvisited, seed)

      while (length(queue) > 0L) {
        node <- queue[1]
        queue <- queue[-1]
        neighbors <- which(adj[node, ] > 0L)
        new_nodes <- intersect(neighbors, unvisited)
        cl[new_nodes] <- current_cluster
        unvisited <- setdiff(unvisited, new_nodes)
        queue <- c(queue, new_nodes)
      }
    }

    if (current_cluster != nclusters) {
      sizes <- tabulate(cl, nbins = current_cluster)
      top_clusters <-
        order(sizes, decreasing = TRUE)[
          seq_len(min(nclusters, current_cluster))
        ]
      remap <- integer(current_cluster)
      remap[top_clusters] <- seq_along(top_clusters)
      cl_new <- remap[cl]
      noise <- cl_new == 0L
      if (any(noise) && any(!noise)) {
        nn <- FNN::get.knnx(
          data = x[!noise, , drop = FALSE],
          query = x[noise, , drop = FALSE],
          k = 1L
        )
        cl_new[noise] <- cl_new[!noise][nn$nn.index[, 1]]
      }
      cl <- cl_new
    }

    data$sp_index <- as.integer(cl)
    data
  }



.assign_spindex_dbscan <-
  function(
    data,
    sp_cols,
    nclusters,
    minPts = NULL,
    search_quantiles = c(0.01, 0.99),
    max_iter = 40L,
    reassign_noise = TRUE
  ) {
    stopifnot(inherits(data, "data.frame"))
    stopifnot(length(sp_cols) >= 2L)
    stopifnot(is.numeric(nclusters), length(nclusters) == 1L, nclusters >= 1L)

    sp_subset <-
      if (data.table::is.data.table(data)) {
        as.data.frame(data)[, sp_cols, drop = TRUE]
      } else {
        data[, sp_cols, drop = FALSE]
      }
    x <- as.matrix(sp_subset)
    storage.mode(x) <- "double"

    if (is.null(minPts)) {
      minPts <- max(4L, floor(log(nrow(x))))
    }
    minPts <- as.integer(minPts)
    k <- max(2L, minPts - 1L)

    kd <- FNN::knn.dist(x, k = k)
    kd_last <- sort(kd[, k], na.last = NA)

    eps_lo <-
      as.numeric(
        stats::quantile(
          kd_last,
          probs = search_quantiles[1],
          names = FALSE,
          type = 8
        )
      )
    eps_hi <-
      as.numeric(
        stats::quantile(
          kd_last,
          probs = search_quantiles[2],
          names = FALSE,
          type = 8
        )
      )

    eps_lo <- max(eps_lo, .Machine$double.eps)
    eps_hi <- max(eps_hi, eps_lo * 1.0001)

    n_clusters_at <- function(eps) {
      fit <- dbscan::dbscan(x, eps = eps, minPts = minPts, borderPoints = TRUE)
      length(setdiff(unique(fit$cluster), 0L))
    }

    n_lo <- n_clusters_at(eps_lo)
    n_hi <- n_clusters_at(eps_hi)

    guard <- 0L
    while ((n_lo < nclusters || n_hi > nclusters) && guard < 25L) {
      if (n_lo < nclusters) eps_lo <- eps_lo * 0.5
      if (n_hi > nclusters) eps_hi <- eps_hi * 1.5
      n_lo <- n_clusters_at(eps_lo)
      n_hi <- n_clusters_at(eps_hi)
      guard <- guard + 1L
    }

    best_fit <- NULL
    best_eps <- NA_real_
    best_diff <- Inf

    for (i in seq_len(max_iter)) {
      eps_mid <- (eps_lo + eps_hi) / 2
      fit <-
        dbscan::dbscan(
          x,
          eps = eps_mid,
          minPts = minPts,
          borderPoints = TRUE
        )
      n_mid <- length(setdiff(unique(fit$cluster), 0L))
      d <- abs(n_mid - nclusters)

      if (d < best_diff) {
        best_fit <- fit
        best_eps <- eps_mid
        best_diff <- d
      }

      if (n_mid == nclusters) {
        best_fit <- fit
        best_eps <- eps_mid
        break
      }

      if (n_mid > nclusters) {
        eps_lo <- eps_mid
      } else {
        eps_hi <- eps_mid
      }
    }

    cl <- best_fit$cluster

    pos <- sort(unique(cl[cl > 0L]))
    if (length(pos) > 0L) {
      remap <- stats::setNames(seq_along(pos), pos)
      out <- integer(length(cl))
      is_pos <- cl > 0L
      out[is_pos] <- unname(remap[as.character(cl[is_pos])])
      cl <- out
    }

    if (reassign_noise && any(cl == 0L) && any(cl > 0L)) {
      nn <- FNN::get.knnx(
        data = x[cl > 0L, , drop = FALSE],
        query = x[cl == 0L, , drop = FALSE],
        k = 1L
      )
      cl[cl == 0L] <- cl[cl > 0L][nn$nn.index[, 1]]
    }

    # Enforce a bounded number of cluster labels for CV folds.
    # DBSCAN can legitimately return many dense components
    # (e.g., repeated coordinates over time), so we normalize
    # to the requested number of spatial clusters.
    positive <- cl[cl > 0L]
    n_pos <- length(unique(positive))
    n_unique_points <- nrow(unique(x))
    target_clusters <- min(as.integer(nclusters), n_unique_points)

    if (n_pos == 0L) {
      km <- stats::kmeans(x, centers = target_clusters, nstart = 10)
      cl <- km$cluster
    } else if (n_pos > target_clusters) {
      sizes <- tabulate(cl, nbins = max(cl))
      top_clusters <- order(sizes, decreasing = TRUE)[seq_len(target_clusters)]
      remap <- integer(max(cl))
      remap[top_clusters] <- seq_along(top_clusters)
      cl_new <- remap[cl]
      dropped <- cl_new == 0L

      if (any(dropped) && any(!dropped)) {
        nn <- FNN::get.knnx(
          data = x[!dropped, , drop = FALSE],
          query = x[dropped, , drop = FALSE],
          k = 1L
        )
        cl_new[dropped] <- cl_new[!dropped][nn$nn.index[, 1]]
      }
      cl <- cl_new
    } else if (n_pos < target_clusters) {
      km <- stats::kmeans(x, centers = target_clusters, nstart = 10)
      cl <- km$cluster
    }

    data$sp_index <- ifelse(cl == 0L, NA_integer_, as.integer(cl))
    attr(data, "dbscan_eps") <- best_eps
    attr(data, "dbscan_minPts") <- minPts
    data
  }


.bin_temporal <-
  function(
    data,
    time_col,
    n_bins
  ) {
    stopifnot(inherits(data, "data.frame"))
    stopifnot(time_col %in% names(data))
    stopifnot(is.numeric(n_bins), length(n_bins) == 1L, n_bins >= 1L)

    time_values <- data[[time_col]]
    if (
      !is.numeric(time_values) && !inherits(time_values, "Date") &&
        !inherits(time_values, "POSIXct")
    ) {
      stop("Time column must be numeric, Date, or POSIXct.")
    }

    breaks <-
      seq(
        min(time_values, na.rm = TRUE),
        max(time_values, na.rm = TRUE),
        length.out = n_bins + 1
      )
    data$.time_bin <-
      cut(
        time_values,
        breaks = breaks,
        include.lowest = TRUE,
        labels = FALSE
      )
    data
  }
