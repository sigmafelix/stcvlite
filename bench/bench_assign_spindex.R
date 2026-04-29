# bench/bench_assign_spindex.R
#
# Benchmarks for .assign_spindex_dbscan / _FNN / _N2R.
#
# NOTE: _FNN and _N2R build an n×n adjacency matrix, so they are O(n²) in
# memory.  At 100 K points that would require ~40 GB RAM.  Those two engines
# are therefore only benchmarked up to N = 5 000.  The dbscan engine is
# benchmarked at all sizes including 100 K+.

pkgload::load_all(quiet = TRUE)  # load package from source without installing
library(bench)
library(ggplot2)

# ── helpers ──────────────────────────────────────────────────────────────────

make_data <- function(n, seed = 42L) {
  set.seed(seed)
  data.frame(
    x = runif(n, 0, 100),
    y = runif(n, 0, 100)
  )
}

# Expose internal functions (they are loaded via pkgload::load_all)
.assign_spindex_dbscan <- .assign_spindex_dbscan
.assign_spindex_FNN    <- .assign_spindex_FNN
.assign_spindex_N2R    <- .assign_spindex_N2R

SP_COLS    <- c("x", "y")
NCLUSTERS  <- 5L

# ── 1. Small-scale comparison: all three engines ──────────────────────────────

sizes_small <- c(500L, 1000L, 2000L, 5000L)

# Check whether N2R works with the installed version before running full bench
n2r_ok <- tryCatch({
  .assign_spindex_N2R(make_data(50L), SP_COLS, NCLUSTERS)
  TRUE
}, error = function(e) {
  message("N2R engine skipped – incompatible installed version: ", conditionMessage(e))
  FALSE
})

results_small <- lapply(sizes_small, function(n) {
  message("Benchmarking n = ", n, " (all engines) ...")
  dat <- make_data(n)

  exprs <- list(
    dbscan = quote(.assign_spindex_dbscan(dat, SP_COLS, NCLUSTERS)),
    FNN    = quote(.assign_spindex_FNN(dat, SP_COLS, NCLUSTERS))
  )
  if (n2r_ok) {
    exprs$N2R <- quote(.assign_spindex_N2R(dat, SP_COLS, NCLUSTERS))
  }

  bm <- do.call(bench::mark, c(exprs, list(iterations = 3L, check = FALSE, memory = TRUE)))
  bm$n <- n
  bm
})

bench_small <- do.call(rbind, results_small)
bench_small$engine <- as.character(bench_small$expression)

print(bench_small[, c("n", "engine", "min", "median", "mem_alloc", "n_itr")])

# ── 2. Large-scale: dbscan only ───────────────────────────────────────────────

sizes_large <- c(10000L, 50000L, 100000L, 150000L)

results_large <- lapply(sizes_large, function(n) {
  message("Benchmarking n = ", n, " (dbscan only) ...")
  dat <- make_data(n)

  bm <- bench::mark(
    dbscan = .assign_spindex_dbscan(dat, SP_COLS, NCLUSTERS),
    iterations = 3L,
    check = FALSE,
    memory = TRUE
  )
  bm$n <- n
  bm
})

bench_large <- do.call(rbind, results_large)
bench_large$engine <- "dbscan"

print(bench_large[, c("n", "engine", "min", "median", "mem_alloc", "n_itr")])

# ── 3. Plots ──────────────────────────────────────────────────────────────────

# 3a. Runtime vs n - small scale, all engines
p_small <- ggplot(bench_small, aes(x = n, y = as.numeric(median), colour = engine)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " s")) +
  labs(
    title  = "assign_spindex runtime - all engines (small n)",
    x      = "n (points)",
    y      = "Median runtime (s)",
    colour = "Engine"
  ) +
  theme_minimal()

# 3b. Runtime vs n - dbscan, large scale
dbscan_small_df <- data.frame(
  n = bench_small$n[bench_small$engine == "dbscan"],
  engine = "dbscan",
  median_sec = as.numeric(bench_small$median[bench_small$engine == "dbscan"]),
  mem_alloc_bytes = as.numeric(bench_small$mem_alloc[bench_small$engine == "dbscan"])
)

dbscan_large_df <- data.frame(
  n = bench_large$n,
  engine = "dbscan",
  median_sec = as.numeric(bench_large$median),
  mem_alloc_bytes = as.numeric(bench_large$mem_alloc)
)

all_dbscan <- rbind(dbscan_small_df, dbscan_large_df)

p_large <- ggplot(all_dbscan, aes(x = n, y = median_sec)) +
  geom_line(colour = "#0072B2", linewidth = 0.8) +
  geom_point(colour = "#0072B2", size = 2) +
  scale_y_continuous(labels = function(x) paste0(round(x, 2), " s")) +
  labs(
    title = "assign_spindex_dbscan runtime - large n",
    x     = "n (points)",
    y     = "Median runtime (s)"
  ) +
  theme_minimal()

# Save plots
ggsave("bench/bench_spindex_small.png", p_small, width = 7, height = 4, dpi = 150)
ggsave("bench/bench_spindex_dbscan_large.png", p_large, width = 7, height = 4, dpi = 150)

# Save benchmark summaries
write.csv(
  bench_small[, c("n", "engine", "min", "median", "mem_alloc", "n_itr")],
  "bench/bench_spindex_small.csv",
  row.names = FALSE
)
write.csv(
  bench_large[, c("n", "engine", "min", "median", "mem_alloc", "n_itr")],
  "bench/bench_spindex_dbscan_large.csv",
  row.names = FALSE
)

message("Done. Plots saved to bench/.")
