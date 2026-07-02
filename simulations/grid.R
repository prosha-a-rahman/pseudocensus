# Factorial design of the simulation study.
#
# The study is organised as a full 3x3x3 core factorial plus one-factor-at-a-
# time sensitivity arms around a fixed anchor cell:
#
#   core (27 cells):  mean_spec x missing_spec x n, at the base configuration
#                     (p = 4, k = 25, gaussian covariates, homoskedastic
#                     errors, response rate 0.7). Crossing n in {500, 2000,
#                     8000} within every mean/missingness pair supports
#                     empirical convergence-rate checks for the bias.
#
#   sensitivity:      anchored at (quadratic, mar, n = 2000) -- a cell where
#                     the working model is misspecified and missingness is
#                     covariate-dependent, so every departure from the base
#                     configuration has a bias mechanism to interact with.
#                       * error_spec  in {heteroskedastic, heavy_tailed}
#                       * covariate_design = skewed
#                       * p = 8
#                       * k in {6, 100}  (near-minimal p + 2, and heavily
#                         overdetermined where local weights fall back to the
#                         minimum-norm solution)

BASE_CONFIG <- list(
  n = 2000L,
  p = 4L,
  k = 25L,
  mean_spec = "quadratic",
  error_spec = "homoskedastic",
  missing_spec = "mar",
  covariate_design = "gaussian",
  rho = 0.5,
  response_rate = 0.7
)

make_cells <- function(block, ...) {
  overrides <- expand.grid(..., stringsAsFactors = FALSE)
  fixed <- BASE_CONFIG[setdiff(names(BASE_CONFIG), names(overrides))]
  cells <- cbind(overrides, fixed, block = block, stringsAsFactors = FALSE)
  cells[c(names(BASE_CONFIG), "block")]
}

build_grid <- function() {
  grid <- rbind(
    make_cells(
      "core",
      mean_spec = c("linear", "quadratic", "rough"),
      missing_spec = c("mcar", "mar", "shift"),
      n = c(500L, 2000L, 8000L)
    ),
    make_cells("sens_error", error_spec = c("heteroskedastic", "heavy_tailed")),
    make_cells("sens_design", covariate_design = "skewed"),
    make_cells("sens_dimension", p = 8L),
    make_cells("sens_neighbours", k = c(6L, 100L))
  )

  grid$cell <- seq_len(nrow(grid))
  grid[c("cell", setdiff(names(grid), "cell"))]
}
