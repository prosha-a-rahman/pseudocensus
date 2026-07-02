# Aggregate raw cell results into bias / variance / RMSE tables and
# convergence-rate diagnostics.
#
# Usage (from the package root):
#   Rscript simulations/summarise.R [--outdir simulations/results]
#
# Outputs, written alongside the cell files:
#   summary.csv        -- per cell x estimator x coefficient metrics
#   summary_norms.csv  -- per cell x estimator aggregate (L2) metrics
#   rate_check.pdf     -- log-log plot of ||bias|| against n for the core grid
#
# Reported metrics:
#   bias_star    -- Monte Carlo bias against the superpopulation projection
#                   parameter beta* (total bias)
#   mc_se        -- Monte Carlo standard error of bias_star
#   bias_census  -- mean paired deviation from the same-replicate census OLS
#                   fit, isolating imputation-induced bias from sampling noise
#   sd, rmse     -- Monte Carlo standard deviation and RMSE against beta*

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  outdir <- file.path("simulations", "results")
  if (length(args) >= 2 && args[1] == "--outdir") outdir <- args[2]
  list(outdir = outdir)
}

summarise_cell <- function(cell) {
  cfg <- cell$config
  estimators <- dimnames(cell$estimates)[[2]]
  coefs <- dimnames(cell$estimates)[[3]]
  census <- cell$estimates[, "census", ]

  rows <- lapply(estimators, function(est) {
    draws <- cell$estimates[, est, ]
    bias_star <- colMeans(draws) - cell$target
    mc_sd <- apply(draws, 2, sd)

    data.frame(
      cell = cfg$cell,
      block = cfg$block,
      mean_spec = cfg$mean_spec,
      error_spec = cfg$error_spec,
      missing_spec = cfg$missing_spec,
      covariate_design = cfg$covariate_design,
      n = cfg$n,
      p = cfg$p,
      k = cfg$k,
      estimator = est,
      coefficient = coefs,
      bias_star = bias_star,
      mc_se = mc_sd / sqrt(nrow(draws)),
      bias_census = colMeans(draws - census),
      sd = mc_sd,
      rmse = sqrt(colMeans(sweep(draws, 2, cell$target)^2)),
      row.names = NULL
    )
  })

  do.call(rbind, rows)
}

l2_norms <- function(summary_df) {
  aggregate(
    cbind(bias_star, bias_census, rmse) ~
      cell + block + mean_spec + error_spec + missing_spec +
      covariate_design + n + p + k + estimator,
    data = summary_df,
    FUN = function(v) sqrt(sum(v^2))
  )
}

plot_rate_check <- function(norms, path) {
  core <- norms[norms$block == "core" & norms$estimator != "census", ]
  if (nrow(core) == 0) return(invisible())

  estimators <- sort(unique(core$estimator))
  palette <- seq_along(estimators)
  names(palette) <- estimators

  pdf(path, width = 10, height = 10)
  op <- par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))
  on.exit({ par(op); dev.off() }, add = TRUE)

  for (mean_spec in c("linear", "quadratic", "rough")) {
    for (missing_spec in c("mcar", "mar", "shift")) {
      panel <- core[core$mean_spec == mean_spec &
                      core$missing_spec == missing_spec, ]
      if (nrow(panel) == 0) {
        plot.new()
        title(main = paste(mean_spec, "/", missing_spec, "(no cells)"))
        next
      }
      plot(
        NA,
        xlim = range(log10(panel$n)),
        ylim = range(log10(pmax(panel$bias_star, 1e-6))),
        xlab = "log10(n)",
        ylab = "log10 ||bias||",
        main = paste(mean_spec, "/", missing_spec)
      )
      for (est in estimators) {
        line <- panel[panel$estimator == est, ]
        line <- line[order(line$n), ]
        lines(log10(line$n), log10(pmax(line$bias_star, 1e-6)),
              col = palette[est], type = "b", pch = 19, cex = 0.6)
      }
      if (mean_spec == "linear" && missing_spec == "mcar") {
        legend("bottomleft", legend = estimators, col = palette,
               lty = 1, pch = 19, cex = 0.55, bty = "n")
      }
    }
  }

  invisible()
}

main <- function() {
  args <- parse_args()
  cell_files <- list.files(args$outdir, pattern = "^cell_\\d+\\.rds$",
                           full.names = TRUE)
  if (length(cell_files) == 0) stop("no cell results found in ", args$outdir)

  summary_df <- do.call(rbind, lapply(cell_files, function(f) {
    summarise_cell(readRDS(f))
  }))
  norms <- l2_norms(summary_df)

  write.csv(summary_df, file.path(args$outdir, "summary.csv"),
            row.names = FALSE)
  write.csv(norms, file.path(args$outdir, "summary_norms.csv"),
            row.names = FALSE)
  plot_rate_check(norms, file.path(args$outdir, "rate_check.pdf"))

  # Console digest: aggregate bias norms at the anchor sample size.
  digest <- norms[norms$n == 2000 & norms$estimator != "census", ]
  digest <- digest[order(digest$cell, digest$bias_star), ]
  cat("\n||bias|| against beta* at n = 2000 (smaller is better):\n\n")
  print(
    digest[c("cell", "block", "mean_spec", "missing_spec", "estimator",
             "bias_star", "bias_census", "rmse")],
    digits = 3,
    row.names = FALSE
  )

  message(sprintf(
    "\nWrote %s, %s and %s",
    file.path(args$outdir, "summary.csv"),
    file.path(args$outdir, "summary_norms.csv"),
    file.path(args$outdir, "rate_check.pdf")
  ))
}

main()
