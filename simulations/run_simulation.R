# Runner for the pseudocensus bias simulation study.
#
# Usage (from the package root):
#   Rscript simulations/run_simulation.R [--reps R] [--workers W]
#                                        [--cells all|core|1,4,9] [--outdir DIR]
#
# Defaults: reps = 1000, workers = all physical cores, cells = all,
# outdir = simulations/results.
#
# Each design cell is written to its own .rds as it completes, so an
# interrupted run can be resumed; existing cell files are skipped.

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- list(
    reps = 1000L,
    workers = max(1L, parallel::detectCores(logical = FALSE)),
    cells = "all",
    outdir = file.path("simulations", "results")
  )

  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(parsed)) stop("unknown argument: ", args[i])
    value <- args[i + 1L]
    parsed[[key]] <- if (key %in% c("reps", "workers")) as.integer(value) else value
    i <- i + 2L
  }

  parsed
}

select_cells <- function(grid, cells) {
  if (identical(cells, "all")) return(grid)
  if (cells %in% unique(grid$block)) return(grid[grid$block == cells, ])
  ids <- as.integer(strsplit(cells, ",")[[1]])
  grid[grid$cell %in% ids, ]
}

main <- function() {
  if (!file.exists("DESCRIPTION")) {
    stop("run from the package root: Rscript simulations/run_simulation.R")
  }

  args <- parse_args()
  dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

  pkgload::load_all(".", quiet = TRUE)
  source(file.path("simulations", "dgp.R"))
  source(file.path("simulations", "replicate.R"))
  source(file.path("simulations", "grid.R"))

  grid <- select_cells(build_grid(), args$cells)
  message(sprintf(
    "Running %d cells x %d replicates on %d workers",
    nrow(grid), args$reps, args$workers
  ))

  # Superpopulation projection targets, one per unique DGP, cached on disk.
  target_path <- file.path(args$outdir, "targets.rds")
  targets <- if (file.exists(target_path)) readRDS(target_path) else list()
  dgp_key <- function(cfg) {
    paste(cfg$mean_spec, cfg$covariate_design, cfg$p, cfg$rho, sep = "|")
  }
  for (i in seq_len(nrow(grid))) {
    cfg <- as.list(grid[i, ])
    key <- dgp_key(cfg)
    if (is.null(targets[[key]])) {
      message("Computing projection target for: ", key)
      targets[[key]] <- projection_target(
        mean_spec = cfg$mean_spec,
        covariate_design = cfg$covariate_design,
        p = cfg$p,
        rho = cfg$rho
      )
      saveRDS(targets, target_path)
    }
  }

  cluster <- NULL
  if (args$workers > 1L) {
    cluster <- parallel::makeCluster(args$workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function() {
      pkgload::load_all(".", quiet = TRUE)
      source(file.path("simulations", "dgp.R"))
      source(file.path("simulations", "replicate.R"))
      NULL
    })
  }

  for (i in seq_len(nrow(grid))) {
    cfg <- as.list(grid[i, ])
    cell_path <- file.path(args$outdir, sprintf("cell_%03d.rds", cfg$cell))
    if (file.exists(cell_path)) {
      message(sprintf("Cell %d already complete, skipping", cfg$cell))
      next
    }

    started <- Sys.time()
    rep_ids <- seq_len(args$reps)
    replicate_list <- if (is.null(cluster)) {
      lapply(rep_ids, function(rep_id) run_replicate(cfg, rep_id))
    } else {
      parallel::parLapplyLB(cluster, rep_ids, function(rep_id, cfg) {
        run_replicate(cfg, rep_id)
      }, cfg = cfg)
    }

    # reps x estimator x coefficient
    estimates <- aperm(simplify2array(replicate_list), c(3, 1, 2))

    saveRDS(
      list(
        config = cfg,
        estimates = estimates,
        target = targets[[dgp_key(cfg)]],
        reps = args$reps
      ),
      cell_path
    )

    message(sprintf(
      "Cell %d/%d (%s, %s, %s, n = %d) done in %.1fs",
      cfg$cell, max(grid$cell), cfg$mean_spec, cfg$missing_spec, cfg$block,
      cfg$n, as.numeric(difftime(Sys.time(), started, units = "secs"))
    ))
  }

  message("All requested cells complete. Summarise with: ",
          "Rscript simulations/summarise.R --outdir ", args$outdir)
}

main()
