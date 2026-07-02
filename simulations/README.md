# Bias simulation study for the pseudocensus OLS estimators

A Monte Carlo study of the finite-sample and asymptotic bias of the estimators
implemented in this package, designed to the reporting standard of a
methodological journal submission: explicit estimands, paired oracle
comparisons, Monte Carlo standard errors on every bias estimate, and empirical
convergence-rate diagnostics.

## Model

Populations of size `n` are drawn from the superpopulation model

    y_i = m(x_i) + e_i,        E[e_i | x_i] = 0,
    x_i ~ F_X  (AR(1) correlation, rho = 0.5),

and each unit's response is observed with propensity `pi(x_i)` (donors) or
missing (recipients). The working model is always linear, so the estimand is
the superpopulation **linear projection parameter**

    beta* = argmin_b E[(y - (1, x')b)^2],

which equals the true coefficients when `m` is linear and remains well defined
under misspecification. `beta*` is evaluated once per data-generating process
by projecting `m(X)` on a noiseless draw of 2 million units.

## Estimators (8 per replicate)

| name              | description                                                    |
|-------------------|----------------------------------------------------------------|
| `census`          | oracle OLS on the fully observed population (reference)        |
| `sample_ols`      | complete-case OLS on donors (`fit_sample_ols`)                 |
| `pseudo_knn`      | pseudocensus OLS, uniform k-NN weights (classical benchmark)   |
| `pseudo_local`    | pseudocensus OLS, `derive_local_debiasing_weights`             |
| `pseudo_global`   | pseudocensus OLS, `derive_global_debiasing_weights`            |
| `debiased_knn`    | `fit_debiased_pseudocensus_ols`, uniform weights               |
| `debiased_local`  | `fit_debiased_pseudocensus_ols`, local weights                 |
| `debiased_global` | `fit_debiased_pseudocensus_ols`, global weights                |

Weights are derived on the design matrix including the intercept column, so
reproducing the recipient covariate imposes the sum-to-one constraint.
Recipient responses are set to `NA` before estimation, guaranteeing that no
estimator can leak unobserved responses.

## Design

**Core factorial (27 cells).** `mean_spec x missing_spec x n`, fully crossed
at the base configuration (p = 4, k = 25, Gaussian covariates, homoskedastic
errors, 70% response rate):

- `mean_spec`: `linear` (correct specification), `quadratic` (smooth
  misspecification: curvature + interaction), `rough` (non-smooth: oscillation
  + kink) — misspecification is what separates complete-case OLS, which is
  biased for `beta*` under covariate-dependent missingness, from local
  imputation, which is not.
- `missing_spec`: `mcar` (benchmark; everything should be unbiased),
  `mar` (logistic propensity in x1, x2, calibrated intercept), `shift`
  (near-complete response below the 0.6 quantile of x1, sparse above —
  recipients concentrate where donors are scarce, stressing nearest-neighbour
  extrapolation at the boundary).
- `n in {500, 2000, 8000}`: three sample sizes per configuration support
  log-log convergence-rate checks — variance-driven error should decay at
  rate n^(-1/2) while genuine bias plateaus.

**Sensitivity arms (6 cells).** One factor at a time around the anchor cell
(`quadratic`, `mar`, n = 2000), chosen because both a misspecification and a
covariate-dependent missingness mechanism are active there:

- errors: heteroskedastic `(0.5 + |x1|) N(0,1)`, heavy-tailed `t_3/sqrt(3)`;
- covariates: standardised-lognormal margins via a Gaussian copula (`skewed`);
- dimension: p = 8;
- neighbourhood size: k = 6 (near-minimal, p + 2) and k = 100 (local weight
  systems underdetermined; minimum-norm solution).

## Metrics

For every cell x estimator x coefficient:

- `bias_star` — Monte Carlo bias against `beta*` (total bias), with `mc_se`,
  its Monte Carlo standard error;
- `bias_census` — mean *paired* deviation from the same-replicate census OLS
  fit; this differences out sampling noise common to both fits and isolates
  the bias induced by imputation itself;
- `sd`, `rmse` — Monte Carlo standard deviation and RMSE against `beta*`.

Aggregate (L2 over coefficients) versions are reported per estimator, and the
core grid produces a 3x3 panel of log-log bias-versus-n plots
(`rate_check.pdf`).

## Hypotheses being tested

1. Under `linear` + any missingness, all estimators are unbiased for `beta*`;
   this is the null calibration check (bias within Monte Carlo error).
2. Under misspecification + `mar`/`shift`, `sample_ols` has non-vanishing bias
   for `beta*` (its bias curve plateaus in n), while the pseudocensus
   estimators' bias shrinks with n as the nearest-neighbour radius contracts.
3. The debiasing correction removes the covariate-reproduction error term:
   `debiased_*` should dominate the corresponding `pseudo_*` on `bias_census`,
   most visibly for `pseudo_knn` (whose uniform weights do not reproduce
   covariates) and under `shift` (boundary extrapolation).
4. Local weights should beat global weights on bias and lose on variance;
   k trades bias (large k, more smoothing) against weight variance (small k).

## Running

From the package root:

```sh
# Full study (33 cells x 1000 replicates)
Rscript simulations/run_simulation.R --reps 1000 --workers 8

# Core factorial only, or specific cells
Rscript simulations/run_simulation.R --cells core
Rscript simulations/run_simulation.R --cells 4,13,22

# Aggregate into summary.csv / summary_norms.csv / rate_check.pdf
Rscript simulations/summarise.R
```

Each cell is saved to `simulations/results/cell_*.rds` as it completes;
re-running skips finished cells, so interrupted runs resume. Every replicate
seeds the RNG deterministically from `(cell, replicate)` — results are
reproducible independently of worker count and scheduling.
