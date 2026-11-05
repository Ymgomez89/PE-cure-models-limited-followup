# Piecewise Exponential Cure Models under Limited Follow-up

R code for parametric cure models (Bernoulli, Poisson and Negative Binomial formulations) with Weibull and piecewise exponential (PE / power PE) baseline distributions, and for a non-cure piecewise exponential proportional hazards comparator.

The code reproduces all results of the study *"Piecewise Exponential Cure Models under Limited Follow-up: An Exploratory Reanalysis of ST-Elevation Myocardial Infarction after Primary Percutaneous Coronary Intervention"* (earlier working title: *"Cure-Model Approaches Using Piecewise Exponential Distributions in STEMI Patients Undergoing Primary PCI"*).

## Authors

- Yolanda M. Gómez (Departamento de Estadística, Universidad del Bío-Bío, Concepción, Chile)
- Allan Burboa (Departamento de Estadística, Universidad del Bío-Bío, Concepción, Chile)

## Overview

Conventional parametric non-cure survival models (e.g., the Weibull model) assume a proper survival distribution, under which all individuals would eventually experience the event. The Kaplan–Meier estimator and the Cox model do not impose this assumption on the observed data (the Cox baseline hazard is unspecified), but neither includes an explicit parameter for individuals who will never experience the event. Cure models introduce such a non-susceptible fraction; its identification, however, requires a follow-up long enough for the survival curve to reach a plateau. In the STEMI cohort analysed here (30-month follow-up), the cure fraction was weakly identified, and a flexible non-cure PE model fitted better than the cure models.

The repository provides R functions (`R/00_functions.R`) to fit:

- Mixture cure models: Bernoulli formulation (`BerPE`, `BerPPE`, `BerW`).
- Non-mixture (promotion time) cure models: Poisson (`PoPE`, `PoPPE`, `PoW`) and Negative Binomial (`NBPE`, `NBPPE`, `NBW`).
- Baseline distributions: Weibull, piecewise exponential (PE) and power piecewise exponential (PPE).
- Non-cure PE proportional hazards model (`PEPH.est2`), with independent verification by a Poisson GLM (`PEPH.glm`).

Main features:

1. Maximum likelihood estimation on the log scale for positive parameters, with Nelder–Mead followed by BFGS refinement; convergence is checked at each stage, and a solution is accepted only if the observed information matrix is positive definite and the gradient is negligible (< 10⁻³).
2. Multiple-start verification of the maximum (`multistart.check`).
3. Standard errors from the observed information matrix (delta method for positive parameters).
4. Individual cure probabilities and the average cure fraction, with simulation-based and profile-likelihood intervals.
5. Residual diagnostics: Cox–Snell residuals and randomized quantile residuals.

## Requirements

```r
install.packages(c("survival", "msm", "pracma", "MASS", "readxl"))
```

## Data

`data/infarto.xlsx`: public dataset of Yang et al. (2017), Dryad Digital Repository, doi:10.5061/dryad.pf56m (464 patients, 118 MACE, 30-month follow-up).

## Software (versions used for the reported results)

- R 4.3.3 (2024-02-29)
- survival 3.5-8, msm 1.7.1, pracma 2.4.4, MASS 7.3-60, readxl 1.4.3

## How to run (reproduction of the manuscript)

Run the scripts from the repository root, in order (e.g., `Rscript R/01_data_preparation.R`). Pre-computed results are provided in `output/` (CSV files), so that they can be compared without rerunning the long analyses (partition scan and bootstrap).

| Script | Output | Manuscript item |
|---|---|---|
| `R/00_functions.R` | Likelihoods and estimation functions (sourced by the other scripts) | Methods (Estimation) |
| `R/01_data_preparation.R` | `output/d12.rds`, `output/data12.rds` | Variable coding |
| `R/02_table1.R` | `output/Table1.csv` | Table 1 |
| `R/03_cox_models.R` | Cox model, Schoenfeld test, sex-stratified model, univariable models, backward elimination | Table 2; Tables S2–S4 |
| `R/04_weibull_cure.R` | Weibull cure models | Table 3; Table 5 |
| `R/05_pe_partition_scan.R` (or `05b`/`05c`) | AIC for L = 1, …, 20 (long run time) | Figure 2 |
| `R/06_pe_cure_models.R` | PE cure models (L = 8), cure fraction, profile likelihood, events / patients at risk / person-time per interval, non-cure PE proportional hazards comparator (estimates, hazard ratios, convergence diagnostics, Poisson-GLM verification), multiple-start verification | Table 4; Table 5; Figure S1; Tables S7–S9 |
| `R/07_sensitivity.R` | Alternative partitions, backward elimination of PoPE(8), bootstrap with per-resample convergence diagnostics (long run time) | Tables S1, S4, S5, S6 |
| `R/07b_bootstrap_chunk.R`, `R/07c_bootstrap_summary.R` | Optional: run the bootstrap of `07_sensitivity.R` in chunks (`Rscript R/07b_bootstrap_chunk.R s b_from b_to`, identical seeds and settings) and assemble/summarize it | Table S6 |
| `R/08_figures.R` | Figures 1–3 and S1 in TIFF (600 dpi, LZW) and PDF; summary of randomized quantile residuals | Figures 1–3, S1; Table S10 |

## Estimation settings

- Parameters estimated on the log scale for hazard rates, α, λ and q + 1.
- Log-likelihoods of the PE cure models are evaluated in a numerically stable form (log νᵢ = `plogis(η, log.p = TRUE)`, `log1p`, `expm1`), so that the evaluation does not fail when νᵢ approaches 0 or 1 (e.g., in bootstrap resamples with monotone likelihood).
- Optimization (function `.fit.optim` in `00_functions.R`), for every start:
  1. Nelder–Mead (`optim`, maxit = 10000, default relative tolerance ≈ 1.5 × 10⁻⁸); the start is discarded unless the convergence code is 0.
  2. BFGS refinement from the Nelder–Mead solution (maxit = 10000, reltol = 10⁻¹²). The BFGS solution is accepted only if **its own convergence code is 0** and it does not increase the negative log-likelihood; otherwise the converged Nelder–Mead solution is retained.
  3. The observed information matrix (numerical Hessian, `pracma::hessian`) must be finite and positive definite (all eigenvalues > 0), and the maximum absolute gradient (`pracma::grad`) must be < 10⁻³; otherwise a new start is drawn.
  The stage of the final solution, both convergence codes, the maximum absolute gradient and the extreme eigenvalues / condition number of the observed information are returned in `fit$diagnostics` (`output/convergence_diagnostics.csv`).
- Initial values: β = 0; λₗ = 1/t̄ (PE cure models), λₗ = (number of events)/(total follow-up) (non-cure PE model) or α = 1, λ = 3t̄ (Weibull models), perturbed by N(0, 0.1²) noise; Poisson estimates are used as initial values for the Negative Binomial models. Up to 200 starts are attempted until a solution satisfying 1–3 is obtained.
- Multiple-start verification (`multistart.check`): each PE model (BerPE, PoPE, NBPE and the non-cure PE model) is refitted from 50 independent random starts (perturbation N(0, 0.5²) on the internal scale; seed 2021) and the maximized log-likelihoods are compared with the reported fit (`output/TableS9_multistart.csv`).
- Weibull cure models (`04_weibull_cure.R`): at most 30 starts per fit. Models without a valid finite solution (diverging coefficients or scale parameter, i.e. likelihood maximized at the boundary of the parameter space) are reported as not identifiable; for them, the supremum of the log-likelihood is obtained from 20 random starts (Nelder–Mead followed by BFGS) and only the log-likelihood, AIC and α are reported (Table 3). With the corrected workflow, this applies to the Poisson and Negative Binomial Weibull models.
- Partition scan (`05_pe_partition_scan.R`): same workflow, at most 30 starts per fit; fits without a valid solution are reported as missing (Figure 2). `05b`/`05c` are optional helpers to run the scan cell by cell.
- Non-cure PE proportional hazards comparator: fitted with the same `.fit.optim` workflow and independently verified by a Poisson GLM on episode-split data with log-exposure offset (`PEPH.glm`), which must reproduce the coefficients, standard errors and log-likelihood.
- Bootstrap (`07_sensitivity.R`): same workflow as the main analysis, with at most 30 starts per resample (`max.tries = 30`); resamples without a valid solution are recorded as NA. Convergence diagnostics (starts used, final stage, maximum absolute gradient, maximum |β|) are stored for each resample (`output/bootstrap_PoPE8.csv`); resamples with max |β| > 10 are flagged as divergent (monotone likelihood: no finite maximum likelihood estimate of the cure component).
- Standard errors: inverse of the observed information matrix at the maximum likelihood estimate; delta method for positive parameters; 95% Wald intervals on the log scale for positive parameters.
- Partitions: empirical quantiles (R `quantile`, type 7) of the observed event times; intervals are closed on the left, [a_{l−1}, a_l).
- Randomized quantile residuals (Dunn and Smyth, 1996): u = F_pop(y) for events and u ~ Uniform(F_pop(y), 1) for censored observations, r = Φ⁻¹(u). **Correction in this version:** the previous code applied `qnorm(runif(u))`, which generated new independent uniforms instead of transforming u; this has been corrected (`rqr()` in `00_functions.R`).

## Random seeds

- Model fitting: `set.seed(2021)` before each fit.
- Average cure fraction interval: 5000 draws, seed 2021.
- Bootstrap: seed `1000 * s + b` (batch s = 1, …, 4; resample b = 1, …, 50).
- Multiple-start verification: seed 2021.
- Randomized quantile residuals: seed 2021 (100 randomizations; Q–Q plot of the median order statistics; median mean, SD and Shapiro–Wilk / Kolmogorov–Smirnov p-values across randomizations).
