# Choosing Covariate Balancing Methods for Causal Inference

*Code for the simulation study and companion Shiny app*

**Repository:** https://github.com/EtiennePeyrot/benchmark_balancing_methods

## Overview

This repo contains code and materials for the article:

> **Choosing Covariate Balancing Methods for Causal Inference: Practical Insights from a Simulation Study**
> Etienne Peyrot, Raphaël Porcher, François Petit

We benchmark weighting approaches including **IPTW**, **standard Covariate Balancing Propensity Score (CBPS-JI and CBPS-OI)**, **Energy Balancing (EB)**, **Kernel Optimal Matching (KOM)**, and **Covariate Balancing by Tailored Loss Functions (TLF)**, paired with **Weighted Least Squares (WLS)** and a **Doubly-Robust (DR)** estimator. We evaluate **ATE** and **ATT** across **36 main simulation scenarios** varying sample size (250, 500, 1000, 2000), treatment prevalence (25%, 50%, 75%), and a complexity factor that jointly increases confounding and reduces overlap. Additional sensitivity analyses repeat the same scenario grid under constant non-null treatment effects.

### TL;DR of findings

* DR generally reduces sensitivity to the choice of weights when the outcome model is reasonable.
* EB and KOM behave similarly under transparent tuning; EB is convenient, KOM requires kernel/penalty choices.
* IPTW and standard CBPS are useful reference propensity-score-based methods.
* IPTW is workable but variance-sensitive as complexity rises.
* TLF often shows low variance but higher bias, leading to sub-nominal CI coverage without outcome modeling.

---

## What’s included

* `display app/` – **R Shiny app** to explore: impact of hyperparameters in the data-generative mechanism, and metrics (bias, variance, RMSE, MAE, CI coverage, variance ratio) across all scenarios.
* `code/` – Scripts/functions to generate data, compute weights, run estimators, and summarize results.
* `results/` – Saved outputs (figures/tables), including main simulation results and non-null treatment-effect sensitivity analyses.

> The data-generating mechanism, scenario grid, and all performance metrics are documented in the paper and supplement. The exact values of `(a0, b0, g)` for each scenario are listed in the **Supplementary Material** (Tables 1–2).

---

## Installation

Tested with **R 4.2.1**.

### Required packages

* [`WeightIt`](https://cran.r-project.org/package=WeightIt) (EB weights, `method = "energy"`)
* [`rootSolve`](https://cran.r-project.org/package=rootSolve) (standard CBPS-JI)
* [`osqp`](https://cran.r-project.org/package=osqp) (QP solver used in KOM implementation)
* `shiny` (for the app)
* Common deps you likely already have: `mvtnorm`, `stats`, `utils`, etc.

```r
install.packages(c("WeightIt", "rootSolve", "osqp", "shiny", "mvtnorm"))
```
