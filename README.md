# Choosing Covariate Balancing Methods for Causal Inference
_Code for the simulation study and companion Shiny app_

**Repository:** https://github.com/EtiennePeyrot/benchmark_balancing_methods

## Overview

This repo contains code and materials for the article:

> **Choosing Covariate Balancing Methods for Causal Inference: Practical Insights from a Simulation Study**  
> Etienne Peyrot, Raphaël Porcher, François Petit

We benchmark four weighting approaches—**IPTW**, **Energy Balancing (EB)**, **Kernel Optimal Matching (KOM)**, and **Covariate Balancing by Tailored Loss Functions (TLF)**—paired with **Weighted Least Squares (WLS)** and a **Doubly-Robust (DR)** estimator. We evaluate **ATE** and **ATT** across **36 scenarios** varying sample size (250, 500, 1000, 2000), treatment prevalence (25%, 50%, 75%), and a complexity factor that jointly increases confounding and reduces overlap. 

### TL;DR of findings
- DR generally reduces sensitivity to the choice of weights when the outcome model is reasonable.  
- EB and KOM behave similarly under transparent tuning; EB is convenient, KOM requires kernel/penalty choices.  
- IPTW is workable but variance-sensitive as complexity rises.  
- TLF often shows low variance but higher bias, leading to sub-nominal CI coverage without outcome modeling. 

---

## What’s included

- `display app/` – **R Shiny app** to explore: impact of hyperparameters in the data-generative mechanism, and metrics (bias, variance, RMSE, MAE, CI coverage) across all scenarios.  
- `code/` – Scripts/functions to generate data, compute weights, run estimators, and summarize results.  
- `results/` – Saved outputs (figures/tables) if you choose to persist them locally.  

> The data-generating mechanism, scenario grid, and all performance metrics are documented in the paper and supplement. The exact values of `(a0, b0, g)` for each scenario are listed in the **Supplementary Material** (Tables 1–2). 

---

## Installation

Tested with **R 4.2.1**.

### Required packages
- [`WeightIt`](https://cran.r-project.org/package=WeightIt) (EB weights, `method = "energy"`)
- [`osqp`](https://cran.r-project.org/package=osqp) (QP solver used in KOM implementation)
- `shiny` (for the app)
- Common deps you likely already have: `mvtnorm`, `stats`, `utils`, etc.

```r
install.packages(c("WeightIt", "osqp", "shiny", "mvtnorm"))
