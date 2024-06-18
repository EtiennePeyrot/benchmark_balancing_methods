This project compare 3 recent (~2020) balancing methods with Inverse Probability of Treatment Weighting (IPTW) ("The central role of the propensity score in observational studies for causal effects" Rosenbaum, Paul R. and Rubin, Donald B.). The recents methods are :
  - kernel optimal matching ("Generalized Optimal Matching Methods for Causal Inference" Nathan Kallus, "Optimal Estimation of Generalized Average Treatment Effects using Kernel Optimal Matching" Nathan Kallus and Michele Santacatterina),
  - energy balancing ("Energy Balancing of Covariate Distributions" Jared D. Huling and Simon Mak),
  - covariates balancing by tailored loss function ("Covariate balancing propensity score by tailored loss functions" Qingyuan Zhao).

To run the app, dowload the folder "display app", open file "app.R" and change the path in the function "setwd" (5th line of the file) to the path of the folder "display app". Run "app.R" with R. This code requires the R package shiny, if not already installed, app.R should automatically install it.
