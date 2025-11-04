## --- write 36 LaTeX tables from metrics.RData --------------------------------
## expected objects in metrics.RData: bias, variance, rmse, CI_coverage
## dims: [n, prop_trt, confdg_lvl, method(4), estimand(2), estimator(2)]

rm(list = ls())

# Path to metrics file
metrics_path = file.choose()
stopifnot(file.exists(metrics_path))
load(metrics_path)

# sanity checks
need = c("bias","variance","rmse","CI_coverage")
missing = setdiff(need, ls())
if (length(missing)) stop("Missing objects in metrics.RData: ", paste(missing, collapse=", "))

# unpack values
dn = dimnames(bias)

# scenario dimensions
Ns          = dn$n                    # c("250","500","1000","2000")
prop_trts   = dn$prop_trt             # c("low","moderate","high")
complex_lvls= dn$confdg_lvl           # c("low","moderate","high")

# method / estimand / estimator dimensions
methods     = dn$method               # should be c("IPTW","EB","KOM","TLF") in your file
estimands   = dn$estimand             # c("ATE","ATT")
estimators  = dn$estimator            # c("WLS","DR")

# enforce method order and estimator order you want in tables
methods_order    = c("IPTW","EB","KOM","TLF")
estimators_order = c("WLS","DR")

# labels
estimator_label = c("WLS" = "WLS",
                     "DR"  = "DR")
method_label    = c("IPTW"="IPTW","EB"="EB","KOM"="KOM","TLF"="TLF")

# column headers (LaTeX-escaped)
col_names = c("Bias","Var","RMSE","MAE","Coverage")

# formatters
fmt3 = function(x) format(round(x, 3), nsmall = 3, trim = TRUE)
fmt2 = function(x) format(round(x, 2), nsmall = 2, trim = TRUE)
fmt_cov = function(x) {                   # CI_coverage is in [0,1]
  p = 100 * as.numeric(x)
  paste0(format(round(p, 1), nsmall = 1, trim = TRUE), "\\%")
}

# output file
out_file = "scenario_tables.tex"
sink(out_file)

# iterate scenarios
for (n in Ns) {
  for (pt in prop_trts) {
    for (cx in complex_lvls) {
      
      # begin table
      cat("\\begin{table}\n")
      cat(sprintf("\\caption{Scenario: $n=%s$, proportion treated: %s; complexity: %s.}\n",
                  n, pt, cx))
      cat(sprintf("\\label{tab:metrics_n%s_pt%s_cx%s}\n", n, pt, cx))
      
      # 2 label columns + 10 metric columns (5 for ATE, 5 for ATT)
      cat("\\centering\n")
      cat("\\small\n")
      cat("\\setlength{\\tabcolsep}{6pt}% tighter columns\n")
      cat("\\begin{tabular}{ll", paste(rep("c", 10), collapse=""), "}\n", sep="")
      cat("\\hline\n")
      
      # header row 1
      cat("\\multirow{2}{*}{Estimator} & \\multirow{2}{*}{Method} &",
          "\\multicolumn{5}{c}{ATE} & \\multicolumn{5}{c}{ATT}\\\\\n", sep = " ")
      # header row 2
      cat(" & & ",
          paste(col_names, collapse=" & "),
          " & ",
          paste(col_names, collapse=" & "),
          "\\\\\n\\hline\n", sep = "")
      
      # rows: WLS block then DR block
      for (est in estimators_order) {
        first_est_row = TRUE
        for (met in methods_order) {
          
          # pull metrics and format
          ATE_bias = fmt3(bias      [n, pt, cx, met, "ATE", est])
          ATE_var  = fmt3(variance  [n, pt, cx, met, "ATE", est])
          ATE_rmse = fmt3(rmse      [n, pt, cx, met, "ATE", est])
          ATE_mae  = fmt3(mae       [n, pt, cx, met, "ATE", est])
          ATE_cov  =     fmt_cov(CI_coverage[n, pt, cx, met, "ATE", est])
          
          ATT_bias = fmt3(bias      [n, pt, cx, met, "ATT", est])
          ATT_var  = fmt3(variance  [n, pt, cx, met, "ATT", est])
          ATT_rmse = fmt3(rmse      [n, pt, cx, met, "ATT", est])
          ATT_mae  = fmt3(mae       [n, pt, cx, met, "ATT", est])
          ATT_cov  =     fmt_cov(CI_coverage[n, pt, cx, met, "ATT", est])
          
          # row
          cat(
            if (first_est_row)
              sprintf("\\multirow{%d}{*}{%s}", length(methods_order), estimator_label[est])
            else "",
            
            sprintf(" & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\\n",
                    method_label[met],
                    ATE_bias, ATE_var, ATE_rmse, ATE_mae, ATE_cov,
                    ATT_bias, ATT_var, ATT_rmse, ATT_mae, ATT_cov),
            sep = ""
          )
          
          first_est_row = FALSE
        }
        # optional: partial rule between WLS and DR blocks
        if (est != tail(estimators_order, 1)) cat("\\cline{1-12}\n")
      }
      
      cat("\\hline\n\\end{tabular}\n")
      cat("\\end{table}\n\n")
    }
  }
}

sink()
cat("Wrote LaTeX tables to ", normalizePath(out_file, mustWork = FALSE), "\n")
