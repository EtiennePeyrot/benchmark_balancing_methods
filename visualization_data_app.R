rm(list=ls())

# ---- Packages ----
library(shiny)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tibble)

# ---- Helpers ----
expit = function(x) exp(x) / (1 + exp(x))
logit = function(p) log(p / (1 - p))

# Fixed coefficients from your DGM
b_vec = c(0.8, -0.25, 0.6, -0.4, -0.8, -0.5, 0.7, 0, 0, 0)
a_vec = c(0.9, -1.08, -2.19, -0.6, 0, 0, 0, 0.71, -0.19, 0.26)

# helper: which covariates are binary vs continuous (matches your DGM)
bin_idx  = c(1, 3, 5, 6, 8, 9)
cont_idx = setdiff(1:10, bin_idx)

# Fixed covariance matrix

# Build Sigma with specified off-diagonals
S = diag(rep(1,10))
S[5,1] = S[1,5] = 0.2
S[8,3] = S[3,8] = 0.2
S[6,2] = S[2,6] = 0.9
S[9,4] = S[4,9] = 0.9

# ---- UI ----
ui = fluidPage(
  titlePanel("Propensity scores & crude ATE sandbox"),
  sidebarLayout(
    sidebarPanel(
      numericInput("nsim", "Sample size (nsim):", value = 100000, min = 1000, step = 5000),
      sliderInput("g", "Confounding strength g:", min = 0, max = 5, value = 2.25, step = 0.05),
      sliderInput("ptarget", "Target proportion treated:", min = 0.05, max = 0.95, value = 0.50, step = 0.05),
      sliderInput("alpha", "Alpha (assess strong overlap):", min = 0.00, max = 0.15, value = 0.01, step = 0.01),
      actionButton("go", "Generate / Update")
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          title = "Propensity & summaries",
          h4("Propensity score density"),
          plotOutput("ps_density", height = 350),
          br(),
          fluidRow(
            column(6,
                   h4("Outcome"),
                   tableOutput("evt_table")
            ),
            column(6,
                   h4("Treatment"),
                   tableOutput("summary_table")
            )
          )
        ),
        tabPanel(
          title = "Covariate overlap",
          h4("Continuous covariates: densities by group"),
          plotOutput("cont_overlap", height = 450),
          br(),
          h4("Binary covariates: stacked proportions by group"),
          plotOutput("bin_overlap", height = 450)
        )
      )
    )
  )
)

# ---- Server ----
server = function(input, output, session) {
  
  # Regenerate data when button pressed or inputs change
  dat = eventReactive(input$go, {
    
    nsim = input$nsim
    g = input$g
    
    # Generate X ~ MVN and discretize selected dims
    X = rmvnorm(nsim, mean = rep(0, 10), sigma = S)
    idx_bin = c(1, 3, 5, 6, 8, 9)
    X[, idx_bin] = (X[, idx_bin] > 0) * 1
    
    colnames(X) = paste0("X", 1:10)
    
    # Linear part for ps
    lin_ps_core = as.numeric(X %*% b_vec + 0.5 * X[,1] * (X[,2]^2))
    
    # calibrate b0 so that mean(ps) ≈ ptarget
    target = input$ptarget
    fn = function(c0) mean(expit(c0 + g * lin_ps_core)) - target
    root = try(uniroot(fn, lower = -15, upper = 15), silent = TRUE)
    b0 = if (inherits(root, "try-error")) 0 else root$root
    
    
    ps = expit(b0 + g * lin_ps_core)
    
    # Treatment A ~ Bernoulli(ps)
    A = rbinom(nsim, 1, ps)
    
    # Randomized B with same marginal prevalence as A
    B = rbinom(nsim, 1, mean(A))
    
    # Outcome generation (no treatment effect)
    lin_y = as.numeric(X %*% a_vec + 0.5 * X[,3] * (X[,4]^2))
    # Calibrate a0 so that mean(Y) ≈ 25%
    fn_a0 = function(a0) mean(expit(a0 + g * lin_y)) - 0.25
    root_a0 = try(uniroot(fn_a0, lower = -15, upper = 15), silent = TRUE)
    a0 = if (inherits(root_a0, "try-error")) -2.22 else root_a0$root
    pY = expit(a0 + g * lin_y)
    Y = rbinom(nsim, 1, pY)
    
    tibble(
      ps = ps,
      A = A,
      B = B,
      Y = Y
    ) %>%
      bind_cols(as_tibble(X)) %>% 
      mutate(b0 = b0, g = g, a0 = a0)
  }, ignoreInit = FALSE)
  
  # ---- Summaries ----
  output$pA = renderText({ sprintf("%.3f", mean(dat()$A)) })
  output$mps = renderText({ sprintf("%.3f", mean(dat()$ps)) })
  output$p_low = renderText({ sprintf("%.3f", mean(dat()$ps < input$alpha)) })
  output$p_high = renderText({ sprintf("%.3f", mean(dat()$ps > (1 - input$alpha))) })
  output$b0 = renderText({ sprintf("%.3f", unique(dat()$b0)) })
  
  output$ateA = renderText({
    d = dat()
    m1 = mean(d$Y[d$A == 1])
    m0 = mean(d$Y[d$A == 0])
    sprintf("%.4f (diff in means)", m1 - m0)
  })
  output$ateB = renderText({
    d = dat()
    m1 = mean(d$Y[d$B == 1])
    m0 = mean(d$Y[d$B == 0])
    sprintf("%.4f (diff in means, randomized B)", m1 - m0)
  })
  
  # ---- Main Plots ----
  output$ps_density = renderPlot({
    d = dat()
    # For speed, sample for plotting if needed
    n_plot = min(nrow(d), 100000)
    idx = sample.int(nrow(d), n_plot)
    dd = d[idx, ]
    dd$Group = ifelse(dd$A,"treated","control")
    
    ggplot(dd, aes(x = ps, group=Group, colour=Group, fill=Group)) +
      annotate("rect", xmin = 0, xmax = input$alpha, ymin = -Inf, ymax = Inf, alpha = 0.15) +
      annotate("rect", xmin = 1 - input$alpha, xmax = 1, ymin = -Inf, ymax = Inf, alpha = 0.15) +
      geom_density(alpha=.25) +
      geom_vline(xintercept = c(input$alpha, 1 - input$alpha), linetype = 2) +
      scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
      labs(title = "Propensity score density",
           subtitle = sprintf("alpha = %.2f; shaded tails are ps < alpha and ps > 1 - alpha", input$alpha),
           x = "Propensity score (ps)", y = "Density") +
      theme_minimal(base_size = 12)
  })
  
  
  
  
  output$evt_table = renderTable({
    d = dat()
    ev = mean(d$Y)
    evA0 = mean(d$Y[d$A == 0])
    evA1 = mean(d$Y[d$A == 1])
    ate = evA1 - evA0
    tibble(
      Metric = c(
        "Event rate: general population",
        "Event rate: A=0",
        "Event rate: A=1",
        "Crude ATE",
        "a0 (calibrated to E(Y)=25%)"
      ),
      Value = c(ev, evA0, evA1, ate,unique(d$a0))
    ) 
  })
  output$summary_table = renderTable({
    d = dat()
    tibble(
      Metric = c(
        "n",
        "Proportion treated (A=1)",
        sprintf("Pr(ps < %.2f[alpha])", input$alpha),
        sprintf("Pr(ps > %.2f[1-alpha])", 1-input$alpha),
        paste0("b0 (calibrated to E(A)=",input$ptarget*100,"%)")
      ),
      Value = c(
        nrow(d),
        sprintf("%.4f", mean(d$A)),
        sprintf("%.4f", mean(d$ps < input$alpha)),
        sprintf("%.4f", mean(d$ps > (1 - input$alpha))),
        sprintf("%.3f", unique(d$b0))
      )
    )
  }, striped = TRUE, bordered = TRUE, digits = 4)
  
  
  # ---- Covariates Overlap Plots ----
  
  output$cont_overlap = renderPlot({
    d = dat()
    
    A_fac = factor(d$A, levels = c(0, 1), labels = c("Control", "Treated"))
    cont_vars = paste0("X", cont_idx)
    
    # long data for continuous vars
    dat_cont_list = lapply(cont_vars, function(v) {
      data.frame(var = v, A = A_fac, value = d[[v]], stringsAsFactors = FALSE)
    })
    dat_cont = do.call(rbind, dat_cont_list)
    
    
    ggplot(dat_cont, aes(x = value, fill = A, colour = A)) +
      geom_density(alpha = 0.35, adjust = 1.2, na.rm = TRUE) +
      facet_wrap(~ var, scales = "free", ncol = 2) +
      labs(title = "Overlap on continuous covariates (pre-weighting)",
           x = NULL, y = "Density", fill = "Group", colour = "Group") +
      theme_minimal()
  })
  
  output$bin_overlap = renderPlot({
    d = dat()
    
    A_fac = factor(d$A, levels = c(0, 1), labels = c("Control", "Treated"))
    bin_vars = paste0("X", bin_idx)
    
    # long data for binary vars
    dat_bin_list = lapply(bin_vars, function(v) {
      data.frame(var = v, A = A_fac, level = factor(d[[v]]), stringsAsFactors = FALSE)
    })
    dat_bin = do.call(rbind, dat_bin_list)
    dat_bin = dat_bin[!is.na(dat_bin$level), , drop = FALSE]
    
    ggplot(dat_bin, aes(x = A, fill = level)) +
      geom_bar(position = "fill", width = 0.7) +
      facet_wrap(~ var, ncol = 3) +
      labs(title = "Binary covariates: distribution by group (pre-weighting)",
           x = NULL, y = "Proportion", fill = "Level") +
      theme_minimal()
  })
}

shinyApp(ui, server)
