# app.R
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

if (!"hp" %in% ls()) load(choose.files())
# hp = hp[,,,,,15:16,,]

# ---- Helpers ---------------------------------------------------------------

# Convert the hp array to a long data.frame (one row per cell & seed)
array_to_df <- function(hp) {
  stopifnot(is.array(hp), !is.null(dimnames(hp)))
  
  df <- as.data.frame(as.table(hp), stringsAsFactors = FALSE)
  names(df) <- c(names(dimnames(hp)), "value")
  
  # Make some useful types/factors with sensible ordering
  dims <- dimnames(hp)
  
  df$estimand   <- factor(df$estimand, levels = dims$estimand)
  df$n          <- factor(df$n, levels = dims$n)
  df$prop_trt   <- factor(df$prop_trt, levels = dims$prop_trt)
  df$confdg_lvl <- factor(df$confdg_lvl, levels = dims$confdg_lvl)
  df$lambda     <- factor(df$lambda, levels = dims$lambda)
  df$sigma      <- factor(df$sigma, levels = dims$sigma)
  df$metric     <- factor(df$metric, levels = dims$metric)
  df$seed       <- as.integer(df$seed)
  
  df
}

# Safe "mode" (first most-frequent value; ignores NA)
stat_mode <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

num_agg_fun <- function(key) {
  switch(key,
         "mean"   = function(x) mean(x, na.rm = TRUE),
         "median" = function(x) median(x, na.rm = TRUE),
         "mode"   = stat_mode,
         "sd"     = function(x) sd(x, na.rm = TRUE))
}

pretty_agg_label <- function(key) {
  switch(key,
         "mean"   = "Mean",
         "median" = "Median",
         "mode"   = "Mode",
         "sd"     = "SD",
         "wins"   = "Wins",
         "Agg")
}

# One function to build the heatmap consistently
# - keeps all factor levels on axes
# - NA tiles are transparent
# - if is_count = TRUE we show integer steps (0..max) in the legend
plot_heat <- function(dat, title, legend_title = "Value",
                      transform_fn = identity, is_count = FALSE) {
  validate(need(nrow(dat) > 0, "No data to display for the current selection."))
  
  p <- ggplot(dat, aes(x = sigma, y = lambda, fill = transform_fn(value))) +
    geom_tile() +
    ggtitle(title) +
    facet_wrap(vars(estimand), nrow = 1, drop = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    coord_flip() +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE)
  
  if (is_count) {
    maxv <- suppressWarnings(max(dat$value, na.rm = TRUE))
    if (!is.finite(maxv)) maxv <- 0
    p + scale_fill_steps(name = legend_title,
                         breaks = 0:maxv,
                         labels = function(x) ifelse(x %in% as.integer(seq(0, maxv, length.out = 4)), x, ""),
                         show.limits = TRUE,
                         na.value = NA)
  } else {
    p + scale_fill_gradient(name = legend_title, na.value = NA)
  }
}

# ---- UI --------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Hyper-parameter heatmaps"),
  sidebarLayout(
    sidebarPanel(
      uiOutput("control_n"),
      uiOutput("control_prop"),
      uiOutput("control_confdg"),
      radioButtons(
        "agg",
        "Aggregate across seeds",
        inline  = TRUE,
        choices = c("Mean" = "mean", "Median" = "median", "Mode (mod)" = "mode",
                    "SD" = "sd", "Wins (best across seeds)" = "wins"),
        selected = "mean"
      ),
      # Cap control for log10(||grad CBSR||)
      checkboxInput("cap_cbsr_enable", "Cap log10(||grad CBSR||)", value = FALSE),
      sliderInput("cap_cbsr_value", "Upper bound (log10 scale):", value = 10, step = .1, min = -1, max = 10)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Max |SMD Pooled|",          plotOutput("plot_smd",      height = "600px")),
        tabPanel("log norm grad CBSR",        plotOutput("plot_cbsr",     height = "600px")),
        tabPanel("Weight Coefficient of Var", plotOutput("plot_wcv",      height = "600px")),
        tabPanel("Effect (raw weights)",        plotOutput("plot_ate_eff",  height = "600px")),
        tabPanel("Effect (normalized weights)", plotOutput("plot_ate_effn", height = "600px")),
        tabPanel("log norm vs lambda", plotOutput("plot_cbsr_lines", height = "600px")),
        tabPanel("SMD vs log norm", plotOutput("plot_smd_vs_cbsr", height = "600px"))
      )
    )
  )
)

# ---- Server ----------------------------------------------------------------

server <- function(input, output, session) {
  # Ensure hp exists
  validate_hp <- reactive({
    validate(need(exists("hp", inherits = TRUE), "Object 'hp' not found in the environment."))
    get("hp", inherits = TRUE)
  })
  
  dims <- reactive({ dimnames(validate_hp()) })
  
  # Dynamic controls based on hp dimnames
  output$control_n <- renderUI({
    selectInput("n", "Sample size (n):", choices = dims()$n, selected = dims()$n[[1]])
  })
  output$control_prop <- renderUI({
    selectInput("prop_trt", "Treatment proportion:", choices = dims()$prop_trt, selected = dims()$prop_trt[[1]])
  })
  output$control_confdg <- renderUI({
    selectInput("confdg_lvl", "Confounding level:", choices = dims()$confdg_lvl, selected = dims()$confdg_lvl[[1]])
  })
  
  # Long data frame from hp
  hp_df <- reactive({ array_to_df(validate_hp()) })
  
  # How many seeds in the current slice (for legend context)
  n_seeds <- reactive({
    req(input$n, input$prop_trt, input$confdg_lvl)
    base <- hp_df() |>
      filter(n == input$n,
             prop_trt == input$prop_trt,
             confdg_lvl == input$confdg_lvl,
             !is.na(value))
    as.integer(dplyr::n_distinct(base$seed))
  })
  
  # Aggregation function selector
  agg_fun <- reactive({
    switch(input$agg,
           "mean"   = function(x) mean(x, na.rm = TRUE),
           "median" = function(x) median(x, na.rm = TRUE),
           "mode"   = stat_mode,
           "sd"     = function(x) sd(x, na.rm = TRUE))
  })
  
  # Aggregate across seeds for selected n/prop_trt/confdg_lvl
  agg_df <- reactive({
    req(input$n, input$prop_trt, input$confdg_lvl)
    base <- hp_df() |>
      filter(n == input$n,
             prop_trt == input$prop_trt,
             confdg_lvl == input$confdg_lvl,
             !is.na(value))
    
    if (input$agg != "wins") {
      out <- base |>
        group_by(estimand, sigma, lambda, metric) |>
        summarise(value = agg_fun()(value), .groups = "drop")
    } else {
      # Wins: one integer winner per seed (lexicographic tie-break: sigma, then lambda)
      sig_levels <- dims()$sigma
      lam_levels <- dims()$lambda
      
      winners <- base |>
        group_by(estimand, metric, seed) |>
        filter(value == min(value, na.rm = TRUE)) |>
        mutate(
          .sigma_i  = as.integer(factor(sigma,  levels = sig_levels)),
          .lambda_i = as.integer(factor(lambda, levels = lam_levels))
        ) |>
        arrange(.sigma_i, .lambda_i) |>
        slice_head(n = 1) |>
        ungroup()
      
      out <- winners |>
        group_by(estimand, sigma, lambda, metric) |>
        summarise(value = as.integer(dplyr::n()), .groups = "drop")
    }
    
    # lock factor orders
    out <- out |>
      mutate(
        lambda   = factor(lambda,   levels = dims()$lambda),
        sigma    = factor(sigma,    levels = dims()$sigma),
        estimand = factor(estimand, levels = dims()$estimand),
        metric   = factor(metric,   levels = dims()$metric)
      )
    
    # ensure FULL grid so rows/cols aren’t dropped
    out <- out |>
      tidyr::complete(
        estimand = factor(dims()$estimand, levels = dims()$estimand),
        metric   = factor(dims()$metric,   levels = dims()$metric),
        sigma    = factor(dims()$sigma,    levels = dims()$sigma),
        lambda   = factor(dims()$lambda,   levels = dims()$lambda)
      )
    
    # for Wins, empty cells mean 0 wins
    if (input$agg == "wins") {
      out <- out |> mutate(value = dplyr::coalesce(value, 0L))
    }
    
    out
  })
  
  # --- NEW: ATE-only aggregation for Effect & Effect norm (always numeric agg; never "wins")
  agg_df_ate_effect <- reactive({
    req(input$n, input$prop_trt, input$confdg_lvl)
    base <- hp_df() |>
      filter(n == input$n,
             prop_trt == input$prop_trt,
             confdg_lvl == input$confdg_lvl,
             metric %in% c("Effect", "Effect norm"),
             !is.na(value))
    
    out <- base |>
      group_by(estimand, sigma, lambda, metric) |>
      summarise(value = agg_fun()(value), .groups = "drop") |>
      mutate(
        lambda   = factor(lambda,   levels = dims()$lambda),
        sigma    = factor(sigma,    levels = dims()$sigma),
        estimand = factor(estimand, levels = dims()$estimand),
        metric   = factor(metric,   levels = dims()$metric)
      ) |>
      tidyr::complete(
        estimand = factor(dims()$estimand, levels = dims()$estimand),         # <-- ATE & ATT
        metric   = factor(c("Effect","Effect norm"), levels = dims()$metric),
        sigma    = factor(dims()$sigma,  levels = dims()$sigma),
        lambda   = factor(dims()$lambda, levels = dims()$lambda)
      )
    
    out
  })  
  # ---- Plots ---------------------------------------------------------------
  
  output$plot_smd <- renderPlot({
    dat <- agg_df() |> filter(metric == "max SMD")
    if (input$agg == "wins") {
      leg <- paste0("Wins (out of ", n_seeds(), ")")
      plot_heat(dat,
                title = "Wins – Max |SMD Pooled|",
                legend_title = leg,
                transform_fn = identity,
                is_count = TRUE)
    } else {
      leg <- paste0(pretty_agg_label(input$agg), " max |SMD|")
      plot_heat(dat,
                title = "Max |SMD Pooled|",
                legend_title = leg)
    }
  })
  
  output$plot_cbsr <- renderPlot({
    dat <- agg_df() |> filter(metric == "grad CBSR")
    if (input$agg == "wins") {
      leg <- paste0("Wins (out of ", n_seeds(), ")")
      plot_heat(dat,
                title = "Wins – grad CBSR",
                legend_title = leg,
                transform_fn = identity,
                is_count = TRUE)
    } else {
      # If capping is enabled, pre-transform to log10 and set values > cap to NA.
      if (isTRUE(input$cap_cbsr_enable)) {
        dat <- dat |>
          mutate(value = suppressWarnings(log10(value)),
                 value = ifelse(!is.finite(value) | value > input$cap_cbsr_value,
                                NA_real_, value))
        leg <- paste0(pretty_agg_label(input$agg),
                      " log10{‖grad CBSR‖} (≤ ", input$cap_cbsr_value, ")")
        plot_heat(dat,
                  title = "log norm grad CBSR (capped)",
                  legend_title = leg,
                  transform_fn = identity,
                  is_count = FALSE)
      } else {
        leg <- paste0(pretty_agg_label(input$agg), " log10{‖grad CBSR‖}")
        plot_heat(dat,
                  title = "log norm grad CBSR",
                  legend_title = leg,
                  transform_fn = function(x) log10(x),
                  is_count = FALSE)
      }
    }
  })
  
  output$plot_wcv <- renderPlot({
    dat <- agg_df() |> filter(metric == "weight CV")
    if (input$agg == "wins") {
      leg <- paste0("Wins (out of ", n_seeds(), ")")
      plot_heat(dat,
                title = "Wins – Weight CV",
                legend_title = leg,
                transform_fn = identity,
                is_count = TRUE)
    } else {
      leg <- paste0(pretty_agg_label(input$agg), " weight CV")
      plot_heat(dat,
                title = "Weight Coefficient of Variation",
                legend_title = leg)
    }
  })
  
  output$plot_ate_eff <- renderPlot({
    dat <- agg_df_ate_effect() |> dplyr::filter(metric == "Effect")
    agg_key <- if (identical(input$agg, "wins")) "mean" else input$agg
    leg <- paste0(pretty_agg_label(agg_key), " effect")
    plot_heat(dat, title = "ATE / ATT – effect (raw weights)", legend_title = leg)
  })
  
  output$plot_ate_effn <- renderPlot({
    dat <- agg_df_ate_effect() |> dplyr::filter(metric == "Effect norm")
    agg_key <- if (identical(input$agg, "wins")) "mean" else input$agg
    leg <- paste0(pretty_agg_label(agg_key), " effect (normalized)")
    plot_heat(dat, title = "ATE / ATT – effect (normalized weights)", legend_title = leg)
  })
  
  output$plot_cbsr_lines <- renderPlot({
    req(input$n, input$prop_trt, input$confdg_lvl)
    
    # If user chose "Wins", fall back to mean for this numeric plot
    agg_key <- if (identical(input$agg, "wins")) "mean" else input$agg
    agg_f   <- num_agg_fun(agg_key)
    
    # Build aggregated data for current slice
    dat <- hp_df() |>
      dplyr::filter(n == input$n,
                    prop_trt == input$prop_trt,
                    confdg_lvl == input$confdg_lvl,
                    metric == "grad CBSR",
                    !is.na(value)) |>
      dplyr::group_by(estimand, sigma, lambda) |>
      dplyr::summarise(value = agg_f(value), .groups = "drop") |>
      # lock factor orders and ensure full grid so all λ are on the x-axis
      dplyr::mutate(
        estimand = factor(estimand, levels = dims()$estimand),
        sigma    = factor(sigma,    levels = dims()$sigma),
        lambda   = factor(lambda,   levels = dims()$lambda)
      ) |>
      tidyr::complete(
        estimand = factor(dims()$estimand, levels = dims()$estimand),
        sigma    = factor(dims()$sigma,    levels = dims()$sigma),
        lambda   = factor(dims()$lambda,   levels = dims()$lambda)
      ) |>
      # transform to log10 and apply optional cap (values above -> NA to break the line)
      dplyr::mutate(
        value = suppressWarnings(log10(value)),
        value = if (isTRUE(input$cap_cbsr_enable))
          ifelse(!is.finite(value) | value > input$cap_cbsr_value, NA_real_, value)
        else value
      )
    
    validate(need(nrow(dat) > 0, "No data to display for the current selection."))
    
    ggplot(dat, aes(x = lambda, y = value, color = sigma, group = sigma)) +
      geom_line(na.rm = TRUE) +
      geom_point(size = 1, na.rm = TRUE) +
      facet_wrap(vars(estimand), nrow = 1, drop = FALSE) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      scale_x_discrete(drop = FALSE) +
      labs(
        title = "log10(||grad CBSR||) vs λ",
        x = "lambda",
        y = "log10(||grad CBSR||)",
        color = "sigma"
      )
  })
  
  output$plot_smd_vs_cbsr <- renderPlot({
    req(input$n, input$prop_trt, input$confdg_lvl)
    
    # Use numeric aggregation; if "wins" is selected, fall back to mean
    agg_key <- if (identical(input$agg, "wins")) "mean" else input$agg
    agg_f <- switch(agg_key,
                    "mean"   = function(x) mean(x, na.rm = TRUE),
                    "median" = function(x) median(x, na.rm = TRUE),
                    "mode"   = stat_mode,
                    "sd"     = function(x) sd(x, na.rm = TRUE))
    
    # Aggregate current slice for both metrics, then pivot wide
    dat <- hp_df() |>
      dplyr::filter(n == input$n,
                    prop_trt == input$prop_trt,
                    confdg_lvl == input$confdg_lvl,
                    metric %in% c("max SMD", "grad CBSR"),
                    !is.na(value)) |>
      dplyr::group_by(estimand, sigma, lambda, metric) |>
      dplyr::summarise(value = agg_f(value), .groups = "drop") |>
      tidyr::pivot_wider(names_from = metric, values_from = value) |>
      # transform CBSR to log10 and optionally cap (remove high points)
      dplyr::mutate(
        log_cbsr = suppressWarnings(log10(`grad CBSR`))
      ) |>
      dplyr::filter(is.finite(log_cbsr)) |>
      # lock factor orders for facets / legend
      dplyr::mutate(
        estimand = factor(estimand, levels = dims()$estimand),
        sigma    = factor(sigma,    levels = dims()$sigma),
        lambda   = factor(lambda,   levels = dims()$lambda)
      )
    
    if (isTRUE(input$cap_cbsr_enable))
      dat = dplyr::filter(dat, log_cbsr <= input$cap_cbsr_value)
    
    validate(need(nrow(dat) > 0, "No points to display for the current selection."))
    
    ggplot(dat, aes(x = `max SMD`, y = log_cbsr, color = sigma)) +
      geom_point(alpha = 0.8) +
      facet_wrap(vars(estimand), nrow = 1, drop = FALSE) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(
        title = "max |SMD| vs log10(||grad CBSR||)",
        x = "max |SMD|",
        y = "log10(||grad CBSR||)",
        color = "sigma"
      )
  })
  
  output$plot_smd_vs_cbsr <- renderPlot({
    req(input$n, input$prop_trt, input$confdg_lvl)
    
    # Use numeric aggregation; if "wins" is selected, fall back to mean
    agg_key <- if (identical(input$agg, "wins")) "mean" else input$agg
    agg_f <- switch(agg_key,
                    "mean"   = function(x) mean(x, na.rm = TRUE),
                    "median" = function(x) median(x, na.rm = TRUE),
                    "mode"   = stat_mode,
                    "sd"     = function(x) sd(x, na.rm = TRUE))
    
    # Aggregate current slice for both metrics, then pivot wide
    dat <- hp_df() |>
      dplyr::filter(n == input$n,
                    prop_trt == input$prop_trt,
                    confdg_lvl == input$confdg_lvl,
                    metric %in% c("max SMD", "grad CBSR"),
                    !is.na(value)) |>
      dplyr::group_by(estimand, sigma, lambda, metric) |>
      dplyr::summarise(value = agg_f(value), .groups = "drop") |>
      tidyr::pivot_wider(names_from = metric, values_from = value) |>
      # transform CBSR to log10 and optionally cap (remove high points)
      dplyr::mutate(
        log_cbsr = suppressWarnings(log10(`grad CBSR`))
      ) |>
      dplyr::filter(is.finite(log_cbsr)) |>
      # lock factor orders for facets / legend
      dplyr::mutate(
        estimand = factor(estimand, levels = dims()$estimand),
        sigma    = factor(sigma,    levels = dims()$sigma),
        lambda   = factor(lambda,   levels = dims()$lambda)
      )
    
    if (isTRUE(input$cap_cbsr_enable))
      dat = dplyr::filter(dat, log_cbsr <= input$cap_cbsr_value)
    
    validate(need(nrow(dat) > 0, "No points to display for the current selection."))
    
    ggplot(dat, aes(x = `max SMD`, y = log_cbsr, color = sigma)) +
      geom_point(alpha = 0.8) +
      facet_wrap(vars(estimand), nrow = 1, drop = FALSE) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(
        title = "max |SMD| vs log10(||grad CBSR||)",
        x = "max |SMD|",
        y = "log10(||grad CBSR||)",
        color = "sigma"
      )
  })
}

shinyApp(ui, server)
