# app.R — minimal deps: shiny + ggplot2
options(shiny.maxRequestSize = 1024^3)

library(shiny)
library(ggplot2)

# --------------------------------------------------------------------------
# Unadjusted SMDs (same for ATE/ATT; indexed by prop_trt & confdg_lvl & variable)
smd_unadjusted <- array(
  data = c(  0.385,  0.363,  0.359,  0.571,  0.520,  0.518,  0.750,  0.669,  0.683,
             -0.260, -0.242, -0.244, -0.354, -0.331, -0.349, -0.426, -0.407, -0.448,
             0.181,  0.182,  0.187,  0.264,  0.260,  0.274,  0.332,  0.328,  0.360,
             -0.247, -0.238, -0.250, -0.343, -0.345, -0.367, -0.441, -0.437, -0.474,
             -0.191, -0.193, -0.202, -0.274, -0.271, -0.295, -0.345, -0.347, -0.391,
             -0.266, -0.255, -0.259, -0.378, -0.362, -0.383, -0.482, -0.461, -0.501,
             0.434,  0.428,  0.439,  0.621,  0.621,  0.658,  0.806,  0.810,  0.884,
             0.022,  0.024,  0.020,  0.033,  0.030,  0.030,  0.043,  0.043,  0.046,
             -0.177, -0.171, -0.180, -0.245, -0.247, -0.262, -0.318, -0.314, -0.340,
             -0.002,  0.003,  0.003, -0.006, -0.001,  0.001,  0.003,  0.002,  0.001),
  dim = c(10, 3, 3),
  dimnames = list(
    variable   = paste0("X", 1:10),
    prop_trt   = c("low", "moderate", "high"),
    confdg_lvl = c("low", "moderate", "high")
  )
)

metric_labs = c("bias" = "Bias",
                "variance" = "Variance",
                "mae" = "MAE",
                "rmse" = "RMSE",
                "CI_coverage" = "CI coverage",
                "CI_width" = "CI width",
                "fail" = "Failure rate",
                "SNegW" = "Sum of negative weights (min)",
                "smd" = "Standardized Mean Difference (SMD)",
                "ess" = "Effective Sample Size (ESS)")

`%||%` <- function(a, b) if (!is.null(a)) a else b

dim_levels <- function(x, nm) {
  if (is.null(dimnames(x))) return(NULL)
  dn <- dimnames(x); if (!is.null(names(dn)) && nm %in% names(dn)) dn[[nm]] else NULL
}

detect_third_dim <- function(arr) {
  if (is.null(dimnames(arr))) return(NULL)
  dnn <- names(dimnames(arr))
  cands <- c("sample_size", "confdg_lvl", "conf_lvl", "alpha")
  m <- cands[cands %in% dnn]
  if (length(m)) m[1] else NULL
}

slice_array <- function(arr, selectors = list()) {
  stopifnot(!is.null(dimnames(arr)))
  dnn <- names(dimnames(arr))
  idx <- rep(list(TRUE), length(dnn)); names(idx) <- dnn
  for (nm in names(selectors)) {
    if (nm %in% dnn && !is.null(selectors[[nm]])) {
      idx[[nm]] <- which(dimnames(arr)[[nm]] == selectors[[nm]])
    }
  }
  arr2 <- do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
  df <- as.data.frame(as.table(arr2), stringsAsFactors = FALSE)
  names(df)[ncol(df)] <- "value"
  
  if (!"method"     %in% names(df)) df$method    <- NA_character_
  if (!"estimand"   %in% names(df)) df$estimand  <- NA_character_
  if (!"estimator"  %in% names(df)) df$estimator <- "—"
  
  keep <- intersect(c("method","estimand","estimator","value"), names(df))
  df <- df[, keep, drop = FALSE]
  
  if ("method" %in% names(df))    df$method    <- factor(df$method,    levels = c("IPTW","KOM","EB","TLF"))
  if ("estimand" %in% names(df))  df$estimand  <- factor(df$estimand,  levels = c("ATE","ATT"))
  if ("estimator" %in% names(df)) df$estimator <- factor(df$estimator, levels = c("WLS","DR","—"))
  df
}

# ---- CLEANED table maker: no weird '— · ATE' / 'WLS · —' / '— · —' columns
wide_by_method <- function(df) {
  orig_has_estimator <- "estimator" %in% names(df) && any(df$estimator != "—")
  orig_has_estimand  <- "estimand"  %in% names(df) && any(df$estimand  != "—")
  
  if (!"estimator" %in% names(df)) df$estimator <- "—"
  if (!"estimand"  %in% names(df)) df$estimand  <- "—"
  
  if (orig_has_estimator && orig_has_estimand) {
    df$colkey <- paste(df$estimator, df$estimand, sep = " · ")
  } else if (!orig_has_estimator && orig_has_estimand) {
    df$colkey <- as.character(df$estimand)
  } else if (orig_has_estimator && !orig_has_estimand) {
    df$colkey <- as.character(df$estimator)
  } else {
    df$colkey <- "value"
  }
  
  mat <- tapply(df$value, list(df$method, df$colkey), function(x) x[1])
  tb  <- as.data.frame.matrix(mat)
  
  if (orig_has_estimator && orig_has_estimand) {
    est_levels   <- c("WLS","DR", setdiff(unique(as.character(df$estimator)), c("WLS","DR","—")))
    estim_levels <- c("ATE","ATT", setdiff(unique(as.character(df$estimand)),  c("ATE","ATT","—")))
    desired <- as.vector(outer(est_levels, estim_levels, paste, sep = " · "))
    desired <- intersect(desired, colnames(tb))
    if (length(desired)) tb <- tb[, desired, drop = FALSE]
  } else if (!orig_has_estimator && orig_has_estimand) {
    desired <- intersect(c("ATE","ATT"), colnames(tb))
    if (length(desired)) tb <- tb[, desired, drop = FALSE]
  } else if (orig_has_estimator && !orig_has_estimand) {
    desired <- intersect(c("WLS","DR"), colnames(tb))
    if (length(desired)) tb <- tb[, desired, drop = FALSE]
  }
  
  row_order <- c("IPTW","KOM","EB","TLF")
  tb <- tb[order(match(rownames(tb), row_order)), , drop = FALSE]
  tb
}

get_unadjusted_smd_vec <- function(prop_val, third_val, third_name, var_levels) {
  arr <- smd_unadjusted
  dn <- dimnames(arr); dnn <- names(dn)
  var_name   <- (c("variables","variable")[c("variables","variable") %in% dnn])[1]
  prop_name  <- "prop_trt"
  third_name_in_arr <- (c(third_name, "confdg_lvl", "sample_size", "conf_lvl", "alpha")[c(third_name, "confdg_lvl", "sample_size", "conf_lvl", "alpha") %in% dnn])[1]
  if (is.null(var_name) || is.null(prop_name %||% NULL) || is.null(third_name_in_arr))
    return(setNames(rep(NA_real_, length(var_levels)), var_levels))
  idx <- rep(list(TRUE), length(dnn)); names(idx) <- dnn
  if (prop_name %in% dnn)  idx[[prop_name]]        <- which(dn[[prop_name]]        == prop_val)
  if (third_name_in_arr %in% dnn) idx[[third_name_in_arr]] <- which(dn[[third_name_in_arr]] == third_val)
  arr2 <- do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
  df <- as.data.frame(as.table(arr2), stringsAsFactors = FALSE)
  names(df)[ncol(df)] <- "value"
  if (!("variables" %in% names(df))) names(df)[names(df) == var_name] <- "variables"
  v <- df$value; names(v) <- df$variables
  v[var_levels]
}

# Colors for methods
method_colors <- c("IPTW" = "#E74C3C", "KOM" = "#F1C40F", "EB" = "#2E86C1", "TLF" = "#229954")
# Shapes for complexity level
shape_map_default <- c("low" = 15, "moderate" = 17, "high" = 19)

# pretty label for estimator facet (already in your file)
est_lab <- function(x) setNames(
  ifelse(x=="WLS","Weighted average estimator",
         ifelse(x=="DR","Augmented estimator","Estimator: —")), x
)

# NEW: clear column strip labels for prop_trt
prop_lab <- function(x) setNames(paste("Prop. treated:", x), x)


# --------------------------------------------------------------------------
# UI

ui <- fluidPage(
  tags$head(tags$title("Balancing methods — metrics explorer")),
  titlePanel("Balancing methods — metrics explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("rdata", "Upload metrics.RData (precomputed metrics only)", accept = ".RData"),
      helpText("Upload the file that contains metrics arrays (bias, variance, rmse, ...). Raw data files are not supported."),
      tags$hr(),
      uiOutput("selectors_ui"),
      selectInput("metric", "Metric",
                  choices = c("Bias" = "bias",
                              "Variance" = "variance",
                              "MAE" = "mae",
                              "RMSE" = "rmse",
                              "CI coverage" = "CI_coverage",
                              "CI width" = "CI_width",
                              "Failure rate" = "fail",
                              "Sum of negative weights (min)" = "SNegW",
                              "Standardized Mean Difference (SMD)" = "smd",
                              "Effective Sample Size (ESS)" = "ess"),
                  selected = "bias"),
      checkboxGroupInput("methods", "Methods",
                         choices = c("IPTW","KOM","EB","TLF"),
                         selected = c("IPTW","KOM","EB","TLF")),
      checkboxGroupInput("estimand", "Estimand",
                         choices = c("ATE","ATT"),
                         selected = c("ATE","ATT")),
      checkboxGroupInput("estimator", "Estimator",
                         choices = c("WLS","DR"),
                         selected = c("WLS","DR"))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Plot",
          div(style = "margin-bottom: 10px; display:flex; gap:10px; align-items:center;",
              uiOutput("dl_btns")   # <- buttons appear here (ATE / ATT / SMD / ESS)
          ),
          uiOutput("plot_ui")
        ),
        tabPanel("Table",
                 helpText("Rows = methods. Columns = estimator · estimand. Rounded to 4 decimals."),
                 tableOutput("metric_table"))
      )
    )
  )
)

# --------------------------------------------------------------------------
# Server

server <- function(input, output, session) {
  store <- reactiveValues(env = NULL, loaded = FALSE, dims = NULL)
  
  # Expected metrics and pretty labels for the selectInput
  metric_choices_all <- c("Bias" = "bias",
                          "Variance" = "variance",
                          "MAE" = "mae",
                          "RMSE" = "rmse",
                          "CI coverage" = "CI_coverage",
                          "CI width" = "CI_width",
                          "Failure rate" = "fail",
                          "Sum of negative weights (min)" = "SNegW",
                          "Standardized Mean Difference (SMD)" = "smd",
                          "Effective Sample Size (ESS)" = "ess")
  
  observeEvent(input$rdata, ignoreInit = TRUE, {
    req(input$rdata$datapath)
    env <- new.env(parent = emptyenv())
    loaded_names <- load(input$rdata$datapath, envir = env)
    
    # Must contain at least one known metric
    has_any_metric <- any(metric_choices_all %in% loaded_names)
    validate(need(has_any_metric, "No metrics found. Please upload a metrics.RData file produced by your pipeline."))
    
    # Update available metrics dropdown to only what exists in the file
    available <- metric_choices_all[metric_choices_all %in% loaded_names]
    updateSelectInput(session, "metric", choices = available, selected = available[1])
    
    # Pick a reference object for levels
    ref_name <- available[1]
    ref_arr  <- get(ref_name, envir = env)
    
    n_levels  <- dim_levels(ref_arr, "n")
    p_levels  <- dim_levels(ref_arr, "prop_trt")
    third_nm  <- detect_third_dim(ref_arr)
    third_lev <- if (!is.null(third_nm)) dim_levels(ref_arr, third_nm) else NULL
    smd_vars  <- if ("smd" %in% loaded_names) dim_levels(get("smd", envir = env), "variables") else NULL
    ess_groups<- if ("ess" %in% loaded_names) dim_levels(get("ess", envir = env), "group") else NULL
    
    store$env   <- env
    store$loaded<- TRUE
    store$dims  <- list(
      n = n_levels, prop_trt = p_levels, third_name = third_nm, third_levels = third_lev,
      smd_vars = smd_vars, ess_groups = ess_groups
    )
  })
  
  # ------------------ Selectors (n / prop_trt / complexity shown when needed) -----------------
  output$selectors_ui <- renderUI({
    req(store$loaded, store$dims)
    d <- store$dims
    cond <- "(['smd','ess'].indexOf(input.metric) >= 0) || (input.tabs === 'Table')"
    
    ui <- tagList(
      conditionalPanel(
        condition = cond,
        selectInput("n", "n", choices = d$n, selected = if (length(d$n)) d$n[1] else NULL),
        selectInput("prop_trt", "prop_trt", choices = d$prop_trt, selected = if (length(d$prop_trt)) d$prop_trt[1] else NULL),
        if (!is.null(d$third_levels)) {
          selectInput("third", d$third_name, choices = d$third_levels, selected = d$third_levels[1])
        }
      ),
      conditionalPanel(
        condition = "input.metric == 'ess'",
        selectInput(
          "ess_group", "ESS group",
          choices = store$dims$ess_groups %||% character(0),
          selected = if (!is.null(store$dims$ess_groups) && length(store$dims$ess_groups)) store$dims$ess_groups[1] else NULL
        )
      )
    )
    ui
  })
  
  # For table and ESS we still slice by selectors
  current_df <- reactive({
    req(store$loaded)
    env <- store$env
    metric <- input$metric
    validate(need(metric %in% ls(env), sprintf("Metric '%s' not found in the uploaded file.", metric)))
    arr <- get(metric, envir = env)
    selectors <- list(n = req(input$n), prop_trt = req(input$prop_trt))
    if (!is.null(store$dims$third_levels)) selectors[[store$dims$third_name]] <- req(input$third)
    if (metric == "ess") selectors[["group"]] <- req(input$ess_group)
    if (metric == "smd") stop("SMD table uses aggregated median |SMD|; handled below.")
    df <- slice_array(arr, selectors)
    if ("method" %in% names(df)    && length(input$methods))   df <- df[df$method %in% input$methods, , drop = FALSE]
    if ("estimand" %in% names(df)  && length(input$estimand))  df <- df[df$estimand %in% input$estimand, , drop = FALSE]
    if ("estimator" %in% names(df) && any(df$estimator != "—") && length(input$estimator)) {
      df <- df[df$estimator %in% input$estimator | df$estimator == "—", , drop = FALSE]
    }
    df
  })
  
  # SMD table: median |SMD| across covariates
  smd_table_df <- reactive({
    req(store$loaded, "smd" %in% ls(store$env))
    env <- store$env
    arr <- get("smd", envir = env)
    dnn <- names(dimnames(arr))
    idx <- rep(list(TRUE), length(dnn)); names(idx) <- dnn
    idx[["n"]]        <- which(dimnames(arr)[["n"]] == req(input$n))
    idx[["prop_trt"]] <- which(dimnames(arr)[["prop_trt"]] == req(input$prop_trt))
    if (!is.null(store$dims$third_levels)) idx[[store$dims$third_name]] <- which(dimnames(arr)[[store$dims$third_name]] == req(input$third))
    arr2 <- do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
    df2 <- as.data.frame(as.table(arr2), stringsAsFactors = FALSE)
    names(df2)[ncol(df2)] <- "smd"
    if (length(input$methods))  df2 <- df2[df2$method %in% input$methods, , drop = FALSE]
    if (length(input$estimand)) df2 <- df2[df2$estimand %in% input$estimand, , drop = FALSE]
    agg <- aggregate(abs(smd) ~ method + estimand, df2, function(z) median(z, na.rm = TRUE))
    names(agg)[names(agg) == "abs(smd)"] <- "value"
    agg$estimator <- "—"
    agg$method    <- factor(agg$method,    levels = c("IPTW","KOM","EB","TLF"))
    agg$estimand  <- factor(agg$estimand,  levels = c("ATE","ATT"))
    agg$estimator <- factor(agg$estimator, levels = c("WLS","DR","—"))
    agg[, c("method","estimand","estimator","value")]
  })
  
  # SMD plot data (signed; add Unadjusted)
  smd_plot_df <- reactive({
    req(store$loaded)
    env <- store$env
    validate(need("smd" %in% ls(env), "SMD metric not found."))
    arr <- get("smd", envir = env)
    dnn <- names(dimnames(arr))
    idx <- rep(list(TRUE), length(dnn)); names(idx) <- dnn
    idx[["n"]]        <- which(dimnames(arr)[["n"]] == req(input$n))
    idx[["prop_trt"]] <- which(dimnames(arr)[["prop_trt"]] == req(input$prop_trt))
    if (!is.null(store$dims$third_levels)) {
      nm <- store$dims$third_name
      idx[[nm]] <- which(dimnames(arr)[[nm]] == req(input$third))
    }
    arr2 <- do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
    df <- as.data.frame(as.table(arr2), stringsAsFactors = FALSE)
    names(df)[ncol(df)] <- "smd"
    if (length(input$methods))  df <- df[df$method %in% input$methods, , drop = FALSE]
    if (length(input$estimand)) df <- df[df$estimand %in% input$estimand, , drop = FALSE]
    var_levels <- dimnames(arr2)$variables
    third_val <- if (!is.null(store$dims$third_levels)) req(input$third) else (dimnames(arr2)[[detect_third_dim(arr2)]] %||% NULL)
    unadj_vec <- get_unadjusted_smd_vec(prop_val = req(input$prop_trt),
                                        third_val = third_val,
                                        third_name = store$dims$third_name,
                                        var_levels = var_levels)
    mk_unadj <- function(est) {
      data.frame(variables = var_levels,
                 method = "Unadjusted",
                 estimand = est,
                 smd = as.numeric(unadj_vec),
                 stringsAsFactors = FALSE)
    }
    estims <- if (length(input$estimand)) input$estimand else dimnames(arr2)$estimand
    unadj  <- do.call(rbind, lapply(estims, mk_unadj))
    out <- rbind(unadj[, c("variables","method","estimand","smd")],
                 df[,     c("variables","method","estimand","smd")])
    out$variables <- factor(out$variables, levels = var_levels)
    out$method    <- factor(out$method, levels = c("Unadjusted","IPTW","KOM","EB","TLF"))
    out
  })
  
  # ------- PLOTS ------------------------------------------------------------
  non_smd_df <- reactive({
    req(store$loaded, input$metric)
    if (input$metric %in% c("smd", "ess")) return(NULL)
    env <- store$env
    arr <- get(input$metric, envir = env)
    dft <- as.data.frame(as.table(arr), stringsAsFactors = FALSE)
    names(dft)[ncol(dft)] <- "value"
    dft$value <- ifelse(abs(dft$value) > 1, NA, dft$value)
    if ("method" %in% names(dft)    && length(input$methods))   dft <- dft[dft$method %in% input$methods, , drop = FALSE]
    if ("estimand" %in% names(dft)  && length(input$estimand))  dft <- dft[dft$estimand %in% input$estimand, , drop = FALSE]
    if ("estimator" %in% names(dft) && any(dft$estimator != "—") && length(input$estimator)) {
      dft <- dft[dft$estimator %in% input$estimator | dft$estimator == "—", , drop = FALSE]
    }
    n_lvls <- dim_levels(arr, "n") %||% unique(dft$n)
    if (!is.null(n_lvls)) {
      n_num <- suppressWarnings(as.numeric(n_lvls))
      if (!anyNA(n_num)) n_lvls <- as.character(sort(n_num))
      dft$n <- factor(dft$n, levels = n_lvls)
    }
    p_lvls <- dim_levels(arr, "prop_trt") %||% unique(dft$prop_trt)
    if (!is.null(p_lvls)) dft$prop_trt <- factor(dft$prop_trt, levels = p_lvls)
    c_lvls <- dim_levels(arr, "confdg_lvl") %||% unique(dft$confdg_lvl)
    if (!is.null(c_lvls)) dft$confdg_lvl <- factor(dft$confdg_lvl, levels = c_lvls)
    if ("method" %in% names(dft)) dft$method <- factor(dft$method, levels = names(method_colors))
    if ("estimator" %in% names(dft)) {
      dft$estimator <- factor(dft$estimator, levels = c("WLS","DR","—"))
    } else {
      dft$estimator <- factor("—", levels = c("WLS","DR","—"))
    }
    dft
  })
  
  shape_map_for <- function(levels_present) {
    base <- setNames(rep(16, length(levels_present)), levels_present)
    look <- intersect(names(shape_map_default), names(base))
    base[look] <- shape_map_default[look]
    base
  }
  
  est_lab <- function(x) setNames(
    ifelse(x=="WLS","WLS estimator",
           ifelse(x=="DR","DR estimator","Estimator: —")), x
  )
  
  make_non_smd_plot <- function(dft_sub, ylab_txt, title_txt,
                                size = 2.8, legend_size = size) {
    validate(need(nrow(dft_sub) > 0, paste("No data for", title_txt)))
    sh_map <- shape_map_for(levels(dft_sub$confdg_lvl))
    ymin <- suppressWarnings(min(dft_sub$value, na.rm = TRUE))
    ymax <- suppressWarnings(max(dft_sub$value, na.rm = TRUE))
    y0   <- if (is.finite(ymin)) min(0, ymin) else 0
    y1   <- if (is.finite(ymax)) ymax * 1.1 else 1
    if (input$metric == "CI_coverage") y1 <- 1
    
    p <- ggplot(dft_sub, aes(x = n, y = value, color = method, shape = confdg_lvl)) +
      { if (input$metric == "CI_coverage")
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") else NULL } +
      geom_point(position = position_dodge(width = 0.6), size = size, na.rm = TRUE) +
      scale_color_manual(values = method_colors, name = "Method", drop = FALSE) +
      scale_shape_manual(values = sh_map, name = "Level of\ncomplexity", drop = FALSE) +
      labs(x = "sample size (n)", y = ylab_txt, title = title_txt) +
      facet_grid(rows = vars(estimator), cols = vars(prop_trt),
                 labeller = labeller(estimator = est_lab, prop_trt = prop_lab)) +
      theme_minimal(base_size = 13) +
      theme(
        panel.background = element_rect(fill = "#f9f9f9", color = NA),
        panel.grid.major.x = element_line(color = "#e2e2e2", size = 0.4),
        panel.grid.major.y = element_line(color = "#e2e2e2", size = 0.4),
        panel.grid.minor = element_blank()
      ) +
      coord_cartesian(ylim = c(y0, y1)) +
      guides(
        color = guide_legend(override.aes = list(size = legend_size)),
        shape = guide_legend(override.aes = list(size = legend_size))
      )
    p
  }
  
  # Dynamic plot area: SMD / ESS / two plots (ATE + ATT)
  output$plot_ui <- renderUI({
    req(input$metric)
    if (input$metric == "smd") {
      plotOutput("metric_plot_smd", height = 520)
    } else if (input$metric == "ess") {
      plotOutput("metric_plot_ess", height = 520)
    } else {
      elems <- list()
      if ("ATE" %in% input$estimand) {
        elems <- c(elems, list(tags$h4("ATE"), plotOutput("metric_plot_ATE", height = 420)))
      }
      if ("ATT" %in% input$estimand) {
        elems <- c(elems, list(tags$h4("ATT"), plotOutput("metric_plot_ATT", height = 420)))
      }
      do.call(tagList, elems)
    }
  })
  
  # SMD plot (signed; colored; no lines)
  output$metric_plot_smd <- renderPlot({
    req(input$metric == "smd")
    df <- smd_plot_df()
    validate(need(nrow(df) > 0, "No SMD data to plot for the current selection."))
    smd_colors <- c("Unadjusted" = "black", method_colors)
    p <- ggplot(df, aes(x = smd, y = variables, color = method)) +
      geom_point(shape = 16, size = 2.8, na.rm = TRUE) +
      scale_color_manual(values = smd_colors, breaks = names(smd_colors), name = "Method") +
      labs(x = "SMD (signed)", y = "Covariate", title = "Standardized Mean Difference by covariate") +
      theme_minimal(base_size = 13) +
      theme(panel.background = element_rect(fill = "#f9f9f9", color = NA))
    if (length(unique(df$estimand)) > 1) p <- p + facet_wrap(~ estimand, ncol = 1)
    p
  })
  
  # ESS plot (bars)
  output$metric_plot_ess <- renderPlot({
    req(input$metric == "ess")
    df <- current_df()
    validate(need(nrow(df) > 0, "No data to plot for the current selection."))
    ggplot(df, aes(x = method, y = value)) +
      geom_col(aes(fill = estimand), position = position_dodge(width = 0.7)) +
      scale_fill_discrete(name = "Estimand") +
      labs(x = "Method", y = "ESS", title = "Effective Sample Size") +
      theme_minimal(base_size = 13) +
      theme(panel.background = element_rect(fill = "#f9f9f9", color = NA))
  })
  
  # ATE plot
  output$metric_plot_ATE <- renderPlot({
    req(input$metric %in% c("bias","variance","mae","rmse","CI_coverage","CI_width","fail","SNegW"),
        "ATE" %in% input$estimand)
    dft <- non_smd_df()
    validate(need(!is.null(dft), "No data to plot."))
    dft_ate <- subset(dft, estimand == "ATE")
    metric_lab = metric_labs[input$metric]
    make_non_smd_plot(dft_ate, ylab_txt = metric_lab, title_txt = paste("Metric:", metric_lab, "— ATE"))
  })
  
  # ATT plot
  output$metric_plot_ATT <- renderPlot({
    req(input$metric %in% c("bias","variance","mae","rmse","CI_coverage","CI_width","fail","SNegW"),
        "ATT" %in% input$estimand)
    dft <- non_smd_df()
    validate(need(!is.null(dft), "No data to plot."))
    dft_att <- subset(dft, estimand == "ATT")
    metric_lab = metric_labs[input$metric]
    make_non_smd_plot(dft_att, ylab_txt = metric_lab, title_txt = paste("Metric:", metric_lab, "— ATT"))
  })
  
  # ------------------------------- TABLE ------------------------------------
  output$metric_table <- renderTable({
    if (input$metric == "smd") {
      df <- smd_table_df()
    } else {
      df <- current_df()
    }
    validate(need(nrow(df) > 0, "No data to display for the current selection."))
    tb <- wide_by_method(df)
    tb[] <- lapply(tb, function(x) if (is.numeric(x)) round(x, 4) else x)
    data.frame(Method = rownames(tb), tb, check.names = FALSE, row.names = NULL)
  }, striped = TRUE, bordered = TRUE, width = "100%", align = "l")
  
  
  # Convert a rendered plotOutput's pixel size to inches for the PDF device
  px_to_in <- function(px, fallback_px) ( (px %||% fallback_px) / 109 )  # assume 96 dpi screen
  
  get_plot_size_in <- function(id, fallback_w = 900, fallback_h = 420) {
    w_px <- session$clientData[[paste0("output_", id, "_width")]]
    h_px <- session$clientData[[paste0("output_", id, "_height")]]
    list(width  = px_to_in(w_px, fallback_w),
         height = px_to_in(h_px, fallback_h))
  }
  
  # Remove title for exported plots (keep everything else identical)
  drop_title <- function(p) p + labs(title = NULL)
  
  output$dl_btns <- renderUI({
    req(input$metric)
    btns <- list()
    
    if (input$metric %in% c("bias","variance","mae","rmse","CI_coverage","CI_width","fail","SNegW")) {
      if ("ATE" %in% input$estimand) btns <- c(btns, list(downloadButton("dl_pdf_ate", "Download ATE (PDF)")))
      if ("ATT" %in% input$estimand) btns <- c(btns, list(downloadButton("dl_pdf_att", "Download ATT (PDF)")))
    } else if (input$metric == "smd") {
      btns <- list(downloadButton("dl_pdf_smd", "Download SMD (PDF)"))
    } else if (input$metric == "ess") {
      btns <- list(downloadButton("dl_pdf_ess", "Download ESS (PDF)"))
    }
    
    do.call(tagList, btns)
  })
  
  # ---- ATE ----
  output$dl_pdf_ate <- downloadHandler(
    filename = function() paste0("plot_", input$metric, "_ATE.pdf"),
    content = function(file) {
      dft <- non_smd_df(); validate(need(!is.null(dft), "No data to plot."))
      dft_ate <- subset(dft, estimand == "ATE")
      metric_lab <- metric_labs[input$metric]
      p <- make_non_smd_plot(
        dft_ate, ylab_txt = metric_lab, title_txt = paste("Metric:", metric_lab, "— ATE"),
        size = 1.5,            # smaller on the plot
        legend_size = 2.8      # normal size in the legends
      )
      
      sz <- get_plot_size_in("metric_plot_ATE", fallback_w = 1100, fallback_h = 420)
      grDevices::cairo_pdf(filename = file, width = sz$width, height = sz$height)  # note: arg is 'file', not 'filename'
      on.exit(grDevices::dev.off(), add = TRUE)
      print(drop_title(p))
    }
  )
  
  # ---- ATT ----
  output$dl_pdf_att <- downloadHandler(
    filename = function() paste0("plot_", input$metric, "_ATT.pdf"),
    content = function(file) {
      dft <- non_smd_df(); validate(need(!is.null(dft), "No data to plot."))
      dft_att <- subset(dft, estimand == "ATT")
      metric_lab <- metric_labs[input$metric]
      p <- make_non_smd_plot(
        dft_att, ylab_txt = metric_lab, title_txt = paste("Metric:", metric_lab, "— ATT"),
        size = 1.5, legend_size = 2.8
      )
      
      sz <- get_plot_size_in("metric_plot_ATT", fallback_w = 1100, fallback_h = 420)
      grDevices::cairo_pdf(filename = file, width = sz$width, height = sz$height)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(drop_title(p))
    }
  )
  
}

shinyApp(ui, server)