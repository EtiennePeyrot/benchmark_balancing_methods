# Shiny metrics explorer for metrics.RData output

library(shiny)
library(ggplot2)
library(grid)

# Fixed display defaults for the user-facing app.
variance_ratio_plot_max_default = 6
missingness_asterisk_threshold_default = 5

# --------------------------------------------------------------------------
# Unadjusted SMDs (same for ATE/ATT; indexed by prop_trt & confdg_lvl & variable)
smd_unadjusted = array(
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
                "variance_ratio" = "Variance ratio",
                "missing" = "Missingness",
                "SNegW" = "Minimum sum of negative weights",
                "smd" = "Standardized Mean Difference (SMD)",
                "ess" = "Effective Sample Size (ESS)")

metric_choices_all = c("Bias" = "bias",
                       "Variance" = "variance",
                       "MAE" = "mae",
                       "RMSE" = "rmse",
                       "CI coverage" = "CI_coverage",
                       "Variance ratio" = "variance_ratio",
                       "Missingness" = "missing",
                       "Sum of negative weights (min)" = "SNegW",
                       "Standardized Mean Difference (SMD)" = "smd",
                       "Effective Sample Size (ESS)" = "ess")

metric_mcse_names = c("bias" = "bias_mcse",
                      "variance" = "variance_mcse",
                      "mae" = "mae_mcse",
                      "rmse" = "rmse_mcse",
                      "CI_coverage" = "CI_coverage_mcse",
                      "variance_ratio" = "variance_ratio_mcse")

missingness_type_labs = c("estimate" = "Point estimate",
                          "standard_error" = "Variance estimate")

default_method_colors = c("IPTW" = "#4AA19E",
                          "KOM" = "#F1C40F",
                          "KOM_none" = "#F1C40F",
                          "KOM_small" = "#F1C40F",
                          "KOM_high" = "#F1C40F",
                          "EB" = "#4153FA",
                          "CBPS-TLF" = "#229954",
                          "CBPS-JI" = "#9D09AB",
                          "CBPS-OI" = "#D35400")

# Shapes for confounding level.
shape_map_default = c("low" = 15, "moderate" = 17, "high" = 19)

`%||%` = function(a, b) if (!is.null(a)) a else b

fmt_num = function(x, digits = 4) {
  out = rep(NA_character_, length(x))
  ok = !is.na(x)
  out[ok] = formatC(x[ok], format = "f", digits = digits)
  out
}

dim_levels = function(x, nm) {
  if (!is.array(x) || is.null(dimnames(x))) return(NULL)
  dn = dimnames(x)
  if (!is.null(names(dn)) && nm %in% names(dn)) dn[[nm]] else NULL
}

has_dim = function(arr, nm) {
  is.array(arr) && !is.null(dimnames(arr)) && nm %in% names(dimnames(arr))
}

detect_third_dim = function(arr) {
  if (!is.array(arr) || is.null(dimnames(arr))) return(NULL)
  dnn = names(dimnames(arr))
  cands = c("confdg_lvl", "sample_size", "conf_lvl", "alpha")
  m = cands[cands %in% dnn]
  if (length(m)) m[1] else NULL
}

array_to_df = function(arr, value_name = "value") {
  df = as.data.frame(as.table(arr), stringsAsFactors = FALSE)
  names(df)[ncol(df)] = value_name
  df
}

slice_array_df = function(arr, selectors = list(), value_name = "value") {
  stopifnot(!is.null(dimnames(arr)))
  dnn = names(dimnames(arr))
  idx = rep(list(TRUE), length(dnn)); names(idx) = dnn
  
  for (nm in names(selectors)) {
    if (nm %in% dnn && !is.null(selectors[[nm]])) {
      idx[[nm]] = which(dimnames(arr)[[nm]] == selectors[[nm]])
    }
  }
  
  arr2 = do.call(`[`, c(list(arr), idx, list(drop = FALSE)))
  array_to_df(arr2, value_name = value_name)
}

find_reference_array = function(env, metric_names = unname(metric_choices_all)) {
  for (nm in metric_names) {
    if (exists(nm, envir = env, inherits = FALSE)) {
      obj = get(nm, envir = env)
      if (is.array(obj)) return(obj)
      if (is.list(obj)) {
        for (elt in obj) if (is.array(elt)) return(elt)
      }
    }
  }
  NULL
}

metric_reference_array = function(env, metric) {
  if (is.null(env) || is.null(metric)) return(NULL)
  
  if (metric == "missing") {
    if (!exists("missing", envir = env, inherits = FALSE)) return(NULL)
    miss = get("missing", envir = env)
    if (!is.list(miss)) return(NULL)
    for (elt in miss) if (is.array(elt)) return(elt)
    return(NULL)
  }
  
  if (!exists(metric, envir = env, inherits = FALSE)) return(NULL)
  obj = get(metric, envir = env)
  if (is.array(obj)) return(obj)
  if (is.list(obj)) {
    for (elt in obj) if (is.array(elt)) return(elt)
  }
  NULL
}

metric_method_names_raw = function(env, metric) {
  arr = metric_reference_array(env, metric)
  dim_levels(arr, "method") %||% character(0)
}

metric_effect_levels_raw = function(env, metric) {
  arr = metric_reference_array(env, metric)
  dim_levels(arr, "effect") %||% character(0)
}

kom_effect_suffixes = function(method_names) {
  out = sub("^KOM_", "", grep("^KOM_", method_names, value = TRUE))
  unique(out[nzchar(out)])
}

method_display_name = function(method_names) {
  out = as.character(method_names)
  out[grepl("^KOM_(none|small|high)$", out)] = "KOM"
  out
}

filter_effect_specific_kom_rows = function(df, arr, effect) {
  if (!"method" %in% names(df)) return(df)
  if (is.null(arr) || has_dim(arr, "effect")) return(df)
  if (is.null(effect) || !nzchar(effect)) return(df)
  
  raw_method = as.character(df$method)
  if (!any(grepl("^KOM_", raw_method))) return(df)
  
  keep = !grepl("^KOM_", raw_method) | raw_method == paste0("KOM_", effect)
  df[keep, , drop = FALSE]
}

global_effect_levels = function(env, metric_names) {
  out = character(0)
  
  for (metric in metric_names) {
    out = c(out, metric_effect_levels_raw(env, metric))
    out = c(out, kom_effect_suffixes(metric_method_names_raw(env, metric)))
  }
  
  preferred = c("none", "small", "high")
  out = unique(out[nzchar(out)])
  c(intersect(preferred, out), setdiff(out, preferred))
}

method_names_for_metric_and_effect = function(env, metric, effect = NULL) {
  methods = metric_method_names_raw(env, metric)
  if (!length(methods)) return(methods)
  
  arr = metric_reference_array(env, metric)
  metric_has_effect_dim = has_dim(arr, "effect")
  
  # If the metric has no effect dimension but has KOM_none/KOM_small/KOM_high,
  # use the selected effect to keep only the corresponding KOM copy.
  if (!metric_has_effect_dim && !is.null(effect) && nzchar(effect) && any(grepl("^KOM_", methods))) {
    keep = !grepl("^KOM_", methods) | methods == paste0("KOM_", effect)
    methods = methods[keep]
  }
  
  # The UI and plots should show the selected KOM copy as KOM.
  unique(method_display_name(methods))
}

assign_method_colors = function(method_names) {
  out = default_method_colors[intersect(names(default_method_colors), method_names)]
  missing_methods = setdiff(method_names, names(out))
  if (length(missing_methods)) {
    fallback = grDevices::rainbow(length(missing_methods), s = 0.65, v = 0.75)
    names(fallback) = missing_methods
    out = c(out, fallback)
  }
  out[method_names]
}

shape_map_for = function(levels_present) {
  levels_present = as.character(levels_present)
  base = setNames(rep(16, length(levels_present)), levels_present)
  look = intersect(names(shape_map_default), names(base))
  base[look] = shape_map_default[look]
  base
}

est_lab = function(x) setNames(
  ifelse(x == "WLS", "WLS estimator",
         ifelse(x == "DR", "DR estimator", "Estimator: —")), x
)

prop_lab = function(x) setNames(paste("Prop. treated:", x), x)

missing_type_lab = function(x) setNames(missingness_type_labs[x] %||% x, x)

safe_filename_part = function(x) {
  x = ifelse(is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x), "all", x)
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

# --------------------------------------------------------------------------
# UI

ui = fluidPage(
  tags$head(tags$title("Balancing methods — metrics explorer")),
  titlePanel("Weighting methods — metrics explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("rdata", "Upload metrics.RData (precomputed metrics only)", accept = ".RData"),
      helpText("Upload the file that contains precomputed metric arrays. Raw simulation files are not supported."),
      tags$hr(),
      selectInput("metric", "Metric", choices = metric_choices_all, selected = "bias"),
      uiOutput("effect_ui"),
      tags$hr(),
      tags$strong("Scenario selectors"),
      helpText("Used by the table. Plots keep n, treatment prevalence, and complexity visible on the axes/facets."),
      uiOutput("selectors_ui"),
      tags$hr(),
      checkboxGroupInput("methods", "Methods", choices = character(0), selected = character(0))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Plot",
          uiOutput("plot_ui")
        ),
        tabPanel(
          "Table",
          helpText("Rows = methods. Columns = estimator · estimand when those dimensions exist."),
          tableOutput("metric_table"),
          uiOutput("missingness_table_ui")
        )
      )
    )
  )
)

# --------------------------------------------------------------------------
# Server

server = function(input, output, session) {
  store = reactiveValues(
    env = NULL,
    loaded = FALSE,
    dims = NULL,
    available = NULL,
    effect_levels = character(0)
  )
  
  fixed_num = function(default) default
  fixed_bool = function(default) default
  
  variance_ratio_plot_max = reactive({
    x = fixed_num(variance_ratio_plot_max_default)
    if (!is.finite(x) || x <= 0) variance_ratio_plot_max_default else x
  })
  
  missingness_asterisk_threshold = reactive({
    x = fixed_num(missingness_asterisk_threshold_default)
    if (!is.finite(x) || x < 0) x = missingness_asterisk_threshold_default
    x / 100
  })
  
  format_pct = function(x) {
    paste0(formatC(100 * x, format = "f", digits = 1), "%")
  }
  
  table_digits = reactive({
    max(0, min(8, round(fixed_num(4))))
  })
  
  plot_theme_settings = reactive({
    list(
      margin_top = fixed_num(5.5),
      margin_right = fixed_num(5.5),
      margin_bottom = fixed_num(5.5),
      margin_left = fixed_num(5.5),
      panel_spacing_x = fixed_num(0.8),
      panel_spacing_y = fixed_num(0.8)
    )
  })
  
  spacing_settings = reactive({
    list(
      sample_size_spacing = fixed_num(2),
      method_spacing = fixed_num(0.25),
      confdg_spacing = fixed_num(0.05)
    )
  })
  
  mcse_settings = reactive({
    list(
      show = fixed_bool(TRUE),
      color = "#333333",
      width = fixed_num(0.5),
      alpha = fixed_num(1)
    )
  })
  
  selected_effect = reactive({
    if (!store$loaded || !length(store$effect_levels)) return(NULL)
    x = input$effect
    if (is.null(x) || !nzchar(x)) return(store$effect_levels[1])
    if (x %in% store$effect_levels) x else store$effect_levels[1]
  })
  
  current_method_names = reactive({
    req(store$loaded, input$metric)
    method_names_for_metric_and_effect(store$env, input$metric, selected_effect())
  })
  
  method_colors_current = reactive({
    req(store$loaded, input$metric)
    assign_method_colors(current_method_names())
  })
  
  output$effect_ui = renderUI({
    if (!store$loaded || !length(store$effect_levels)) return(NULL)
    selectInput(
      "effect",
      "Effect",
      choices = store$effect_levels,
      selected = store$effect_levels[1]
    )
  })
  
  observeEvent(input$rdata, ignoreInit = TRUE, {
    req(input$rdata$datapath)
    
    env = new.env(parent = emptyenv())
    loaded_names = load(input$rdata$datapath, envir = env)
    
    available = metric_choices_all[unname(metric_choices_all) %in% loaded_names]
    validate(need(length(available) > 0, "No supported metrics found. Please upload the current metrics.RData file."))
    
    ref_arr = find_reference_array(env, unname(available))
    validate(need(!is.null(ref_arr), "No metric array with dimensions was found in the uploaded file."))
    
    n_levels = dim_levels(ref_arr, "n")
    p_levels = dim_levels(ref_arr, "prop_trt")
    third_nm = detect_third_dim(ref_arr)
    third_lev = if (!is.null(third_nm)) dim_levels(ref_arr, third_nm) else NULL
    smd_vars = if (exists("smd", envir = env, inherits = FALSE)) dim_levels(get("smd", envir = env), "variables") else NULL
    ess_groups = if (exists("ess", envir = env, inherits = FALSE)) dim_levels(get("ess", envir = env), "group") else NULL
    effect_levels = global_effect_levels(env, unname(available))
    
    store$env = env
    store$loaded = TRUE
    store$available = available
    store$effect_levels = effect_levels
    store$dims = list(
      n = n_levels,
      prop_trt = p_levels,
      third_name = third_nm,
      third_levels = third_lev,
      smd_vars = smd_vars,
      ess_groups = ess_groups
    )
    
    first_metric = unname(available[1])
    first_methods = method_names_for_metric_and_effect(env, first_metric, if (length(effect_levels)) effect_levels[1] else NULL)
    
    validate(need(length(first_methods) > 0, "No method dimension was found in the selected metric."))
    
    updateSelectInput(session, "metric", choices = available, selected = first_metric)
    if (length(effect_levels)) updateSelectInput(session, "effect", choices = effect_levels, selected = effect_levels[1])
    
    updateCheckboxGroupInput(
      session,
      "methods",
      choices = first_methods,
      selected = first_methods
    )
  })
  
  observeEvent(list(input$metric, input$effect), {
    req(store$loaded, input$metric)
    
    methods = current_method_names()
    old_selected = isolate(input$methods)
    selected = intersect(old_selected, methods)
    if (!length(selected)) selected = methods
    
    updateCheckboxGroupInput(
      session,
      "methods",
      choices = methods,
      selected = selected
    )
  }, ignoreInit = FALSE)
  
  output$selectors_ui = renderUI({
    req(store$loaded, store$dims)
    d = store$dims
    
    tagList(
      selectInput("n", "Sample size", choices = d$n, selected = if (length(d$n)) d$n[1] else NULL),
      selectInput("prop_trt", "Treatment prevalence", choices = d$prop_trt, selected = if (length(d$prop_trt)) d$prop_trt[1] else NULL),
      if (!is.null(d$third_levels)) {
        selectInput("third", "Level of complexity", choices = d$third_levels, selected = d$third_levels[1])
      }
    )
  })
  
  factorize_common = function(df, arr = NULL) {
    if (!"method" %in% names(df)) df$method = NA_character_
    if (!"effect" %in% names(df)) df$effect = "—"
    if (!"estimand" %in% names(df)) df$estimand = "—"
    if (!"estimator" %in% names(df)) df$estimator = "—"
    if (!"confdg_lvl" %in% names(df)) df$confdg_lvl = "—"
    if (!"prop_trt" %in% names(df)) df$prop_trt = "—"
    if (!"n" %in% names(df)) df$n = "—"
    
    # For non-effect metrics, the raw data may contain KOM_none, KOM_small,
    # and KOM_high because KOM tuning depends on the DGM effect. The selected
    # effect chooses the raw row, but the displayed method remains KOM.
    df = filter_effect_specific_kom_rows(df, arr, selected_effect())
    df$method = method_display_name(df$method)
    
    n_lvls = if (!is.null(arr)) dim_levels(arr, "n") else NULL
    n_lvls = n_lvls %||% store$dims$n %||% unique(as.character(df$n))
    
    p_lvls = if (!is.null(arr)) dim_levels(arr, "prop_trt") else NULL
    p_lvls = p_lvls %||% store$dims$prop_trt %||% unique(as.character(df$prop_trt))
    
    c_lvls = if (!is.null(arr)) dim_levels(arr, "confdg_lvl") else NULL
    c_lvls = c_lvls %||% store$dims$third_levels %||% unique(as.character(df$confdg_lvl))
    
    effect_lvls = if (!is.null(arr)) dim_levels(arr, "effect") else NULL
    effect_lvls = effect_lvls %||% store$effect_levels %||% unique(as.character(df$effect))
    effect_lvls = unique(c(effect_lvls, "—"))
    
    method_lvls = current_method_names()
    if (!length(method_lvls)) {
      raw_method_lvls = if (!is.null(arr)) dim_levels(arr, "method") else NULL
      method_lvls = unique(method_display_name(raw_method_lvls %||% as.character(df$method)))
    }
    
    df$n = factor(df$n, levels = n_lvls)
    df$prop_trt = factor(df$prop_trt, levels = p_lvls)
    df$confdg_lvl = factor(df$confdg_lvl, levels = c_lvls)
    df$effect = factor(df$effect, levels = effect_lvls)
    df$method = factor(df$method, levels = method_lvls)
    df$estimand = factor(df$estimand, levels = c("ATE", "ATT", "—"))
    df$estimator = factor(df$estimator, levels = c("WLS", "DR", "—"))
    
    df
  }
  
  filter_common = function(df) {
    if ("method" %in% names(df)) {
      available_methods = current_method_names()
      if (length(available_methods)) {
        df = df[as.character(df$method) %in% available_methods, , drop = FALSE]
      }
      if (length(input$methods)) {
        df = df[as.character(df$method) %in% input$methods, , drop = FALSE]
      }
    }
    df
  }
  
  add_x_plot_coordinates = function(df) {
    sp = spacing_settings()
    n_levels = store$dims$n %||% levels(df$n) %||% unique(as.character(df$n))
    method_names = current_method_names()
    if (!length(method_names)) method_names = levels(df$method) %||% unique(as.character(df$method))
    confdg_levels = store$dims$third_levels %||% levels(df$confdg_lvl) %||% unique(as.character(df$confdg_lvl))
    
    df$n_index = as.numeric(factor(as.character(df$n), levels = n_levels))
    df$method_index = as.numeric(factor(as.character(df$method), levels = method_names))
    df$confdg_index = as.numeric(factor(as.character(df$confdg_lvl), levels = confdg_levels))
    
    method_center = mean(seq_along(method_names))
    confdg_center = mean(seq_along(confdg_levels))
    
    df$method_centered = df$method_index - method_center
    df$confdg_centered = df$confdg_index - confdg_center
    df$x_plot = df$n_index * sp$sample_size_spacing +
      df$method_centered * sp$method_spacing +
      df$confdg_centered * sp$confdg_spacing
    
    df
  }
  
  axis_breaks = function() {
    sp = spacing_settings()
    n_levels = store$dims$n
    list(
      breaks = seq_along(n_levels) * sp$sample_size_spacing,
      labels = n_levels
    )
  }
  
  effect_selector_for_metric = function(metric) {
    arr = metric_reference_array(store$env, metric)
    if (!is.null(arr) && has_dim(arr, "effect")) {
      return(selected_effect())
    }
    NULL
  }
  
  selected_table_selectors = reactive({
    req(input$n, input$prop_trt)
    
    selectors = list(
      n = input$n,
      prop_trt = input$prop_trt
    )
    
    if (!is.null(store$dims$third_levels)) {
      selectors[[store$dims$third_name]] = req(input$third)
    }
    
    eff = effect_selector_for_metric(input$metric)
    if (!is.null(eff)) selectors[["effect"]] = eff
    
    selectors
  })
  
  selected_plot_selectors = reactive({
    selectors = list()
    eff = effect_selector_for_metric(input$metric)
    if (!is.null(eff)) selectors[["effect"]] = eff
    selectors
  })
  
  metric_has_mcse = function(metric) {
    metric %in% names(metric_mcse_names) && exists(metric_mcse_names[[metric]], envir = store$env, inherits = FALSE)
  }
  
  metric_label_current = reactive({
    unname(metric_labs[input$metric])
  })
  
  non_smd_df = reactive({
    req(store$loaded, input$metric)
    validate(need(input$metric %in% c("bias", "variance", "mae", "rmse", "CI_coverage", "variance_ratio", "SNegW"),
                  "This plot is handled separately."))
    
    env = store$env
    validate(need(exists(input$metric, envir = env, inherits = FALSE), sprintf("Metric '%s' was not found.", input$metric)))
    
    arr = get(input$metric, envir = env)
    dft = slice_array_df(arr, selected_plot_selectors(), value_name = "value")
    
    if (metric_has_mcse(input$metric)) {
      mcse_arr = get(metric_mcse_names[[input$metric]], envir = env)
      mcse_df = slice_array_df(mcse_arr, selected_plot_selectors(), value_name = "mcse")
      by_cols = intersect(names(dft), names(mcse_df))
      dft = merge(dft, mcse_df, by = by_cols, all.x = TRUE, sort = FALSE)
      dft$lower = dft$value - 1.96 * dft$mcse
      dft$upper = dft$value + 1.96 * dft$mcse
    } else {
      dft$mcse = NA_real_
      dft$lower = NA_real_
      dft$upper = NA_real_
    }
    
    dft = factorize_common(dft, arr)
    dft = filter_common(dft)
    add_x_plot_coordinates(dft)
  })
  
  missingness_keys_for_metric = function(metric) {
    if (!exists("metric_missingness", envir = store$env, inherits = FALSE)) return(character(0))
    mm = get("metric_missingness", envir = store$env)
    if (!is.list(mm) || !(metric %in% names(mm))) return(character(0))
    unique(as.character(mm[[metric]]))
  }
  
  missingness_df_for_keys = function(keys, selectors = list()) {
    req(store$loaded)
    env = store$env
    validate(need(exists("missing", envir = env, inherits = FALSE), "Missingness object was not found."))
    
    miss = get("missing", envir = env)
    keys = intersect(keys, names(miss))
    validate(need(length(keys) > 0, "No matching missingness quantity was found."))
    
    out = list()
    for (key in keys) {
      arr = miss[[key]]
      df = if (length(selectors)) slice_array_df(arr, selectors, value_name = "value") else array_to_df(arr, value_name = "value")
      df$missing_key = key
      df$missing_type = missingness_type_labs[key] %||% key
      out[[key]] = df
    }
    
    df = do.call(rbind, out)
    ref_arr = miss[[keys[1]]]
    df = factorize_common(df, ref_arr)
    df = filter_common(df)
    df$missing_type = factor(df$missing_type, levels = unname(missingness_type_labs[keys]))
    df
  }
  
  missing_metric_df = reactive({
    req(store$loaded, input$metric == "missing")
    df = missingness_df_for_keys(c("estimate", "standard_error"), selectors = selected_plot_selectors())
    add_x_plot_coordinates(df)
  })
  
  missingness_scope_text = function(keys) {
    if (all(c("estimate", "standard_error") %in% keys)) {
      "point estimate and/or variance estimate"
    } else if ("standard_error" %in% keys) {
      "variance estimate"
    } else {
      "point estimate"
    }
  }
  
  missingness_df_for_plot_context = function(metric, dft_context = NULL) {
    keys = missingness_keys_for_metric(metric)
    if (!length(keys) || metric == "missing") return(NULL)
    
    df = missingness_df_for_keys(keys, selectors = selected_plot_selectors())
    if (is.null(df) || !nrow(df) || is.null(dft_context) || !nrow(dft_context)) return(df)
    
    context_cols = intersect(
      c("n", "prop_trt", "confdg_lvl", "effect", "method", "estimand", "estimator"),
      intersect(names(df), names(dft_context))
    )
    
    for (nm in context_cols) {
      vals = unique(as.character(dft_context[[nm]]))
      vals = vals[!is.na(vals)]
      if (length(vals)) df = df[as.character(df[[nm]]) %in% vals, , drop = FALSE]
    }
    
    df
  }
  
  methods_with_missingness_above_threshold = function(metric, dft_context = NULL) {
    df = missingness_df_for_plot_context(metric, dft_context)
    if (is.null(df) || !nrow(df)) return(character(0))
    
    threshold = missingness_asterisk_threshold()
    
    max_missing = tapply(
      df$value,
      as.character(df$method),
      function(x) {
        if (all(is.na(x))) return(NA_real_)
        max(x, na.rm = TRUE)
      }
    )
    
    flagged = !is.na(max_missing) & max_missing > threshold
    names(max_missing)[flagged]
  }
  
  method_legend_labels = function(metric, dft_context = NULL) {
    labs = setNames(current_method_names(), current_method_names())
    
    flagged = methods_with_missingness_above_threshold(metric, dft_context)
    flagged = intersect(flagged, names(labs))
    
    if (length(flagged)) labs[flagged] = paste0(labs[flagged], "*")
    
    labs
  }
  
  missingness_note_text_for_plot = function(metric, dft_context = NULL) {
    if (metric == "missing") return("")
    
    keys = missingness_keys_for_metric(metric)
    if (!length(keys)) return("")
    
    flagged = methods_with_missingness_above_threshold(metric, dft_context)
    if (!length(flagged)) return("")
    
    paste0(
      "* Methods marked with an asterisk have more than ",
      format_pct(missingness_asterisk_threshold()),
      " missing ",
      missingness_scope_text(keys),
      " in the displayed selection. ",
      "The corresponding metric was computed after excluding replications with missing values. ",
      "See the Missingness plot for the corresponding proportions."
    )
  }
  
  caption_for_metric = function(metric) {
    if (metric == "CI_coverage") {
      base = "The dashed line indicates nominal 95% coverage. The shaded band indicates the range around nominal coverage corresponding to ±1.96 Monte Carlo standard errors with 5,000 replications. Error bars represent ±1.96 Monte Carlo standard errors of the observed coverage."
    } else if (metric == "variance_ratio") {
      base = "The dashed line corresponds to equality between the average method-estimated variance and the empirical Monte Carlo variance. Error bars represent ±1.96 Monte Carlo standard errors."
    } else if (metric == "missing") {
      base = "Points show the proportion of Monte Carlo replications excluded because the corresponding quantity was missing."
    } else if (metric %in% names(metric_mcse_names)) {
      base = "Points show Monte Carlo estimates of the performance metric. Error bars represent ±1.96 Monte Carlo standard errors."
    } else {
      base = "Points show Monte Carlo estimates of the performance metric."
    }
    base
  }
  
  plot_note_text_for_plot = function(metric, dft_context = NULL) {
    parts = character(0)
    
    base_note = caption_for_metric(metric)
    if (nzchar(base_note)) parts = c(parts, base_note)
    
    if (metric == "variance_ratio") {
      parts = c(
        parts,
        paste0(
          "The y-axis is shown on a log2 scale from ",
          formatC(1 / variance_ratio_plot_max(), format = "fg", digits = 3),
          " to ",
          formatC(variance_ratio_plot_max(), format = "fg", digits = 3),
          "."
        )
      )
    }
    
    miss_note = missingness_note_text_for_plot(metric, dft_context)
    if (nzchar(miss_note)) parts = c(parts, miss_note)
    
    paste(parts, collapse = "\n")
  }
  
  base_plot_theme = function() {
    ts = plot_theme_settings()
    theme_minimal(base_size = 13) +
      theme(
        panel.background = element_rect(fill = "#ffffff", color = NA),
        panel.grid.major.x = element_line(color = "#e2e2e2", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "#e2e2e2", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        plot.margin = margin(ts$margin_top, ts$margin_right, ts$margin_bottom, ts$margin_left),
        panel.spacing.x = unit(ts$panel_spacing_x, "lines"),
        panel.spacing.y = unit(ts$panel_spacing_y, "lines"),
        plot.caption = element_text(hjust = 0, size = 9)
      )
  }
  
  make_non_smd_plot = function(dft_sub, ylab_txt, title_txt, metric,
                               size = 2.8, legend_size = size, show_title = TRUE) {
    validate(need(nrow(dft_sub) > 0, paste("No data for", title_txt)))
    
    dft_plot = dft_sub
    dft_plot$lower_plot = dft_plot$lower
    dft_plot$upper_plot = dft_plot$upper
    
    sh_map = shape_map_for(levels(dft_plot$confdg_lvl))
    ax = axis_breaks()
    mc = mcse_settings()
    
    ymin = suppressWarnings(min(dft_plot$value, dft_plot$lower_plot, na.rm = TRUE))
    ymax = suppressWarnings(max(dft_plot$value, dft_plot$upper_plot, na.rm = TRUE))
    
    y0 = if (metric %in% c("variance", "mae", "rmse", "CI_coverage", "variance_ratio")) {
      if (metric == "variance_ratio" && is.finite(ymin) && ymin < 0) ymin else 0
    } else if (metric == "bias") {
      if (is.finite(ymin)) min(0, ymin) else 0
    } else {
      if (is.finite(ymin)) ymin else 0
    }
    
    y1 = if (is.finite(ymax)) ymax * 1.08 else 1
    
    if (metric == "CI_coverage") {
      y0 = 0.50
      y1 = 1
    }
    
    if (metric == "variance_ratio") {
      vr_max = variance_ratio_plot_max()
      vr_min = 1 / vr_max
      y0 = vr_min
      y1 = vr_max
      
      dft_plot$value[dft_plot$value <= 0] = NA_real_
      dft_plot$lower_plot = ifelse(is.finite(dft_plot$lower_plot), pmax(dft_plot$lower_plot, vr_min), NA_real_)
      dft_plot$upper_plot = ifelse(is.finite(dft_plot$upper_plot), dft_plot$upper_plot, NA_real_)
    }
    
    p = ggplot(dft_plot, aes(x = x_plot, y = value, color = method, shape = confdg_lvl))
    
    if (metric == "CI_coverage") {
      lower_nominal = 0.95 - 1.96 * sqrt(0.95 * 0.05 / 5000)
      upper_nominal = 0.95 + 1.96 * sqrt(0.95 * 0.05 / 5000)
      p = p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = lower_nominal, ymax = upper_nominal, alpha = 0.15)
      p = p + geom_hline(yintercept = 0.95, linetype = "dashed", color = "black")
    }
    
    if (metric == "bias") {
      p = p + geom_hline(yintercept = 0, color = "black")
    }
    
    if (metric == "variance_ratio") {
      p = p + geom_hline(yintercept = 1, linetype = "dashed", color = "black")
    }
    
    method_cols = method_colors_current()
    method_labs = method_legend_labels(metric, dft_plot)
    
    # Draw points first, then MCSE intervals on top so short 95% MC confidence intervals remain visible.
    p = p + geom_point(size = size, na.rm = TRUE)
    
    if (metric_has_mcse(metric) && isTRUE(mc$show)) {
      p = p + geom_linerange(
        aes(ymin = lower_plot, ymax = upper_plot),
        linewidth = mc$width,
        color = mc$color,
        alpha = mc$alpha,
        na.rm = TRUE,
        show.legend = FALSE
      )
    }
    
    p = p +
      scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
      scale_color_manual(values = method_cols, labels = method_labs[names(method_cols)], name = "Method", drop = FALSE) +
      scale_shape_manual(values = sh_map, name = "Level of\ncomplexity") +
      labs(
        x = "sample size (n)",
        y = ylab_txt,
        title = if (isTRUE(show_title)) title_txt else NULL,
        caption = NULL
      ) +
      base_plot_theme() +
      guides(
        color = guide_legend(override.aes = list(size = legend_size)),
        shape = guide_legend(override.aes = list(size = legend_size))
      )
    
    if (any(as.character(dft_plot$estimator) != "—", na.rm = TRUE)) {
      p = p + facet_grid(
        rows = vars(estimator),
        cols = vars(prop_trt),
        labeller = labeller(estimator = est_lab, prop_trt = prop_lab)
      )
    } else {
      p = p + facet_grid(
        cols = vars(prop_trt),
        labeller = labeller(prop_trt = prop_lab)
      )
    }
    
    if (metric == "variance_ratio") {
      p = p + scale_y_continuous(transform = "log2") + coord_cartesian(ylim = c(y0, y1))
    } else {
      p = p + coord_cartesian(ylim = c(y0, y1))
    }
    
    p
  }
  
  smd_plot_df = reactive({
    req(store$loaded, input$metric == "smd")
    env = store$env
    validate(need(exists("smd", envir = env, inherits = FALSE), "SMD metric not found."))
    
    arr = get("smd", envir = env)
    df = array_to_df(arr, value_name = "smd")
    df = factorize_common(df, arr)
    df = filter_common(df)
    validate(need(nrow(df) > 0, "No SMD data for the current selection."))
    
    group_cols = intersect(c("n", "prop_trt", "confdg_lvl", "method", "estimand"), names(df))
    df$smd_abs = abs(df$smd)
    
    out = aggregate(
      df$smd_abs,
      df[, group_cols, drop = FALSE],
      function(z) median(z, na.rm = TRUE)
    )
    
    names(out)[ncol(out)] = "value"
    out = factorize_common(out, arr)
    add_x_plot_coordinates(out)
  })
  
  make_smd_plot = function(df, estimand_keep, size = 2.8, legend_size = size, show_title = TRUE) {
    df = df[as.character(df$estimand) == estimand_keep, , drop = FALSE]
    validate(need(nrow(df) > 0, paste("No SMD data to plot for", estimand_keep)))
    
    sh_map = shape_map_for(levels(df$confdg_lvl))
    ax = axis_breaks()
    
    ggplot(df, aes(x = x_plot, y = value, color = method, shape = confdg_lvl)) +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
      geom_point(size = size, na.rm = TRUE) +
      scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
      scale_color_manual(values = method_colors_current(), name = "Method", drop = FALSE) +
      scale_shape_manual(values = sh_map, name = "Level of\ncomplexity") +
      labs(
        x = "sample size (n)",
        y = "Median absolute SMD",
        title = if (isTRUE(show_title)) paste("Standardized Mean Difference —", estimand_keep) else NULL,
        caption = NULL
      ) +
      facet_grid(cols = vars(prop_trt), labeller = labeller(prop_trt = prop_lab)) +
      base_plot_theme() +
      coord_cartesian(ylim = c(0, NA)) +
      guides(
        color = guide_legend(override.aes = list(size = legend_size)),
        shape = guide_legend(override.aes = list(size = legend_size))
      )
  }
  
  ess_plot_df = reactive({
    req(store$loaded, input$metric == "ess")
    env = store$env
    validate(need(exists("ess", envir = env, inherits = FALSE), "ESS metric not found."))
    
    arr = get("ess", envir = env)
    df = array_to_df(arr, value_name = "value")
    
    if ("group" %in% names(df) && !is.null(store$dims$ess_groups) && length(store$dims$ess_groups)) {
      df = df[df$group == store$dims$ess_groups[1], , drop = FALSE]
    }
    
    df = factorize_common(df, arr)
    df = filter_common(df)
    add_x_plot_coordinates(df)
  })
  
  make_ess_plot = function(df, estimand_keep, size = 2.8, legend_size = size, show_title = TRUE) {
    df = df[as.character(df$estimand) == estimand_keep, , drop = FALSE]
    validate(need(nrow(df) > 0, paste("No ESS data to plot for", estimand_keep)))
    
    sh_map = shape_map_for(levels(df$confdg_lvl))
    ax = axis_breaks()
    
    ggplot(df, aes(x = x_plot, y = value, color = method, shape = confdg_lvl)) +
      geom_point(size = size, na.rm = TRUE) +
      scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
      scale_color_manual(values = method_colors_current(), name = "Method", drop = FALSE) +
      scale_shape_manual(values = sh_map, name = "Level of\ncomplexity") +
      labs(
        x = "sample size (n)",
        y = "ESS",
        title = if (isTRUE(show_title)) paste("Effective Sample Size —", estimand_keep) else NULL,
        caption = NULL
      ) +
      facet_grid(cols = vars(prop_trt), labeller = labeller(prop_trt = prop_lab)) +
      base_plot_theme() +
      coord_cartesian(ylim = c(0, NA)) +
      guides(
        color = guide_legend(override.aes = list(size = legend_size)),
        shape = guide_legend(override.aes = list(size = legend_size))
      )
  }
  
  make_missing_plot = function(df, estimand_keep, size = 2.8, legend_size = size, show_title = TRUE) {
    df = df[as.character(df$estimand) == estimand_keep, , drop = FALSE]
    validate(need(nrow(df) > 0, paste("No missingness data for", estimand_keep)))
    
    sh_map = shape_map_for(levels(df$confdg_lvl))
    ax = axis_breaks()
    
    ggplot(df, aes(x = x_plot, y = value, color = method, shape = confdg_lvl)) +
      geom_hline(yintercept = 0, color = "black") +
      geom_point(size = size, na.rm = TRUE) +
      scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
      scale_color_manual(values = method_colors_current(), name = "Method", drop = FALSE) +
      scale_shape_manual(values = sh_map, name = "Level of\ncomplexity") +
      labs(
        x = "sample size (n)",
        y = "Proportion missing",
        title = if (isTRUE(show_title)) paste("Missingness —", estimand_keep) else NULL,
        caption = NULL
      ) +
      facet_grid(
        rows = vars(missing_type, estimator),
        cols = vars(prop_trt),
        labeller = labeller(estimator = est_lab, prop_trt = prop_lab)
      ) +
      base_plot_theme() +
      coord_cartesian(ylim = c(0, 1)) +
      guides(
        color = guide_legend(override.aes = list(size = legend_size)),
        shape = guide_legend(override.aes = list(size = legend_size))
      )
  }
  
  plot_for_estimand = function(estimand_keep, show_title = TRUE) {
    req(store$loaded, input$metric)
    metric = input$metric
    
    if (metric == "smd") {
      return(make_smd_plot(smd_plot_df(), estimand_keep, show_title = show_title))
    }
    
    if (metric == "ess") {
      return(make_ess_plot(ess_plot_df(), estimand_keep, show_title = show_title))
    }
    
    if (metric == "missing") {
      return(make_missing_plot(missing_metric_df(), estimand_keep, show_title = show_title))
    }
    
    dft = non_smd_df()
    dft = dft[as.character(dft$estimand) == estimand_keep, , drop = FALSE]
    lab = metric_label_current()
    make_non_smd_plot(
      dft,
      ylab_txt = lab,
      title_txt = paste("Metric:", lab, "—", estimand_keep),
      metric = metric,
      show_title = show_title
    )
  }
  
  plot_note_for_estimand = function(estimand_keep) {
    if (input$metric %in% c("smd", "ess")) return(caption_for_metric(input$metric))
    if (input$metric == "missing") return(caption_for_metric("missing"))
    
    dft = non_smd_df()
    dft = dft[as.character(dft$estimand) == estimand_keep, , drop = FALSE]
    plot_note_text_for_plot(input$metric, dft)
  }
  
  download_pdf_filename = function(estimand_keep) {
    paste0(
      safe_filename_part(selected_effect() %||% "all"), "-",
      safe_filename_part(input$metric), "-",
      safe_filename_part(estimand_keep),
      ".pdf"
    )
  }
  
  download_plot_width = function() 11
  download_plot_height = function() 5.8
  
  output$plot_ui = renderUI({
    req(input$metric)
    
    lab = unname(metric_labs[input$metric])
    
    tagList(
      tags$h4(paste("Metric:", lab, "— ATE")),
      plotOutput("metric_plot_ATE", height = 450),
      downloadButton("download_ATE", "Download ATE PDF"),
      textOutput("plot_note_ATE"),
      tags$hr(),
      tags$h4(paste("Metric:", lab, "— ATT")),
      plotOutput("metric_plot_ATT", height = 450),
      downloadButton("download_ATT", "Download ATT PDF"),
      textOutput("plot_note_ATT")
    )
  })
  
  output$metric_plot_ATE = renderPlot({
    plot_for_estimand("ATE", show_title = FALSE)
  })
  
  output$metric_plot_ATT = renderPlot({
    plot_for_estimand("ATT", show_title = FALSE)
  })
  
  output$plot_note_ATE = renderText({
    plot_note_for_estimand("ATE")
  })
  
  output$plot_note_ATT = renderText({
    plot_note_for_estimand("ATT")
  })
  
  output$download_ATE = downloadHandler(
    filename = function() download_pdf_filename("ATE"),
    content = function(file) {
      ggsave(filename = file, plot = plot_for_estimand("ATE", show_title = FALSE), width = download_plot_width(), height = download_plot_height(), units = "in")
    }
  )
  
  output$download_ATT = downloadHandler(
    filename = function() download_pdf_filename("ATT"),
    content = function(file) {
      ggsave(filename = file, plot = plot_for_estimand("ATT", show_title = FALSE), width = download_plot_width(), height = download_plot_height(), units = "in")
    }
  )
  
  smd_table_df = reactive({
    req(store$loaded, exists("smd", envir = store$env, inherits = FALSE))
    arr = get("smd", envir = store$env)
    selectors = selected_table_selectors()
    df2 = slice_array_df(arr, selectors, value_name = "smd")
    df2 = factorize_common(df2, arr)
    df2 = filter_common(df2)
    validate(need(nrow(df2) > 0, "No SMD data for the current selection."))
    
    agg = aggregate(abs(smd) ~ method + estimand, df2, function(z) median(z, na.rm = TRUE))
    names(agg)[names(agg) == "abs(smd)"] = "value"
    agg$estimator = "—"
    factorize_common(agg, arr)[, c("method", "estimand", "estimator", "value")]
  })
  
  current_df = reactive({
    req(store$loaded, input$metric)
    env = store$env
    validate(need(exists(input$metric, envir = env, inherits = FALSE), sprintf("Metric '%s' not found.", input$metric)))
    
    arr = get(input$metric, envir = env)
    selectors = selected_table_selectors()
    
    if (input$metric == "ess" && !is.null(store$dims$ess_groups) && length(store$dims$ess_groups)) {
      selectors[["group"]] = store$dims$ess_groups[1]
    }
    
    df = slice_array_df(arr, selectors, value_name = "value")
    df = factorize_common(df, arr)
    filter_common(df)
  })
  
  wide_by_method = function(df, value_col = "display", method_order = current_method_names(), append_mcse_label = FALSE) {
    orig_has_estimator = "estimator" %in% names(df) && any(as.character(df$estimator) != "—")
    orig_has_estimand = "estimand" %in% names(df) && any(as.character(df$estimand) != "—")
    
    if (!"estimator" %in% names(df)) df$estimator = "—"
    if (!"estimand" %in% names(df)) df$estimand = "—"
    
    if (orig_has_estimator && orig_has_estimand) {
      df$colkey = paste(df$estimator, df$estimand, sep = " · ")
    } else if (!orig_has_estimator && orig_has_estimand) {
      df$colkey = as.character(df$estimand)
    } else if (orig_has_estimator && !orig_has_estimand) {
      df$colkey = as.character(df$estimator)
    } else {
      df$colkey = "value"
    }
    
    mat = tapply(df[[value_col]], list(as.character(df$method), df$colkey), function(x) x[1])
    tb = as.data.frame.matrix(mat, stringsAsFactors = FALSE)
    
    if (orig_has_estimator && orig_has_estimand) {
      est_levels = c("WLS", "DR", setdiff(unique(as.character(df$estimator)), c("WLS", "DR", "—")))
      estim_levels = c("ATE", "ATT", setdiff(unique(as.character(df$estimand)), c("ATE", "ATT", "—")))
      desired = as.vector(outer(est_levels, estim_levels, paste, sep = " · "))
      desired = intersect(desired, colnames(tb))
      if (length(desired)) tb = tb[, desired, drop = FALSE]
    } else if (!orig_has_estimator && orig_has_estimand) {
      desired = intersect(c("ATE", "ATT"), colnames(tb))
      if (length(desired)) tb = tb[, desired, drop = FALSE]
    } else if (orig_has_estimator && !orig_has_estimand) {
      desired = intersect(c("WLS", "DR"), colnames(tb))
      if (length(desired)) tb = tb[, desired, drop = FALSE]
    }
    
    if (append_mcse_label && ncol(tb)) colnames(tb) = paste0(colnames(tb), " value (MCSE)")
    tb = tb[order(match(rownames(tb), method_order)), , drop = FALSE]
    tb
  }
  
  add_display_value = function(df, include_mcse = FALSE) {
    digits = table_digits()
    if (include_mcse && "mcse" %in% names(df)) {
      val = fmt_num(df$value, digits)
      mcse = fmt_num(df$mcse, digits)
      df$display = ifelse(is.na(df$mcse), val, paste0(val, " (", mcse, ")"))
    } else {
      df$display = fmt_num(df$value, digits)
    }
    df
  }
  
  metric_table_df = reactive({
    req(store$loaded, input$metric)
    if (input$metric == "smd") return(smd_table_df())
    if (input$metric == "missing") return(NULL)
    
    df = current_df()
    
    if (metric_has_mcse(input$metric)) {
      mcse_arr = get(metric_mcse_names[[input$metric]], envir = store$env)
      mcse_df = slice_array_df(mcse_arr, selected_table_selectors(), value_name = "mcse")
      mcse_df = factorize_common(mcse_df, mcse_arr)
      mcse_df = filter_common(mcse_df)
      by_cols = intersect(names(df), names(mcse_df))
      by_cols = setdiff(by_cols, c("value", "mcse"))
      df = merge(df, mcse_df, by = by_cols, all.x = TRUE, sort = FALSE)
    }
    
    df
  })
  
  missing_table_for_metric = reactive({
    req(store$loaded, input$metric)
    
    keys = if (input$metric == "missing") c("estimate", "standard_error") else missingness_keys_for_metric(input$metric)
    if (!length(keys)) return(NULL)
    
    df = missingness_df_for_keys(keys, selectors = selected_table_selectors())
    if (!nrow(df) || (input$metric != "missing" && !any(df$value > 0, na.rm = TRUE))) return(NULL)
    
    keep = c("method", "estimand", "estimator", "missing_key", "value")
    df = df[, intersect(keep, names(df)), drop = FALSE]
    names(df)[names(df) == "value"] = "missing_value"
    
    wide = reshape(
      df,
      idvar = c("method", "estimand", "estimator"),
      timevar = "missing_key",
      direction = "wide"
    )
    
    names(wide) = sub("missing_value.estimate", "missing_estimate", names(wide), fixed = TRUE)
    names(wide) = sub("missing_value.standard_error", "missing_variance_estimate", names(wide), fixed = TRUE)
    
    if ("missing_estimate" %in% names(wide)) wide$missing_estimate = fmt_num(wide$missing_estimate, table_digits())
    if ("missing_variance_estimate" %in% names(wide)) wide$missing_variance_estimate = fmt_num(wide$missing_variance_estimate, table_digits())
    
    wide[order(match(as.character(wide$method), current_method_names())), , drop = FALSE]
  })
  
  output$metric_table = renderTable({
    if (input$metric == "missing") {
      mt = missing_table_for_metric()
      validate(need(!is.null(mt) && nrow(mt) > 0, "No missingness data to display for the current selection."))
      return(mt)
    }
    
    df = metric_table_df()
    validate(need(!is.null(df) && nrow(df) > 0, "No data to display for the current selection."))
    
    include_mcse = metric_has_mcse(input$metric)
    df = add_display_value(df, include_mcse = include_mcse)
    tb = wide_by_method(df, append_mcse_label = include_mcse)
    
    data.frame(Method = rownames(tb), tb, check.names = FALSE, row.names = NULL)
  }, striped = TRUE, bordered = TRUE, width = "100%", align = "l")
  
  output$missingness_table_ui = renderUI({
    req(store$loaded, input$metric)
    if (input$metric == "missing") return(NULL)
    
    mt = missing_table_for_metric()
    if (is.null(mt) || !nrow(mt)) return(NULL)
    
    tagList(tags$h4("Missingness for this selection"), tableOutput("missingness_table"))
  })
  
  output$missingness_table = renderTable({
    mt = missing_table_for_metric()
    validate(need(!is.null(mt) && nrow(mt) > 0, "No missingness for this selection."))
    mt
  }, striped = TRUE, bordered = TRUE, width = "100%", align = "l")
}

shinyApp(ui, server)
