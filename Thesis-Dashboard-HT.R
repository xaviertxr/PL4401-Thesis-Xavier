# ============================================================
#  R SHINY App: Network Recovery Simulation Dashboard
#  National University of Singapore
#  PL4401 Honours Thesis (HT)
#  Title of HT: Evaluating Imputation Strategies for Missing Data in Network Analysis
#  
#  Author: Tan Xuan Rong Xavier
#  Supervisor: Dr Sacha Epskamp
#
#  PURPOSE:
#  This dashboard visualises the results of a simulation study
#  comparing 22 network estimation methods under varying levels
#  of missingness (0–90%) and sample sizes (N = 250, 500, 750).
#
#  The "true" network was generated using genGGM() with 8 nodes,
#  20% edge density, and 80% positive edges. Data were drawn from
#  this true GGM and estimation methods were applied to recover it.
#  Performance was measured using sensitivity, specificity, F1,
#  MCC, and Jaccard index across 200 replications per condition.
#
#  DATA INPUTS (must be in environment before running):
#    - all_batches         : raw results (one row per rep)
#    - all_batches_summary : summarised results (mean/SD per condition)
#
#  TABS:
#    1. PERFORMANCE   — line plots of recovery metrics vs missingness/N
#    2. EDGE COUNTS   — how many edges each method recovers vs truth
#    3. FP vs FN      — false positives and false negatives separately
#    4. TOTAL ERRORS  — combined FP + FN classification errors
#    5. FAILURES      — % of reps where the method threw an error
#    6. HEATMAP       — all estimators × missingness in one grid view
# ============================================================

library(shiny)
library(ggplot2)
library(plotly)   
library(dplyr)
library(tidyr)
library(DT)      

all_batches         <- readRDS("all_batches.rds")
all_batches_summary <- readRDS("all_batches_summary.rds")

# ── Colour palette ───────────────────────────────────────────────────────────
# One fixed colour per estimator, shared across ALL tabs so the same
# estimator is always the same colour regardless of tab
estimator_levels <- sort(unique(all_batches$estimator))
n_est <- length(estimator_levels)
pal <- setNames(
  colorRampPalette(c("#003f5c","#2f4b7c","#665191","#a05195",
                     "#d45087","#f95d6a","#ff7c43","#ffa600",
                     "#7bc8a4","#4cc9f0","#4361ee","#3a0ca3"))(n_est),
  estimator_levels
)

# ── Estimator family groupings ───────────────────────────────────────────────
# Used for the "Quick-select family" buttons in the sidebar,
family_map <- tibble(estimator = estimator_levels) %>%
  mutate(family = case_when(
    startsWith(estimator, "ggm_prune")    ~ "GGM Prune",
    startsWith(estimator, "EBICglasso")   ~ "EBICglasso",
    startsWith(estimator, "ggmModSelect") ~ "ggmModSelect",
    startsWith(estimator, "EM_")          ~ "EM",
    TRUE                                  ~ "Other"
  ))


# ════════════════════════════════════════════════════════════════════════════
#  UI
#  Defines the visual layout: sidebar (controls) + main panel (plots/tables)
# ════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  
  # ── Styling (dark GitHub-inspired theme) ──────────────────────────────────
  tags$head(
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=DM+Sans:wght@300;400;600&display=swap"),
    tags$style(HTML("
      * { box-sizing: border-box; }
      body {
        background: #0d1117;
        color: #e6edf3;
        font-family: 'DM Sans', sans-serif;
        font-size: 14px;
      }
      .navbar-title {
        font-family: 'DM Mono', monospace;
        font-size: 1.05rem;
        letter-spacing: 0.08em;
        color: #58a6ff;
        padding: 18px 24px 6px;
      }
      .subtitle {
        font-size: 0.78rem;
        color: #8b949e;
        padding: 0 24px 14px;
        font-family: 'DM Mono', monospace;
        letter-spacing: 0.04em;
      }
      .well {
        background: #161b22 !important;
        border: 1px solid #30363d !important;
        border-radius: 8px !important;
        padding: 16px !important;
      }
      label { color: #8b949e !important; font-size: 0.78rem; letter-spacing: 0.05em; }
      .selectize-input, .selectize-dropdown {
        background: #21262d !important;
        border: 1px solid #30363d !important;
        color: #e6edf3 !important;
        font-size: 0.82rem;
      }
      .selectize-dropdown-content .option { color: #e6edf3 !important; }
      .selectize-dropdown-content .option:hover { background: #30363d !important; }
      .btn { border-radius: 6px; font-size: 0.8rem; }
      .nav-tabs { border-bottom: 1px solid #30363d !important; }
      .nav-tabs > li > a {
        color: #8b949e !important;
        background: transparent !important;
        border: none !important;
        font-family: 'DM Mono', monospace;
        font-size: 0.78rem;
        letter-spacing: 0.06em;
        padding: 8px 14px;
      }
      .nav-tabs > li.active > a {
        color: #58a6ff !important;
        border-bottom: 2px solid #58a6ff !important;
        background: transparent !important;
      }
      .tab-content { padding-top: 16px; }
      .dataTables_wrapper { color: #e6edf3 !important; font-size: 0.72rem; }
      .dataTables_wrapper td, .dataTables_wrapper th { padding: 4px 8px !important; font-size: 0.72rem; }
      table.dataTable thead th {
        background: #161b22 !important;
        color: #58a6ff !important;
        border-bottom: 1px solid #30363d !important;
      }
      table.dataTable tbody tr { background: #0d1117 !important; color: #e6edf3 !important; }
      table.dataTable tbody tr:hover { background: #161b22 !important; }
      .dataTables_filter input, .dataTables_length select {
        background: #21262d !important;
        color: #e6edf3 !important;
        border: 1px solid #30363d !important;
        border-radius: 4px;
      }
      .dataTables_info, .dataTables_paginate { color: #8b949e !important; }
      .paginate_button { color: #8b949e !important; }
      .paginate_button.current { color: #58a6ff !important; font-weight: 600; }
      hr { border-color: #30363d; }
      .shiny-plot-output { border-radius: 8px; overflow: hidden; }
      .stat-box {
        background: #161b22;
        border: 1px solid #30363d;
        border-radius: 8px;
        padding: 12px 16px;
        margin-bottom: 10px;
        text-align: center;
      }
      .stat-box .val {
        font-family: 'DM Mono', monospace;
        font-size: 1.6rem;
        color: #58a6ff;
        font-weight: 500;
      }
      .stat-box .lbl { font-size: 0.72rem; color: #8b949e; margin-top: 2px; }
      .family-btn {
        margin: 2px;
        font-size: 0.72rem;
        padding: 3px 8px;
        background: #21262d;
        border: 1px solid #30363d;
        color: #8b949e;
        cursor: pointer;
        border-radius: 4px;
      }
    "))
  ),
  
  # ── Title bar ─────────────────────────────────────────────────────────────
  div(class = "navbar-title", "● NETWORK RECOVERY — SIMULATION DASHBOARD"),
  div(class = "subtitle", paste0(
    "n_nodes = 8  ·  reps = 200  ·  estimators = ", n_est,
    "  ·  missingness = 0–0.9  ·  n_people ∈ {250, 500, 750}"
  )),
  
  sidebarLayout(
    
    # ── SIDEBAR ─────────────────────────────────────────────────────────────
    sidebarPanel(width = 3,
                 
                 # QUICK-SELECT FAMILY BUTTONS
                 # Clicking e.g. "EBICglasso" instantly loads all EBICglasso estimators
                 # into the estimator selector below. "All" loads everything.
                 # These use JavaScript to fire a Shiny input event (family_select),
                 # which is caught by observeEvent() in the server.
                 tags$div(
                   tags$p(style = "color:#8b949e;font-size:0.72rem;margin-bottom:4px;letter-spacing:0.05em;",
                          "QUICK-SELECT FAMILY"),
                   tags$div(id = "family_btns",
                            lapply(unique(family_map$family), function(fam) {
                              tags$button(fam, class = "family-btn",
                                          onclick = sprintf(
                                            "Shiny.setInputValue('family_select', '%s', {priority: 'event'})", fam))
                            })
                   ),
                   tags$button("All", class = "family-btn",
                               style = "color:#58a6ff;border-color:#58a6ff;",
                               onclick = "Shiny.setInputValue('family_select', 'ALL', {priority: 'event'})")
                 ),
                 
                 hr(),
                 
                 # ESTIMATOR MULTI-SELECT
                 # Can select any subset of the 22 estimators.
                 # The × button on each tag removes that estimator.
                 selectizeInput("estimators", "ESTIMATOR(S):",
                                choices  = estimator_levels,
                                multiple = TRUE,
                                selected = c("EBICglasso_listwise", "EBICglasso_fiml",
                                             "EBICglasso_miceMI", "ggm_prune", "EM_EBIC"),
                                options  = list(plugins = list("remove_button"))),
                 
                 # METRIC SELECTOR
                 # Switches the y-axis on the Performance tab.
                 # sensitivity = TP/(TP+FN): did we find the true edges?
                 # specificity = TN/(TN+FP): did we correctly exclude false edges?
                 # f1          = harmonic mean of precision and sensitivity
                 # mcc         = Matthews Correlation Coefficient (best single metric)
                 # jaccard     = TP/(TP+FP+FN): overlap between true and estimated edges
                 selectInput("metric", "METRIC:",
                             choices  = c("sensitivity","specificity","f1","mcc","jaccard"),
                             selected = "f1"),
                 
                 # X-AXIS TOGGLE
                 # "Missingness": x = 0 to 0.9, one plot per selected N (below)
                 # "Sample size": x = 250/500/750, faceted by missingness level
                 radioButtons("xaxis", "X-AXIS:",
                              choices  = c("Missingness" = "missing",
                                           "Sample size" = "n_people"),
                              selected = "missing",
                              inline   = TRUE),
                 
                 # DISPLAY OPTIONS
                 checkboxInput("show_points", "Show points", TRUE),
                 # SD ribbon shows mean ± 1 SD across the 200 reps, to gives a sense of
                 # how variable the estimates are, not just the average performance
                 checkboxInput("show_ribbon", "Show ± 1 SD ribbon", FALSE),
                 
                 # SAMPLE SIZE SELECTOR
                 # Only active when X-axis = Missingness.
                 # Shows one clean plot for the chosen N instead of three cluttered facets.
                 # When X-axis = Sample size this is ignored (all N shown on x-axis).
                 radioButtons("n_people_sel", "SAMPLE SIZE (N):",
                              choices  = c("250" = 250, "500" = 500, "750" = 750),
                              selected = 500,
                              inline   = TRUE),
                 
                 hr(),
                 
                 # MINI STAT BOXES
                 uiOutput("miniStats")
    ),
    
    # ── MAIN PANEL ──────────────────────────────────────────────────────────
    mainPanel(width = 9,
              tabsetPanel(id = "tabs",
                          
                          # TAB 1: PERFORMANCE
                          # Primary results tab. Y-axis = chosen metric (e.g. F1).
                          # Use the plotly toolbar 
                          tabPanel("PERFORMANCE",
                                   plotlyOutput("metricPlot", height = "480px"),
                                   br(),
                                   # Sortable table of the same data
                                   DTOutput("metricTable")),
                          
                          # TAB 2: EDGE COUNTS
                          # Shows mean number of edges estimated vs the true edge count (~6).
                          # Three reference lines:
                          #   Blue dotted  = empty network (0 edges) — method predicts nothing
                          #   Purple dotted = complete network (28 edges) — method over-connects
                          # Key diagnostic: if a method hits 0 edges at high missingness,
                          # its high specificity in Tab 1 is artificial (see Tab 3).
                          tabPanel("EDGE COUNTS",
                                   plotlyOutput("edgePlot", height = "480px"),
                                   br(),),
                          
                          # TAB 3: FP vs FN (SEPARATE)
                          # Solid lines = False Positives (edges estimated that don't exist)
                          # Dashed lines = False Negatives (true edges that were missed)
                          # Key insight: methods that "collapse" to an empty network at high
                          # missingness will show FP = 0 (good-looking!) but FN = ~6 (all
                          # true edges missed). This exposes artificial specificity.
                          tabPanel("FP vs FN",
                                   checkboxGroupInput("fpfn_show", label = NULL,
                                                      choices  = c("False Positives (FP)" = "FP",
                                                                   "False Negatives (FN)" = "FN"),
                                                      selected = c("FP", "FN"), inline = TRUE),
                                   plotlyOutput("fpfnPlot", height = "450px"),
                                   br(),),
                          
                          # TAB 4: TOTAL ERRORS (FP + FN combined)
                          # Simpler summary of overall classification error.
                          tabPanel("TOTAL ERRORS",
                                   plotlyOutput("errorPlot", height = "480px"),
                                   br(),),
                          
                          # TAB 5: METHOD FAILURES
                          # % of the 200 reps where the estimator threw an R error
                          tabPanel("FAILURES",
                                   plotlyOutput("failurePlot", height = "480px"),
                                   br(),),
                          
                          # TAB 6: HEATMAP
                          # Bird's-eye summary: estimators (rows) × missingness (columns).
                          tabPanel("HEATMAP",
                                   fluidRow(
                                     column(4,
                                            selectInput("hm_metric", "Metric for heatmap:",
                                                        choices  = c("sensitivity","specificity","f1","mcc","jaccard"),
                                                        selected = "f1"),
                                            selectInput("hm_n", "Sample size:",
                                                        choices  = c(250, 500, 750),
                                                        selected = 500)
                                     )
                                   ),
                                   plotOutput("heatmapPlot", height = "600px"))
              )
    )
  )
)


# ════════════════════════════════════════════════════════════════════════════
#  SERVER
#  Reactive logic: filters data in response to UI inputs, renders plots/tables
# ════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  # ── Family quick-select ────────────────────────────────────────────────────
  observeEvent(input$family_select, {
    sel <- if (input$family_select == "ALL") {
      estimator_levels
    } else {
      family_map %>% filter(family == input$family_select) %>% pull(estimator)
    }
    updateSelectizeInput(session, "estimators", selected = sel)
  })
  
  # ── Reactive data filters ──────────────────────────────────────────────────
  # filt_summary: used by Tab 1 (performance metrics)
  #   - filters all_batches_summary by chosen metric, estimators, and N
  #   - when x-axis = missingness, only one N is shown (cleaner plot)
  filt_summary <- reactive({
    d <- all_batches_summary %>% filter(metric == input$metric)
    if (length(input$estimators) > 0) d <- filter(d, estimator %in% input$estimators)
    if (input$xaxis == "missing") d <- filter(d, n_people == as.numeric(input$n_people_sel))
    d
  })
  
  # filt_raw: used by Tabs 2–5 (edge counts, FP/FN, errors, failures)
  filt_raw <- reactive({
    d <- all_batches
    if (length(input$estimators) > 0) d <- filter(d, estimator %in% input$estimators)
    if (input$xaxis == "missing") d <- filter(d, n_people == as.numeric(input$n_people_sel))
    d
  })
  
  # ── Mini stat boxes ────────────────────────────────────────────────────────
  # Shows number of estimators currently selected.
  output$miniStats <- renderUI({
    if (length(input$estimators) == 0) return(NULL)
    tagList(
      tags$div(class = "stat-box",
               tags$div(class = "val", length(input$estimators)),
               tags$div(class = "lbl", "ESTIMATORS SELECTED")
      )
    )
  })
  
  # ── Shared ggplot dark theme ───────────────────────────────────────────────
  # Applied to every plot for visual consistency.
  # Uses DM Sans + DM Mono (Google Fonts loaded in the UI head).
  dark_theme <- function() {
    theme_minimal(base_size = 12, base_family = "DM Sans") +
      theme(
        plot.background   = element_rect(fill = "#0d1117", color = NA),
        panel.background  = element_rect(fill = "#161b22", color = NA),
        panel.grid.major  = element_line(color = "#21262d"),
        panel.grid.minor  = element_blank(),
        axis.text         = element_text(color = "#8b949e", size = 9),
        axis.title        = element_text(color = "#8b949e", size = 10),
        strip.text        = element_text(color = "#58a6ff", size = 9, family = "DM Mono"),
        strip.background  = element_rect(fill = "#161b22", color = "#30363d"),
        legend.background = element_rect(fill = "#161b22", color = NA),
        legend.key        = element_rect(fill = "#161b22", color = NA),
        legend.text       = element_text(color = "#8b949e", size = 8),
        legend.title      = element_text(color = "#58a6ff", size = 9, family = "DM Mono"),
        plot.title        = element_text(color = "#e6edf3", size = 13, family = "DM Mono", face = "bold"),
        plot.subtitle     = element_text(color = "#8b949e", size = 9, family = "DM Mono"),
        legend.position   = "bottom"
      )
  }
  
  # ── X-axis / faceting helper ───────────────────────────────────────────────
  # When x = missingness → single panel (N already filtered upstream)
  # When x = sample size → facet by missingness level (one panel per level)
  apply_xaxis <- function(p, xaxis) {
    if (xaxis == "missing") {
      p + scale_x_continuous(breaks = seq(0, 0.9, 0.1))
    } else {
      p + scale_x_continuous(breaks = c(250, 500, 750)) +
        facet_wrap(~missing, labeller = label_both)
    }
  }
  
  # ── Fix duplicate legend entries from ggplotly ────────────────────────────
  # When ggplot has facets, ggplotly creates one trace per facet panel and
  # names them "(estimator_name,1)", "(estimator_name,2)" etc.
  # This helper strips the parentheses and numeric suffix, then hides
  # duplicate legend entries so each estimator appears only once.
  fix_legend <- function(pl) {
    seen <- c()
    for (i in seq_along(pl$x$data)) {
      nm <- pl$x$data[[i]]$name
      if (!is.null(nm)) {
        clean <- gsub("^\\(|,\\d+\\)$|,\\d+$", "", nm)
        pl$x$data[[i]]$name <- clean
        if (clean %in% seen) {
          pl$x$data[[i]]$showlegend <- FALSE
        } else {
          seen <- c(seen, clean)
          pl$x$data[[i]]$showlegend <- TRUE
        }
      }
    }
    pl
  }
  
  
  # ════════════════════════════════════════════════════════════════════════
  #  TAB 1 — PERFORMANCE METRICS
  #  Y-axis: mean metric value (sensitivity / specificity / F1 / MCC / Jaccard)
  #  X-axis: missingness (0–0.9) or sample size (250/500/750)
  #  Each line = one estimator, averaged over 200 reps
  #  Optional SD ribbon shows variability across reps
  # ════════════════════════════════════════════════════════════════════════
  output$metricPlot <- renderPlotly({
    dat <- filt_summary()
    validate(need(nrow(dat) > 0, "No data — select at least one estimator."))
    
    p <- ggplot(dat, aes(.data[[input$xaxis]], mean_value,
                         color = estimator, group = estimator,
                         # Custom hover tooltip (IMPortant)
                         text = paste0("Estimator: ", estimator,
                                       "<br>", ifelse(input$xaxis=="missing","Missing","N"),
                                       ": ", .data[[input$xaxis]],
                                       "<br>Mean ", input$metric, ": ", round(mean_value, 3),
                                       "<br>SD: ", round(sd_value, 3),
                                       "<br>N people: ", n_people))) +
      # SD ribbon (optional): shaded band of mean ± 1 SD, clamped to [0,1]
      { if (input$show_ribbon)
        geom_ribbon(aes(ymin = pmax(0, mean_value - sd_value),
                        ymax = pmin(1, mean_value + sd_value),
                        fill = estimator), alpha = 0.12, color = NA)
        else NULL } +
      geom_line(linewidth = 1) +
      { if (input$show_points) geom_point(size = 2.5) else NULL } +
      scale_color_manual(values = pal) +
      scale_fill_manual(values  = pal) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(x     = ifelse(input$xaxis == "missing", "Missingness", "Sample size (N)"),
           y     = paste("Mean", toupper(input$metric)),
           color = "Estimator", fill = "Estimator",
           title = paste0("PERFORMANCE — ", toupper(input$metric),
                          ifelse(input$xaxis=="missing",
                                 paste0("  (N = ", input$n_people_sel, ")"), ""))) +
      dark_theme() +
      guides(color = guide_legend(nrow = 3))
    
    ggplotly(apply_xaxis(p, input$xaxis), tooltip = "text") %>%
      layout(paper_bgcolor = "#0d1117", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3"),
             legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE, displaylogo = FALSE) %>%  # always show toolbar
      fix_legend()
  })
  
  # Sortable/searchable table of the same summarised data
  output$metricTable <- renderDT({
    filt_summary() %>%
      select(estimator, n_people, missing, mean_value, sd_value, n_reps) %>%
      mutate(across(c(mean_value, sd_value), ~round(., 4))) %>%
      arrange(desc(mean_value)) %>%
      datatable(rownames = FALSE,
                options  = list(pageLength = 15, dom = "ftp", fontSize = "12px"),
                colnames = c("Estimator","N","Missing","Mean","SD","Reps"))
  })
  
  
  # ════════════════════════════════════════════════════════════════════════
  #  TAB 2 — EDGE COUNT DIAGNOSTICS
  #  Shows mean number of edges in the estimated network vs the true count.
  #  True network has ~6 edges (8 nodes, 20% density =  8*7/2 * 0.2 ≈ 5.6).
  #  Max possible edges = 28 (complete graph for 8 nodes).
  #
  #  Reference lines:
  #    Blue  = 0  = empty network (method predicts no edges at all)
  #    Purple = 28 = complete network (method connects everything)
  #
  #  Why this matters: a method that collapses to 0 edges will have
  #  perfect specificity (no FP) but terrible sensitivity (all FN missed).
  #  This tab helps diagnose that pathological behaviour.
  # ════════════════════════════════════════════════════════════════════════
  output$edgePlot <- renderPlotly({
    dat <- filt_raw() %>%
      group_by(estimator, n_people, missing) %>%
      summarise(est = mean(n_estimated_edges, na.rm=TRUE),
                tru = mean(n_true_edges, na.rm=TRUE), .groups="drop")
    validate(need(nrow(dat) > 0, "No data."))
    
    true_val <- mean(dat$tru)
    lx <- ifelse(input$xaxis == "missing", 0.02, 255)
    
    p <- ggplot(dat, aes(.data[[input$xaxis]], est,
                         color = estimator, group = estimator,
                         text = paste0("Estimator: ", estimator,
                                       "<br>", ifelse(input$xaxis=="missing","Missing","N"),
                                       ": ", .data[[input$xaxis]],
                                       "<br>Mean edges: ", round(est, 2)))) +
      geom_hline(yintercept = 0,        linetype="dotted", color="#4cc9f0", linewidth=.6, show.legend=FALSE) +
      geom_hline(yintercept = 28,       linetype="dotted", color="#a05195", linewidth=.6, show.legend=FALSE) +
      geom_hline(yintercept = true_val, linetype="dashed", color="#f95d6a", linewidth=.8, show.legend=FALSE) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      annotate("text", x=lx, y=1.8,          label="Empty (0)",     color="#4cc9f0", hjust=0, size=3) +
      annotate("text", x=lx, y=26.5,         label="Complete (28)", color="#a05195", hjust=0, size=3) +
      annotate("text", x=lx, y=true_val+1.5, label=paste0("True (~", round(true_val),")"),
               color="#f95d6a", hjust=0, size=3) +
      scale_color_manual(values = pal) +
      scale_y_continuous(limits=c(0,30), breaks=seq(0,28,4)) +
      labs(x     = ifelse(input$xaxis=="missing","Missingness","Sample size (N)"),
           y     = "Mean estimated edges",
           color = "Estimator",
           title = paste0("EDGE COUNT DIAGNOSTICS",
                          ifelse(input$xaxis=="missing",
                                 paste0("  (N = ", input$n_people_sel, ")"), ""))) +
      dark_theme() +
      guides(color = guide_legend(nrow = 3))
    
    ggplotly(apply_xaxis(p, input$xaxis), tooltip = "text") %>%
      layout(paper_bgcolor = "#0d1117", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3"),
             legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE, displaylogo = FALSE) %>%
      fix_legend()
  })
  
  
  
  # ════════════════════════════════════════════════════════════════════════
  #  TAB 3 — FALSE POSITIVES vs FALSE NEGATIVES (SEPARATE)
  #  Solid lines  = FP: edges estimated that do NOT exist in the true network
  #  Dashed lines = FN: true edges that the method FAILED to recover
  #
  #  Key pattern to watch for (network collapse):
  #    At high missingness, some methods (e.g. EM_EBIC, mean imputation)
  #    predict an empty network. This makes FP = 0 (looks good!) but
  #    FN = ~6 (all true edges missed). Cross-reference with Tab 2
  #    to confirm: if edge count = 0, the specificity in Tab 1 is spurious.
  # ════════════════════════════════════════════════════════════════════════
  output$fpfnPlot <- renderPlotly({
    dat <- filt_raw() %>%
      group_by(estimator, n_people, missing) %>%
      summarise(FP = mean(FP, na.rm=TRUE),
                FN = mean(FN, na.rm=TRUE), .groups="drop")
    validate(need(nrow(dat) > 0, "No data."))
    
    # Two separate geom_line layers (FP solid, FN dashed) with the same colour
    # this is to keep the legend clean (one entry per estimator, not two).
    show_fp <- "FP" %in% input$fpfn_show
    show_fn <- "FN" %in% input$fpfn_show
    validate(need(show_fp | show_fn, "Select at least FP or FN above."))
    
    p <- ggplot(dat, aes(x = .data[[input$xaxis]], color = estimator, group = estimator)) +
      { if (show_fp) geom_line(aes(y = FP,
                                   text = paste0("Estimator: ", estimator,
                                                 "<br>Type: FP (False Positive)",
                                                 "<br>", ifelse(input$xaxis=="missing","Missing","N"),
                                                 ": ", .data[[input$xaxis]],
                                                 "<br>Count: ", round(FP, 2))),
                               linewidth = 1, linetype = "solid") else NULL } +
      { if (show_fn) geom_line(aes(y = FN,
                                   text = paste0("Estimator: ", estimator,
                                                 "<br>Type: FN (False Negative)",
                                                 "<br>", ifelse(input$xaxis=="missing","Missing","N"),
                                                 ": ", .data[[input$xaxis]],
                                                 "<br>Count: ", round(FN, 2))),
                               linewidth = 1, linetype = "dashed") else NULL } +
      { if (show_fp) geom_point(aes(y = FP), size = 2) else NULL } +
      { if (show_fn) geom_point(aes(y = FN), size = 2, shape = 1) else NULL } +
      scale_color_manual(values = pal) +
      labs(x        = ifelse(input$xaxis=="missing","Missingness","Sample size (N)"),
           y        = "Mean count",
           color    = "Estimator",
           title    = paste0("FALSE POSITIVES vs FALSE NEGATIVES",
                             ifelse(input$xaxis=="missing",
                                    paste0("  (N = ", input$n_people_sel, ")"), "")),
           subtitle = paste0(
             if(show_fp) "Solid + filled = FP" else "",
             if(show_fp & show_fn) "  ·  " else "",
             if(show_fn) "Dashed + open = FN" else "")) +
      dark_theme() +
      guides(color = guide_legend(nrow = 3))
    
    ggplotly(apply_xaxis(p, input$xaxis), tooltip = "text") %>%
      layout(paper_bgcolor = "#0d1117", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3"),
             legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE, displaylogo = FALSE) %>%
      fix_legend()
  })
  
  
  # ════════════════════════════════════════════════════════════════════════
  #  TAB 4 — TOTAL CLASSIFICATION ERRORS (FP + FN combined)
  #  Simpler summary of overall edge recovery error.
  # ════════════════════════════════════════════════════════════════════════
  output$errorPlot <- renderPlotly({
    dat <- filt_raw() %>%
      group_by(estimator, n_people, missing) %>%
      summarise(errors = mean(FP+FN, na.rm=TRUE), .groups="drop")
    validate(need(nrow(dat) > 0, "No data."))
    
    p <- ggplot(dat, aes(.data[[input$xaxis]], errors,
                         color=estimator, group=estimator,
                         text = paste0("Estimator: ", estimator,
                                       "<br>", ifelse(input$xaxis=="missing","Missing","N"),
                                       ": ", .data[[input$xaxis]],
                                       "<br>FP+FN: ", round(errors, 2)))) +
      geom_line(linewidth=1) +
      geom_point(size=2.5) +
      scale_color_manual(values=pal) +
      labs(x     = ifelse(input$xaxis=="missing","Missingness","Sample size (N)"),
           y     = "Mean FP + FN",
           color = "Estimator",
           title = paste0("TOTAL CLASSIFICATION ERRORS (FP + FN)",
                          ifelse(input$xaxis=="missing",
                                 paste0("  (N = ", input$n_people_sel, ")"), ""))) +
      dark_theme() +
      guides(color=guide_legend(nrow=3))
    
    ggplotly(apply_xaxis(p, input$xaxis), tooltip = "text") %>%
      layout(paper_bgcolor = "#0d1117", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3"),
             legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE, displaylogo = FALSE) %>%
      fix_legend()
  })
  
  
  
  # ════════════════════════════════════════════════════════════════════════
  #  TAB 5 — METHOD FAILURE RATES
  #  % of the 200 reps where the method threw an R error (error = TRUE
  #  in the parSim output) and returned NA instead of a network.
  # ════════════════════════════════════════════════════════════════════════
  output$failurePlot <- renderPlotly({
    dat <- all_batches
    if (!"error" %in% names(dat)) dat$error <- FALSE
    if (length(input$estimators)>0) dat <- filter(dat, estimator %in% input$estimators)
    if (input$xaxis == "missing") dat <- filter(dat, n_people == as.numeric(input$n_people_sel))
    
    dat <- dat %>%
      group_by(estimator, n_people, missing) %>%
      summarise(failure_rate=mean(error,na.rm=TRUE)*100, .groups="drop")
    validate(need(nrow(dat)>0,"No data."))
    
    p <- ggplot(dat, aes(.data[[input$xaxis]], failure_rate,
                         color=estimator, group=estimator,
                         text = paste0("Estimator: ", estimator,
                                       "<br>", ifelse(input$xaxis=="missing","Missing","N"),
                                       ": ", .data[[input$xaxis]],
                                       "<br>Failure rate: ", round(failure_rate, 1), "%"))) +
      geom_line(linewidth=1) +
      geom_point(size=2.5) +
      scale_color_manual(values=pal) +
      scale_y_continuous(limits=c(0,100)) +
      labs(x        = ifelse(input$xaxis=="missing","Missingness","Sample size (N)"),
           y        = "Failure rate (%)",
           color    = "Estimator",
           title    = paste0("METHOD FAILURE RATES",
                             ifelse(input$xaxis=="missing",
                                    paste0("  (N = ", input$n_people_sel, ")"), "")),
           subtitle = "% of reps where the method threw an error and returned NA") +
      dark_theme() +
      guides(color=guide_legend(nrow=3))
    
    ggplotly(apply_xaxis(p, input$xaxis), tooltip = "text") %>%
      layout(paper_bgcolor = "#0d1117", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3"),
             legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE, displaylogo = FALSE) %>%
      fix_legend()
  })
  
  
  
  # ════════════════════════════════════════════════════════════════════════
  #  TAB 6 — HEATMAP
  #  Rows    = estimators (all selected in sidebar)
  #  Columns = missingness levels (0 to 0.9)
  #  Fill    = mean metric value (colour gradient: dark = low, bright = high)
  #  Numbers = exact mean value printed in each cell
  #
  #  Has its own metric and N dropdowns (top of tab), independent of sidebar.
  #  Useful as a summary figure 
  #  Note: uses renderPlot (not plotly) so the plotly toolbar isn't available,
  # ════════════════════════════════════════════════════════════════════════
  output$heatmapPlot <- renderPlot({
    dat <- all_batches_summary %>%
      filter(metric == input$hm_metric, n_people == as.numeric(input$hm_n)) %>%
      { if(length(input$estimators)>0) filter(., estimator %in% input$estimators) else . }
    validate(need(nrow(dat)>0,"No data."))
    
    ggplot(dat, aes(x=factor(missing), y=estimator, fill=mean_value)) +
      geom_tile(color="#0d1117", linewidth=0.5) +
      geom_text(aes(label=round(mean_value,2)), color="white", size=3, family="DM Mono") +
      scale_fill_gradientn(
        colors = c("#161b22","#003f5c","#2f4b7c","#58a6ff","#4cc9f0","#7bc8a4"),
        limits = c(0,1),
        name   = toupper(input$hm_metric)) +
      labs(x        = "Missingness",
           y        = NULL,
           title    = paste("HEATMAP —", toupper(input$hm_metric)),
           subtitle = paste("N =", input$hm_n,
                            " · Rows = estimators · Columns = missingness level")) +
      dark_theme() +
      theme(legend.position  = "right",
            axis.text.y      = element_text(family="DM Mono", size=8))
  }, bg="#0d1117")
  
}

# ── Launch ───────────────────────────────────────────────────────────────────
shinyApp(ui, server)
