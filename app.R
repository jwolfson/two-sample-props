# =============================================================================
#  Two Proportions: Estimates, Estimators & Sampling Distributions
#  An interactive Shiny app for teaching inference on p1 - p2
#
#  Tabs
#   1. Parameter vs. estimate   - one sample, what is fixed and what is random
#   2. Sampling distribution    - repeated sampling of the estimator p̂1 - p̂2
#   3. Confidence intervals     - repeated intervals and long-run coverage
#   4. Hypothesis test          - pooled z-test, permutation null, error rates
#   5. Guide                    - notation and suggested classroom activities
#
#  Packages: shiny, bslib, ggplot2
# =============================================================================

library(shiny)
library(bslib)
library(ggplot2)

# ---- Palette -----------------------------------------------------------------
col_g1    <- "#0F766E"  # group 1 (teal)
col_g2    <- "#C2410C"  # group 2 (burnt orange)
col_fail  <- "#E5E7EB"  # failures in waffle
col_param <- "#111827"  # the parameter (truth)
col_hit   <- "#2563EB"  # interval captured parameter
col_miss  <- "#DC2626"  # interval missed / rejection region
col_null  <- "#6B7280"  # null distribution
col_sim   <- "#93C5FD"  # simulated histogram fill

MAX_STORE <- 20000      # max simulated samples kept in memory per tab

# ---- Helpers -----------------------------------------------------------------
fmt  <- function(x, d = 3) formatC(x, format = "f", digits = d)
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, "< 0.0001", fmt(p, 4)))
pct  <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")

true_se <- function(p1, p2, n1, n2) sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)

draw_samples <- function(k, p1, p2, n1, n2) {
  data.frame(x1 = rbinom(k, n1, p1), x2 = rbinom(k, n2, p2))
}

# Confidence interval for p1 - p2 from counts (vectorised)
ci_two_prop <- function(x1, x2, n1, n2, conf = 0.95, method = "wald") {
  zstar <- qnorm(1 - (1 - conf) / 2)
  if (method == "ac") {               # Agresti-Caffo: add 1 success + 1 failure per group
    pt1 <- (x1 + 1) / (n1 + 2); pt2 <- (x2 + 1) / (n2 + 2)
    centre <- pt1 - pt2
    se <- sqrt(pt1 * (1 - pt1) / (n1 + 2) + pt2 * (1 - pt2) / (n2 + 2))
  } else {                            # Wald: plug the sample proportions into the SE
    ph1 <- x1 / n1; ph2 <- x2 / n2
    centre <- ph1 - ph2
    se <- sqrt(ph1 * (1 - ph1) / n1 + ph2 * (1 - ph2) / n2)
  }
  data.frame(est = x1 / n1 - x2 / n2, centre = centre, se = se, zstar = zstar,
             lo = pmax(-1, centre - zstar * se), hi = pmin(1, centre + zstar * se))
}

# Pooled two-proportion z statistic (vectorised)
z_pooled <- function(x1, x2, n1, n2) {
  pp  <- (x1 + x2) / (n1 + n2)
  se0 <- sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
  d   <- x1 / n1 - x2 / n2
  z   <- ifelse(se0 > 0, d / se0, 0)   # all successes or all failures: no evidence
  list(d = d, pp = pp, se0 = se0, z = z)
}

p_value <- function(z, alt) {
  switch(alt,
         two.sided = 2 * pnorm(-abs(z)),
         greater   = pnorm(z, lower.tail = FALSE),
         less      = pnorm(z))
}

alt_symbol <- function(alt) switch(alt, two.sided = "\\neq", greater = ">", less = "<")

approx_power <- function(p1, p2, n1, n2, alpha, alt) {
  pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
  se0  <- sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
  se1  <- true_se(p1, p2, n1, n2)
  del  <- p1 - p2
  switch(alt,
    two.sided = { c <- qnorm(1 - alpha / 2) * se0
                  pnorm((-c - del) / se1) + pnorm((c - del) / se1, lower.tail = FALSE) },
    greater   = pnorm((qnorm(1 - alpha) * se0 - del) / se1, lower.tail = FALSE),
    less      = pnorm((-qnorm(1 - alpha) * se0 - del) / se1))
}

# Histogram bin width that respects the lattice of possible values when n1 == n2
nice_binwidth <- function(se, n1, n2, step = if (n1 == n2) 1 / n1 else NULL) {
  if (!is.null(step)) {
    k <- max(1, round(se / 3.5 / step))
    if (k %% 2 == 0) k <- k + 1       # odd multiple -> bin edges fall between lattice points
    list(width = k * step, center = 0)
  } else {
    list(width = max(se / 4, 1e-4), center = 0)
  }
}

waffle_df <- function(x, n, label, ncol = 20) {
  i <- seq_len(n) - 1
  data.frame(group = label, col = i %% ncol, row = -(i %/% ncol), success = i < x)
}

base_theme <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "#4B5563", size = 12),
        legend.position = "top")

key_ideas <- function(...) {
  card(class = "border-0 bg-light",
       card_body(tags$h6(icon("lightbulb"), " Key ideas", class = "fw-bold mb-2"),
                 tags$ul(class = "mb-0", ...)))
}

btn <- function(id, label, cls = "btn-outline-primary") {
  actionButton(id, label, class = paste("btn-sm", cls))
}

header_row <- function(title, ...) {
  card_header(div(class = "d-flex flex-wrap align-items-center gap-2",
                  span(class = "fw-semibold me-auto", title), ...))
}

presets <- list(
  equal    = list(p1 = 0.50, p2 = 0.50, n1 = 100, n2 = 100),
  moderate = list(p1 = 0.60, p2 = 0.45, n1 = 100, n2 = 100),
  big_n    = list(p1 = 0.60, p2 = 0.45, n1 = 1000, n2 = 1000),
  rare     = list(p1 = 0.10, p2 = 0.04, n1 = 30,  n2 = 30),
  unequal  = list(p1 = 0.35, p2 = 0.25, n1 = 250, n2 = 50)
)

# =============================================================================
#  UI
# =============================================================================
ui <- page_navbar(
  title = "Two Proportions: Estimates & Estimators",
  fillable = FALSE,
  theme = bs_theme(version = 5, primary = "#1D4ED8", font_scale = 0.95),
  header = tagList(
    withMathJax(),
    tags$style(HTML("
      .bslib-value-box .value-box-value { font-size: 1.6rem; }
      .worked { font-size: 0.95rem; line-height: 1.7; }
      .worked .step { margin-bottom: .35rem; }
      .callout { border-left: 4px solid #1D4ED8; background: #EFF6FF; padding: .6rem .8rem; border-radius: .25rem; }
      .callout-warn { border-left-color: #C2410C; background: #FFF7ED; }
      .sidebar-title { font-weight: 700; }
      html body nav.navbar, html body nav.navbar.navbar-default { background-color: #1E3A8A !important; --bslib-navbar-default-bg: #1E3A8A; }
      html body nav.navbar .navbar-brand, html body nav.navbar .nav-link { color: rgba(255,255,255,.8) !important; }
      html body nav.navbar .nav-link.active, html body nav.navbar .nav-link:hover, html body nav.navbar .navbar-brand { color: #fff !important; }
      .navbar .nav-link.active { border-bottom: 3px solid #93C5FD; }
    "))
  ),
  sidebar = sidebar(
    width = 300, title = "The population (the truth)",
    selectInput("preset", "Load a scenario",
      c("Choose a scenario…" = "",
        "No difference: p₁ = p₂ = 0.50, n = 100"  = "equal",
        "Moderate difference, n = 100"            = "moderate",
        "Same difference, n = 1000"               = "big_n",
        "Rare outcome, small samples"             = "rare",
        "Unequal sample sizes"                    = "unequal")),
    layout_columns(col_widths = c(6, 6),
      textInput("g1", "Group 1 name", "Group 1"),
      textInput("g2", "Group 2 name", "Group 2")),
    sliderInput("p1", HTML("True proportion in group 1, \\(p_1\\)"), 0.01, 0.99, 0.60, 0.01),
    sliderInput("p2", HTML("True proportion in group 2, \\(p_2\\)"), 0.01, 0.99, 0.45, 0.01),
    layout_columns(col_widths = c(6, 6),
      numericInput("n1", HTML("\\(n_1\\)"), 100, min = 5, max = 5000, step = 5),
      numericInput("n2", HTML("\\(n_2\\)"), 100, min = 5, max = 5000, step = 5)),
    uiOutput("param_summary"),
    tags$small(class = "text-muted",
      "In real research these values are unknown. Here you play the role of nature, ",
      "so you can see how estimates behave. Changing them clears every simulation.")
  ),

  # ---------------------------------------------------------------- Tab 1 ----
  nav_panel("1 · Parameter vs. estimate", icon = icon("bullseye"),
    layout_columns(col_widths = c(4, 4, 4), fill = FALSE,
      value_box(title = "Parameter (fixed, unknown in practice)",
                value = uiOutput("vb_param"), showcase = icon("lock"),
                theme = "dark", p("The true difference in the populations")),
      value_box(title = "Estimator (a rule — random before sampling)",
                value = HTML("\\(\\hat p_1 - \\hat p_2\\)"), showcase = icon("shuffle"),
                theme = "primary", p("Apply it to any sample to get an estimate")),
      value_box(title = "Estimate (one number from this sample)",
                value = textOutput("vb_est"), showcase = icon("hashtag"),
                theme = "secondary", uiOutput("vb_est_detail"))
    ),
    layout_columns(col_widths = c(7, 5),
      card(full_screen = TRUE,
        header_row("The sample you just drew",
                   btn("t1_draw", "Draw a new sample", "btn-primary")),
        plotOutput("t1_waffle", height = 330),
        uiOutput("t1_text")),
      card(full_screen = TRUE,
        header_row("Every estimate you have drawn so far",
                   btn("t1_reset", "Clear", "btn-outline-secondary")),
        plotOutput("t1_strip", height = 330),
        uiOutput("t1_strip_text"))
    ),
    key_ideas(
      tags$li(HTML("The <b>parameter</b> \\(p_1 - p_2\\) is a fixed number describing the populations. It does not change when you click the button.")),
      tags$li(HTML("The <b>estimator</b> \\(\\hat p_1 - \\hat p_2\\) is a <i>random variable</i>: a recipe whose output depends on which individuals happen to be sampled.")),
      tags$li(HTML("An <b>estimate</b> is the single number the estimator produces for one particular sample. Different samples → different estimates.")),
      tags$li("Click the button several times. Do the estimates cluster around the parameter? How far off can a single estimate be?")
    )
  ),

  # ---------------------------------------------------------------- Tab 2 ----
  nav_panel("2 · Sampling distribution", icon = icon("chart-column"),
    layout_columns(col_widths = c(8, 4),
      card(full_screen = TRUE,
        header_row(HTML("Sampling distribution of \\(\\hat p_1 - \\hat p_2\\)"),
                   btn("t2_1", "+1 sample"), btn("t2_10", "+10"),
                   btn("t2_100", "+100"), btn("t2_1000", "+1,000", "btn-primary"),
                   btn("t2_reset", "Reset", "btn-outline-secondary")),
        plotOutput("t2_hist", height = 420),
        layout_columns(col_widths = c(4, 4, 4), fill = FALSE,
          checkboxInput("t2_curve", "Show theoretical normal curve", TRUE),
          checkboxInput("t2_param", "Show the parameter", TRUE),
          checkboxInput("t2_band", "Shade ±2 SE around parameter", FALSE))),
      card(
        card_header("Simulation vs. theory"),
        tableOutput("t2_table"),
        uiOutput("t2_formula"),
        uiOutput("t2_conditions"))
    ),
    key_ideas(
      tags$li(HTML("Each bar is built from many <i>estimates</i>; the whole histogram is a picture of the <i>estimator</i>'s behavior.")),
      tags$li(HTML("<b>Center:</b> the estimates average out to the parameter — \\(\\hat p_1 - \\hat p_2\\) is <b>unbiased</b> for \\(p_1 - p_2\\).")),
      tags$li(HTML("<b>Spread:</b> the standard deviation of the sampling distribution is the <b>standard error</b>. Try quadrupling both \\(n\\)'s — the SE is cut in half.")),
      tags$li(HTML("<b>Shape:</b> when every group has at least 10 successes and 10 failures expected, the shape is close to normal. Try the “rare outcome” scenario to see when it is not.")),
      tags$li(HTML("About 95% of estimates land within 2 SE of the parameter. Flip that around and you have the idea behind a confidence interval."))
    )
  ),

  # ---------------------------------------------------------------- Tab 3 ----
  nav_panel("3 · Confidence intervals", icon = icon("arrows-left-right"),
    layout_columns(col_widths = c(3, 9),
      card(
        card_header("Interval settings"),
        sliderInput("conf", "Confidence level", 0.50, 0.99, 0.95, 0.01),
        radioButtons("ci_method", "Method",
          c("Wald (plug in \\(\\hat p\\)'s)" = "wald",
            "Agresti–Caffo (add 1 success & 1 failure per group)" = "ac")),
        tags$hr(),
        div(class = "d-grid gap-2",
            btn("t3_1", "Draw 1 sample & interval"),
            btn("t3_20", "Draw 20"),
            btn("t3_100", "Draw 100", "btn-primary"),
            btn("t3_reset", "Reset", "btn-outline-secondary")),
        tags$small(class = "text-muted mt-2",
          "Changing the confidence level or method recomputes the intervals for the same samples.")),
      div(
        layout_columns(col_widths = c(4, 4, 4), fill = FALSE,
          value_box("Intervals drawn", textOutput("t3_n"), theme = "light"),
          value_box("Captured the parameter", textOutput("t3_cov"), theme = "primary"),
          value_box("Nominal confidence level", textOutput("t3_nom"), theme = "light")),
        layout_columns(col_widths = c(7, 5),
          card(full_screen = TRUE,
            card_header("The most recent 100 intervals (newest on top)"),
            plotOutput("t3_plot", height = 470)),
          card(card_header("How the newest interval was built"),
               uiOutput("t3_worked")))
      )
    ),
    key_ideas(
      tags$li(HTML("The interval is <b>estimate ± (critical value) × (estimated SE)</b>. Both the center and the width come from the sample, so the <i>interval</i> is random — the parameter is not.")),
      tags$li(HTML("The true SE uses \\(p_1, p_2\\), which we don't know, so we plug in the estimates \\(\\hat p_1, \\hat p_2\\). This <i>estimated</i> SE changes from sample to sample too.")),
      tags$li(HTML("“95% confidence” describes the <b>method</b>: in the long run about 95% of intervals built this way capture the parameter. Any single interval either does or doesn't.")),
      tags$li(HTML("Try the “rare outcome” scenario with the Wald method: coverage falls short of the nominal level. Switch to Agresti–Caffo and compare."))
    )
  ),

  # ---------------------------------------------------------------- Tab 4 ----
  nav_panel("4 · Hypothesis test", icon = icon("scale-balanced"),
    layout_columns(col_widths = c(3, 9),
      card(
        card_header("Test settings"),
        withMathJax(HTML("<p class='mb-1'>\\(H_0: p_1 - p_2 = 0\\)</p>")),
        selectInput("alt", HTML("Alternative \\(H_A\\)"),
          c("p₁ − p₂ ≠ 0 (two-sided)" = "two.sided",
            "p₁ − p₂ > 0" = "greater",
            "p₁ − p₂ < 0" = "less")),
        sliderInput("alpha", HTML("Significance level \\(\\alpha\\)"), 0.01, 0.20, 0.05, 0.01),
        tags$small(class = "text-muted",
          "These settings apply to both the single test and the long-run simulation.")),
      navset_card_tab(
        full_screen = TRUE,
        nav_panel("Test one sample",
          layout_columns(col_widths = c(4, 8),
            div(
              radioButtons("t4_src", "Data",
                c("Draw a sample from the population" = "pop",
                  "Enter my own counts" = "own")),
              conditionalPanel("input.t4_src == 'pop'",
                btn("t4_draw", "Draw a new sample", "btn-primary"),
                checkboxInput("t4_truth", "Reveal the true difference", FALSE)),
              conditionalPanel("input.t4_src == 'own'",
                layout_columns(col_widths = c(6, 6),
                  numericInput("o_x1", HTML("Successes \\(x_1\\)"), 42, min = 0),
                  numericInput("o_n1", HTML("\\(n_1\\)"), 80, min = 1),
                  numericInput("o_x2", HTML("Successes \\(x_2\\)"), 30, min = 0),
                  numericInput("o_n2", HTML("\\(n_2\\)"), 85, min = 1))),
              tags$hr(),
              checkboxInput("t4_perm", "Overlay simulated null distribution (2,000 random re-shuffles of group labels)", FALSE),
              uiOutput("t4_worked")),
            div(plotOutput("t4_plot", height = 420),
                uiOutput("t4_plot_caption")))
        ),
        nav_panel("Long run: error rates & power",
          layout_columns(col_widths = c(8, 4),
            div(
              div(class = "d-flex flex-wrap gap-2 mb-2",
                  btn("t5_100", "Run 100 tests"),
                  btn("t5_1000", "Run 1,000 tests", "btn-primary"),
                  btn("t5_reset", "Reset", "btn-outline-secondary")),
              plotOutput("t5_plot", height = 400)),
            uiOutput("t5_summary"))
        )
      )
    ),
    key_ideas(
      tags$li(HTML("The test asks: <i>if</i> \\(H_0\\) were true, how unusual would an estimate like ours be? So we need the sampling distribution of \\(\\hat p_1 - \\hat p_2\\) <b>assuming \\(H_0\\)</b>.")),
      tags$li(HTML("Under \\(H_0\\) the two groups share one proportion, so we estimate it with the <b>pooled</b> proportion \\(\\hat p = (x_1 + x_2)/(n_1 + n_2)\\). That is why the test's SE differs from the CI's SE.")),
      tags$li(HTML("The p-value is a probability about the <i>estimator</i> (over repeated samples), computed using the one <i>estimate</i> we observed. It is not the probability that \\(H_0\\) is true.")),
      tags$li(HTML("In the long-run tab, set \\(p_1 = p_2\\): the rejection rate settles near \\(\\alpha\\) (Type I errors). Make them differ: the rejection rate is the <b>power</b>."))
    )
  ),

  # ---------------------------------------------------------------- Guide ----
  nav_panel("Guide", icon = icon("book-open"),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Notation"),
        withMathJax(HTML("
<table class='table table-sm'>
<thead><tr><th>Symbol</th><th>Meaning</th><th>Type</th></tr></thead>
<tbody>
<tr><td>\\(p_1, p_2\\)</td><td>True proportions of successes in populations 1 and 2</td><td>Parameters (fixed)</td></tr>
<tr><td>\\(p_1 - p_2\\)</td><td>True difference</td><td>Parameter (fixed)</td></tr>
<tr><td>\\(\\hat p_1 = X_1/n_1\\)</td><td>Sample proportion in group 1</td><td>Estimator (random)</td></tr>
<tr><td>\\(\\hat p_1 - \\hat p_2\\)</td><td>Difference in sample proportions</td><td>Estimator (random)</td></tr>
<tr><td>\\(0.53 - 0.41 = 0.12\\)</td><td>The value computed from one sample</td><td>Estimate (a number)</td></tr>
<tr><td>\\(SE = \\sqrt{\\tfrac{p_1(1-p_1)}{n_1} + \\tfrac{p_2(1-p_2)}{n_2}}\\)</td><td>Standard deviation of the sampling distribution</td><td>Depends on unknown parameters</td></tr>
<tr><td>\\(\\widehat{SE} = \\sqrt{\\tfrac{\\hat p_1(1-\\hat p_1)}{n_1} + \\tfrac{\\hat p_2(1-\\hat p_2)}{n_2}}\\)</td><td>Estimated SE used for confidence intervals</td><td>Estimator / estimate</td></tr>
<tr><td>\\(SE_0 = \\sqrt{\\hat p(1-\\hat p)\\left(\\tfrac{1}{n_1}+\\tfrac{1}{n_2}\\right)}\\)</td><td>SE under \\(H_0\\), using the pooled \\(\\hat p\\)</td><td>Used for the z-test</td></tr>
</tbody></table>"))),
      card(card_header("Suggested classroom activities"),
        tags$ol(
          tags$li(tags$b("Fixed vs. random (Tab 1). "),
                  "Have each student draw one sample and report their estimate. Collect them on the board, then compare with the app's strip plot. Ask: which quantity was the same for everyone?"),
          tags$li(tags$b("What controls the spread? (Tab 2). "),
                  "Simulate 1,000 samples, note the SE, then double n₁ and n₂. Predict the new SE before simulating. Repeat with p₁, p₂ near 0.5 vs. near 0.05."),
          tags$li(tags$b("What does 95% mean? (Tab 3). "),
                  "Draw 100 intervals. Count the red ones. Then slide the confidence level to 80% and 99% without drawing new samples — which intervals change color?"),
          tags$li(tags$b("When methods break (Tab 3). "),
                  "Load the “rare outcome” scenario and run several hundred Wald intervals. Is the coverage 95%? Switch to Agresti–Caffo."),
          tags$li(tags$b("Two SEs (Tab 4). "),
                  "For one sample, compare the pooled SE used in the test with the unpooled SE used in the interval. Why are they different?"),
          tags$li(tags$b("Errors and power (Tab 4, long run). "),
                  "With p₁ = p₂, run 1,000 tests and compare the rejection rate with α. Then create a small true difference and investigate how power changes with n.")
        ),
        tags$p(class = "text-muted small mb-0",
          "Theory reference: the normal approximation for p̂₁ − p̂₂ is reasonable when n₁p₁, n₁(1−p₁), n₂p₂, n₂(1−p₂) are all at least about 10.")
      )
    )
  )
)

# =============================================================================
#  Server
# =============================================================================
server <- function(input, output, session) {

  # ---- Population --------------------------------------------------------
  observeEvent(input$preset, {
    pr <- presets[[input$preset]]
    if (is.null(pr)) return()
    updateSliderInput(session, "p1", value = pr$p1)
    updateSliderInput(session, "p2", value = pr$p2)
    updateNumericInput(session, "n1", value = pr$n1)
    updateNumericInput(session, "n2", value = pr$n2)
  }, ignoreInit = TRUE)

  pop <- reactive({
    n1 <- suppressWarnings(as.integer(input$n1)); n2 <- suppressWarnings(as.integer(input$n2))
    validate(need(!is.na(n1) && n1 >= 2 && n1 <= 5000, "n₁ must be between 2 and 5,000."),
             need(!is.na(n2) && n2 >= 2 && n2 <= 5000, "n₂ must be between 2 and 5,000."))
    g1 <- trimws(input$g1); g2 <- trimws(input$g2)
    if (g1 == "") g1 <- "Group 1"
    if (g2 == "") g2 <- "Group 2"
    if (g1 == g2) g2 <- paste(g2, "(2)")
    list(p1 = input$p1, p2 = input$p2, n1 = n1, n2 = n2,
         diff = input$p1 - input$p2,
         se = true_se(input$p1, input$p2, n1, n2), g1 = g1, g2 = g2)
  })

  output$param_summary <- renderUI({
    P <- pop()
    withMathJax(div(class = "callout small",
      HTML(sprintf("Parameter: \\(p_1 - p_2 = %s\\)<br>True SE of \\(\\hat p_1 - \\hat p_2\\): %s",
                   fmt(P$diff, 2), fmt(P$se, 4)))))
  })

  rv <- reactiveValues(t1 = NULL, t2 = NULL, t3 = NULL, t4 = NULL, t5 = NULL)

  # Any change to the population clears every simulation
  observeEvent(list(input$p1, input$p2, input$n1, input$n2), {
    req(pop())
    P <- pop()
    s <- draw_samples(1, P$p1, P$p2, P$n1, P$n2)
    s$jit <- runif(1)
    rv$t1 <- s
    rv$t2 <- NULL; rv$t3 <- NULL; rv$t5 <- NULL
    rv$t4 <- draw_samples(1, P$p1, P$p2, P$n1, P$n2)
  }, priority = 100)   # run before outputs so nothing renders stale samples with new n's

  add_rows <- function(old, new) {
    out <- if (is.null(old)) new else rbind(old, new)
    if (nrow(out) > MAX_STORE) out <- out[(nrow(out) - MAX_STORE + 1):nrow(out), , drop = FALSE]
    out
  }

  # =========================================================== Tab 1 =====
  observeEvent(input$t1_draw, {
    P <- pop()
    s <- draw_samples(1, P$p1, P$p2, P$n1, P$n2); s$jit <- runif(1)
    rv$t1 <- add_rows(rv$t1, s)
  })
  observeEvent(input$t1_reset, { rv$t1 <- rv$t1[nrow(rv$t1), , drop = FALSE] })

  t1_last <- reactive({ req(rv$t1); rv$t1[nrow(rv$t1), ] })

  output$vb_param <- renderUI(withMathJax(sprintf("\\(%s\\)", fmt(pop()$diff, 2))))
  output$vb_est <- renderText({
    s <- t1_last(); P <- pop()
    fmt(s$x1 / P$n1 - s$x2 / P$n2, 3)
  })
  output$vb_est_detail <- renderUI({
    s <- t1_last(); P <- pop()
    p(sprintf("%s − %s  (sample #%d)", fmt(s$x1 / P$n1, 3), fmt(s$x2 / P$n2, 3), nrow(rv$t1)))
  })

  output$t1_waffle <- renderPlot({
    s <- t1_last(); P <- pop()
    lab1 <- sprintf("%s: %d of %d successes, p̂₁ = %s", P$g1, s$x1, P$n1, fmt(s$x1 / P$n1))
    lab2 <- sprintf("%s: %d of %d successes, p̂₂ = %s", P$g2, s$x2, P$n2, fmt(s$x2 / P$n2))
    if (max(P$n1, P$n2) <= 600) {
      nc <- if (max(P$n1, P$n2) <= 100) 10 else if (max(P$n1, P$n2) <= 300) 20 else 30
      d <- rbind(waffle_df(s$x1, P$n1, lab1, nc), waffle_df(s$x2, P$n2, lab2, nc))
      d$group <- factor(d$group, levels = c(lab1, lab2))
      d$fill <- ifelse(!d$success, "Failure",
                       ifelse(d$group == lab1, "Success (group 1)", "Success (group 2)"))
      ggplot(d, aes(col, row, fill = fill)) +
        geom_tile(colour = "white", linewidth = 0.6) +
        facet_wrap(~group, ncol = 1, scales = "free_y") +
        scale_fill_manual(values = c("Failure" = col_fail, "Success (group 1)" = col_g1,
                                     "Success (group 2)" = col_g2), name = NULL) +
        coord_cartesian(expand = FALSE) +
        labs(subtitle = "Each square is one sampled individual (successes listed first)") +
        theme_void(base_size = 14) +
        theme(strip.text = element_text(face = "bold", hjust = 0, size = 12,
                                        margin = margin(4, 0, 4, 0)),
              legend.position = "top", plot.subtitle = element_text(colour = "#4B5563"))
    } else {
      d <- data.frame(group = factor(c(P$g1, P$g2), levels = c(P$g1, P$g2)),
                      phat = c(s$x1 / P$n1, s$x2 / P$n2), p = c(P$p1, P$p2))
      ggplot(d, aes(group, phat, fill = group)) +
        geom_col(width = 0.55) +
        geom_errorbar(aes(ymin = p, ymax = p), width = 0.75, linewidth = 1.1,
                      colour = col_param, linetype = "22") +
        geom_text(aes(label = paste0("p̂ = ", fmt(phat))), vjust = -0.6, size = 5) +
        scale_fill_manual(values = c(col_g1, col_g2), guide = "none") +
        scale_y_continuous(limits = c(0, 1.08), labels = scales::percent) +
        labs(x = NULL, y = "Proportion of successes",
             subtitle = "Bars = sample proportions; dashed lines = true proportions") +
        base_theme
    }
  })

  output$t1_text <- renderUI({
    s <- t1_last(); P <- pop()
    est <- s$x1 / P$n1 - s$x2 / P$n2
    withMathJax(div(class = "worked",
      HTML(sprintf("Estimate: \\(\\hat p_1 - \\hat p_2 = \\frac{%d}{%d} - \\frac{%d}{%d} = %s - %s = %s.\\)
                   &nbsp; This estimate is off from the parameter by \\(%s\\).",
                   s$x1, P$n1, s$x2, P$n2, fmt(s$x1 / P$n1), fmt(s$x2 / P$n2), fmt(est),
                   fmt(est - P$diff)))))
  })

  output$t1_strip <- renderPlot({
    req(rv$t1); P <- pop(); d <- rv$t1
    d$est <- d$x1 / P$n1 - d$x2 / P$n2
    d$latest <- seq_len(nrow(d)) == nrow(d)
    lim <- P$diff + c(-1, 1) * max(4.2 * P$se, abs(d$est - P$diff) + 0.01)
    ggplot(d, aes(est, jit)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = NA) +
      geom_vline(xintercept = P$diff, colour = col_param, linewidth = 1.1) +
      annotate("label", x = P$diff, y = 1.12, label = sprintf("parameter = %s", fmt(P$diff, 2)),
               size = 4, label.size = 0, fill = "white") +
      geom_point(data = d[!d$latest, ], colour = col_hit, alpha = 0.55, size = 3) +
      geom_point(data = d[d$latest, ], colour = col_miss, size = 5.5, shape = 18) +
      scale_y_continuous(limits = c(-0.05, 1.18), breaks = NULL) +
      scale_x_continuous(limits = lim) +
      labs(x = expression(hat(p)[1] - hat(p)[2]), y = NULL,
           subtitle = "Red diamond = newest estimate; blue dots = earlier ones") +
      base_theme + theme(panel.grid.major.y = element_blank())
  })

  output$t1_strip_text <- renderUI({
    req(rv$t1); P <- pop(); d <- rv$t1
    if (nrow(d) < 2) return(p(class = "text-muted small", "Draw more samples to watch the estimates pile up."))
    est <- d$x1 / P$n1 - d$x2 / P$n2
    p(class = "small",
      sprintf("%d estimates so far: smallest %s, largest %s, average %s. The parameter never moved.",
              length(est), fmt(min(est)), fmt(max(est)), fmt(mean(est))),
      tags$br(), "Want hundreds at once? Head to Tab 2.")
  })

  # =========================================================== Tab 2 =====
  t2_add <- function(k) {
    P <- pop()
    rv$t2 <- add_rows(rv$t2, draw_samples(k, P$p1, P$p2, P$n1, P$n2))
  }
  observeEvent(input$t2_1,    t2_add(1))
  observeEvent(input$t2_10,   t2_add(10))
  observeEvent(input$t2_100,  t2_add(100))
  observeEvent(input$t2_1000, t2_add(1000))
  observeEvent(input$t2_reset, rv$t2 <- NULL)

  output$t2_hist <- renderPlot({
    P <- pop()
    lim <- P$diff + c(-1, 1) * 4.5 * P$se
    base <- ggplot() +
      labs(x = expression("Estimate  " * hat(p)[1] - hat(p)[2]), y = "Density") +
      base_theme
    if (isTRUE(input$t2_band))
      base <- base + annotate("rect", xmin = P$diff - 2 * P$se, xmax = P$diff + 2 * P$se,
                              ymin = -Inf, ymax = Inf, fill = "#FDE68A", alpha = 0.45)
    if (is.null(rv$t2)) {
      base <- base +
        annotate("text", x = P$diff, y = 0.5, size = 5, colour = "#6B7280",
                 label = "Click “+1,000” to draw many samples") +
        scale_y_continuous(limits = c(0, 1), breaks = NULL)
    } else {
      d <- rv$t2; d$est <- d$x1 / P$n1 - d$x2 / P$n2
      lim <- range(c(lim, d$est))
      bw <- nice_binwidth(P$se, P$n1, P$n2)
      last <- d$est[nrow(d)]
      base <- base +
        geom_histogram(data = d, aes(est, after_stat(density)), binwidth = bw$width,
                       center = bw$center, fill = col_sim, colour = "white", linewidth = 0.2) +
        annotate("point", x = last, y = 0, shape = 17, size = 5, colour = col_miss) +
        annotate("text", x = last, y = 0, vjust = 2, size = 3.8, colour = col_miss,
                 label = "latest") +
        labs(title = sprintf("%s simulated estimates", format(nrow(d), big.mark = ",")))
    }
    if (isTRUE(input$t2_curve))
      base <- base + stat_function(fun = dnorm, args = list(mean = P$diff, sd = P$se),
                                   colour = col_param, linewidth = 1, n = 400, xlim = lim)
    if (isTRUE(input$t2_param))
      base <- base + geom_vline(xintercept = P$diff, colour = col_param,
                                linetype = "22", linewidth = 1)
    base + coord_cartesian(xlim = lim, clip = "off")
  })

  output$t2_table <- renderTable({
    P <- pop()
    if (is.null(rv$t2)) {
      sim_mean <- sim_sd <- within2 <- "—"; k <- 0
    } else {
      est <- rv$t2$x1 / P$n1 - rv$t2$x2 / P$n2; k <- length(est)
      sim_mean <- fmt(mean(est), 4)
      sim_sd   <- if (k > 1) fmt(sd(est), 4) else "—"
      within2  <- pct(mean(abs(est - P$diff) <= 2 * P$se))
    }
    data.frame(
      ` ` = c("Number of samples", "Center (mean)", "Spread (SD = SE)", "Within 2 SE of parameter"),
      Simulation = c(format(k, big.mark = ","), sim_mean, sim_sd, within2),
      Theory = c("∞", fmt(P$diff, 4), fmt(P$se, 4), pct(pnorm(2) - pnorm(-2))),
      check.names = FALSE)
  }, striped = TRUE, align = "lrr", width = "100%")

  output$t2_formula <- renderUI({
    P <- pop()
    withMathJax(div(class = "worked small",
      HTML(sprintf("Theory says the sampling distribution is approximately
        $$\\hat p_1 - \\hat p_2 \\sim N\\!\\left(p_1 - p_2,\\ \\sqrt{\\tfrac{p_1(1-p_1)}{n_1} + \\tfrac{p_2(1-p_2)}{n_2}}\\right)$$
        $$SE = \\sqrt{\\tfrac{%s(%s)}{%d} + \\tfrac{%s(%s)}{%d}} = %s$$",
        fmt(P$p1, 2), fmt(1 - P$p1, 2), P$n1, fmt(P$p2, 2), fmt(1 - P$p2, 2), P$n2, fmt(P$se, 4)))))
  })

  output$t2_conditions <- renderUI({
    P <- pop()
    v <- c(P$n1 * P$p1, P$n1 * (1 - P$p1), P$n2 * P$p2, P$n2 * (1 - P$p2))
    ok <- all(v >= 10)
    div(class = paste("callout small", if (!ok) "callout-warn"),
      tags$b(if (ok) "Normal approximation: conditions met" else "Normal approximation: conditions NOT met"),
      tags$br(),
      sprintf("Expected counts: %s, %s (group 1); %s, %s (group 2). Each should be ≥ 10.",
              fmt(v[1], 1), fmt(v[2], 1), fmt(v[3], 1), fmt(v[4], 1)),
      if (!ok) tagList(tags$br(), "Look for skew and gaps in the histogram compared with the curve."))
  })

  # =========================================================== Tab 3 =====
  t3_add <- function(k) {
    P <- pop()
    rv$t3 <- add_rows(rv$t3, draw_samples(k, P$p1, P$p2, P$n1, P$n2))
  }
  observeEvent(input$t3_1,   t3_add(1))
  observeEvent(input$t3_20,  t3_add(20))
  observeEvent(input$t3_100, t3_add(100))
  observeEvent(input$t3_reset, rv$t3 <- NULL)

  t3_ci <- reactive({
    req(rv$t3); P <- pop()
    ci <- ci_two_prop(rv$t3$x1, rv$t3$x2, P$n1, P$n2, input$conf, input$ci_method)
    ci$hit <- ci$lo <= P$diff & ci$hi >= P$diff
    ci$id <- seq_len(nrow(ci))
    ci
  })

  output$t3_n   <- renderText(if (is.null(rv$t3)) "0" else format(nrow(rv$t3), big.mark = ","))
  output$t3_cov <- renderText(if (is.null(rv$t3)) "—" else
                     sprintf("%s  (%d)", pct(mean(t3_ci()$hit)), sum(t3_ci()$hit)))
  output$t3_nom <- renderText(pct(input$conf, 0))

  output$t3_plot <- renderPlot({
    P <- pop()
    if (is.null(rv$t3)) {
      return(ggplot() + annotate("text", x = 0, y = 0, size = 5, colour = "#6B7280",
                                 label = "Draw some samples to build intervals") + theme_void())
    }
    ci <- t3_ci(); ci <- ci[ci$id > max(ci$id) - 100, ]
    ci$status <- factor(ifelse(ci$hit, "Captured the parameter", "Missed the parameter"),
                        levels = c("Captured the parameter", "Missed the parameter"))
    ggplot(ci, aes(y = id)) +
      geom_vline(xintercept = P$diff, colour = col_param, linewidth = 1) +
      geom_segment(aes(x = lo, xend = hi, yend = id, colour = status),
                   linewidth = if (nrow(ci) > 60) 0.9 else 1.4) +
      geom_point(aes(x = est, colour = status), size = if (nrow(ci) > 60) 1.2 else 2) +
      annotate("label", x = P$diff, y = max(ci$id) + 0.5, vjust = -0.3,
               label = sprintf("parameter = %s", fmt(P$diff, 2)), label.size = 0, size = 4) +
      scale_colour_manual(values = c("Captured the parameter" = col_hit,
                                     "Missed the parameter" = col_miss), name = NULL, drop = FALSE) +
      scale_y_continuous(expand = expansion(add = c(1, 4))) +
      labs(x = expression(p[1] - p[2]), y = "Sample number") +
      coord_cartesian(clip = "off") + base_theme
  })

  output$t3_worked <- renderUI({
    P <- pop()
    if (is.null(rv$t3)) return(p(class = "text-muted", "No intervals yet."))
    ci <- t3_ci(); k <- nrow(ci); s <- rv$t3[k, ]; r <- ci[k, ]
    ph1 <- s$x1 / P$n1; ph2 <- s$x2 / P$n2
    lvl <- pct(input$conf, 0)
    if (input$ci_method == "wald") {
      se_line <- sprintf("\\(\\widehat{SE} = \\sqrt{\\tfrac{%s(%s)}{%d} + \\tfrac{%s(%s)}{%d}} = %s\\)",
                         fmt(ph1), fmt(1 - ph1), P$n1, fmt(ph2), fmt(1 - ph2), P$n2, fmt(r$se, 4))
      centre_line <- sprintf("\\(\\hat p_1 - \\hat p_2 = %s - %s = %s\\)", fmt(ph1), fmt(ph2), fmt(r$est))
    } else {
      pt1 <- (s$x1 + 1) / (P$n1 + 2); pt2 <- (s$x2 + 1) / (P$n2 + 2)
      centre_line <- sprintf("Adjusted: \\(\\tilde p_1 = \\tfrac{%d+1}{%d+2} = %s,\\ \\tilde p_2 = \\tfrac{%d+1}{%d+2} = %s\\), centre \\(= %s\\)",
                             s$x1, P$n1, fmt(pt1), s$x2, P$n2, fmt(pt2), fmt(r$centre))
      se_line <- sprintf("\\(\\widehat{SE} = \\sqrt{\\tfrac{\\tilde p_1(1-\\tilde p_1)}{n_1+2} + \\tfrac{\\tilde p_2(1-\\tilde p_2)}{n_2+2}} = %s\\)",
                         fmt(r$se, 4))
    }
    withMathJax(div(class = "worked",
      div(class = "step", HTML(sprintf("<b>Sample #%d:</b> \\(x_1 = %d\\) of \\(%d\\), \\(x_2 = %d\\) of \\(%d\\)",
                                       k, s$x1, P$n1, s$x2, P$n2))),
      div(class = "step", HTML(paste("<b>1. Estimate:</b>", centre_line))),
      div(class = "step", HTML(paste("<b>2. Estimated SE:</b>", se_line))),
      div(class = "step", HTML(sprintf("<b>3. Critical value</b> for %s: \\(z^* = %s\\)", lvl, fmt(r$zstar)))),
      div(class = "step", HTML(sprintf("<b>4. Interval:</b> \\(%s \\pm %s \\times %s = (%s,\\ %s)\\)",
                                       fmt(r$centre), fmt(r$zstar), fmt(r$se, 4), fmt(r$lo), fmt(r$hi)))),
      div(class = paste("callout small mt-2", if (!r$hit) "callout-warn"),
        HTML(sprintf("This interval %s the parameter \\(p_1 - p_2 = %s\\).<br>
                      Its estimated SE (%s) vs. the true SE you can't normally see (%s).",
                     if (r$hit) "<b>captured</b>" else "<b>missed</b>", fmt(P$diff, 2),
                     fmt(r$se, 4), fmt(P$se, 4))))
    ))
  })

  # =========================================================== Tab 4 =====
  observeEvent(input$t4_draw, {
    P <- pop(); rv$t4 <- draw_samples(1, P$p1, P$p2, P$n1, P$n2)
  })

  t4_data <- reactive({
    if (input$t4_src == "own") {
      x1 <- input$o_x1; n1 <- input$o_n1; x2 <- input$o_x2; n2 <- input$o_n2
      validate(need(all(!is.na(c(x1, n1, x2, n2))), "Enter all four counts."),
               need(x1 >= 0 && x2 >= 0 && n1 >= 1 && n2 >= 1, "Counts must be non-negative and n ≥ 1."),
               need(x1 <= n1 && x2 <= n2, "Successes cannot exceed the sample size."))
      list(x1 = round(x1), x2 = round(x2), n1 = round(n1), n2 = round(n2), src = "own")
    } else {
      req(rv$t4); P <- pop()
      list(x1 = rv$t4$x1, x2 = rv$t4$x2, n1 = P$n1, n2 = P$n2, src = "pop")
    }
  })

  t4_test <- reactive({
    D <- t4_data()
    zt <- z_pooled(D$x1, D$x2, D$n1, D$n2)
    pv <- if (zt$se0 > 0) p_value(zt$z, input$alt) else 1
    c(D, zt, list(pv = pv, ph1 = D$x1 / D$n1, ph2 = D$x2 / D$n2,
                  se_unpooled = sqrt(D$x1 / D$n1 * (1 - D$x1 / D$n1) / D$n1 +
                                     D$x2 / D$n2 * (1 - D$x2 / D$n2) / D$n2)))
  })

  t4_perm <- reactive({
    req(isTRUE(input$t4_perm)); T <- t4_test()
    X <- T$x1 + T$x2; N <- T$n1 + T$n2
    x1s <- rhyper(2000, X, N - X, T$n1)
    ds <- x1s / T$n1 - (X - x1s) / T$n2
    tol <- 1e-10
    pv <- switch(input$alt,
                 two.sided = mean(abs(ds) >= abs(T$d) - tol),
                 greater   = mean(ds >= T$d - tol),
                 less      = mean(ds <= T$d + tol))
    list(ds = ds, pv = pv)
  })

  output$t4_worked <- renderUI({
    T <- t4_test(); a <- input$alpha
    reject <- T$pv <= a
    sym <- alt_symbol(input$alt)
    pv_txt <- if (input$alt == "two.sided") sprintf("2P(Z \\ge |%s|)", fmt(T$z, 2)) else
              if (input$alt == "greater") sprintf("P(Z \\ge %s)", fmt(T$z, 2)) else
              sprintf("P(Z \\le %s)", fmt(T$z, 2))
    withMathJax(div(class = "worked",
      div(class = "step", HTML(sprintf("\\(H_0: p_1 - p_2 = 0 \\quad H_A: p_1 - p_2 %s 0\\)", sym))),
      div(class = "step", HTML(sprintf("<b>Estimate:</b> \\(\\hat p_1 - \\hat p_2 = \\tfrac{%d}{%d} - \\tfrac{%d}{%d} = %s\\)",
                                       T$x1, T$n1, T$x2, T$n2, fmt(T$d)))),
      div(class = "step", HTML(sprintf("<b>Pooled:</b> \\(\\hat p = \\tfrac{%d + %d}{%d + %d} = %s\\)",
                                       T$x1, T$x2, T$n1, T$n2, fmt(T$pp)))),
      div(class = "step", HTML(sprintf("<b>Null SE:</b> \\(SE_0 = \\sqrt{%s(%s)\\left(\\tfrac{1}{%d} + \\tfrac{1}{%d}\\right)} = %s\\)",
                                       fmt(T$pp), fmt(1 - T$pp), T$n1, T$n2, fmt(T$se0, 4)))),
      div(class = "step", HTML(sprintf("<b>Test statistic:</b> \\(z = \\dfrac{%s - 0}{%s} = %s\\)",
                                       fmt(T$d), fmt(T$se0, 4), fmt(T$z, 2)))),
      div(class = "step", HTML(sprintf("<b>p-value:</b> \\(%s = %s\\)", pv_txt, fmtp(T$pv)))),
      if (isTRUE(input$t4_perm))
        div(class = "step", HTML(sprintf("<b>Simulated p-value</b> (re-shuffling): %s", fmtp(t4_perm()$pv)))),
      div(class = paste("callout small", if (reject) "callout-warn"),
        HTML(sprintf("<b>Decision at \\(\\alpha = %s\\):</b> %s<br>
                      For comparison, the unpooled SE used for a confidence interval would be %s.",
                     fmt(a, 2),
                     if (reject) "p-value ≤ α, so reject \\(H_0\\)." else "p-value > α, so do not reject \\(H_0\\).",
                     fmt(T$se_unpooled, 4))))
    ))
  })

  output$t4_plot <- renderPlot({
    T <- t4_test()
    validate(need(T$se0 > 0, "Every observation is a success (or every one a failure), so the null distribution has no spread and no test is possible."))
    span <- max(4.5 * T$se0, abs(T$d) * 1.15)
    showtruth <- T$src == "pop" && isTRUE(input$t4_truth)
    if (showtruth) span <- max(span, abs(pop()$diff) * 1.15)
    xs <- seq(-span, span, length.out = 600)
    dens <- data.frame(x = xs, y = dnorm(xs, 0, T$se0))
    shade <- switch(input$alt,
      two.sided = dens[abs(dens$x) >= abs(T$d), ],
      greater   = dens[dens$x >= T$d, ],
      less      = dens[dens$x <= T$d, ])
    if (input$alt == "two.sided") {
      shade$side <- ifelse(shade$x < 0, "L", "R")
    } else shade$side <- "A"
    ymax <- max(dens$y)
    if (isTRUE(input$t4_perm)) {
      pm <- t4_perm()
      # with equal n's a re-shuffled difference (2*x1 - X)/n moves in steps of 2/n
      bw <- nice_binwidth(T$se0, T$n1, T$n2, step = if (T$n1 == T$n2) 2 / T$n1 else NULL)
      ymax <- max(ymax, max(table(round((pm$ds - T$d) / bw$width))) / (length(pm$ds) * bw$width))
    }
    g <- ggplot(dens, aes(x, y))
    if (isTRUE(input$t4_perm)) {
      g <- g + geom_histogram(data = data.frame(x = pm$ds), aes(x, after_stat(density)),
                              inherit.aes = FALSE, binwidth = bw$width, center = T$d,  # observed d sits on the lattice
                              fill = col_sim, colour = "white", linewidth = 0.2, alpha = 0.9)
    }
    g <- g +
      geom_area(data = shade, aes(group = side), fill = col_miss, alpha = 0.45) +
      geom_line(colour = col_null, linewidth = 1.1) +
      geom_vline(xintercept = T$d, colour = col_param, linewidth = 1.2) +
      annotate("label", x = T$d, y = ymax * 1.08, label.size = 0, size = 4,
               label = sprintf("our estimate = %s", fmt(T$d))) +
      geom_vline(xintercept = 0, colour = col_null, linetype = "22")
    if (showtruth)
      g <- g + geom_vline(xintercept = pop()$diff, colour = col_g1, linewidth = 1.1, linetype = "41") +
        annotate("label", x = pop()$diff, y = ymax * 0.6, label.size = 0, size = 3.8,
                 colour = col_g1, label = sprintf("true p₁ − p₂ = %s\n(unknown in practice)", fmt(pop()$diff, 2)))
    zax <- pretty(c(-span, span) / T$se0)
    g + scale_x_continuous(name = expression("Difference  " * hat(p)[1] - hat(p)[2]),
                           sec.axis = sec_axis(~ . / T$se0, name = "z-score", breaks = zax)) +
      scale_y_continuous(name = "Density", limits = c(0, ymax * 1.15)) +
      labs(title = expression("Null distribution of " * hat(p)[1] - hat(p)[2] * "  (if " * H[0] * " is true)"),
           subtitle = sprintf("Normal curve centred at 0 with SE₀ = %s.  Red area = p-value = %s",
                              fmt(T$se0, 4), fmtp(T$pv))) +
      base_theme
  })

  output$t4_plot_caption <- renderUI({
    p(class = "small text-muted",
      "The grey curve shows where estimates would fall if the two population proportions were equal. ",
      "The black line is the estimate from our sample. ",
      if (isTRUE(input$t4_perm)) "The blue histogram comes from randomly re-assigning the observed successes to the two groups, which is another way to create a world where H₀ is true.")
  })

  # ----- Long run --------------------------------------------------------
  t5_add <- function(k) {
    P <- pop(); s <- draw_samples(k, P$p1, P$p2, P$n1, P$n2)
    zt <- z_pooled(s$x1, s$x2, P$n1, P$n2)
    rv$t5 <- add_rows(rv$t5, data.frame(z = zt$z, se0 = zt$se0))
  }
  observeEvent(input$t5_100,  t5_add(100))
  observeEvent(input$t5_1000, t5_add(1000))
  observeEvent(input$t5_reset, rv$t5 <- NULL)

  t5_reject <- reactive({
    req(rv$t5)
    pv <- ifelse(rv$t5$se0 > 0, p_value(rv$t5$z, input$alt), 1)
    pv <= input$alpha
  })

  output$t5_plot <- renderPlot({
    a <- input$alpha; P <- pop()
    crit <- switch(input$alt, two.sided = qnorm(1 - a / 2), qnorm(1 - a))
    zs <- if (is.null(rv$t5)) numeric(0) else rv$t5$z
    lim <- range(c(-4, 4, zs))
    reg <- switch(input$alt,
      two.sided = data.frame(xmin = c(-Inf, crit), xmax = c(-crit, Inf)),
      greater   = data.frame(xmin = crit, xmax = Inf),
      less      = data.frame(xmin = -Inf, xmax = -crit))
    g <- ggplot() +
      geom_rect(data = reg, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                fill = col_miss, alpha = 0.12)
    if (length(zs)) {
      g <- g + geom_histogram(data = data.frame(z = zs), aes(z, after_stat(density)),
                              binwidth = 0.25, center = 0, fill = col_sim,
                              colour = "white", linewidth = 0.2) +
        labs(title = sprintf("z-statistics from %s simulated studies", format(length(zs), big.mark = ",")))
    } else {
      g <- g + annotate("text", x = 0, y = 0.25, size = 5, colour = "#6B7280",
                        label = "Run some tests to see their z-statistics")
    }
    g + stat_function(fun = dnorm, colour = col_null, linewidth = 1, n = 400, xlim = lim) +
      geom_vline(xintercept = switch(input$alt, two.sided = c(-crit, crit), greater = crit, less = -crit),
                 colour = col_miss, linetype = "22") +
      labs(x = "z-statistic", y = "Density",
           subtitle = "Grey curve: distribution of z if H₀ is true.  Red zone: reject H₀.") +
      coord_cartesian(xlim = lim) + base_theme
  })

  output$t5_summary <- renderUI({
    P <- pop(); a <- input$alpha
    h0_true <- isTRUE(all.equal(P$p1, P$p2))
    pw <- approx_power(P$p1, P$p2, P$n1, P$n2, a, input$alt)
    k <- if (is.null(rv$t5)) 0 else nrow(rv$t5)
    rate <- if (k) mean(t5_reject()) else NA
    truth <- if (h0_true)
      div(class = "callout small", tags$b("H₀ is TRUE in this population"), tags$br(),
          sprintf("p₁ = p₂ = %s. Every rejection is a Type I error. The rejection rate should settle near α = %s.",
                  fmt(P$p1, 2), fmt(a, 2)))
    else
      div(class = "callout callout-warn small", tags$b("H₀ is FALSE in this population"), tags$br(),
          sprintf("p₁ − p₂ = %s. Every rejection is a correct decision; the rejection rate estimates the power. Every non-rejection is a Type II error.",
                  fmt(P$diff, 2)),
          if ((input$alt == "greater" && P$diff < 0) || (input$alt == "less" && P$diff > 0))
            tagList(tags$br(), tags$i("Note: the true difference points the opposite way from Hₐ, so this one-sided test will almost never reject.")))
    tagList(
      truth,
      value_box(title = "Tests run", value = format(k, big.mark = ","), theme = "light"),
      value_box(title = if (h0_true) "Rejection rate (Type I error rate)" else "Rejection rate (estimated power)",
                value = if (k) pct(rate) else "—",
                theme = if (h0_true) "secondary" else "primary",
                p(if (h0_true) sprintf("Target: α = %s", pct(a, 0))
                  else sprintf("Normal-theory power ≈ %s", pct(pw)))),
      if (!h0_true && k)
        p(class = "small", sprintf("Type II error rate ≈ %s (= 1 − power).", pct(1 - rate)))
    )
  })
}

shinyApp(ui, server)
