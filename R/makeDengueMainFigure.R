library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggplotify)
library(grid)

# ============================================================
# Helpers
# ============================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

panel_tag_theme <- theme(
  plot.tag = element_text(
    size = 28,
    face = "bold",
    family = "Helvetica"
  ),
  plot.tag.position = c(0, 1)
)

dengue_theme <- theme(
  text = element_text(size = 25, family = "Helvetica"),
  legend.position = "bottom"
)

.decomp_colorbar <- function(p) {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)

  sc <- p$scales$get_scales("fill")
  if (is.null(sc) || !inherits(sc, "ScaleContinuous")) return(p)

  p + guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = grid::unit(2.8, "in"),
      barheight = grid::unit(0.20, "in")
    )
  )
}

.clean_term_plot <- function(
    p,
    text_size = 25,
    keep_legend = FALSE,
    aspect_ratio = 1
) {
  p <- p +
    labs(title = NULL) +
    theme(
      text = element_text(size = text_size, family = "Helvetica"),
      legend.position = if (keep_legend) "bottom" else "none"#,
      #aspect.ratio = aspect_ratio
    )

  if (keep_legend) .decomp_colorbar(p) else p
}

# plotUTerms() returns the fitted SSANOVA terms.  We identify the month,
# IOBW, and interaction plots by their names rather than relying on order.
.extract_dengue_terms <- function(
    fit,
    data,
    component,
    month_var,
    iobw_var,
    facet = NULL,
    fill = "typical",
    plotUTerms_args = list()
) {
  term_obj <- do.call(
    mySvdx::plotUTerms,
    c(
      list(
        fit = fit,
        data = data,
        component = component,
        facet = facet,
        fill = fill
      ),
      plotUTerms_args
    )
  )

  if (inherits(term_obj, "ggplot")) {
    stop(
      "plotUTerms() returned a single ggplot. Expected separate month, IOBW, ",
      "and interaction terms for the dengue model.",
      call. = FALSE
    )
  }

  term_plots <- Filter(function(x) inherits(x, "ggplot"), term_obj)

  if (length(term_plots) < 3L) {
    stop(
      "Expected at least three fitted term plots from plotUTerms(); got ",
      length(term_plots), ".",
      call. = FALSE
    )
  }

  nm <- names(term_plots)
  if (is.null(nm)) nm <- rep("", length(term_plots))

  # Exact names first.
  month_idx <- which(nm == month_var)
  iobw_idx  <- which(nm == iobw_var)

  # Fall back to simple string matching.
  if (length(month_idx) != 1L) {
    month_idx <- which(
      grepl(month_var, nm, fixed = TRUE) &
        !grepl(iobw_var, nm, fixed = TRUE)
    )
  }

  if (length(iobw_idx) != 1L) {
    iobw_idx <- which(
      grepl(iobw_var, nm, fixed = TRUE) &
        !grepl(month_var, nm, fixed = TRUE)
    )
  }

  int_idx <- which(
    grepl(month_var, nm, fixed = TRUE) &
      grepl(iobw_var, nm, fixed = TRUE)
  )

  # If the interaction name is not informative, use the remaining 2-D/fill plot.
  if (length(int_idx) != 1L) {
    candidate <- setdiff(
      seq_along(term_plots),
      unique(c(month_idx, iobw_idx))
    )

    has_cont_fill <- vapply(
      term_plots[candidate],
      function(p) {
        sc <- p$scales$get_scales("fill")
        !is.null(sc) && inherits(sc, "ScaleContinuous")
      },
      logical(1)
    )

    int_idx <- candidate[has_cont_fill]
  }

  if (length(month_idx) != 1L ||
      length(iobw_idx) != 1L ||
      length(int_idx) != 1L) {
    stop(
      "Could not uniquely identify month, IOBW, and interaction terms.\n",
      "Available plotUTerms() names: ",
      paste(ifelse(nzchar(nm), nm, "<unnamed>"), collapse = ", "),
      call. = FALSE
    )
  }

  list(
    month = term_plots[[month_idx]],
    iobw = term_plots[[iobw_idx]],
    interaction = term_plots[[int_idx]]
  )
}



# If example_countries = NULL, reproduce the original behavior:
# choose one illustrative country for each requested component from the
# largest absolute SSANOVA-SVD loading, forcing the two choices to differ.
.auto_example_countries <- function(rw_single, components = c(1, 2)) {
  svdx_result <- rw_single$benchmark_results$svdx

  if (is.null(svdx_result) || is.null(svdx_result$V)) {
    stop(
      "Could not find rw_single$benchmark_results$svdx$V for automatic ",
      "country selection. Supply example_countries explicitly.",
      call. = FALSE
    )
  }

  V <- as.matrix(svdx_result$V)

  if (ncol(V) < max(components)) {
    stop("Requested components are not available in stored SSANOVA-SVD V.")
  }

  feature_names <- rownames(V) %||% rownames(rw_single$Y)
  if (is.null(feature_names)) {
    stop("Could not determine country/feature names.")
  }

  chosen <- character(0)

  for (cc in components) {
    ord <- order(abs(V[, cc]), decreasing = TRUE)
    cand <- feature_names[ord]
    cand <- cand[!cand %in% chosen]

    if (!length(cand)) {
      stop("Could not choose distinct illustrative countries.")
    }

    chosen <- c(chosen, cand[[1L]])
  }

  chosen
}


# Build a long data frame for user-supplied representative countries.
# Dengue incidence is standardized within country for display only.
.country_long_data <- function(
    rw_single,
    countries,
    date_var = "date",
    iobw_var = "iobw_lag_1_3"
) {
  Y <- as.matrix(rw_single$Y)
  data <- as.data.frame(rw_single$data)

  feature_names <- rownames(Y)
  if (is.null(feature_names)) {
    stop("Could not determine country names from rownames(rw_single$Y).")
  }

  missing_countries <- setdiff(countries, feature_names)
  if (length(missing_countries)) {
    stop(
      "Country/countries not found in rownames(rw_single$Y): ",
      paste(missing_countries, collapse = ", ")
    )
  }

  missing_vars <- setdiff(c(date_var, iobw_var), names(data))
  if (length(missing_vars)) {
    stop(
      "Variable(s) not found in rw_single$data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  out <- lapply(countries, function(country) {
    y <- as.numeric(Y[country, ])
    yz <- as.numeric(scale(y))

    data.frame(
      country = country,
      date = data[[date_var]],
      iobw = data[[iobw_var]],
      value = yz,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out) |>
    mutate(country = factor(country, levels = countries))
}


# C: Compare representative-country dengue trajectories over calendar time.
# Default is a single overlaid plot; use time_layout = "facet" if the two
# trajectories are clearer as separate panels.
.country_time_plot <- function(
    d,
    time_label = "Time",
    time_layout = c("overlay", "facet"),
    text_size = 25
) {
  time_layout <- match.arg(time_layout)

  p <- ggplot(
    d,
    aes(
      x = date,
      y = value,
      group = country,
      color = country
    )
  ) +
    geom_line(linewidth = 1.0, alpha = 0.9) +
    geom_point(size = 2.4, alpha = 0.9) +
    labs(
      x = time_label,
      y = "Standardized dengue incidence",
      color = NULL,
      title = NULL
    ) + theme_bw() +
    theme(
      text = element_text(size = text_size, family = "Helvetica"),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      #legend.background = element_rect(fill = scales::alpha("white", 0.70), colour = NA),
      legend.key = element_blank(), 
      axis.title.y = element_text(vjust=-2)
    ) 

  if (identical(time_layout, "facet")) {
    p <- p +
      facet_wrap(~country, nrow = 1) +
      guides(color = "none")  #+ 
      #theme(aspect.ratio = 1)
  }

  p
}


# D: Country-specific dengue-IOBW relationships.  These are deliberately
# shown as separate facets with LOESS smooths and standard-error ribbons.
.country_iobw_plot <- function(
    d,
    iobw_label = "IOBW",
    text_size = 25
) {
  ggplot(d, aes(x = iobw, y = value)) +
    geom_point(size = 2.6, alpha = 0.8) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      linewidth = 1.1
    ) +
    facet_wrap(~country, nrow = 1) +
    labs(
      x = iobw_label,
      y = "Standardized dengue incidence",
      title = NULL
    ) + theme_bw() +
    theme(
      text = element_text(size = text_size, family = "Helvetica"),
      legend.position = "none"
    ) #+ theme(aspect.ratio=1)
}


# Compact expression heatmap for countries/features selected by the largest
# absolute loading in one component.
.loading_heatmap_grob <- function(
    fit,
    rw_single,
    component,
    top_n = 10L,
    order_by = NULL,
    annotate = NULL
) {
  ht <- mySvdx::plotVHeatmap(
    object = fit,
    Y = as.matrix(rw_single$Y),
    data = as.data.frame(rw_single$data),
    zdata = if (!is.null(rw_single$zdata)) as.data.frame(rw_single$zdata) else NULL,
    component = component,
    select_by = "loading",
    direction = "absolute",
    n_features = top_n,
    order_by = order_by,
    annotate = annotate,
    fdata = rw_single$fdata
  )
  ht@row_names_param$gp <- grid::gpar(fontsize = 18)
  ht@column_names_param$gp <- grid::gpar(fontsize = 14)
  ht@heatmap_param$heatmap_legend_param$labels_gp <-
    grid::gpar(fontsize = 15)
  ht@heatmap_param$heatmap_legend_param$title_gp <-
    grid::gpar(fontsize = 15)
  ht@heatmap_param$heatmap_legend_param$direction <- "horizontal"
  
  ht@heatmap_param$show_heatmap_legend <- FALSE
  
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legends = TRUE
    )
  )
}


# ============================================================
# Main figure
# ============================================================

make_Dengue_main_figure <- function(
    rw_single,
    output_dir,
    filename = "DengueMainFigure.pdf",
    components = c(1, 2),
    month_var = "month",
    iobw_var = "iobw_lag_1_3",
    month_label = "Month",
    iobw_label = "IOBW",
    date_var = "date",
    time_label = "Time",
    example_countries = NULL,
    time_layout = c("overlay", "facet"),
    heatmap_top_n = 10L,
    heatmap_order_by = NULL,
    heatmap_annotate = NULL,
    include_heatmap = TRUE,
    plotUTerms_args = list(),
    width = 20,
    height = 25,
    text_size = 25
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (length(components) != 2L) {
    stop("This main-figure maker expects exactly two SSANOVA-SVD components.")
  }

  fit <- rw_single$raw_fits$svdx$fit
  if (is.null(fit)) {
    stop("Could not find rw_single$raw_fits$svdx$fit.")
  }

  data <- as.data.frame(rw_single$data)

  # ============================================================
  # A/B: SSANOVA decomposition of Components 1 and 2
  # ============================================================

  terms1 <- .extract_dengue_terms(
    fit = fit,
    data = data,
    component = components[[1L]],
    month_var = month_var,
    iobw_var = iobw_var,
    plotUTerms_args = plotUTerms_args
  )

  terms2 <- .extract_dengue_terms(
    fit = fit,
    data = data,
    component = components[[2L]],
    month_var = month_var,
    iobw_var = iobw_var,
    plotUTerms_args = plotUTerms_args
  )

  pA_month <- .clean_term_plot(
    terms1$month, text_size,
    keep_legend = FALSE#,
    #aspect_ratio = 1
  ) +
    labs(x = month_label, tag = "A") +
    panel_tag_theme

  pA_iobw <- .clean_term_plot(
    terms1$iobw, text_size,
    keep_legend = FALSE#,
    #aspect_ratio = 1
  ) +
    labs(x = iobw_label)

  pA_int <- .clean_term_plot(
    terms1$interaction, text_size,
    keep_legend = TRUE#,
    #aspect_ratio = 1
  ) +
    labs(x = month_label, y = iobw_label)

  pA <- wrap_plots(
    pA_month, pA_iobw, pA_int,
    nrow = 1,
    widths = c(1, 1, 1.20)
  ) +
    plot_annotation(
      title = paste0("Component ", components[[1L]])
    ) &
    theme(
      plot.title = element_text(
        size = text_size,
        face = "bold",
        hjust = 0.5,
        family = "Helvetica"
      )
    )


  pB_month <- .clean_term_plot(
    terms2$month, text_size,
    keep_legend = FALSE#,
    #aspect_ratio = 1
  ) +
    labs(x = month_label, tag = "B") +
    panel_tag_theme

  pB_iobw <- .clean_term_plot(
    terms2$iobw, text_size,
    keep_legend = FALSE#,
    #aspect_ratio = 1
  ) +
    labs(x = iobw_label)

  pB_int <- .clean_term_plot(
    terms2$interaction, text_size,
    keep_legend = TRUE#,
    #aspect_ratio = 1
  ) +
    labs(x = month_label, y = iobw_label)

  pB <- wrap_plots(
    pB_month, pB_iobw, pB_int,
    nrow = 1,
    widths = c(1, 1, 1.20)
  ) +
    plot_annotation(
      title = paste0("Component ", components[[2L]])
    ) &
    theme(
      plot.title = element_text(
        size = text_size,
        face = "bold",
        hjust = 0.5,
        family = "Helvetica"
      )
    )


  # ============================================================
  # C/D: representative countries supplied by the user
  # ============================================================

  time_layout <- match.arg(time_layout)

  if (is.null(example_countries)) {
    example_countries <- .auto_example_countries(
      rw_single = rw_single,
      components = components
    )
  }

  if (length(example_countries) != 2L) {
    stop("example_countries must be NULL or contain exactly two country names.")
  }

  country_data <- .country_long_data(
    rw_single = rw_single,
    countries = example_countries,
    date_var = date_var,
    iobw_var = iobw_var
  )

  pC <- .country_time_plot(
    d = country_data,
    time_label = time_label,
    time_layout = time_layout,
    text_size = text_size
  ) +
    labs(tag = "C") +
    panel_tag_theme

  pD <- .country_iobw_plot(
    d = country_data,
    iobw_label = iobw_label,
    text_size = text_size
  ) +
    labs(tag = "D") +
    panel_tag_theme

  country_row <- wrap_plots(
    pC, pD,
    nrow = 1,
    widths = c(1.10, 1)
  )


  # ============================================================
  # E: compact loading-based heatmaps
  # ============================================================

  if (isTRUE(include_heatmap)) {
    h1 <- .loading_heatmap_grob(
      fit = fit,
      rw_single = rw_single,
      component = components[[1L]],
      top_n = heatmap_top_n,
      order_by = heatmap_order_by,
      annotate = heatmap_annotate
    )

    h2 <- .loading_heatmap_grob(
      fit = fit,
      rw_single = rw_single,
      component = components[[2L]],
      top_n = heatmap_top_n,
      order_by = heatmap_order_by,
      annotate = heatmap_annotate
    )

    pE1 <- ggplotify::as.ggplot(h1) +
      labs(tag = "E") +
      panel_tag_theme

    pE2 <- ggplotify::as.ggplot(h2)

    pE <- wrap_plots(
      pE1, pE2, plot_spacer(),
      nrow = 1,
      widths = c(1, 1, .08)
    ) +
      plot_annotation(
        title = "Countries with the largest absolute loadings"
      ) &
      theme(
        plot.title = element_text(
          size = text_size,
          face = "bold",
          hjust = 0.5,
          family = "Helvetica"
        )
      )
  }


  # ============================================================
  # Assembly
  # ============================================================

  rows <- list(pA, pB, country_row) 
  heights <- c(1.0, 1.0, 0.85)

  if (isTRUE(include_heatmap)) {
    rows <- c(rows, list(pE))
    heights <- c(heights, 1.10)
  }

  fig_dengue <- wrap_plots(
    rows,
    ncol = 1,
    heights = heights
  ) 

  output_file <- file.path(output_dir, filename)

  ggsave(
    filename = output_file,
    plot = fig_dengue,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  message(
    "Representative countries: ",
    paste(example_countries, collapse = ", "),
    "; time layout = ", time_layout
  )

  invisible(fig_dengue)
}
