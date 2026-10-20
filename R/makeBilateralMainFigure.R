library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)

# ============================================================
# Helpers
# ============================================================

incr_repel_labels <- function(p, size = 5) {
  for (i in seq_along(p$layers)) {
    geom_class <- class(p$layers[[i]]$geom)[1]
    if (geom_class %in% c("GeomTextRepel", "GeomLabelRepel")) {
      p$layers[[i]]$aes_params$size <- size
    }
  }
  p
}

clean_paper_plot <- function(
    p,
    .title = NULL,
    .textSize = 25,
    .repel_size = 5,
    aspectRatio = NULL
) {
  if (is.null(.title)) {
    p <- p + theme(plot.title = element_blank())
  } else {
    p <- p + ggtitle(.title)
  }

  p <- incr_repel_labels(p, size = .repel_size)
  p <- p + theme(text = element_text(size = .textSize))

  if (!is.null(aspectRatio)) {
    p <- p + theme(aspect.ratio = aspectRatio)
  }

  p
}


set_line_color <- function(p, color) {
  for (i in seq_along(p$layers)) {
    geom_class <- class(p$layers[[i]]$geom)[1]
    if (geom_class %in% c("GeomLine", "GeomPath", "GeomSmooth")) {
      p$layers[[i]]$aes_params$colour <- color
    }
  }
  p
}

increase_plot_geom_sizes <- function(
    p,
    line_delta = 0,
    point_delta = 0
) {
  for (i in seq_along(p$layers)) {
    geom_class <- class(p$layers[[i]]$geom)[1]

    if (geom_class %in% c("GeomLine", "GeomPath", "GeomSmooth") &&
        line_delta != 0) {
      current <- p$layers[[i]]$aes_params$linewidth
      if (is.null(current)) current <- p$layers[[i]]$geom$default_aes$linewidth
      if (is.null(current) || !is.numeric(current)) current <- 0.5
      p$layers[[i]]$aes_params$linewidth <- current + line_delta
    }

    if (geom_class == "GeomPoint" && point_delta != 0) {
      current <- p$layers[[i]]$aes_params$size
      if (is.null(current)) current <- p$layers[[i]]$geom$default_aes$size
      if (is.null(current) || !is.numeric(current)) current <- 1.5
      p$layers[[i]]$aes_params$size <- current + point_delta
    }
  }

  p
}

clean_enrichment_label <- function(x, width = 35) {
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  stringr::str_to_sentence(
    stringr::str_trunc(x, width = width, side = "right")
  )
}

panel_tag_theme <- theme(
  plot.tag = element_text(
    size = 28,
    face = "bold",
    family = "Helvetica"
  ),
  plot.tag.position = c(0, 1)
)

make_bilateral_main_figure <- function(
    rw_single,
    output_dir="/tmp/bilateral",
    filename = "BilateralMainFigure.pdf",
    width = 30,
    height = 20
) {
  
  if (dir.exists(output_dir)) {
    message("Output directory with the same name exists. Overwriting. ")
    unlink(output_dir, recursive = TRUE) # Deletes the folder and its contents
  }
  dir.create(output_dir, recursive = TRUE)
  

  # ============================================================
  # Extract inference
  # ============================================================

  inf <- mySvdxBenchmark:::.benchmark_report_extract_svdx_inf(rw_single)

  comp1_name <- mySvdx:::.report_get_inference_component_name(inf, 1)
  comp2_name <- mySvdx:::.report_get_inference_component_name(inf, 2)

  # ============================================================
  # A: reconstructed processes
  # ============================================================

  fit_svdx <- rw_single$raw_fits$svdx$fit

  gt <- theme(
    text = element_text(size = 25),
    legend.position = "none",
    aspect.ratio = 1
  ) + labs(title = NULL)

  pA1 <- plotU(
    fit_svdx,
    data = rw_single$data,
    component = 1,
    x = "Timepoint",
    group = "Group", 
    n_grid = 200,
    observed = FALSE
  ) +
    ylab("Component Function 1") +
    labs(title = NULL, tag = "A") +
    gt +
    panel_tag_theme

  pA2 <- plotU(
    fit_svdx,
    data = rw_single$data,
    component = 2,
    x = "Timepoint",
    group = "Group",
    n_grid = 200,
    observed = FALSE
  ) +
    ylab("Component Function 2") +
    labs(title = NULL) +
    gt

  # Use the same SVDX color for both reconstructed component functions.
  # SVDX_COMPONENT_COLOR <- "#FF64B0"
  # 
  # pA1 <- set_line_color(
  #   pA1,
  #   SVDX_COMPONENT_COLOR
  # )
  # 
  # pA2 <- set_line_color(
  #   pA2,
  #   SVDX_COMPONENT_COLOR
  # )

  pA1 <- increase_plot_geom_sizes(
    pA1,
    line_delta = 1.0
  )

  pA2 <- increase_plot_geom_sizes(
    pA2,
    line_delta = 1.0
  )

  # ============================================================
  # B: conditional-inference volcano plots
  # ============================================================

  pB1 <- clean_paper_plot(
    mySvdx::plot_component_volcano(
      inf,
      component = comp1_name,
      fdata = rw_single$fdata,
      label_top_n = 10
    )
  ) +
    labs(title = NULL, tag = "B") +
    gt +
    panel_tag_theme

  pB2 <- clean_paper_plot(
    mySvdx::plot_component_volcano(
      inf,
      component = comp2_name,
      fdata = rw_single$fdata,
      label_top_n = 10
    )
  ) +
    labs(title = NULL) +
    gt

  # ============================================================
  # C: representative feature profiles
  # ============================================================

  ftrs <- c( "Rspo3", "Six4", "Cd3d", "Itk", "Slamf1",
             "Tgtp1", "Podxl2", "Frem2", "Lrrn2", "Apol9b"
             )
  pC <- mySvdx::plotFeatureProfiles(
    inf,
    features = ftrs,
    Y = rw_single$Y,
    data = rw_single$data,
    source = "passed",
    x = "Timepoint",
    group = "Group",
    inference_scope = "component",
    component = c(1,1,1,1,1,
                  2,2,2,2,2
                  ),
    show_significance = F,
    select_by = "q",
    fdr = 0.1,
    top_n = 5,
    color = "Group",
    facet = "feature",
    summarize = "mean",
    show_se = TRUE,
    scale = "zscore",
    line = TRUE,
    points = TRUE,
    fdata = rw_single$fdata
  ) ; 
  
  feature_levels <- vapply(
    ftrs,
    function(ff) {
      unique(pC$data$feature_label[pC$data$feature == ff])[1]
    },
    character(1)
  )
  
  
  pC$data$feature_label <- factor(
    pC$data$feature_label,
    levels = feature_levels
    )

  pC <- pC +
    facet_wrap(~feature_label, nrow = 2, scales = "free_y") +
    labs(tag = "C") +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      text = element_text(size = 25)
    ) +
    panel_tag_theme

  # ============================================================
  # D: enrichment
  # ============================================================

  enr <- mySvdxBenchmark:::.benchmark_report_extract_svdx_inference_enrichment(
    rw_single
  )

  enr_obj <- mySvdxBenchmark:::.benchmark_report_enrichment_dotplot_object(
    enr,
    context = "bilateral inference enrichment",
    verbose = TRUE
  )

  enr_df <- as.data.frame(enr_obj$results)

  e1 <- enr_df |>
    dplyr::filter(
      scope == "component",
      component == "IRCF_1"
    )

  e2 <- enr_df |>
    dplyr::filter(
      scope == "component",
      component == "IRCF_2"
    )

  pD1 <- mySvdx::plot_enrichment_dotplot(
    structure(
      list(results = e1),
      class = c("svdxEnrichment", "list")
    ),
    top_n = 10
  )

  pD2 <- mySvdx::plot_enrichment_dotplot(
    structure(
      list(results = e2),
      class = c("svdxEnrichment", "list")
    ),
    top_n = 10
  )

  pD1 <- pD1 +
    scale_y_discrete(
      labels = function(x) clean_enrichment_label(x, width = 35)
    ) +
    labs(title = "Component 1", tag = "D") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 25, hjust = 0.5),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 25),
      text = element_text(size = 25, family = "Helvetica")
    ) +
    panel_tag_theme

  pD2 <- pD2 +
    scale_y_discrete(
      labels = function(x) clean_enrichment_label(x, width = 35)
    ) +
    labs(title = "Component 2") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 25, hjust = 0.5),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 25),
      text = element_text(size = 25, family = "Helvetica")
    )

  # ============================================================
  # Comparator plots
  # ============================================================

  pltU12Benchmark <- function(method) {
    pF1 <- plot_real_world_U(
      rw_single,
      method = method,
      component = 1,
      x = "Timepoint",
      color = "Group",
      connect = TRUE,
      line_group = "Group"
    ) +
      labs(title = NULL) +
      gt

    if (ncol(rw_single$benchmark_results[[method]]$U) > 1) {
      pF2 <- plot_real_world_U(
        rw_single,
        method = method,
        component = 2,
        x = "Timepoint",
        color = "Group",
        connect = TRUE,
        line_group = "Group"
      ) +
        labs(title = NULL) +
        gt

      pF1 + pF2
    } else {
      pF1
    }
  }

  methods <- names(rw_single$benchmark_results)[-c(1:2)]
  names(methods) <- methods

  uplots <- lapply(
    methods,
    function(x) pltU12Benchmark(x)
  )

  pca_idx <- which(tolower(names(uplots)) == "pca")
  nmf_idx <- which(tolower(names(uplots)) == "nmf")

  if (length(pca_idx) != 1) {
    stop(
      "Could not uniquely identify PCA in uplots. Available names: ",
      paste(names(uplots), collapse = ", ")
    )
  }

  if (length(nmf_idx) != 1) {
    stop(
      "Could not uniquely identify NMF in uplots. Available names: ",
      paste(names(uplots), collapse = ", ")
    )
  }

  pE <- uplots[[pca_idx]]
  pF <- uplots[[nmf_idx]]

  # Flatten these patchworks so the bottom row has four equal plots.
  pE1 <- increase_plot_geom_sizes(
    pE[[1]],
    line_delta = 0.5,
    point_delta = 0.5
  ) +
    labs(tag = "E") +
    panel_tag_theme

  pE2 <- increase_plot_geom_sizes(
    pE[[2]],
    line_delta = 0.5,
    point_delta = 0.5
  )

  pF1 <- increase_plot_geom_sizes(
    pF[[1]],
    line_delta = 0.5,
    point_delta = 0.5
  ) +
    labs(tag = "F") +
    panel_tag_theme

  pF2 <- increase_plot_geom_sizes(
    pF[[2]],
    line_delta = 0.5,
    point_delta = 0.5
  )

  # ============================================================
  # Assembly
  # ============================================================

  # Four individual plots in one layout: no artificial A/B gap.
  top_row <- wrap_plots(
    pA1, pA2, pB1, pB2,
    ncol = 4
  )

  # C is deliberately narrower than the two enrichment plots.
  middle_core <- wrap_plots(
    pC,
    pD1,
    #pD2, #no enrichment for second
    ncol = 2,
    widths = c(2.2, 1)
  )

  # Narrow and center the middle row to compensate for GO labels
  # and enrichment legends.
  middle_row <- wrap_plots(
    plot_spacer(),
    middle_core,
    plot_spacer(),
    ncol = 3,
    widths = c(0.02, 0.96, 0.02)
  )

  # Four individual plots in one layout: no artificial PCA/NMF gap.
  bottom_row <- wrap_plots(
    pE1, pE2, pF1, pF2,
    ncol = 4
  )

  fig_bilateral <- wrap_plots(
    top_row,
    middle_row,
    bottom_row,
    ncol = 1
  )


  output_file <- file.path(
    output_dir,
    filename
  )

  ggsave(
    filename = output_file,
    plot = fig_bilateral,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  invisible(fig_bilateral)
}

