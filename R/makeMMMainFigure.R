library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(grImport)

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

make_MM_main_figure <- function(
    rw_single,
    Y,
    output_dir,
    filename = "MMMainFigure.pdf",
    width = 30,
    height = 20
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # ============================================================
  # Extract inference
  # ============================================================

  inf <- mySvdxBenchmark:::.benchmark_report_extract_svdx_inf(rw_single)

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
    x = "group",group = "group", 
    observed = FALSE
  ); #plots as points, we need to remove that, and go to geom col 
  pA1$layers = list()
  pA1 <- pA1 + geom_col( aes(fill=group) ) +
    ylab("Component Function 1") +
    labs(title = NULL, tag = "A") +
    gt +
    panel_tag_theme

  # ============================================================
  # B: heatmap , no inference here due to small sample size
  # ============================================================

  tmpfn <- function(direction){
    
    # --- positive 
    feature_names <- rownames(Y)
    idx <- mySvdx:::.select_features_from_loadings(fit_svdx, component = 1, 
                                                   direction = direction, 
                                                   n_features = 100)
    features <- feature_names[idx]
    features <- features[!grepl("^IG", features, ignore.case = TRUE)][1:25]
    pB1 <- plotVHeatmap(
      fit_svdx,
      data = rw_single$data,
      component = 1,
      order_by = "group", Y = Y, 
      features = features
    ) 
    
    library(grid)
    library(ggplotify)
    grid::grid.grabExpr(
      ComplexHeatmap::draw(
        pB1,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legends = TRUE
      )
    )
    
  }
  
  pB1_grob <- tmpfn("positive")
  pB1 <- ggplotify::as.ggplot(pB1_grob) + labs(tag="B") +panel_tag_theme
  pB2_grob <- tmpfn("negative")
  pB2 <- ggplotify::as.ggplot(pB2_grob) 
  pB <- wrap_plots( pB1, pB2 ) + 
    gt

  # ============================================================
  # C: representative feature profiles
  # ============================================================

  ftrs <- c("MLLT3", "XAF1", "MX1",  #interferon block & MLLT3 added
            "SLC2A3", "BNIP3",  "CD27"
            )
  
  pC <- mySvdx::plotFeatureProfiles(
    inf,
    features = ftrs,
    Y = rw_single$Y,
    data = rw_single$data,
    source = "passed",
    x = "group",
    inference_scope = "component",
    component = comp1_name,
    select_by = "q",
    fdr = 0.1,
    top_n = 5,
    color = "feature",
    facet = "feature",
    summarize = "mean",
    show_se = TRUE,
    scale = "zscore",
    line = TRUE,
    points = TRUE,
    fdata = rw_single$fdata
  );  
  pC$data$feature_label <- factor(
    pC$data$feature_label,
    levels = ftrs
  )
  pC <- pC  +
    facet_wrap(~feature_label, ncol = 3) +
    labs(tag = "C") +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      text = element_text(size = 25)
    ) +
    panel_tag_theme

  # ============================================================
  # D: staining
  # ============================================================

  img <- pdftools::pdf_render_page(
    "./data/experiments/stainingMM.pdf",
    page = 1,
    dpi = 300
  )
  img <- magick::image_read(img)
  img <- magick::image_trim(img)
  pD <- ggplotify::as.ggplot( img ) +
    labs(tag="D") + 
    panel_tag_theme

  # ============================================================
  # Comparator plots
  # ============================================================

  pltU12Benchmark <- function(method) {
    pF1 <- plot_real_world_U(
      rw_single,
      method = method,
      component = 1,
      x = "group",
      color = "group",
      connect = TRUE,
      line_group = "group"
    ) +
      labs(title = NULL) +
      gt

    if (ncol(rw_single$benchmark_results[[method]]$U) > 1) {
      pF2 <- plot_real_world_U(
        rw_single,
        method = method,
        component = 2,
        x = "group",
        color = "group",
        connect = TRUE,
        line_group = "group"
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
    point_delta = 1.25
  ) +
    labs(tag = "E") +
    panel_tag_theme

  pE2 <- increase_plot_geom_sizes(
    pE[[2]],
    line_delta = 0.5,
    point_delta = 1.25
  )

  pF1 <- increase_plot_geom_sizes(
    pF[[1]],
    line_delta = 0.5,
    point_delta = 1.25
  ) +
    labs(tag = "F") +
    panel_tag_theme

  pF2 <- increase_plot_geom_sizes(
    pF[[2]],
    line_delta = 0.5,
    point_delta = 1.25
  )

  # ============================================================
  # Assembly
  # ============================================================

  # Four individual plots in one layout: no artificial A/B gap.
  top_row <- wrap_plots(
    pA1, pB1, pB2,
    ncol = 3, 
    widths = c(.85, 1.15, 1.15)
  )

  # C is deliberately narrower than the two enrichment plots.
  middle_core <- wrap_plots(
    pC,
    pD,
    ncol = 2,
    widths = c(0.85, 2.30)
  )

  # Narrow and center the middle row to compensate for GO labels
  # and enrichment legends.
  middle_row <- wrap_plots(
    plot_spacer(),
    middle_core,
    plot_spacer(),
    ncol = 3,
    widths = c(0.04, 0.92, 0.04)
  )

  # Four individual plots in one layout: no artificial PCA/NMF gap.
  bottom_row <- wrap_plots(
    pE1, pE2, pF1, pF2,
    ncol = 4
  )

  fig_mm <- wrap_plots(
    top_row,
    middle_core,
    bottom_row,
    ncol = 1#, 
    #heights = c(.85, 1.15, 1.05)
  )


  output_file <- file.path(
    output_dir,
    filename
  )

  ggsave(
    filename = output_file,
    plot = fig_mm,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  invisible(fig_mm)
}
