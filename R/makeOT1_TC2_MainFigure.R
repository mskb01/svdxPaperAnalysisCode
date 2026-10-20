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

make_OT1_TC2_main_figure <- function(
    rw_single,
    Y,
    xdata,
    fdata,
    output_dir="/tmp/Malaria",
    filename = "MalariaMainFigure.pdf",
    width = 30,
    height = 20
) {
  
  
  if (dir.exists(output_dir)) {
    message('Overwriting directory ')
    unlink(output_dir, recursive=T)
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
    x = "day",
    group = c("treatment"),
    #color=c("treatment", "tissue"), 
    #line_group = c("treatment", "tissue"),
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
    x = "day",
    group = c("treatment"),
    n_grid = 200,
    observed = FALSE
  ) +
    ylab("Component Function 2") +
    labs(title = NULL) +
    gt +
    theme(legend.position = "inside", 
          legend.position.inside = c(0.32, 0.98),
          legend.justification = c(1, 1)
    )
  
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
  # C: representative feature profiles
  # ============================================================
  
  ftrs <- c( "Arg1", "Msr1","Lyz2", "C1qc", #myeloid assoc with comp 1
             "Prdm1", "Cx3cr1", "Entpd1", "Zfp683"  #T cells differentiation assoc with  comp 2
  )
  
  pC <- mySvdx::plotFeatureProfiles(
    inf,
    features = ftrs,
    Y = rw_single$Y,
    data = rw_single$data,
    source = "passed",
    x = "day",
    group = "treatment",
    inference_scope = "component",
    show_significance = F,
    component = comp1_name,
    select_by = "q",
    fdr = 0.1,
    top_n = 5,
    color = "treatment",
    facet = "feature",
    summarize = "mean",
    show_se = TRUE,
    scale = "zscore",
    line = TRUE,
    points = TRUE,
    fdata = rw_single$fdata
  ) ;
  pC$data$feature_label <- factor(
    pC$data$feature_label,
    levels = ftrs
  )
  pC <- pC +
    facet_wrap(~feature_label, nrow = 2, scales = "free_y") +
    labs(tag = "C") +
    theme(
      aspect.ratio = 1,
      legend.position = "none", #same colors as pA2, so no need to replot
      text = element_text(size = 25)
    ) +
    panel_tag_theme
  
  # pC <- pC + theme(
  #   strip.text = element_text(size = 15),
  #   axis.text.x = element_text(size = 15),
  #   axis.text.y = element_text(size = 15),
  #   #axis.title.x = element_text(size = 11),
  #   #axis.title.y = element_text(size = 11)
  # )
  
  # ============================================================
  # Comparator plots
  # ============================================================
  
  pltU12Benchmark <- function(method) {
    pF1 <- plot_real_world_U(
      rw_single,
      method = method,
      component = 1,
      x = "day",
      color = "treatment",
      connect = TRUE,
      line_group = "treatment"
    ) +
      labs(title = NULL) +
      gt
    
    if (ncol(rw_single$benchmark_results[[method]]$U) > 1) {
      pF2 <- plot_real_world_U(
        rw_single,
        method = method,
        component = 2,
        x = "day",
        color = "treatment",
        connect = TRUE,
        line_group = "treatment"
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
  #nmf_idx <- pca_idx
  
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
    labs(tag = "D") +
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
    labs(tag = "E") +
    panel_tag_theme
  
  pF2 <- increase_plot_geom_sizes(
    pF[[2]],
    line_delta = 0.5,
    point_delta = 1.25
  )
  
  competitor_theme <- theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      vjust = 1#,
      #size = 8
    ),
    #axis.title.x = element_text(size = 10),
    #axis.title.y = element_text(size = 10)
  )
  
  pE1 <- pE1 + competitor_theme
  pE2 <- pE2 + competitor_theme
  pF1 <- pF1 + competitor_theme
  pF2 <- pF2 + competitor_theme
  
  
  #colors for other
  pb <- ggplot_build(pE1)
  d <- pb$data[[1]]
  group_cols <- unique(
    data.frame(
      treatment = pE1$data$treatment,
      colour = d$colour
    )
  )
  group_cols <- setNames(group_cols$colour, group_cols$treatment)
  
  
  
  
  # ============================================================
  # D: heatmaps
  # ============================================================
  
  tmpfn <- function(direction, component = 1, numberOfFtrs=20) {
    
    feature_names <- rownames(Y)
    feature_names_symbols <- fdata$symbol
    names(feature_names_symbols) <- feature_names
    
    idx <- mySvdx:::.select_features_from_inference(
      inf,
      component = component,
      select_by = "stat",
      direction = direction,
      n_features = numberOfFtrs
    )
    
    features <- feature_names_symbols[idx]
    # features <- features[
    #   !grepl("^PBANKA", features, ignore.case = TRUE)
    # ][1:numberOfFtrs]
    
    pB1 <- plotVHeatmap(
      inf,
      data = xdata,
      component = component,
      inference_scope = "component",
      fdr = 0.1,
      order_by = c("treatment", "day"),
      Y = Y,
      features = names(features)
    )
    
    #replace displayed feature IDs with symbols
    old <- pB1@row_names_param$labels
    new <- unname(features[old])
    new[is.na(new)] <- old[is.na(new)]
    
    pB1@row_names_param$labels <- new
    pB1@row_names_param$anno@var_env$value <- new
    
    # suppress ALL legends inside individual heatmaps
    pB1@top_annotation@anno_list[[1]]@show_legend <- FALSE
    pB1@heatmap_param$show_heatmap_legend <- FALSE
    
    pB1@row_names_param$gp <- grid::gpar(fontsize = 20)
    
    pB1@top_annotation@anno_list[[1]]@name_param$gp <-
      grid::gpar(fontsize = 20)
    pB1@top_annotation@anno_list[[2]]@name_param$gp <-
      grid::gpar(fontsize = 20)
    
    #pB1@top_annotation@anno_size <- grid::unit(5, "mm")
    #pB1@top_annotation@height    <- grid::unit(5, "mm")
    n_ann <- length(pB1@top_annotation@anno_list)
    pB1@top_annotation@anno_size <-
      grid::unit(rep(5, n_ann), "mm")
    
    ann <- pB1@top_annotation@anno_list[[1]]
    hm_levs <- ann@color_mapping@levels
    ann@color_mapping <- ComplexHeatmap::ColorMapping(
      name   = ann@color_mapping@name,
      colors = group_cols[hm_levs]
    )
    pB1@top_annotation@anno_list[[1]] <- ann
    
    # pB1@top_annotation@anno_list[[1]]@color_mapping@colors <-
    #   group_cols[
    #     names(pB1@top_annotation@anno_list[[1]]@color_mapping@colors)
    #   ]
    
    # grab as grob
    grid::grid.grabExpr(
      ComplexHeatmap::draw(
        pB1,
        show_heatmap_legend = FALSE,
        show_annotation_legend = FALSE
      )
    )
  }
  
  
  #NOTE WE are switching off legends for z scores for space, which in principle will vary
  #.  across plots. So the interpretation should be : row wise z scores are plotted. 
  #  High is red, low is blue. 
  # CF1
  pB11 <- ggplotify::as.ggplot(
    tmpfn("positive", component = 1) 
  ) +
    labs(tag = "B") + 
    panel_tag_theme
  
  pB12 <- ggplotify::as.ggplot(
    tmpfn("negative", component = 1)
  )
  
  # CF7: ookinete / midgut invasion
  pB71 <- ggplotify::as.ggplot(
    tmpfn("positive", component = 2)
  )
  
  pB72 <- ggplotify::as.ggplot(
    tmpfn("negative", component = 2)
  )
  
  #for group legend
  tmp <- plotVHeatmap(
    fit_svdx,
    data = xdata,
    component = 1,
    order_by = c( "treatment", "day"),
    Y = Y,
    features = rownames(Y)[1:2]
  )
  group_cols <- tmp@top_annotation@anno_list[[1]]@color_mapping@colors
  lgd_group <- ComplexHeatmap::Legend(
    title = "treatment",
    labels = names(group_cols),
    legend_gp = grid::gpar(fill = group_cols)
  )
  group_legend_grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(lgd_group)
  )
  p_group_legend <- ggplotify::as.ggplot(group_legend_grob)
  

  lefty <- pB11/pB71
  righty <- pB12/(pB72)
  
  pB_heatmaps <- wrap_plots(
    lefty, righty,
    #pB71, pB72,
    #pB91, pB92,
    nrow = 1, widths = c(1,1),
    byrow = FALSE
  )
  pB <- pB_heatmaps #| p_group_legend +
  #plot_layout(widths = c(12, 0.5))
  
  pB <- pB + gt + theme(aspect.ratio=1)
  
  
  
  
  # ============================================================
  # Assembly
  # ============================================================
  
  left_top <- (pA1 | pA2) / pC
  
  top_row <- wrap_plots(
    
    left_top,
    
    pB,
    
    ncol = 2,
    
    widths = c(2.1, 2.9)
    
  )
  
  bottom_row <- wrap_plots(
    
    pE1, pE2,
    
    pF1, pF2,
    
    #p_group_legend,
    
    ncol = 4,
    
    widths = c(1, 1, 1, 1)
    
  )
  
  fig_ottc2 <- wrap_plots(
    
    top_row,
    
    bottom_row,
    
    ncol = 1,
    
    heights = c(3.7, 1.5)
    
  ) + 
    plot_annotation(
      theme = theme(
        plot.margin = margin(t = 10, r = 10, b = 7, l = 5)
      )
    )
  
  # fig_ottc2 <- wrap_plots(
  #   top_row,
  #   #middle_row,
  #   bottom_row,
  #   ncol = 1,
  #   # heights = c(4.2, #2, 
  #   #             1)
  #   heights = c(3.5, 
  #               1.7)
  # ) + 
  #   plot_annotation(
  #     theme = theme(
  #       plot.margin = margin(t = 10, r = 7, b = 7, l = 5)
  #     )
  #   )
  
  
  output_file <- file.path(
    output_dir,
    filename
  )
  
  ggsave(
    filename = output_file,
    plot = fig_ottc2,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}
  
  