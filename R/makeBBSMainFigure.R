library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(rnaturalearth)

# ============================================================
# BBS main figure
#
# A-B: Component 1 — continental east/west process + inference
# C-D: Component 3 — localized South Texas process + inference
# E-G: Component 5 — spatially discontinuous process,
#      selected SSANOVA decomposition, + inference
#
# Comparator methods remain in the BBS supplement.
# ============================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

.get_ns <- function(pkg, name) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.", call. = FALSE)
  }
  getFromNamespace(name, pkg)
}

panel_tag_theme <- theme(
  plot.tag = element_text(
    size = 28,
    face = "bold",
    family = "Helvetica"
  ),
  plot.tag.position = c(0, 1)
)

paper_theme <- theme(
  text = element_text(size = 25, family = "Helvetica"),
  legend.text = element_text(size = 20, family = "Helvetica"),
  legend.title = element_text(size = 20, family = "Helvetica")
)

# ------------------------------------------------------------
# BBS map helpers
# ------------------------------------------------------------

set_ggrepel_text_size <- function(p, size = 5) {
  for (i in seq_along(p$layers)) {
    geom_class <- class(p$layers[[i]]$geom)[1]
    
    if (geom_class %in% c("GeomTextRepel", "GeomLabelRepel")) {
      p$layers[[i]]$aes_params$size <- size
    }
  }
  
  p
}


.bbs_prepare_map_data <- function(
    scores,
    data,
    grid_id_col,
    lat_col,
    lon_col,
    year_col
) {
  if (length(scores) != nrow(data)) {
    stop(
      "Component score length (", length(scores),
      ") does not match nrow(rw_single$data) (", nrow(data), ").",
      call. = FALSE
    )
  }

  df <- as.data.frame(data)
  df$.bbs_component_score <- as.numeric(scores)

  keep <- is.finite(df[[lat_col]]) &
    is.finite(df[[lon_col]]) &
    is.finite(df$.bbs_component_score) &
    !is.na(df[[grid_id_col]]) &
    !is.na(df[[year_col]])

  df <- df[keep, , drop = FALSE]

  df |>
    dplyr::group_by(
      .data[[grid_id_col]],
      .data[[lat_col]],
      .data[[lon_col]]
    ) |>
    dplyr::summarize(
      component_score = mean(.bbs_component_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    as.data.frame()
}

.bbs_svdx_U <- function(fit) {
  U <- tryCatch(fit$fitAtOptCV$umat, error = function(e) NULL)
  if (is.null(U)) U <- fit$U %||% NULL
  if (is.null(U)) U <- fit$u %||% NULL

  if (is.null(U)) {
    stop(
      "Could not find fitted SSANOVA-SVD component scores in ",
      "fit$fitAtOptCV$umat, fit$U, or fit$u.",
      call. = FALSE
    )
  }

  as.matrix(U)
}

.bbs_component_map <- function(
    scores,
    data,
    component,
    north_america,
    grid_id_col = "grid_id",
    lat_col = "lat_center",
    lon_col = "lon_center",
    year_col = "year",
    point_size = 1.4,
    map_xlim = c(-170, -50),
    map_ylim = c(20, 75),
    map_lon_breaks = c(-160, -120, -80, -60),
    tag = NULL
) {
  df <- .bbs_prepare_map_data(
    scores = scores,
    data = data,
    grid_id_col = grid_id_col,
    lat_col = lat_col,
    lon_col = lon_col,
    year_col = year_col
  )

  lim <- max(abs(df$component_score), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1

  ggplot() +
    geom_sf(
      data = north_america,
      fill = "grey95",
      color = "grey70",
      linewidth = 0.2
    ) +
    geom_point(
      data = df,
      aes(
        x = .data[[lon_col]],
        y = .data[[lat_col]],
        color = component_score
      ),
      size = point_size,
      alpha = 0.85
    ) +
    coord_sf(
      xlim = map_xlim,
      ylim = map_ylim,
      expand = FALSE
    ) +
    scale_x_continuous(breaks = map_lon_breaks) +
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-lim, lim),
      name = "Component score"
    ) +
    guides(
      color = guide_colorbar(
        direction = "horizontal",
        title.position = "left",
        barwidth = grid::unit(2.8, "in"),
        barheight = grid::unit(0.22, "in")
      )
    ) +
    labs(
      title = paste0("Component ", component),
      x = NULL,
      y = NULL,
      tag = tag
    ) +
    theme_minimal(base_size = 25) +
    theme(
      text = element_text(size = 25, family = "Helvetica"),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "pt"),
      legend.margin = margin(t = -60, r = 0, b = 0, l = 0, unit = "pt"),
      legend.title = element_text(size = 20),
      plot.title = element_text(size = 25, hjust = 0.5)#,
      #plot.margin = margin(t = 10, r = 10, b = 12, l = 10, unit = "pt")
    ) +
    panel_tag_theme
  
  
  
}

# ------------------------------------------------------------
# Inference helper
# ------------------------------------------------------------

.clean_volcano <- function(p, tag = NULL) {
  
  p <- set_ggrepel_text_size(p, size = 7)
  
  # Shorten package x-axis label if necessary.
  xlab <- p$labels$x
  if (is.character(xlab) && length(xlab) == 1L && !is.na(xlab)) {
    xlab <- sub("^Estimated effect size", "Effect size", xlab)
    p <- p + labs(x = xlab)
  }

  # Give repel labels a little more room at the top.
  built <- tryCatch(ggplot_build(p), error = function(e) NULL)
  if (!is.null(built) && length(built$data)) {
    yvals <- unlist(lapply(
      built$data,
      function(d) if (!is.null(d$y)) d$y else numeric(0)
    ), use.names = FALSE)
    yvals <- yvals[is.finite(yvals)]
    if (length(yvals) && max(yvals) > 0) {
      p <- p + coord_cartesian(
        ylim = c(NA_real_, 1.18 * max(yvals)),
        clip = "off"
      )
    }
  }

  p +
    labs(title = NULL, tag = tag) +
    paper_theme +
    theme(
      #aspect.ratio = 1#,
      #plot.margin = margin(t = 15, r = 15, b = 10, l = 10, unit = "pt")
    ) +
    panel_tag_theme
}

# ------------------------------------------------------------
# C5 decomposition helper
# ------------------------------------------------------------

.clean_term_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "_", x))
}

.select_c5_terms <- function(term_plots) {
  if (inherits(term_plots, "ggplot")) term_plots <- list(term_plots)
  term_plots <- Filter(function(x) inherits(x, "ggplot"), term_plots)

  nm <- names(term_plots)
  if (is.null(nm)) nm <- rep("", length(term_plots))
  key <- .clean_term_name(nm)

  # We deliberately show only the latitude, longitude, and latitude×longitude
  # terms in the BBS main figure. Elevation-containing terms remain in the
  # supplement.
  is_lat_main <- grepl("lat", key) & !grepl("lon", key) & !grepl("elev", key)
  is_lon_main <- grepl("lon", key) & !grepl("lat", key) & !grepl("elev", key)
  is_lat_lon  <- grepl("lat", key) & grepl("lon", key) & !grepl("elev", key)

  idx <- c(
    which(is_lat_main)[1],
    which(is_lon_main)[1],
    which(is_lat_lon)[1]
  )
  idx <- idx[is.finite(idx)]

  if (length(idx) != 3L) {
    stop(
      "Could not uniquely identify latitude, longitude, and latitude×longitude ",
      "plots from plotUTerms(). Available term names: ",
      paste(ifelse(nzchar(nm), nm, "<unnamed>"), collapse = ", "),
      call. = FALSE
    )
  }

  out <- term_plots[idx]
  names(out) <- c("Latitude", "Longitude", "Latitude × longitude")
  out
}

.clean_decomposition_plot <- function(p, title) {
  p +
    labs(title = title) +
    theme(
      text = element_text(size = 22, family = "Helvetica"),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 17),
      legend.title = element_text(size = 17),
      plot.margin = margin(t = 10, r = 8, b = 10, l = 8, unit = "pt")
    )
}

# ============================================================
# Main maker
# ============================================================

make_bbs_main_figure <- function(
    rw_single,
    output_dir,
    filename = "BBSMainFigure.pdf",
    width = 25,
    height = 15,
    grid_id_col = "grid_id",
    lat_col = "lat_center",
    lon_col = "lon_center",
    year_col = "year",
    map_xlim = c(-170, -50),
    map_ylim = c(20, 75),
    map_lon_breaks = c(-160, -120, -80, -60),
    map_point_size = 1.4,
    inference_label_top_n = 10L,
    svdx_report_args = list()
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  fit <- mySvdxBenchmark:::.benchmark_report_extract_svdx_fit(rw_single)
  inf <- mySvdxBenchmark:::.benchmark_report_extract_svdx_inf(rw_single)

  if (is.null(fit)) stop("No stored SSANOVA-SVD fit found.", call. = FALSE)
  if (is.null(inf)) stop("No stored SSANOVA-SVD inference found.", call. = FALSE)

  data <- as.data.frame(rw_single$data)
  U <- .bbs_svdx_U(fit)

  if (nrow(U) != nrow(data)) {
    stop("Stored SVDX U does not match rw_single$data.", call. = FALSE)
  }
  if (ncol(U) < 5L) {
    stop("BBS main figure requires at least five fitted components.", call. = FALSE)
  }

  north_america <- rnaturalearth::ne_countries(
    continent = "North America",
    returnclass = "sf"
  )

  get_comp_name <- mySvdx:::.report_get_inference_component_name

  # ============================================================
  # A-B: Component 1
  # ============================================================

  pA <- .bbs_component_map(
    scores = U[, 1],
    data = data,
    component = 1,
    north_america = north_america,
    grid_id_col = grid_id_col,
    lat_col = lat_col,
    lon_col = lon_col,
    year_col = year_col,
    point_size = map_point_size,
    map_xlim = map_xlim,
    map_ylim = map_ylim,
    map_lon_breaks = map_lon_breaks,
    tag = "A"
  )

  pB <- mySvdx::plot_component_volcano(
    inf,
    component = get_comp_name(inf, 1),
    fdata = rw_single$fdata,
    label_top_n = inference_label_top_n
  )
  pB <- .clean_volcano(pB)#, tag = "B")

  # ============================================================
  # C-D: Component 3
  # ============================================================

  pC <- .bbs_component_map(
    scores = U[, 3],
    data = data,
    component = 3,
    north_america = north_america,
    grid_id_col = grid_id_col,
    lat_col = lat_col,
    lon_col = lon_col,
    year_col = year_col,
    point_size = map_point_size,
    map_xlim = map_xlim,
    map_ylim = map_ylim,
    map_lon_breaks = map_lon_breaks,
    #tag = "C"
    tag = "B"
  )

  pD <- mySvdx::plot_component_volcano(
    inf,
    component = get_comp_name(inf, 3),
    fdata = rw_single$fdata,
    label_top_n = inference_label_top_n
  )
  pD <- .clean_volcano(pD)#, tag = "D")

  # ============================================================
  # E-G: Component 5
  # ============================================================

  pE <- .bbs_component_map(
    scores = U[, 5],
    data = data,
    component = 5,
    north_america = north_america,
    grid_id_col = grid_id_col,
    lat_col = lat_col,
    lon_col = lon_col,
    year_col = year_col,
    point_size = map_point_size,
    map_xlim = map_xlim,
    map_ylim = map_ylim,
    map_lon_breaks = map_lon_breaks,
    #tag = "E"
    tag = "C"
  )

  # Reuse the exact reporting argument path used by the BBS supplement.
  svdx_plot_args <- .get_ns(
    "mySvdxBenchmark",
    ".benchmark_report_svdx_plot_args"
  )

  term_obj <- do.call(
    mySvdx::plotUTerms,
    c(
      list(
        fit = fit,
        data = data,
        component = 5,
        facet = NULL,
        fill = "typical"
      ),
      svdx_plot_args(svdx_report_args, "plotUTerms")
    )
  )

  c5_terms <- .select_c5_terms(term_obj)

  if(FALSE){
    
    pF1 <- .clean_decomposition_plot(c5_terms[[1]], "Latitude")
    pF2 <- .clean_decomposition_plot(c5_terms[[2]], "Longitude")
    pF3 <- .clean_decomposition_plot(c5_terms[[3]] + coord_flip(), "Longitude × latitude")
    
    # F is one conceptual panel containing the three SSANOVA terms.
    pF <- wrap_plots(
      pF1, pF2, pF3,
      nrow = 1,
      guides = "collect"
    ) +
      plot_annotation(tag_levels = NULL) &
      theme(legend.position = "bottom")
    
    # Put the F tag on the decomposition block rather than each term.
    pF <- wrap_elements(full = pF) +
      labs(tag = "F") +
      panel_tag_theme
    
  }
  

  pG <- mySvdx::plot_component_volcano(
    inf,
    component = get_comp_name(inf, 5),
    fdata = rw_single$fdata,
    label_top_n = inference_label_top_n
  )
  pG <- .clean_volcano(pG)#, tag = "G")

  # ============================================================
  # Assembly
  # ============================================================

  row1 <- wrap_plots(
    pA, pB,
    nrow = 1,
    widths = c(1, 1)
  )

  row2 <- wrap_plots(
    pC, pD,
    nrow = 1,
    widths = c(1, 1)
  )

  # Give the three-term decomposition extra horizontal space.
  row3 <- wrap_plots(
    pE, #pF, 
    pG,
    nrow = 1,
    widths = c(1, #2.15, 
               1)
  )

  # fig_bbs <- wrap_plots(
  #   row1,
  #   row2,
  #   row3,
  #   ncol = 1,
  #   heights = c(1, 1, 1.08)
  # )
  
  # fig_bbs <- wrap_plots(
  #   pA, pC, pE,
  #   pB, pD, pG, 
  #   nrow=2, 
  #   widths=c(1,1,1)
  # )
  
  
  
  # fig_bbs <- wrap_plots(
  #   pA/pB, pC/pD, pE/pG,
  #   nrow=1,
  #   widths=c(1,1,1)
  # )
  
  pAB <- wrap_plots( pA, pB, nrow=2, heights = c(1,1) )
  pCD <- wrap_plots( pC, pD, nrow=2, heights = c(1,1) )
  pEG <- wrap_plots( pE, pG, nrow=2, heights = c(1,1) )
  fig_bbs <- wrap_plots(
    pAB, pCD, pEG,
    nrow=1,
    widths=c(1,1,1)
  )
  
  
  

  output_file <- file.path(output_dir, filename)

  ggsave(
    filename = output_file,
    plot = fig_bbs,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  invisible(fig_bbs)
}
