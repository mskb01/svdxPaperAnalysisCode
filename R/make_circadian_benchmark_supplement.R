# Designed against the current svdxBenchmark::runSingleBenchmarkReports()
# and svdx reporting APIs. Dataset-specific plotting metadata are passed
# through the same arguments used by runSingleBenchmarkReports(); nothing is
# inferred from hard-coded column names.

`%||%` <- function(x, y) if (!is.null(x)) x else y

.supp_get_ns <- function(pkg, name) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.", call. = FALSE)
  }
  getFromNamespace(name, pkg)
}

.supp_require <- function() {
  pkgs <- c("ggplot2", "patchwork", "mySvdx", "mySvdxBenchmark")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

.supp_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}

.supp_component_pairs <- function(components) {
  if (!length(components)) return(list())
  split(components, ceiling(seq_along(components) / 2))
}

.supp_chunks <- function(x, n) {
  if (!length(x)) return(list())
  split(x, ceiling(seq_along(x) / n))
}

.supp_theme <- function(text_size = 25) {
  ggplot2::theme(
    text = ggplot2::element_text(size = text_size),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 0.75 * text_size),
    legend.title = ggplot2::element_text(size = 0.75 * text_size),
    aspect.ratio = 1,
    plot.tag = ggplot2::element_text(size = 0.75 * text_size, face = "bold")
  )
}

.supp_clean_plot <- function(p, text_size = 25, keep_title = FALSE) {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)
  if (!isTRUE(keep_title)) p <- p + ggplot2::labs(title = NULL)
  p + .supp_theme(text_size = text_size)
}

# Use Cairo PDF when available so standard Helvetica text is embedded for
# downstream vector editing (e.g. Illustrator) without changing requested
# text sizes. Fall back to the base PDF device if Cairo is unavailable.
.supp_pdf_device <- function(filename, width, height, ...) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(
      filename = filename, width = width, height = height,
      family = "Helvetica", ...
    )
  } else {
    grDevices::pdf(
      file = filename, width = width, height = height,
      family = "Helvetica", ...
    )
  }
}

.supp_open_pdf <- function(file, width, height) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(file, width = width, height = height, family = "Helvetica")
  } else {
    grDevices::pdf(file, width = width, height = height, family = "Helvetica")
  }
}

# A continuous fill legend in a decomposition heatmap needs an explicitly
# horizontal color bar. Keep the established legend typography unchanged;
# only allocate enough physical bar width for its tick labels.
.supp_decomposition_colorbar <- function(p) {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)
  sc <- p$scales$get_scales("fill")
  if (is.null(sc) || !inherits(sc, "ScaleContinuous")) return(p)
  p + ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = grid::unit(3.5, "in"),
      barheight = grid::unit(0.22, "in")
    )
  )
}

.supp_is_discrete_x <- function(data, x) {
  if (is.null(x) || !x %in% colnames(data)) return(FALSE)
  z <- data[[x]]
  is.factor(z) || is.character(z) || is.logical(z)
}

.supp_discrete_axis_theme <- function(text_size = 25, hide_legend = FALSE) {
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 0.55 * text_size
    ),
    plot.margin = ggplot2::margin(t = 12, r = 12, b = 14, l = 22, unit = "pt"),
    legend.position = if (isTRUE(hide_legend)) "none" else "bottom"
  )
}

# Keep small legends unchanged. For larger displayed grouping interactions,
# explicitly cap the number of legend columns so ggplot does not choose a
# horizontally over-wide layout. The rule is based on the number of observed
# group combinations, not on the number of grouping variables.
.supp_group_combination_count <- function(data, group) {
  if (is.null(group) || !length(group) || !all(group %in% names(data))) return(0L)
  z <- data[group]
  keep <- stats::complete.cases(z)
  if (!any(keep)) return(0L)
  if (length(group) == 1L) {
    return(length(unique(z[[1L]][keep])))
  }
  length(unique(do.call(
    interaction,
    c(unname(z[keep, , drop = FALSE]), list(drop = TRUE, sep = " - "))
  )))
}

# For SSANOVA-SVD, an unordered/ordered discrete x represents fitted
# component-function values, not sample-level scores. Display these fitted
# effects as columns about zero. Comparator plots remain sample-level points.
.supp_svdx_discrete_columns <- function(p, x, group = NULL) {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)
  p$layers <- list()
  if (!is.null(group) && group %in% names(p$data)) {
    p <- p + ggplot2::geom_col(
      ggplot2::aes(fill = .data[[group]]),
      width = 0.72
    )
  } else {
    p <- p + ggplot2::geom_col(width = 0.72)
  }
  p
}

.supp_shorten_volcano_xlab <- function(p) {
  if (is.null(p) || !inherits(p, "ggplot")) return(p)
  xlab <- p$labels$x
  if (is.character(xlab) && length(xlab) == 1L && !is.na(xlab)) {
    xlab <- sub("^Estimated effect size", "Effect size", xlab)
    p <- p + ggplot2::labs(x = xlab)
  }
  p
}

.supp_panel_caption <- function(components) {
  components <- as.integer(components)
  if (length(components) == 1L) return(paste0("Component ", components, "."))
  tags <- LETTERS[seq_along(components)]
  paste0("(", tags, ") Component ", components, ".", collapse = " ")
}

.supp_method_result <- function(rw_single, method) {
  getter <- .supp_get_ns("mySvdxBenchmark", ".rwbench_get_method_result")
  tryCatch(getter(rw_single, method), error = function(e) NULL)
}

.supp_method_display_name <- function(rw_single, method) {
  method <- as.character(method)
  method_lower <- tolower(method)

  # Canonical display labels for benchmark method identifiers.
  if (identical(method_lower, "pca")) return("PCA")
  if (identical(method_lower, "nmf")) return("NMF")
  if (identical(method_lower, "mixomics_pls")) return("PLS")
  if (identical(method_lower, "mefisto")) return("MEFISTO")
  if (grepl("^mefisto_duplicate_view", method_lower)) return("MEFISTO: duplicate views")

  res <- .supp_method_result(rw_single, method)
  hp_name <- res$extra$hyperparameter_name %||% NULL
  hp <- res$extra$hyperparameter %||% NULL

  if (startsWith(method, "sPCA_keepX_")) {
    value <- if (!is.null(hp) && length(hp)) hp[[1L]] else sub("^sPCA_keepX_", "", method)
    return(paste0("sPCA: keepX=", format(value, trim = TRUE, scientific = FALSE)))
  }
  if (startsWith(method, "sNMF_lambda_")) {
    if (!is.null(hp) && length(hp)) {
      value <- format(hp[[1L]], trim = TRUE, scientific = FALSE)
    } else {
      value <- sub("^sNMF_lambda_", "", method)
      value <- sub("p", ".", value, fixed = TRUE)
    }
    return(paste0("sNMF: lambda=", value))
  }

  if (!is.null(hp_name) && length(hp_name) && !is.na(hp_name[[1L]]) && nzchar(hp_name[[1L]]) &&
      !is.null(hp) && length(hp)) {
    return(paste0(method, ": ", hp_name[[1L]], "=", format(hp[[1L]], trim = TRUE, scientific = FALSE)))
  }
  method
}


.supp_add_tags <- function(p, n_panels) {
  if (n_panels <= 1L) return(p)
  p + patchwork::plot_annotation(tag_levels = "A")
}

.supp_save_gg <- function(p, file, n_panels = 1L,
                          panel_width = 5, panel_height = 5.5,
                          text_size = 25) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (is.null(p)) return(NULL)
  p <- .supp_add_tags(p, n_panels = n_panels)
  ggplot2::ggsave(
    filename = file,
    plot = p,
    width = panel_width * if (n_panels > 1L) 2 else 1,
    height = panel_height,
    units = "in",
    device = .supp_pdf_device
  )
  normalizePath(file, mustWork = FALSE)
}

.supp_save_heatmap <- function(ht, file, width = 10, height = 8) {
  if (is.null(ht)) return(NULL)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to save inference heatmaps.", call. = FALSE)
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  .supp_open_pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)
  normalizePath(file, mustWork = FALSE)
}

.supp_write_tsv_allow_empty <- function(df, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    as.data.frame(df), file = file, sep = "\t",
    row.names = FALSE, quote = FALSE, na = ""
  )
  normalizePath(file, mustWork = FALSE)
}



# Export complete feature-wise conditional-inference results.  These are audit
# tables, not significance-filtered tables: the FDR decision is recorded in
# explicit columns so every tested feature remains available downstream.
.supp_export_feature_inference_tables <- function(
    inf,
    components,
    tables_dir,
    get_inf_component_name,
    fdr = 0.1
) {
  if (is.null(inf)) return(character(0))
  if (!is.numeric(fdr) || length(fdr) != 1L || !is.finite(fdr) || fdr <= 0 || fdr > 1) {
    stop("fdr must be a single number in (0, 1].", call. = FALSE)
  }

  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  saved <- character(0)
  component_tables <- list()

  fr <- inf$feature_results
  if (!is.null(fr) && is.data.frame(fr) && nrow(fr)) {
    fr <- as.data.frame(fr)
    if (!all(c("feature", "component", "qval") %in% names(fr))) {
      stop("Inference feature_results is missing one or more required columns: feature, component, qval.", call. = FALSE)
    }

    for (cc in components) {
      comp_name <- get_inf_component_name(inf, cc)
      d <- fr[fr$component == comp_name, , drop = FALSE]
      d$significant_fdr <- !is.na(d$qval) & is.finite(d$qval) & d$qval <= fdr
      d$fdr_threshold <- fdr
      component_tables[[as.character(cc)]] <- d
      file <- file.path(tables_dir, sprintf("svdx_inference_component_%02d.tsv", cc))
      saved <- c(saved, .supp_write_tsv_allow_empty(d, file))
    }

    all_components <- if (length(component_tables)) {
      do.call(rbind, component_tables)
    } else {
      fr[0, , drop = FALSE]
    }
    rownames(all_components) <- NULL
    saved <- c(saved, .supp_write_tsv_allow_empty(
      all_components,
      file.path(tables_dir, "svdx_inference_component_all.tsv")
    ))
  }

  ur <- inf$union_results
  if (!is.null(ur) && is.data.frame(ur) && nrow(ur)) {
    ur <- as.data.frame(ur)

    d_union <- ur
    d_union$significant_fdr <- if ("union_q" %in% names(d_union)) {
      !is.na(d_union$union_q) & is.finite(d_union$union_q) & d_union$union_q <= fdr
    } else {
      rep(NA, nrow(d_union))
    }
    d_union$fdr_threshold <- fdr
    saved <- c(saved, .supp_write_tsv_allow_empty(
      d_union,
      file.path(tables_dir, "svdx_inference_union.tsv")
    ))

    d_global <- ur
    d_global$significant_fdr <- if ("global_q" %in% names(d_global)) {
      !is.na(d_global$global_q) & is.finite(d_global$global_q) & d_global$global_q <= fdr
    } else {
      rep(NA, nrow(d_global))
    }
    d_global$fdr_threshold <- fdr
    saved <- c(saved, .supp_write_tsv_allow_empty(
      d_global,
      file.path(tables_dir, "svdx_inference_global.tsv")
    ))
  }

  saved
}

.supp_escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

# Escape ordinary caption text while rendering <= thresholds in LaTeX math.
.supp_caption_tex <- function(x) {
  x <- gsub("<=\\s*([0-9]+(?:\\.[0-9]+)?)", "SUPPLEQMARK\\1SUPPENDMARK", x, perl = TRUE)
  x <- .supp_escape_tex(x)
  gsub("SUPPLEQMARK([0-9]+(?:\\.[0-9]+)?)SUPPENDMARK", "$\\\\leq \\1$", x, perl = TRUE)
}

.supp_manifest_add <- function(manifest, figure_number, section, method,
                               components = NA_character_, plot_type,
                               file, page = 1L, caption) {
  row <- data.frame(
    figure_number = as.integer(figure_number),
    section = as.character(section),
    method = as.character(method),
    components = as.character(components),
    plot_type = as.character(plot_type),
    page = as.integer(page),
    file = as.character(file),
    caption = as.character(caption),
    stringsAsFactors = FALSE
  )
  rbind(manifest, row)
}

.supp_label_slug <- function(x) {
  x <- tolower(as.character(x))
  x[is.na(x)] <- ""
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("^-+|-+$", "", x)
  x
}

.supp_assign_labels <- function(manifest, figure_label_prefix) {
  if (!nrow(manifest)) {
    manifest$label <- character(0)
    return(manifest)
  }
  prefix <- as.character(figure_label_prefix)[1L]
  if (is.na(prefix) || !nzchar(prefix)) stop("figure_label_prefix must be a single non-empty value.", call. = FALSE)
  if (!grepl("^[A-Za-z0-9_-]+$", prefix)) stop("figure_label_prefix may contain only letters, numbers, underscores, and hyphens.", call. = FALSE)
  suffix <- vapply(seq_len(nrow(manifest)), function(i) {
    parts <- c(manifest$method[i], manifest$plot_type[i], manifest$components[i])
    parts <- .supp_label_slug(parts)
    parts <- parts[nzchar(parts)]
    if (!length(parts)) "figure" else paste(parts, collapse = "-")
  }, character(1))
  suffix <- make.unique(suffix, sep = "-")
  manifest$label <- paste0("fig:", prefix, "-", suffix)
  manifest
}

.supp_write_tex <- function(manifest, tex_file, note_number, dataset_title) {
  if (!nrow(manifest)) return(NULL)
  dir.create(dirname(tex_file), recursive = TRUE, showWarnings = FALSE)

  rel <- function(path) {
    p <- normalizePath(path, mustWork = FALSE)
    root <- normalizePath(dirname(tex_file), mustWork = FALSE)
    prefix <- paste0(root, .Platform$file.sep)
    if (startsWith(p, prefix)) substring(p, nchar(prefix) + 1L) else p
  }

  lines <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage[margin=0.65in]{geometry}",
    "\\usepackage{graphicx}",
    "\\usepackage{caption}",
    "\\usepackage{float}",
    "\\usepackage{placeins}",
    "\\setlength{\\parindent}{0pt}",
    "\\renewcommand{\\figurename}{Supplementary Figure}",
    "\\captionsetup[figure]{labelfont=bf,labelsep=period}",
    "\\begin{document}",
    paste0("\\setcounter{figure}{", min(manifest$figure_number) - 1L, "}"),
    paste0("\\section*{Supplementary Note ", .supp_escape_tex(note_number), "}"),
    paste0("This note presents the supplementary figures for ", .supp_escape_tex(dataset_title), "."),
    ""
  )

  for (i in seq_len(nrow(manifest))) {
    r <- manifest[i, , drop = FALSE]
    cap <- .supp_caption_tex(r$caption)
    lines <- c(
      lines,
      "\\begin{figure}[p]",
      "\\centering",
      paste0("\\includegraphics[width=\\textwidth,height=0.82\\textheight,keepaspectratio]{", rel(r$file), "}"),
      paste0("\\caption{", cap, "}"),
      paste0("\\label{", r$label, "}"),
      "\\end{figure}",
      "\\clearpage"
    )
  }

  lines <- c(lines, "\\end{document}")
  writeLines(lines, tex_file)
  normalizePath(tex_file, mustWork = FALSE)
}

.supp_compile_tex <- function(tex_file, verbose = TRUE) {
  latexmk <- Sys.which("latexmk")
  if (!nzchar(latexmk)) {
    if (isTRUE(verbose)) message("latexmk not found; wrote TeX but did not compile it.")
    return(NULL)
  }
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(dirname(tex_file))
  status <- system2(
    latexmk,
    args = c("-pdf", "-interaction=nonstopmode", "-halt-on-error", basename(tex_file)),
    stdout = if (isTRUE(verbose)) "" else FALSE,
    stderr = if (isTRUE(verbose)) "" else FALSE
  )
  if (!identical(status, 0L)) {
    warning("latexmk returned non-zero status: ", status, call. = FALSE)
    return(NULL)
  }
  pdf <- sub("\\.tex$", ".pdf", tex_file, ignore.case = TRUE)
  if (file.exists(pdf)) normalizePath(pdf) else NULL
}

#' Generate supplementary figures from one real-world benchmark object.
#'
#' @param rw_single A realWorldMethodBenchmark returned by
#'   run_real_world_single_benchmark().
#' @param outdir Destination directory supplied by the caller.
#' @param note_number Supplementary Note number displayed on the title page.
#' @param dataset_title Dataset title displayed on the title page.
#' @param start_figure First Supplementary Figure number.
#' @param methods Comparator methods. NULL uses every benchmark_results method
#'   other than svdx and svdx_adjusted, if present.
#' @param components SSANOVA-SVD components to include. NULL means all fitted
#'   components.
#' @param u_plot_specs,xaxis,color,group,line_group,facet Same role as in
#'   runSingleBenchmarkReports().
#' @param features_top_n Number of features shown in each inference or loading heatmap.
#' @param loading_heatmaps If TRUE, add component-specific heatmaps based on the largest absolute SSANOVA-SVD loadings.
#' @param loading_enrichment If TRUE, add component-specific loading-based enrichment plots when stored loading enrichment is available.
#' @param enrichment_top_n Maximum enrichment terms displayed per component.
#' @param enrichment_fdr Strict adjusted-P threshold; default padj < 0.1.
#' @param volcano_labels_top_n Number of volcano labels.
#' @param volcano_label_features Optional explicit volcano labels.
#' @param svdx_report_args Named plot-specific argument lists, following the
#'   runSingleBenchmarkReports() convention.
#' @param text_size Source-PDF base text size.
#' @param compile Compile the generated TeX with latexmk when available. Defaults to FALSE.
#' @param verbose Print progress messages.
#'
#' @return Invisibly returns paths, manifest, and next figure number.
make_circadian_benchmark_supplement <- function(
    rw_single,
    outdir,
    note_number,
    dataset_title,
    start_figure = 1L,
    methods = NULL,
    components = NULL,
    u_plot_specs = NULL,
    color = NULL,
    group = NULL,
    line_group = NULL,
    xaxis = NULL,
    facet = NULL,
    features_top_n = 50L,
    loading_heatmaps = FALSE,
    loading_enrichment = FALSE,
    enrichment_top_n = 20L,
    enrichment_fdr = 0.1,
    volcano_labels_top_n = 10L,
    volcano_label_features = NULL,
    svdx_report_args = list(),
    text_size = 25,
    compile = FALSE,
    verbose = TRUE
) {
  .supp_require()

  if (!inherits(rw_single, "realWorldMethodBenchmark")) {
    stop("rw_single must inherit from 'realWorldMethodBenchmark'.", call. = FALSE)
  }
  if (missing(outdir) || !length(outdir) || !nzchar(outdir)) {
    stop("outdir must be supplied.", call. = FALSE)
  }
  if (missing(note_number) || length(note_number) != 1L || is.na(note_number) || !nzchar(as.character(note_number))) {
    stop("note_number must be supplied as a single non-empty value.", call. = FALSE)
  }
  if (missing(dataset_title) || length(dataset_title) != 1L || is.na(dataset_title) || !nzchar(as.character(dataset_title))) {
    stop("dataset_title must be supplied as a single non-empty value.", call. = FALSE)
  }
  if (!is.numeric(enrichment_fdr) || length(enrichment_fdr) != 1L ||
      !is.finite(enrichment_fdr) || enrichment_fdr <= 0 || enrichment_fdr > 1) {
    stop("enrichment_fdr must be a single number in (0, 1].", call. = FALSE)
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fig_dir <- file.path(outdir, "figures")
  svdx_dir <- file.path(fig_dir, "svdx")
  comparator_dir <- file.path(fig_dir, "comparators")
  enr_dir <- file.path(outdir, "tables", "enrichment")
  inf_table_dir <- file.path(outdir, "tables", "inference")
  dir.create(svdx_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(comparator_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(enr_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(inf_table_dir, recursive = TRUE, showWarnings = FALSE)

  # Exact helpers used by runSingleBenchmarkReports().
  prepare_u_specs <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_prepare_u_plot_specs")
  svdx_u_specs <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_svdx_u_plot_specs")
  spec_value <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_spec_value")
  plotU_call_args <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_plotU_call_args")
  add_connected_u <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_add_connected_u_layer")
  svdx_plot_args <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_svdx_plot_args")
  extract_fit <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_extract_svdx_fit")
  extract_inf <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_extract_svdx_inf")
  extract_enr <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_extract_svdx_inference_enrichment")
  extract_loading_enr <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_extract_svdx_loading_enrichment")
  loading_component_names <- .supp_get_ns("mySvdxBenchmark", ".benchmark_loading_component_names")
  normalize_enr <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_enrichment_dotplot_object")
  component_count <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_component_count_from_fit")
  u_plot_one <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_u_plot_one")
  report_methods <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_methods")
  filter_method_components <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_filter_components_for_method")
  svdx_model_vars_fn <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_svdx_model_vars")
  svdx_process_vars_fn <- .supp_get_ns("mySvdxBenchmark", ".benchmark_report_svdx_process_vars")

  plotU_fn <- .supp_get_ns("mySvdx", "plotU")
  plotUTerms_fn <- .supp_get_ns("mySvdx", "plotUTerms")
  plotVHeatmap_fn <- .supp_get_ns("mySvdx", "plotVHeatmap")
  plot_component_volcano_fn <- .supp_get_ns("mySvdx", "plot_component_volcano")
  plot_enrichment_dotplot_fn <- .supp_get_ns("mySvdx", "plot_enrichment_dotplot")
  get_inf_component_name <- .supp_get_ns("mySvdx", ".report_get_inference_component_name")

  fit <- extract_fit(rw_single)
  if (is.null(fit)) stop("No stored SSANOVA-SVD fit found in rw_single.", call. = FALSE)

  inf <- extract_inf(rw_single)

  # Complete feature-wise conditional-inference flat files at the paper-wide
  # 10% FDR threshold.  Rows are not filtered; significance is annotated.
  .supp_export_feature_inference_tables(
    inf = inf,
    components = components,
    tables_dir = inf_table_dir,
    get_inf_component_name = get_inf_component_name,
    fdr = 0.1
  )

  data <- as.data.frame(rw_single$data)
  Y <- as.matrix(rw_single$Y)
  zdata <- if (!is.null(rw_single$zdata)) as.data.frame(rw_single$zdata) else NULL

  if (ncol(Y) != nrow(data)) {
    stop("ncol(rw_single$Y) must equal nrow(rw_single$data).", call. = FALSE)
  }
  if (!is.null(zdata) && nrow(zdata) != ncol(Y)) {
    stop("nrow(rw_single$zdata) must equal ncol(rw_single$Y).", call. = FALSE)
  }

  prepared_specs <- prepare_u_specs(
    data = data,
    xformula = rw_single$xformula,
    u_plot_specs = u_plot_specs,
    xaxis = xaxis
  )
  xaxes <- vapply(
    prepared_specs,
    function(sp) spec_value(sp, "var", spec_value(sp, "x", NA_character_)),
    character(1)
  )
  u_specs <- svdx_u_specs(data, xaxes = xaxes, u_plot_specs = prepared_specs)

  # SSANOVA-SVD model variables determine prediction; plotting variables are
  # display coordinates and may also come from the supplied design data.
  # plotU() decides whether to evaluate on a model-space grid or at supplied
  # design points.
  svdx_model_vars <- svdx_model_vars_fn(fit, fallback_formula = rw_single$xformula)
  process_facet <- facet

  K_eff <- component_count(fit, default = 0L)
  if (!is.finite(K_eff) || K_eff < 1L) stop("Could not determine fitted SSANOVA-SVD component count.", call. = FALSE)
  if (is.null(components)) {
    components <- seq_len(K_eff)
  } else {
    components <- unique(suppressWarnings(as.integer(components)))
    components <- components[!is.na(components) & components >= 1L & components <= K_eff]
    if (!length(components)) stop("No requested SSANOVA-SVD components are available.", call. = FALSE)
  }

  manifest <- data.frame(
    figure_number = integer(0), section = character(0), method = character(0),
    components = character(0), plot_type = character(0), page = integer(0),
    file = character(0), caption = character(0), stringsAsFactors = FALSE
  )
  fig_no <- as.integer(start_figure)

  # ------------------------------------------------------------------------
  # 1. SSANOVA-SVD component functions: paired, for each report plot spec.
  #
  # Circadian-specific handling for the time plot only:
  # the fitted SSANOVA-SVD process is a function of time, while tissue_type is
  # not part of that fitted process model. We therefore evaluate the fitted
  # time term once on the dense plotUTerms() grid and replicate that same
  # fitted curve across the observed tissue_type levels only so tissue_type can
  # be used as a plotting colour/group. This does NOT create tissue-specific
  # fitted functions or tissue-specific predictions.
  # ------------------------------------------------------------------------
  for (spec_i in seq_along(u_specs)) {
    spec <- u_specs[[spec_i]]
    xa <- spec_value(spec, "var", spec_value(spec, "x", NULL))
    if (is.null(xa) || !xa %in% colnames(data)) next
    spec_group <- spec_value(spec, "group", group)
    spec_line_group <- spec_value(spec, "line_group", line_group)
    # plotU() owns its default geometry. connect=TRUE is retained only as an
    # explicit request to add the benchmark helper's connected-line layer.
    spec_connect <- "connect" %in% names(spec) && isTRUE(spec$connect)

    for (pair in .supp_component_pairs(components)) {
      plots <- lapply(pair, function(cc) {
        if (identical(as.character(xa), "time") && "tissue_type" %in% colnames(data)) {
          p_term <- plotUTerms_fn(
            fit = fit,
            data = data,
            component = cc,
            terms = "time",
            fill = "typical",
            n_grid = 200
          )
          if (is.list(p_term) && !inherits(p_term, "ggplot")) {
            p_term <- p_term[["time"]]
          }
          if (is.null(p_term) || !inherits(p_term, "ggplot")) {
            stop("plotUTerms() did not return a ggplot for the Circadian time term.", call. = FALSE)
          }

          curve_df <- as.data.frame(p_term$data)
          if (!all(c("time", "contribution") %in% colnames(curve_df))) {
            stop("plotUTerms() output is missing the expected time/contribution columns.", call. = FALSE)
          }
          curve_df <- curve_df[order(curve_df$time), , drop = FALSE]

          tissue_values <- unique(data$tissue_type)
          tissue_values <- tissue_values[!is.na(tissue_values)]
          plot_df <- do.call(rbind, lapply(tissue_values, function(tt) {
            d <- curve_df
            d$tissue_type <- tt
            d
          }))
          rownames(plot_df) <- NULL
          if (is.factor(data$tissue_type)) {
            plot_df$tissue_type <- factor(plot_df$tissue_type, levels = levels(data$tissue_type))
          }

          p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(
              x = .data$time,
              y = .data$contribution,
              colour = .data$tissue_type,
              group = .data$tissue_type
            )
          ) +
            ggplot2::geom_line(linewidth = 1) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank()
            )
        } else {
          p <- do.call(
            plotU_fn,
            plotU_call_args(
              fit = fit,
              data = data,
              component = cc,
              x = xa,
              group = spec_group,
              facet = process_facet,
              fill = "typical",
              xformula = rw_single$xformula,
              plot_args = svdx_plot_args(svdx_report_args, "plotU")
            )
          )
          p <- add_connected_u(
            p, connect = spec_connect,
            group = spec_group, line_group = spec_line_group
          )
        }
        discrete_x <- .supp_is_discrete_x(data, xa)
        if (isTRUE(discrete_x)) {
          p <- .supp_svdx_discrete_columns(p, x = xa, group = spec_group)
        }
        p <- .supp_clean_plot(p, text_size = text_size, keep_title = FALSE) +
          ggplot2::ylab(paste0("Component Function ", cc))
        if (isTRUE(discrete_x)) {
          # If the grouping variable is the x variable, its legend repeats the
          # category labels and only steals panel area.
          hide_legend <- !is.null(spec_group) && identical(as.character(spec_group), as.character(xa))
          p <- p + ggplot2::labs(x = NULL) +
            .supp_discrete_axis_theme(text_size = text_size, hide_legend = hide_legend)
        }
        p
      })
      discrete_x <- .supp_is_discrete_x(data, xa)
      p <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
        ggplot2::theme(legend.position = if (isTRUE(discrete_x) && !is.null(spec_group) && identical(as.character(spec_group), as.character(xa))) "none" else "bottom")
      stem <- sprintf("svdx_components_%s_%s", .supp_safe_name(xa), paste(sprintf("%02d", pair), collapse = "_"))
      file <- file.path(svdx_dir, paste0(stem, ".pdf"))
      if (isTRUE(discrete_x)) {
        file <- .supp_save_gg(
          p, file, n_panels = length(pair),
          panel_width = 6.5, panel_height = 7,
          text_size = text_size
        )
      } else {
        # Preserve the existing 5-inch panel width for small legends. When a
        # displayed grouping produces many observed combinations, give the
        # collected bottom legend more horizontal canvas rather than changing
        # its typography or guide layout.
        n_group_comb <- .supp_group_combination_count(data, spec_group)
        n_pair_panels <- length(pair)
        if (n_group_comb <= 4L) {
          # Preserve the previous sizing exactly for ordinary legends.
          panel_width <- 5
        } else {
          # Large legends need a minimum total figure width even when the last
          # component is alone on its page. This keeps a 6-combination legend
          # at 12 inches for both a 2-panel pair and a 1-panel remainder.
          min_total_width <- if (n_group_comb <= 8L) 12 else 14
          panel_width <- min_total_width / if (n_pair_panels > 1L) 2 else 1
        }
        file <- .supp_save_gg(
          p, file, n_panels = n_pair_panels,
          panel_width = panel_width, text_size = text_size
        )
      }
      manifest <- .supp_manifest_add(
        manifest, fig_no, "SSANOVA-SVD", "svdx",
        components = paste(pair, collapse = ","), plot_type = "component_functions",
        file = file,
        caption = paste0("SSANOVA-SVD reconstructed component functions.", .supp_panel_caption(pair))
      )
      fig_no <- fig_no + 1L
    }
  }

  # ------------------------------------------------------------------------
  # 2. Component inference volcanoes: paired.
  # ------------------------------------------------------------------------
  if (!is.null(inf)) {
    for (pair in .supp_component_pairs(components)) {
      plots <- lapply(pair, function(cc) {
        comp_name <- get_inf_component_name(inf, cc)
        p <- plot_component_volcano_fn(
          inf,
          component = comp_name,
          fdata = rw_single$fdata,
          label_top_n = volcano_labels_top_n,
          label_features = volcano_label_features
        )
        p <- .supp_clean_plot(p, text_size = text_size, keep_title = FALSE)
        .supp_shorten_volcano_xlab(p)
      })
      p <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
      file <- file.path(
        svdx_dir,
        paste0("svdx_volcanoes_", paste(sprintf("%02d", pair), collapse = "_"), ".pdf")
      )
      file <- .supp_save_gg(p, file, n_panels = length(pair), text_size = text_size)
      manifest <- .supp_manifest_add(
        manifest, fig_no, "SSANOVA-SVD", "svdx",
        components = paste(pair, collapse = ","), plot_type = "inference_volcano",
        file = file,
        caption = paste0("SSANOVA-SVD based conditional inference for feature participation. ", .supp_panel_caption(pair))
      )
      fig_no <- fig_no + 1L
    }

    # ----------------------------------------------------------------------
    # 3. Inference heatmaps: exactly one component per figure/page.
    # ----------------------------------------------------------------------
    for (cc in components) {
      comp_name <- get_inf_component_name(inf, cc)
      ht <- plotVHeatmap_fn(
        object = inf,
        Y = Y,
        data = data,
        zdata = zdata,
        inference_scope = "component",
        component = comp_name,
        select_by = "q",
        fdr = 0.1,
        n_features = features_top_n,
        order_by = svdx_report_args$heatmap_order_by %||% NULL,
        annotate = svdx_report_args$heatmap_annotate %||% NULL,
        fdata = rw_single$fdata
      )
      if (is.null(ht)) next
      file <- file.path(svdx_dir, sprintf("svdx_inference_heatmap_component_%02d.pdf", cc))
      file <- .supp_save_heatmap(ht, file)
      manifest <- .supp_manifest_add(
        manifest, fig_no, "SSANOVA-SVD", "svdx",
        components = as.character(cc), plot_type = "inference_heatmap",
        file = file,
        caption = paste0("Heatmap of features inferred to be associated with component ", cc, " by SSANOVA-SVD.")
      )
      fig_no <- fig_no + 1L
    }
  } else if (isTRUE(verbose)) {
    message("No stored SSANOVA-SVD inference object; volcanoes and inference heatmaps skipped.")
  }

  # ------------------------------------------------------------------------
  # 4. Inference enrichment: flat files for every component; plot only rows
  #    with padj < enrichment_fdr, one component per figure.
  # ------------------------------------------------------------------------
  enr_raw <- extract_enr(rw_single)
  enr_obj <- if (!is.null(enr_raw)) {
    normalize_enr(enr_raw, context = "supplement inference enrichment", verbose = verbose)
  } else NULL

  if (!is.null(inf) && !is.null(enr_obj) && is.data.frame(enr_obj$results)) {
    enr_df <- as.data.frame(enr_obj$results)
    if (!all(c("scope", "component", "padj") %in% names(enr_df))) {
      stop("Normalized inference enrichment is missing one or more required columns: scope, component, padj.", call. = FALSE)
    }

    component_name_map <- setNames(
      vapply(components, function(cc) get_inf_component_name(inf, cc), character(1)),
      as.character(components)
    )

    component_pass <- enr_df[0, , drop = FALSE]
    pass_by_component <- vector("list", length(components))
    names(pass_by_component) <- as.character(components)

    fdr_tag <- format(enrichment_fdr, trim = TRUE, scientific = FALSE)

    for (cc in components) {
      comp_name <- component_name_map[[as.character(cc)]]
      d <- enr_df[
        enr_df$scope == "component" &
          enr_df$component == comp_name &
          is.finite(enr_df$padj) &
          enr_df$padj < enrichment_fdr,
        , drop = FALSE
      ]
      pass_by_component[[as.character(cc)]] <- d
      component_pass <- rbind(component_pass, d)
      .supp_write_tsv_allow_empty(
        d,
        file.path(enr_dir, sprintf("svdx_inference_enrichment_component_%02d_fdr%s.tsv", cc, fdr_tag))
      )
    }

    # Explicit aggregate of component-level enrichment.  Do not use an
    # unsuffixed filename: that is ambiguous with union/global inference.
    .supp_write_tsv_allow_empty(
      component_pass,
      file.path(enr_dir, paste0("svdx_inference_enrichment_component_all_fdr", fdr_tag, ".tsv"))
    )

    # Export the non-component inference scopes separately when they are
    # present in the stored enrichment results.  These are not plotted in the
    # supplement, but retaining them as flat files makes the inferential scope
    # explicit and preserves the complete analysis output.
    for (scope_name in c("union", "global")) {
      d_scope <- enr_df[
        enr_df$scope == scope_name &
          is.finite(enr_df$padj) &
          enr_df$padj < enrichment_fdr,
        , drop = FALSE
      ]
      .supp_write_tsv_allow_empty(
        d_scope,
        file.path(
          enr_dir,
          paste0("svdx_inference_enrichment_", scope_name, "_fdr", fdr_tag, ".tsv")
        )
      )
    }

    # Complete FDR-filtered enrichment table, retaining the scope column so
    # component, union and global results remain distinguishable.
    all_scope_pass <- enr_df[
      is.finite(enr_df$padj) & enr_df$padj < enrichment_fdr,
      , drop = FALSE
    ]
    .supp_write_tsv_allow_empty(
      all_scope_pass,
      file.path(enr_dir, paste0("svdx_inference_enrichment_all_scopes_fdr", fdr_tag, ".tsv"))
    )

    for (cc in components) {
      d <- pass_by_component[[as.character(cc)]]
      if (is.null(d) || !nrow(d)) next
      obj <- structure(list(results = d), class = c("svdxEnrichment", "list"))
      # Keep the package's default enrichment styling and place exactly one
      # component enrichment plot on each supplementary figure page.
      p <- plot_enrichment_dotplot_fn(obj, top_n = enrichment_top_n) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = 0.75 * text_size),
          legend.title = ggplot2::element_text(size = 0.75 * text_size)
        )
      file <- file.path(svdx_dir, sprintf("svdx_enrichment_component_%02d.pdf", cc))
      file <- .supp_save_gg(p, file, n_panels = 1L, panel_width = 8, panel_height = 6, text_size = text_size)
      manifest <- .supp_manifest_add(
        manifest, fig_no, "SSANOVA-SVD", "svdx",
        components = as.character(cc), plot_type = "inference_enrichment",
        file = file,
        caption = paste0("Enrichment among features inferred to be associated with SSANOVA-SVD component ", cc,
                         "; only terms with adjusted P <= ", enrichment_fdr, " are shown.")
      )
      fig_no <- fig_no + 1L
    }
  } else if (isTRUE(verbose)) {
    message("No usable SSANOVA-SVD inference enrichment (or no inference object); inference enrichment tables/plots skipped.")
  }

  # ------------------------------------------------------------------------
  # 5. Optional loading-based heatmaps and enrichment. These are deliberately
  #    separate from conditional-inference displays and may be requested even
  #    when no inference object exists.
  # ------------------------------------------------------------------------
  if (isTRUE(loading_heatmaps)) {
    for (cc in components) {
      ht <- plotVHeatmap_fn(
        object = fit,
        Y = Y,
        data = data,
        zdata = zdata,
        component = cc,
        select_by = "loading",
        direction = "absolute",
        n_features = features_top_n,
        order_by = svdx_report_args$heatmap_order_by %||% NULL,
        annotate = svdx_report_args$heatmap_annotate %||% NULL,
        fdata = rw_single$fdata
      )
      if (is.null(ht)) next
      file <- file.path(svdx_dir, sprintf("svdx_loading_heatmap_component_%02d.pdf", cc))
      file <- .supp_save_heatmap(ht, file)
      manifest <- .supp_manifest_add(
        manifest, fig_no, "SSANOVA-SVD", "svdx",
        components = as.character(cc), plot_type = "loading_heatmap",
        file = file,
        caption = paste0("Heatmap of features with the largest absolute loadings for SSANOVA-SVD component ", cc, ".")
      )
      fig_no <- fig_no + 1L
    }
  }

  if (isTRUE(loading_enrichment)) {
    loading_enr_raw <- extract_loading_enr(rw_single)
    loading_enr_obj <- if (!is.null(loading_enr_raw)) {
      normalize_enr(loading_enr_raw, context = "supplement loading enrichment", verbose = verbose)
    } else NULL

    if (!is.null(loading_enr_obj) && is.data.frame(loading_enr_obj$results)) {
      loading_enr_df <- as.data.frame(loading_enr_obj$results)
      if (!all(c("scope", "component", "padj") %in% names(loading_enr_df))) {
        stop("Normalized loading enrichment is missing one or more required columns: scope, component, padj.", call. = FALSE)
      }

      # Use the exact component names generated by the benchmark loading-
      # enrichment machinery from the stored SSANOVA-SVD benchmark V matrix.
      svdx_result <- .supp_method_result(rw_single, "svdx")
      if (is.null(svdx_result) || is.null(svdx_result$V)) {
        stop("loading_enrichment=TRUE but the stored SSANOVA-SVD benchmark result has no V matrix.", call. = FALSE)
      }
      loading_names <- loading_component_names(svdx_result$V)
      if (length(loading_names) < max(components)) {
        stop("Stored SSANOVA-SVD loading component names do not cover the requested components.", call. = FALSE)
      }

      fdr_tag <- format(enrichment_fdr, trim = TRUE, scientific = FALSE)
      loading_component_pass <- loading_enr_df[0, , drop = FALSE]

      for (cc in components) {
        comp_name <- loading_names[[cc]]
        d <- loading_enr_df[
          loading_enr_df$scope == "component" &
            loading_enr_df$component == comp_name &
            is.finite(loading_enr_df$padj) &
            loading_enr_df$padj < enrichment_fdr,
          , drop = FALSE
        ]
        loading_component_pass <- rbind(loading_component_pass, d)
        .supp_write_tsv_allow_empty(
          d,
          file.path(enr_dir, sprintf("svdx_loading_enrichment_component_%02d_fdr%s.tsv", cc, fdr_tag))
        )

        if (!nrow(d)) next
        obj <- structure(list(results = d), class = c("svdxEnrichment", "list"))
        p <- plot_enrichment_dotplot_fn(obj, top_n = enrichment_top_n) +
          ggplot2::theme(
            legend.text = ggplot2::element_text(size = 0.75 * text_size),
            legend.title = ggplot2::element_text(size = 0.75 * text_size)
          )
        file <- file.path(svdx_dir, sprintf("svdx_loading_enrichment_component_%02d.pdf", cc))
        file <- .supp_save_gg(p, file, n_panels = 1L, panel_width = 8, panel_height = 6, text_size = text_size)
        manifest <- .supp_manifest_add(
          manifest, fig_no, "SSANOVA-SVD", "svdx",
          components = as.character(cc), plot_type = "loading_enrichment",
          file = file,
          caption = paste0(
            "Enrichment based on features with the largest absolute loadings for SSANOVA-SVD component ", cc,
            "; only terms with adjusted P <= ", enrichment_fdr, " are shown."
          )
        )
        fig_no <- fig_no + 1L
      }

      .supp_write_tsv_allow_empty(
        loading_component_pass,
        file.path(enr_dir, paste0("svdx_loading_enrichment_component_all_fdr", fdr_tag, ".tsv"))
      )
    } else if (isTRUE(verbose)) {
      message("loading_enrichment=TRUE, but no normalized SSANOVA-SVD loading enrichment is stored; loading enrichment skipped.")
    }
  }

  # ------------------------------------------------------------------------
  # 6. SSANOVA term decompositions. plotUTerms() determines the actual fitted
  #    terms. To keep text readable, at most two term plots are placed on a
  #    page. The package's current diverging interaction scale is used as-is.
  # ------------------------------------------------------------------------
  for (cc in components) {
    term_obj <- do.call(
      plotUTerms_fn,
      c(
        list(fit = fit, data = data, component = cc, facet = process_facet, fill = "typical"),
        svdx_plot_args(svdx_report_args, "plotUTerms")
      )
    )
    if (is.null(term_obj)) next
    term_plots <- if (inherits(term_obj, "ggplot")) list(term_obj) else term_obj
    term_plots <- Filter(function(x) inherits(x, "ggplot"), term_plots)
    if (!length(term_plots)) next

    pages <- .supp_chunks(seq_along(term_plots), 2L)
    term_names <- names(term_plots)
    if (is.null(term_names)) term_names <- rep("", length(term_plots))

    for (page_i in seq_along(pages)) {
      idx <- pages[[page_i]]
      plots <- lapply(term_plots[idx], function(x) {
        x <- .supp_clean_plot(x, text_size = text_size, keep_title = FALSE)
        .supp_decomposition_colorbar(x)
      })
      has_continuous_fill <- any(vapply(
        plots,
        function(x) {
          sc <- x$scales$get_scales("fill")
          !is.null(sc) && inherits(sc, "ScaleContinuous")
        },
        logical(1)
      ))
      p <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
      file <- file.path(
        svdx_dir,
        sprintf("svdx_decomposition_component_%02d_page_%02d.pdf", cc, page_i)
      )
      file <- .supp_save_gg(
        p, file, n_panels = length(idx),
        panel_height = if (isTRUE(has_continuous_fill)) 6.25 else 5.5,
        text_size = text_size
      )
      panel_terms <- term_names[idx]
      panel_terms <- as.character(panel_terms)
      panel_detail <- ""
      if (length(idx) > 1L && all(!is.na(panel_terms) & nzchar(panel_terms))) {
        tags <- LETTERS[seq_along(panel_terms)]
        panel_detail <- paste0(" ", paste0("(", tags, ") ", panel_terms, ".", collapse = " "))
      }
      manifest <- .supp_manifest_add(
        manifest, fig_no, "SSANOVA-SVD", "svdx",
        components = as.character(cc), plot_type = "decomposition",
        file = file, page = page_i,
        caption = paste0("Decomposition of reconstructed SSANOVA-SVD component function ", cc, ".", panel_detail)
      )
      fig_no <- fig_no + 1L
    }
  }

  # ------------------------------------------------------------------------
  # 7. Comparator component galleries. Every returned U component is shown.
  #    Four component plots per page; one collected legend below.
  # ------------------------------------------------------------------------
  available_methods <- report_methods(rw_single, integrated = FALSE)
  if (is.null(methods)) {
    methods_use <- setdiff(available_methods, c("svdx", "svdx_adjusted"))
  } else {
    methods_use <- intersect(as.character(methods), available_methods)
  }

  if (isTRUE(verbose)) {
    message("Comparator methods available: ", paste(available_methods, collapse = ", "))
    message("Comparator methods selected: ", paste(methods_use, collapse = ", "))
  }

  for (method in methods_use) {
    method_display <- .supp_method_display_name(rw_single, method)
    # Use the same component-discovery helper as runSingleBenchmarkReports().
    # Do not assume a particular benchmark_results substructure here.
    method_components <- filter_method_components(
      rw_single,
      method = method,
      components = NULL,
      integrated = FALSE,
      dataset = NULL
    )
    if (!length(method_components)) {
      if (isTRUE(verbose)) message("Skipping comparator '", method, "': no plottable U components.")
      next
    }

    method_dir <- file.path(comparator_dir, .supp_safe_name(method))
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

    for (spec_i in seq_along(prepared_specs)) {
      spec <- prepared_specs[[spec_i]]
      xa <- spec_value(spec, "var", spec_value(spec, "x", NULL))
      if (is.null(xa) || !xa %in% colnames(data)) {
        if (isTRUE(verbose)) message("Skipping comparator spec for '", method, "': invalid x variable.")
        next
      }

      component_pages <- .supp_chunks(method_components, 4L)
      for (page_i in seq_along(component_pages)) {
        cc_page <- component_pages[[page_i]]
        plots <- lapply(cc_page, function(cc) {
          p <- tryCatch(
            u_plot_one(
              rw_single,
              method = method,
              component = cc,
              spec = spec,
              integrated = FALSE,
              color = color,
              group = group,
              line_group = line_group,
              verbose = verbose
            ),
            error = function(e) {
              if (isTRUE(verbose)) {
                message("Could not plot comparator '", method,
                        "' component ", cc, " against '", xa, "': ", e$message)
              }
              NULL
            }
          )
          p <- .supp_clean_plot(p, text_size = text_size, keep_title = FALSE)
          if (!is.null(p) && inherits(p, "ggplot")) {
            p <- p + ggplot2::ylab(paste0("Component ", cc))
            if (.supp_is_discrete_x(data, xa)) {
              hide_legend <- (!is.null(color) && identical(as.character(color), as.character(xa))) ||
                (!is.null(group) && identical(as.character(group), as.character(xa)))
              p <- p + ggplot2::labs(x = NULL) +
                .supp_discrete_axis_theme(text_size = text_size, hide_legend = hide_legend)
            }
          }
          p
        })
        keep <- !vapply(plots, is.null, logical(1))
        plots <- plots[keep]
        plotted_components <- cc_page[keep]
        if (!length(plots)) next

        discrete_x <- .supp_is_discrete_x(data, xa)
        redundant_legend <- isTRUE(discrete_x) && (
          (!is.null(color) && identical(as.character(color), as.character(xa))) ||
          (!is.null(group) && identical(as.character(group), as.character(xa)))
        )

        if (length(plots) == 1L) {
          # A single returned component should use the full page width rather
          # than being placed into a nominal two-column patchwork layout.
          p_out <- plots[[1]] +
            ggplot2::labs(title = method_display) +
            ggplot2::theme(
              legend.position = if (isTRUE(redundant_legend)) "none" else "bottom",
              plot.title = ggplot2::element_text(hjust = 0.5)
            )
          save_width <- if (isTRUE(discrete_x)) 12 else 10
          save_height <- if (isTRUE(discrete_x)) 7 else 6
        } else {
          p <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect") &
            ggplot2::theme(legend.position = if (isTRUE(redundant_legend)) "none" else "bottom")
          p_out <- p + patchwork::plot_annotation(
            title = method_display,
            tag_levels = "A"
          )
          save_width <- if (isTRUE(discrete_x)) 13 else 10
          save_height <- if (isTRUE(discrete_x)) 12 else 10
        }

        file <- file.path(
          method_dir,
          sprintf("%s_%s_components_page_%02d.pdf",
                  .supp_safe_name(method), .supp_safe_name(xa), page_i)
        )
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(
          filename = file,
          plot = p_out,
          width = save_width,
          height = save_height,
          units = "in",
          device = grDevices::pdf
        )
        file <- normalizePath(file, mustWork = FALSE)
        manifest <- .supp_manifest_add(
          manifest,
          fig_no,
          "Comparator components",
          method,
          components = paste(plotted_components, collapse = ","),
          plot_type = "component_gallery",
          file = file,
          page = page_i,
          caption = paste0("Components for ", method_display, ". ", .supp_panel_caption(plotted_components))
        )
        fig_no <- fig_no + 1L
      }
    }
  }

  manifest <- .supp_assign_labels(manifest, "supp_circadian")

  manifest_file <- file.path(outdir, "supplement_manifest.tsv")
  .supp_write_tsv_allow_empty(manifest, manifest_file)

  tex_file <- file.path(outdir, "supplement_figures.tex")
  tex_file <- .supp_write_tex(manifest, tex_file, note_number = note_number, dataset_title = dataset_title)
  pdf_file <- if (isTRUE(compile) && !is.null(tex_file)) .supp_compile_tex(tex_file, verbose = verbose) else NULL

  out <- list(
    outdir = normalizePath(outdir, mustWork = FALSE),
    start_figure = as.integer(start_figure),
    end_figure = if (nrow(manifest)) max(manifest$figure_number) else as.integer(start_figure) - 1L,
    next_figure = fig_no,
    manifest = manifest,
    manifest_file = manifest_file,
    tex_file = tex_file,
    pdf_file = pdf_file,
    note_number = as.character(note_number),
    dataset_title = as.character(dataset_title)
  )
  invisible(out)
}
