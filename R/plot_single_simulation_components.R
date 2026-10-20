# plot_single_simulation_components.R
# Visualize each complete reconstructed SVDX component over group x time x severity.

plot_single_simulation_components <- function(
    oneRun,
    fit = NULL,
    components = NULL,
    severity_probs = c(0.10, 0.50, 0.90),
    n_grid = 200
) {
  stopifnot(length(severity_probs) == 3L)

  if (is.null(fit)) {
    candidate <- oneRun$fits$factorization$svdx
    if (inherits(candidate, "SSANOVASVD")) {
      fit <- candidate
    } else if (!is.null(candidate$fit) && inherits(candidate$fit, "SSANOVASVD")) {
      fit <- candidate$fit
    } else {
      stop("Could not identify retained SVDX factorization fit; pass fit explicitly.")
    }
  }

  dat <- oneRun$sim$meta
  train_dat <- dat[oneRun$split$train, , drop = FALSE]
  if (!all(c("group", "t", "severity") %in% names(train_dat))) {
    stop("Expected group, t, and severity in oneRun$sim$meta.")
  }

  group_values <- if (is.factor(train_dat$group)) levels(droplevels(train_dat$group)) else sort(unique(train_dat$group))
  sev_values <- as.numeric(stats::quantile(train_dat$severity, probs = severity_probs, na.rm = TRUE, names = FALSE))
  sev_labels <- c("Low severity", "Median severity", "High severity")
  t_grid <- seq(min(train_dat$t, na.rm = TRUE), max(train_dat$t, na.rm = TRUE), length.out = n_grid)

  grid <- expand.grid(
    t = t_grid,
    group = group_values,
    severity_index = seq_along(sev_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$severity <- sev_values[grid$severity_index]
  grid$severity_level <- factor(sev_labels[grid$severity_index], levels = sev_labels)
  grid$severity_index <- NULL
  if (is.factor(train_dat$group)) grid$group <- factor(grid$group, levels = levels(train_dat$group))

  # SSANOVASVD stores the fitted component matrices under fitAtOptCV.
  # Columns of umat / vmat are the reconstructed components.
  if (is.null(fit$fitAtOptCV$umat) || is.null(fit$fitAtOptCV$vmat)) {
    stop("Expected SSANOVASVD fit$fitAtOptCV$umat and fit$fitAtOptCV$vmat.")
  }

  K <- ncol(fit$fitAtOptCV$vmat)

  if (ncol(fit$fitAtOptCV$umat) != K) {
    stop("fit$fitAtOptCV$umat and fit$fitAtOptCV$vmat have inconsistent component counts.")
  }

  if (is.null(components)) components <- seq_len(K)

  components <- as.integer(components)
  if (!length(components) || anyNA(components) ||
      any(components < 1L | components > K)) {
    stop("components must be valid indices in 1,...,K.")
  }

  # predict.SSANOVASVD() returns the COMPLETE predicted U matrix:
  # one row per newdata row and one column per reconstructed component.
  # Select the requested component column explicitly.
  predict_component <- function(k) {
    nd <- grid

    Uhat <- as.matrix(
      stats::predict(
        fit,
        newdata = nd
      )
    )

    if (nrow(Uhat) != nrow(nd)) {
      stop(
        "predict(fit, newdata) returned ", nrow(Uhat),
        " rows for ", nrow(nd), " prediction rows."
      )
    }

    if (k > ncol(Uhat)) {
      stop(
        "Requested component ", k,
        " but predict(fit, newdata) returned only ",
        ncol(Uhat), " component columns."
      )
    }

    nd$fitted <- Uhat[, k]
    nd$component <- paste0("CF", k)
    nd
  }

  pred_data <- do.call(rbind, lapply(components, predict_component))
  pred_data$component <- factor(pred_data$component, levels = paste0("CF", components))

  make_component_plot <- function(k) {
    dd <- pred_data[pred_data$component == paste0("CF", k), , drop = FALSE]
    ggplot2::ggplot(dd, ggplot2::aes(t, fitted, color = group, group = group)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::facet_wrap(~severity_level, nrow = 1) +
      ggplot2::labs(
        title = paste0("CF", k), x = "Time",
        y = "Reconstructed component function", color = "Group"
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "right")
  }

  plots <- lapply(components, make_component_plot)
  names(plots) <- paste0("CF", components)
  combined <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = 1) +
      patchwork::plot_annotation(title = "Reconstructed SVDX component functions")
  }

  invisible(list(
    combined = combined,
    plots = plots,
    data = pred_data,
    severity_values = stats::setNames(sev_values, sev_labels),
    fit = fit
  ))
}

# Example:
# source("plot_single_simulation_components.R")
# out <- plot_single_simulation_components(oneRun, components = 1:4)
# out$combined
# out$plots$CF1
# out$severity_values
#
# ggplot2::ggsave(
#   "single_simulation_reconstructed_components.pdf",
#   out$combined,
#   device = grDevices::cairo_pdf,
#   width = 14, height = 14, units = "in"
# )
