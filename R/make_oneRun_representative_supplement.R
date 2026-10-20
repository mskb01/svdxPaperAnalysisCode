# make_oneRun_representative_supplement.R
#
# Build the representative-simulation supplementary figures from a retained
# oneRun object
#
# The script includes the representative-simulation
# plotting and recovery helpers and combines their outputs into a coherent
# supplement:
#
#   Figure 1  True simulated biological processes
#   Figure 2  Complete reconstructed SVDX component functions (CF1-CF4)
#   Figure 3  Recovery diagnostics:
#               A-D: each true process reconstructed from the full CF basis
#               E:   held-out truth$biology prediction
#   Figures 4+ Term decomposition, one figure per recovered CF
#
# It also writes compact TSVs containing the quantitative recovery results.
#
#
#
# Example:
#
# source("make_oneRun_representative_supplement.R")
#
# supp <- make_oneRun_representative_supplement(
#   oneRun,
#   out_dir = "./results/simulations/representative_supplement",
#   note_number = "X",
#   dataset_title = "representative simulation",
#   start_figure = 1,
#   components = 1:4,
#   save = TRUE,
#   compile = FALSE
# )
#
# supp$figures$true_processes
# supp$figures$reconstructed_components
# supp$figures$recovery
# supp$figures$term_decomposition$CF1
#
# supp$process_recovery$summary
# supp$heldout_biology$summary




# =============================================================================
# representative-simulation helpers
# =============================================================================

# plot_true_simulation_processes.R
#
# Plot the four TRUE biological processes from a single benchmark simulation.
# Processes are identified by truth$module_type, never by module number, so
# this remains correct when process-to-module assignments change across runs.

plot_true_simulation_processes <- function(
    oneRun,
    severity_probs = c(1/3, 2/3),
    base_size = 13
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required.")

  truth <- oneRun$sim$truth
  meta  <- oneRun$sim$meta

  if (ncol(truth$module_signal) != nrow(meta)) {
    stop("ncol(truth$module_signal) must equal nrow(oneRun$sim$meta).")
  }
  if (length(truth$module_type) != nrow(truth$module_signal)) {
    stop("module_type and module_signal have inconsistent module counts.")
  }

  get_process <- function(type) {
    idx <- which(truth$module_type == type)
    if (length(idx) != 1L) {
      stop("Expected exactly one module of type '", type,
           "'; found ", length(idx), ".")
    }

    out <- meta
    out$.process <- as.numeric(truth$module_signal[idx, ])
    out$.module  <- idx
    out
  }

  d_group  <- get_process("group_mean")
  d_shared <- get_process("shared_oscillation")
  d_sev    <- get_process("severity_modulated_oscillation")
  d_phase  <- get_process("group_phase_oscillation")

  gt <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  # A: group-specific mean process
  group_summary <- stats::aggregate(
    .process ~ group,
    data = d_group,
    FUN = mean
  )

  p_group <- ggplot2::ggplot(
    group_summary,
    ggplot2::aes(x = factor(group), y = .process, group = 1)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = "Group",
      x = "Group",
      y = "True process"
    ) +
    gt

  # B: temporal oscillation shared across groups
  # Repeated t values can occur across groups; collapse exact/near-exact t
  # values only for a clean line representation of the process.
  shared_order <- order(d_shared$t)
  d_shared <- d_shared[shared_order, , drop = FALSE]

  p_shared <- ggplot2::ggplot(
    d_shared,
    ggplot2::aes(x = t, y = .process)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Oscillation: shared",
      x = "Time",
      y = "True process"
    ) +
    gt

  # C: severity-modulated temporal oscillation.
  # Keep the exact realized process values; severity bins are only a visual
  # device to reveal the changing oscillation amplitude.
  sev_breaks <- as.numeric(stats::quantile(
    d_sev$severity,
    probs = severity_probs,
    na.rm = TRUE,
    names = FALSE
  ))

  d_sev$.severity_level <- cut(
    d_sev$severity,
    breaks = c(-Inf, sev_breaks, Inf),
    labels = c("Low", "Middle", "High"),
    include.lowest = TRUE
  )

  p_sev <- ggplot2::ggplot(
    d_sev,
    ggplot2::aes(
      x = t,
      y = .process,
      color = .severity_level
    )
  ) +
    ggplot2::geom_point(alpha = 0.35, size = 1) +
    ggplot2::geom_smooth(
      se = FALSE,
      linewidth = 1,
      method = "loess",
      formula = y ~ x
    ) +
    ggplot2::labs(
      title = "Oscillation: severity-amplified",
      x = "Time",
      y = "True process",
      color = "Severity"
    ) +
    gt

  # D: temporal oscillation with group-specific phase shifts
  d_phase <- d_phase[order(d_phase$group, d_phase$t), , drop = FALSE]

  p_phase <- ggplot2::ggplot(
    d_phase,
    ggplot2::aes(
      x = t,
      y = .process,
      color = factor(group),
      group = factor(group)
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Oscillation: group-specific phases",
      x = "Time",
      y = "True process",
      color = "Group"
    ) +
    gt

  plots <- list(
    group_mean = p_group,
    shared_oscillation = p_shared,
    severity_modulated_oscillation = p_sev,
    group_phase_oscillation = p_phase
  )

  combined <- patchwork::wrap_plots(plots, ncol = 2) +
    patchwork::plot_annotation(
      title = "Simulated biological processes"
    )

  invisible(list(
    combined = combined,
    plots = plots,
    module_mapping = data.frame(
      module = seq_along(truth$module_type),
      process = truth$module_type
    )
  ))
}


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


# check_true_process_recovery.R
#
# Quantify whether the four simulated biological processes are contained in
# the subspace spanned by the reconstructed SVDX component functions.
#
# The key quantity is process-wise R^2 after projecting each TRUE process onto
# the full recovered CF basis. This is invariant to rotations/mixing of the CFs.
#
# Usage:
#   source("check_true_process_recovery.R")
#   rec <- check_true_process_recovery(oneRun)
#   rec$summary
#   rec$coefficients
#   rec$overall_variance_explained
#   rec$plots$combined

check_true_process_recovery <- function(
    oneRun,
    fit = NULL,
    components = 1:4,
    make_plots = TRUE
) {
  truth <- oneRun$sim$truth

  if (is.null(fit)) {
    fit <- oneRun$fits$factorization$svdx

    # Allow for benchmark wrappers, but prefer the SSANOVASVD object directly.
    if (is.null(fit$fitAtOptCV$umat) &&
        !is.null(fit$fit) &&
        !is.null(fit$fit$fitAtOptCV$umat)) {
      fit <- fit$fit
    }
  }

  if (is.null(fit$fitAtOptCV$umat)) {
    stop("Expected an SSANOVASVD object with fit$fitAtOptCV$umat.")
  }

  U <- as.matrix(fit$fitAtOptCV$umat)

  components <- as.integer(components)
  if (anyNA(components) ||
      any(components < 1L) ||
      any(components > ncol(U))) {
    stop("Invalid component indices.")
  }

  U <- U[, components, drop = FALSE]
  colnames(U) <- paste0("CF", components)

  P_all <- t(as.matrix(truth$module_signal))
  colnames(P_all) <- truth$module_type

  # Align truth to the rows represented by umat.
  train_idx <- oneRun$split$train

  if (!is.null(train_idx) && nrow(U) == length(train_idx)) {
    P <- P_all[train_idx, , drop = FALSE]
    sample_index <- train_idx
    sample_scope <- "training"
  } else if (nrow(U) == nrow(P_all)) {
    P <- P_all
    sample_index <- seq_len(nrow(P_all))
    sample_scope <- "all"
  } else {
    stop(
      "Could not align fit$fitAtOptCV$umat (", nrow(U),
      " rows) with module_signal (", nrow(P_all),
      " samples) or the training split (",
      if (is.null(train_idx)) "NULL" else length(train_idx), " samples)."
    )
  }

  # Center both spaces. This tests whether variation in each true process lies
  # in the recovered CF subspace, independent of arbitrary mean offsets.
  Uc <- scale(U, center = TRUE, scale = FALSE)
  Pc <- scale(P, center = TRUE, scale = FALSE)

  # Orthogonal projection onto span(Uc). QR is numerically stable and makes
  # no assumption that individual CFs correspond one-to-one with processes.
  qrU <- qr(Uc)
  rankU <- qrU$rank
  Q <- qr.Q(qrU)[, seq_len(rankU), drop = FALSE]

  P_hat_c <- Q %*% (t(Q) %*% Pc)

  sst <- colSums(Pc^2)
  sse <- colSums((Pc - P_hat_c)^2)

  r2 <- 1 - sse / sst
  r2[sst == 0] <- NA_real_

  # Recover original process means for fitted-value plots.
  process_means <- colMeans(P)
  P_hat <- sweep(P_hat_c, 2, process_means, "+")

  # Coefficients expressing each centered true process as a combination of CFs.
  # Columns correspond to true processes; rows correspond to CFs.
  coef_mat <- qr.coef(qr(Uc), Pc)
  rownames(coef_mat) <- colnames(U)
  colnames(coef_mat) <- colnames(P)

  # Overall fraction of centered true-process variation explained by the
  # recovered CF subspace.
  overall_variance_explained <-
    1 - sum((Pc - P_hat_c)^2) / sum(Pc^2)

  summary <- data.frame(
    process = colnames(P),
    R2_from_all_CFs = as.numeric(r2),
    truth_SD = apply(P, 2, stats::sd),
    reconstruction_RMSE = sqrt(colMeans((P - P_hat)^2)),
    stringsAsFactors = FALSE
  )

  fitted_long <- do.call(
    rbind,
    lapply(seq_len(ncol(P)), function(j) {
      data.frame(
        sample_index = sample_index,
        process = colnames(P)[j],
        truth = P[, j],
        reconstructed = P_hat[, j],
        residual = P[, j] - P_hat[, j],
        stringsAsFactors = FALSE
      )
    })
  )

  plots <- NULL

  if (isTRUE(make_plots)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not installed; skipping plots.")
    } else {
      p_list <- lapply(colnames(P), function(proc) {
        dd <- fitted_long[fitted_long$process == proc, , drop = FALSE]
        rr <- summary$R2_from_all_CFs[summary$process == proc]

        ggplot2::ggplot(
          dd,
          ggplot2::aes(x = truth, y = reconstructed)
        ) +
          ggplot2::geom_point(alpha = 0.45, size = 1.4) +
          ggplot2::geom_abline(
            slope = 1,
            intercept = 0,
            linetype = 2,
            linewidth = 0.7
          ) +
          ggplot2::coord_equal() +
          ggplot2::labs(
            title = proc,
            subtitle = sprintf("R² = %.3f", rr),
            x = "True process",
            y = "Reconstructed from all CFs"
          ) +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
          )
      })

      names(p_list) <- colnames(P)

      combined <- NULL
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- patchwork::wrap_plots(p_list, ncol = 2) +
          patchwork::plot_annotation(
            title = "Recovery of simulated processes from the full CF basis"
          )
      }

      plots <- list(
        individual = p_list,
        combined = combined
      )
    }
  }

  invisible(list(
    summary = summary,
    coefficients = coef_mat,
    overall_variance_explained = overall_variance_explained,
    recovered_rank = rankU,
    sample_scope = sample_scope,
    truth = P,
    reconstructed = P_hat,
    fitted_long = fitted_long,
    plots = plots
  ))
}


# check_heldout_biology_recovery.R
#
# Quantify recovery of the TRUE biological matrix in completely held-out
# samples. For Gaussian benchmark runs, truth$biology has already been
# centered/scaled onto the same feature-wise standardized scale used for fit.
#
# Biological prediction only:
#   Y_bio,test_hat = V %*% t(U_test_hat)
#
# The nuisance term Zeta Z^T is deliberately excluded, and comparison is
# against truth$biology rather than noisy Y.

check_heldout_biology_recovery <- function(
    oneRun,
    fit = NULL,
    components = NULL,
    make_plot = TRUE
) {
  if (is.null(oneRun$sim$truth$biology)) stop("oneRun$sim$truth$biology is missing.")
  if (is.null(oneRun$split$test)) stop("oneRun$split$test is missing.")
  if (is.null(oneRun$sim$meta)) stop("oneRun$sim$meta is missing.")

  if (is.null(fit)) {
    fit <- oneRun$fits$factorization$svdx
    if (is.null(fit$fitAtOptCV$vmat) &&
        !is.null(fit$fit) &&
        !is.null(fit$fit$fitAtOptCV$vmat)) {
      fit <- fit$fit
    }
  }

  if (is.null(fit$fitAtOptCV$vmat)) {
    stop("Expected an SSANOVASVD object with fit$fitAtOptCV$vmat.")
  }

  V <- as.matrix(fit$fitAtOptCV$vmat)
  K <- ncol(V)
  if (is.null(components)) components <- seq_len(K)
  components <- as.integer(components)

  if (anyNA(components) || any(components < 1L) || any(components > K)) {
    stop("Invalid component indices.")
  }

  test_idx <- oneRun$split$test
  x_test <- oneRun$sim$meta[test_idx, , drop = FALSE]

  # predict.SSANOVASVD returns all reconstructed component functions.
  U_test <- as.matrix(stats::predict(fit, newdata = as.data.frame(x_test)))

  if (nrow(U_test) != length(test_idx)) {
    stop("Predicted U rows do not match the test split.")
  }
  if (ncol(U_test) < max(components)) {
    stop("Predicted U has fewer columns than requested components.")
  }

  U_test <- U_test[, components, drop = FALSE]
  V_use <- V[, components, drop = FALSE]

  biology_hat <- V_use %*% t(U_test)
  biology_truth <- as.matrix(oneRun$sim$truth$biology[, test_idx, drop = FALSE])

  if (!identical(dim(biology_hat), dim(biology_truth))) {
    stop("Predicted and true biological matrices have different dimensions.")
  }

  resid <- biology_truth - biology_hat
  sse <- sum(resid^2, na.rm = TRUE)
  truth_mean <- mean(biology_truth, na.rm = TRUE)
  sst <- sum((biology_truth - truth_mean)^2, na.rm = TRUE)
  r2 <- if (is.finite(sst) && sst > 0) 1 - sse / sst else NA_real_

  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  truth_sd <- stats::sd(as.numeric(biology_truth), na.rm = TRUE)
  correlation <- suppressWarnings(stats::cor(
    as.numeric(biology_truth),
    as.numeric(biology_hat),
    use = "complete.obs"
  ))
  relative_frobenius_error <- sqrt(sse) / sqrt(sum(biology_truth^2, na.rm = TRUE))

  signal_idx <- oneRun$sim$truth$signal_features
  signal_idx <- signal_idx[
    is.finite(signal_idx) & signal_idx >= 1L & signal_idx <= nrow(biology_truth)
  ]

  signal_r2 <- NA_real_
  signal_rmse <- NA_real_
  if (length(signal_idx) > 0L) {
    bt_s <- biology_truth[signal_idx, , drop = FALSE]
    bh_s <- biology_hat[signal_idx, , drop = FALSE]
    sse_s <- sum((bt_s - bh_s)^2, na.rm = TRUE)
    mean_s <- mean(bt_s, na.rm = TRUE)
    sst_s <- sum((bt_s - mean_s)^2, na.rm = TRUE)
    signal_r2 <- if (is.finite(sst_s) && sst_s > 0) 1 - sse_s / sst_s else NA_real_
    signal_rmse <- sqrt(mean((bt_s - bh_s)^2, na.rm = TRUE))
  }

  summary <- data.frame(
    n_test = length(test_idx),
    n_features = nrow(biology_truth),
    n_components = length(components),
    heldout_biology_R2 = r2,
    heldout_biology_correlation = correlation,
    heldout_biology_RMSE = rmse,
    truth_biology_SD = truth_sd,
    relative_frobenius_error = relative_frobenius_error,
    signal_features_R2 = signal_r2,
    signal_features_RMSE = signal_rmse,
    stringsAsFactors = FALSE
  )

  p <- NULL
  if (isTRUE(make_plot)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 is not installed; skipping plot.")
    } else {
      dd <- data.frame(
        truth = as.numeric(biology_truth),
        predicted = as.numeric(biology_hat)
      )
      max_plot_points <- 20000L
      if (nrow(dd) > max_plot_points) {
        keep <- unique(round(seq(1, nrow(dd), length.out = max_plot_points)))
        dd <- dd[keep, , drop = FALSE]
      }

      p <- ggplot2::ggplot(dd, ggplot2::aes(x = truth, y = predicted)) +
        ggplot2::geom_point(alpha = 0.15, size = 0.7) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.7) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = "Held-out biological-signal recovery",
          subtitle = sprintf("Entire held-out truth matrix: R² = %.4f", r2),
          x = "True biological signal",
          y = "Predicted biological signal"
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    }
  }

  invisible(list(
    summary = summary,
    truth = biology_truth,
    predicted = biology_hat,
    residual = resid,
    U_test = U_test,
    V = V_use,
    plot = p
  ))
}

# =============================================================================
# Supplement assembly utilities
# =============================================================================

.clean_process_label <- function(x) {
  out <- c(
    group_mean = "Group",
    shared_oscillation = "Oscillation: shared",
    severity_modulated_oscillation = "Oscillation: severity-amplified",
    group_phase_oscillation = "Oscillation: group-specific phases"
  )

  ans <- unname(out[x])
  ans[is.na(ans)] <- x[is.na(ans)]
  ans
}


.wrap_rep_title <- function(x, width = 28) {
  paste(strwrap(x, width = width), collapse = "\n")
}


.rep_write_tsv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  normalizePath(file, mustWork = FALSE)
}


.rep_escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}


.rep_manifest_add <- function(
    manifest,
    figure_number,
    section,
    components = NA_character_,
    plot_type,
    file,
    caption
) {
  row <- data.frame(
    figure_number = as.integer(figure_number),
    section = as.character(section),
    method = "svdx",
    components = as.character(components),
    plot_type = as.character(plot_type),
    page = 1L,
    file = as.character(file),
    caption = as.character(caption),
    stringsAsFactors = FALSE
  )
  rbind(manifest, row)
}


.rep_label_slug <- function(x) {
  x <- tolower(as.character(x))
  x[is.na(x)] <- ""
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("^-+|-+$", "", x)
  x
}

.rep_assign_labels <- function(manifest, figure_label_prefix = "supp_onerun") {
  if (!nrow(manifest)) {
    manifest$label <- character(0)
    return(manifest)
  }
  suffix <- vapply(seq_len(nrow(manifest)), function(i) {
    parts <- .rep_label_slug(c(manifest$plot_type[i], manifest$components[i]))
    parts <- parts[nzchar(parts)]
    if (!length(parts)) "figure" else paste(parts, collapse = "-")
  }, character(1))
  suffix <- make.unique(suffix, sep = "-")
  manifest$label <- paste0("fig:", figure_label_prefix, "-", suffix)
  manifest
}

.rep_write_tex <- function(
    manifest,
    tex_file,
    note_number,
    dataset_title
) {
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
    paste0(
      "\\section*{Supplementary Note ",
      .rep_escape_tex(note_number),
      "}"
    ),
    paste0(
      "This note presents the supplementary figures for ",
      .rep_escape_tex(dataset_title),
      "."
    ),
    ""
  )

  for (i in seq_len(nrow(manifest))) {
    r <- manifest[i, , drop = FALSE]
    cap <- .rep_escape_tex(r$caption)

    lines <- c(
      lines,
      "\\begin{figure}[p]",
      "\\centering",
      paste0(
        "\\includegraphics[width=\\textwidth,height=0.82\\textheight,keepaspectratio]{",
        rel(r$file),
        "}"
      ),
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


.rep_compile_tex <- function(tex_file, verbose = TRUE) {
  latexmk <- Sys.which("latexmk")
  if (!nzchar(latexmk)) {
    if (isTRUE(verbose)) {
      message("latexmk not found; wrote TeX but did not compile it.")
    }
    return(NULL)
  }

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(dirname(tex_file))

  status <- system2(
    latexmk,
    args = c(
      "-pdf",
      "-interaction=nonstopmode",
      "-halt-on-error",
      basename(tex_file)
    ),
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


.rep_theme <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 2)
    )
}


make_oneRun_representative_supplement <- function(
    oneRun,
    out_dir = "./results/simulations/representative_supplement",
    note_number = "X",
    dataset_title = "representative simulation",
    start_figure = 1L,
    components = 1:4,
    severity_probs = c(0.10, 0.50, 0.90),
    n_grid = 200,
    base_size = 13,
    save = TRUE,
    compile = FALSE,
    verbose = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.", call. = FALSE)
  }
  if (!length(out_dir) || !nzchar(out_dir)) {
    stop("out_dir must be supplied.", call. = FALSE)
  }
  if (length(note_number) != 1L || is.na(note_number) ||
      !nzchar(as.character(note_number))) {
    stop("note_number must be a single non-empty value.", call. = FALSE)
  }
  if (length(dataset_title) != 1L || is.na(dataset_title) ||
      !nzchar(as.character(dataset_title))) {
    stop("dataset_title must be a single non-empty value.", call. = FALSE)
  }
  start_figure <- as.integer(start_figure)
  if (length(start_figure) != 1L || is.na(start_figure) || start_figure < 1L) {
    stop("start_figure must be a positive integer.", call. = FALSE)
  }

  components <- as.integer(components)

  # Resolve the retained SSANOVA-SVD fit once.
  fit <- oneRun$fits$factorization$svdx
  if (is.null(fit$fitAtOptCV$vmat) &&
      !is.null(fit$fit) &&
      !is.null(fit$fit$fitAtOptCV$vmat)) {
    fit <- fit$fit
  }
  if (is.null(fit$fitAtOptCV$vmat) || is.null(fit$fitAtOptCV$umat)) {
    stop(
      "Could not resolve retained SSANOVA-SVD fit with fitAtOptCV$umat/vmat.",
      call. = FALSE
    )
  }

  if (anyNA(components) ||
      any(components < 1L) ||
      any(components > ncol(fit$fitAtOptCV$vmat))) {
    stop("Invalid requested components.", call. = FALSE)
  }

  train_dat <- oneRun$sim$meta[oneRun$split$train, , drop = FALSE]

  # ---------------------------------------------------------------------------
  # 1. TRUE generating processes
  # ---------------------------------------------------------------------------

  truth_out <- plot_true_simulation_processes(
    oneRun = oneRun,
    base_size = base_size
  )

  truth_plots <- truth_out$plots
  truth_names <- names(truth_plots)

  truth_plots <- lapply(seq_along(truth_plots), function(i) {
    p <- truth_plots[[i]]
    key <- truth_names[i]

    p +
      ggplot2::labs(
        title = .clean_process_label(key),
        tag = LETTERS[i]
      ) +
      ggplot2::theme(
        plot.tag = ggplot2::element_text(
          face = "bold",
          size = base_size + 2
        )
      )
  })
  names(truth_plots) <- truth_names

  # Replace the group-mean line/point display with filled bars, matching the
  # group colors used elsewhere.  Suppress its redundant legend because the
  # group-specific phase panel already identifies the group colors.
  truth <- oneRun$sim$truth
  meta <- oneRun$sim$meta
  group_idx <- which(truth$module_type == "group_mean")
  if (length(group_idx) != 1L) {
    stop("Expected exactly one group_mean module.", call. = FALSE)
  }

  group_df <- data.frame(
    group = factor(meta$group),
    process = as.numeric(truth$module_signal[group_idx, ])
  )
  group_df <- stats::aggregate(process ~ group, data = group_df, FUN = mean)

  truth_plots$group_mean <- ggplot2::ggplot(
    group_df,
    ggplot2::aes(x = group, y = process, fill = group)
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::labs(
      title = "Group",
      x = "Group",
      y = "True process",
      tag = "A"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = base_size + 2
      )
    )

  fig_true <- patchwork::wrap_plots(truth_plots, ncol = 2)

  # ---------------------------------------------------------------------------
  # 2. Complete reconstructed CF functions
  # ---------------------------------------------------------------------------

  component_out <- plot_single_simulation_components(
    oneRun = oneRun,
    fit = fit,
    components = components,
    severity_probs = severity_probs,
    n_grid = n_grid
  )

  component_plots <- component_out$plots

  component_plots <- lapply(seq_along(component_plots), function(i) {
    k <- components[i]
    component_plots[[i]] +
      ggplot2::labs(
        title = NULL,
        y = paste0("Component Function ", k),
        tag = LETTERS[i]
      ) +
      ggplot2::theme(
        text = ggplot2::element_text(size = base_size + 5),
        axis.title = ggplot2::element_text(size = base_size + 5),
        axis.text = ggplot2::element_text(size = base_size + 3),
        strip.text = ggplot2::element_text(
          size = base_size + 4,
          face = "bold"
        ),
        legend.title = ggplot2::element_text(size = base_size + 4),
        legend.text = ggplot2::element_text(size = base_size + 3),
        aspect.ratio = 1,
        plot.tag = ggplot2::element_text(
          face = "bold",
          size = base_size + 5
        )
      )
  })
  names(component_plots) <- paste0("CF", components)

  fig_components <- patchwork::wrap_plots(
    component_plots,
    ncol = 1,
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )

  # ---------------------------------------------------------------------------
  # 3. Process-space recovery: P approximately U A
  # ---------------------------------------------------------------------------

  process_out <- check_true_process_recovery(
    oneRun = oneRun,
    fit = fit,
    components = components,
    make_plots = TRUE
  )

  process_plots <- process_out$plots$individual

  # Put the generating-process names into reader-facing form and preserve the
  # process-wise R^2 already computed by the helper.
  process_plots <- lapply(seq_along(process_plots), function(i) {
    key <- names(process_plots)[i]
    rr <- process_out$summary$R2_from_all_CFs[
      match(key, process_out$summary$process)
    ]

    process_plots[[i]] +
      ggplot2::labs(
        title = .wrap_rep_title(.clean_process_label(key), width = 27),
        subtitle = sprintf("R² = %.4f", rr),
        tag = LETTERS[i]
      ) +
      .rep_theme(base_size = base_size + 4) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = base_size + 4),
        axis.text = ggplot2::element_text(size = base_size + 2),
        plot.title = ggplot2::element_text(size = base_size + 5, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = base_size + 3),
        plot.tag = ggplot2::element_text(size = base_size + 5, face = "bold")
      )
  })
  names(process_plots) <- names(process_out$plots$individual)

  process_grid <- patchwork::wrap_plots(
    process_plots,
    ncol = 2
  )

  # ---------------------------------------------------------------------------
  # 4. Held-out recovery of the TRUE biological matrix
  # ---------------------------------------------------------------------------

  heldout_out <- check_heldout_biology_recovery(
    oneRun = oneRun,
    fit = fit,
    components = components,
    make_plot = TRUE
  )

  hs <- heldout_out$summary

  heldout_plot <- heldout_out$plot +
    ggplot2::labs(
      title = .wrap_rep_title(
        "Held-out biological-signal recovery",
        width = 27
      ),
      subtitle = sprintf(
        "All features: R² = %.4f\nSignal features: R² = %.4f",
        hs$heldout_biology_R2,
        hs$signal_features_R2
      ),
      tag = "E"
    ) +
    .rep_theme(base_size = base_size + 4) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = base_size + 4),
      axis.text = ggplot2::element_text(size = base_size + 2),
      plot.title = ggplot2::element_text(size = base_size + 5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size + 3),
      plot.tag = ggplot2::element_text(size = base_size + 5, face = "bold")
    )

  # A-D show P ~= U A; E tests the full low-rank biological reconstruction on
  # completely unseen samples.  Give the four process panels more width than
  # the single held-out matrix panel.
  fig_recovery <- (
    process_grid |
      heldout_plot
  ) +
    patchwork::plot_layout(widths = c(2, 1))

  # ---------------------------------------------------------------------------
  # 5. Term decomposition, one figure per recovered CF.
  #
  # One figure per CF is intentional: each decomposition already contains five
  # panels, and forcing all four CFs into one page makes the interaction
  # surfaces too small to interpret.
  # ---------------------------------------------------------------------------

  make_term_decomposition <- function(k) {
    enforce_min_effect_ylim <- function(p, min_abs = 1) {
      built <- ggplot2::ggplot_build(p)

      vals <- unlist(
        lapply(
          built$data,
          function(d) {
            cols <- intersect(c("y", "ymin", "ymax", "yend"), names(d))
            if (!length(cols)) return(numeric(0))
            unlist(d[cols], use.names = FALSE)
          }
        ),
        use.names = FALSE
      )
      vals <- vals[is.finite(vals)]

      if (!length(vals)) {
        lim <- c(-min_abs, min_abs)
      } else {
        lim <- c(
          min(-min_abs, min(vals)),
          max(min_abs, max(vals))
        )
      }

      p + ggplot2::coord_cartesian(ylim = lim)
    }

    term_names <- c(
      "group",
      "t",
      "severity",
      "group:t",
      "t:severity"
    )
    term_titles <- c(
      "Group",
      "Time",
      "Severity",
      "Group × time",
      "Time × severity"
    )

    pp <- lapply(seq_along(term_names), function(i) {
      p <- mySvdx::plotUTerms(
        fit = fit,
        data = train_dat,
        component = k,
        terms = term_names[i],
        fill = "typical",
        include_intercept = FALSE
      )

      if (is.list(p) && !inherits(p, "ggplot") && length(p) == 1L) {
        p <- p[[1L]]
      }

      # For the categorical group main effect, use filled bars with the same
      # default group palette as the group-specific interaction curves.  The
      # bar legend is suppressed because the group × time panel supplies the
      # shared group legend.
      if (identical(term_names[i], "group")) {
        d <- p$data
        p <- ggplot2::ggplot(
          d,
          ggplot2::aes(
            x = factor(group),
            y = contribution,
            fill = factor(group)
          )
        ) +
          ggplot2::geom_col(show.legend = FALSE) +
          ggplot2::theme_bw() +
          ggplot2::labs(
            x = "group",
            y = "Fitted effect function"
          ) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(
              angle = 90,
              hjust = 1,
              vjust = 0.5
            )
          )
      }

      p <- p + ggplot2::labs(title = term_titles[i])

      # Fix the two actual guides at the panel level.  This prevents ggplot2
      # from switching them to horizontal simply because the collected guide
      # box is later placed in the unused sixth cell.
      if (identical(term_names[i], "group:t")) {
        p <- p +
          ggplot2::guides(
            colour = ggplot2::guide_legend(
              direction = "vertical",
              ncol = 1
            ),
            color = ggplot2::guide_legend(
              direction = "vertical",
              ncol = 1
            )
          )
      }

      if (identical(term_names[i], "t:severity")) {
        p <- p +
          ggplot2::guides(
            fill = ggplot2::guide_colourbar(
              direction = "vertical",
              title.position = "top",
              title.hjust = 0.5,
              barwidth = grid::unit(4, "mm"),
              barheight = grid::unit(34, "mm")
            )
          )
      }

      if (!identical(term_names[i], "t:severity")) {
        p <- enforce_min_effect_ylim(p, min_abs = 1)
      }

      p
    })

    patchwork::wrap_plots(
      c(pp, list(patchwork::guide_area())),
      ncol = 3,
      guides = "collect"
    ) +
      patchwork::plot_annotation(
        title = paste0("CF", k, ": Term decomposition")
      ) &
      ggplot2::theme(
        text = ggplot2::element_text(size = base_size),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "right",
        legend.box = "horizontal",
        legend.box.just = "left"
      )
  }

  term_plots <- lapply(components, make_term_decomposition)
  names(term_plots) <- paste0("CF", components)

  # ---------------------------------------------------------------------------
  # Save figures, quantitative tables, manifest, and LaTeX supplementary note.
  # Figure order is fixed to mirror the scientific narrative:
  #   1. true processes
  #   2. reconstructed component functions
  #   3. quantitative recovery
  #   4+. term decompositions, CF1 through CFK
  # ---------------------------------------------------------------------------

  files <- list()
  manifest <- data.frame(
    figure_number = integer(0),
    section = character(0),
    method = character(0),
    components = character(0),
    plot_type = character(0),
    page = integer(0),
    file = character(0),
    caption = character(0),
    stringsAsFactors = FALSE
  )
  fig_no <- start_figure
  manifest_file <- NULL
  tex_file <- NULL
  pdf_file <- NULL

  if (isTRUE(save)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fig_dir <- file.path(out_dir, "figures")
    table_dir <- file.path(out_dir, "tables")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

    save_one <- function(plot, stem, width, height) {
      file <- file.path(
        fig_dir,
        sprintf(
          "supplementary_figure_%03d_%s.pdf",
          fig_no,
          stem
        )
      )
      ggplot2::ggsave(
        filename = file,
        plot = plot,
        device = grDevices::cairo_pdf,
        width = width,
        height = height,
        units = "in"
      )
      normalizePath(file, mustWork = FALSE)
    }

    # Supplementary Figure 1: true generating processes.
    files$true_processes <- save_one(
      fig_true,
      "true_processes",
      width = 12,
      height = 8
    )
    manifest <- .rep_manifest_add(
      manifest,
      fig_no,
      section = "Representative simulation",
      components = NA_character_,
      plot_type = "true_processes",
      file = files$true_processes,
      caption = paste0(
        "True biological processes used in the representative simulation. ",
        "(A) Group-specific mean shifts. ",
        "(B) An oscillation shared across groups. ",
        "(C) A severity-amplified oscillation whose amplitude varies continuously ",
        "with sample-specific severity. For visualization, samples are divided ",
        "into low, middle, and high severity strata and smoothed within each ",
        "stratum; these three curves summarize a continuous severity effect and ",
        "do not represent three discrete simulated severity states. ",
        "(D) Oscillations with group-specific phase shifts."
      )
    )
    fig_no <- fig_no + 1L

    # Supplementary Figure 2: reconstructed CFs.
    files$reconstructed_components <- save_one(
      fig_components,
      "component_functions",
      width = 14,
      height = 14
    )
    manifest <- .rep_manifest_add(
      manifest,
      fig_no,
      section = "Representative simulation",
      components = paste(components, collapse = ","),
      plot_type = "component_functions",
      file = files$reconstructed_components,
      caption = paste0(
        "Reconstructed SSANOVA-SVD component functions for the representative ",
        "simulation. Component functions are shown across the simulated design ",
        "for low, median, and high severity values."
      )
    )
    fig_no <- fig_no + 1L

    # Supplementary Figure 3: process-space and held-out biological recovery.
    files$recovery <- save_one(
      fig_recovery,
      "recovery",
      width = 12,
      height = 10
    )
    manifest <- .rep_manifest_add(
      manifest,
      fig_no,
      section = "Representative simulation",
      components = paste(components, collapse = ","),
      plot_type = "recovery",
      file = files$recovery,
      caption = paste0(
        "Recovery of the simulated biological structure. ",
        "Panels A-D compare each true generating process with its reconstruction ",
        "from the full learned component-function basis; the recovered component ",
        "space explains ",
        sprintf("%.4f", process_out$overall_variance_explained),
        " of the combined true-process variation. Panel E compares the true and ",
        "predicted biological signal in held-out samples."
      )
    )
    fig_no <- fig_no + 1L

    # Supplementary Figures 4+: one decomposition page per CF.
    files$term_decomposition <- vector("list", length(components))
    names(files$term_decomposition) <- paste0("CF", components)

    for (k in components) {
      file <- save_one(
        term_plots[[paste0("CF", k)]],
        sprintf("CF%d_term_decomposition", k),
        width = 12,
        height = 8
      )
      files$term_decomposition[[paste0("CF", k)]] <- file

      manifest <- .rep_manifest_add(
        manifest,
        fig_no,
        section = "Representative simulation",
        components = as.character(k),
        plot_type = "term_decomposition",
        file = file,
        caption = paste0(
          "Term decomposition of reconstructed component function ",
          k,
          " into group, time, severity, group-by-time, and time-by-severity ",
          "contributions."
        )
      )
      fig_no <- fig_no + 1L
    }

    # Compact quantitative tables retained for audit and later writing.
    process_summary <- process_out$summary
    process_summary$process_label <- .clean_process_label(
      process_summary$process
    )
    process_summary$overall_process_space_R2 <-
      process_out$overall_variance_explained

    files$process_recovery_table <- .rep_write_tsv(
      process_summary,
      file.path(table_dir, "representative_process_recovery.tsv")
    )

    coef_df <- data.frame(
      CF = rownames(process_out$coefficients),
      process_out$coefficients,
      check.names = FALSE,
      row.names = NULL
    )
    files$process_mixing_coefficients <- .rep_write_tsv(
      coef_df,
      file.path(table_dir, "representative_process_mixing_coefficients.tsv")
    )

    files$heldout_biology_table <- .rep_write_tsv(
      heldout_out$summary,
      file.path(table_dir, "representative_heldout_biology_recovery.tsv")
    )

    files$module_mapping <- .rep_write_tsv(
      truth_out$module_mapping,
      file.path(table_dir, "representative_module_mapping.tsv")
    )

    manifest <- .rep_assign_labels(manifest, "supp_onerun")

    manifest_file <- .rep_write_tsv(
      manifest,
      file.path(out_dir, "supplement_manifest.tsv")
    )

    tex_file <- .rep_write_tex(
      manifest,
      file.path(out_dir, "supplement_figures.tex"),
      note_number = note_number,
      dataset_title = dataset_title
    )

    pdf_file <- if (isTRUE(compile) && !is.null(tex_file)) {
      .rep_compile_tex(tex_file, verbose = verbose)
    } else {
      NULL
    }

    files$manifest <- manifest_file
    files$tex <- tex_file
    files$compiled_supplement <- pdf_file
  }

  invisible(list(
    outdir = normalizePath(out_dir, mustWork = FALSE),
    start_figure = start_figure,
    end_figure = if (nrow(manifest)) {
      max(manifest$figure_number)
    } else {
      start_figure - 1L
    },
    next_figure = fig_no,
    note_number = as.character(note_number),
    dataset_title = as.character(dataset_title),
    manifest = manifest,
    manifest_file = manifest_file,
    tex_file = tex_file,
    pdf_file = pdf_file,
    figures = list(
      true_processes = fig_true,
      reconstructed_components = fig_components,
      recovery = fig_recovery,
      term_decomposition = term_plots
    ),
    truth = truth_out,
    reconstructed_components = component_out,
    process_recovery = process_out,
    heldout_biology = heldout_out,
    files = files
  ))
}
