# make_single_simulation_processes_main.R
#
# Reconstruct the four NET biological process structures used in the
# single-dataset simulation from one fitted SVDX component.
#
# Simulation processes:
#   1. group-specific mean shift              -> group
#   2. shared oscillation                     -> t
#   3. group-specific phase oscillation       -> t + group:t
#   4. severity-modulated oscillation         -> t + severity + t:severity
#
# Usage:
#   p_main <- make_single_simulation_processes_main(
#     fit = fit,
#     dat = dat,
#     component = 1
#   )
#   p_main
#
# Optionally:
#   ggsave("single_simulation_processes_main.pdf",
#          p_main, width = 12, height = 8)

library(ggplot2)
library(patchwork)

# ---- small helpers -----------------------------------------------------------

.pick_component <- function(x, component = 1) {
  x <- as.matrix(x)
  if (ncol(x) < component) {
    stop("Requested component exceeds the number of predicted components.",
         call. = FALSE)
  }
  as.numeric(x[, component])
}

.typical_value <- function(x) {
  if (is.factor(x)) {
    lev <- names(which.max(table(x)))
    return(factor(lev, levels = levels(x)))
  }

  if (is.character(x)) {
    return(names(which.max(table(x))))
  }

  if (is.logical(x)) {
    return(names(which.max(table(x))) == "TRUE")
  }

  stats::median(x, na.rm = TRUE)
}

.restore_class <- function(x, template) {
  if (is.factor(template)) {
    return(factor(x, levels = levels(template)))
  }
  x
}

.predict_terms <- function(fit, newdata, include, component) {
  .pick_component(
    predict(
      fit,
      newdata = newdata,
      include = include
    ),
    component = component
  )
}

# ---- main figure -------------------------------------------------------------

make_single_simulation_processes_main <- function(
    fit,
    dat,
    component = 1,
    n_grid = 200,
    severity_probs = c(0.10, 0.50, 0.90)
) {

  req <- c("group", "t", "severity")
  miss <- setdiff(req, names(dat))
  if (length(miss)) {
    stop("dat is missing required variables: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  # Typical values used only for variables outside the process being displayed.
  group_typical <- .typical_value(dat$group)
  severity_typical <- .typical_value(dat$severity)

  t_grid <- seq(
    min(dat$t, na.rm = TRUE),
    max(dat$t, na.rm = TRUE),
    length.out = n_grid
  )

  # ---------------------------------------------------------------------------
  # 1. Group-specific mean shift: f_g(group)
  # ---------------------------------------------------------------------------

  group_vals <- if (is.factor(dat$group)) levels(dat$group) else sort(unique(dat$group))

  nd_group <- data.frame(
    group = group_vals,
    t = stats::median(dat$t, na.rm = TRUE),
    severity = severity_typical
  )
  nd_group$group <- .restore_class(nd_group$group, dat$group)

  nd_group$process <- .predict_terms(
    fit = fit,
    newdata = nd_group,
    include = "group",
    component = component
  )

  p_group <- ggplot(nd_group, aes(x = group, y = process)) +
    geom_col(width = 0.75) +
    labs(
      title = "Group-specific mean",
      x = "Group",
      y = "Reconstructed process"
    ) +
    theme_bw()

  # ---------------------------------------------------------------------------
  # 2. Shared oscillation: f_t(t)
  # ---------------------------------------------------------------------------

  nd_time <- data.frame(
    t = t_grid,
    severity = severity_typical
  )
  nd_time$group <- rep(group_typical, nrow(nd_time))
  nd_time$group <- .restore_class(nd_time$group, dat$group)

  nd_time$process <- .predict_terms(
    fit = fit,
    newdata = nd_time,
    include = "t",
    component = component
  )

  p_time <- ggplot(nd_time, aes(x = t, y = process)) +
    geom_line(linewidth = 1.05) +
    labs(
      title = "Shared oscillation",
      x = "Time",
      y = "Reconstructed process"
    ) +
    theme_bw()

  # ---------------------------------------------------------------------------
  # 3. Group-specific phase oscillation:
  #       f_t(t) + f_{group:t}(group, t)
  #
  # Severity is held at its typical value.
  # ---------------------------------------------------------------------------

  nd_phase <- expand.grid(
    t = t_grid,
    group = group_vals,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  nd_phase$group <- .restore_class(nd_phase$group, dat$group)
  nd_phase$severity <- severity_typical

  nd_phase$process <- .predict_terms(
    fit = fit,
    newdata = nd_phase,
    include = c("t", "group", "group:t"),
    component = component
  )

  p_phase <- ggplot(
    nd_phase,
    aes(x = t, y = process, color = group, group = group)
  ) +
    geom_line(linewidth = 1.0) +
    labs(
      title = "Group-specific phase",
      x = "Time",
      y = "Reconstructed process",
      color = "Group"
    ) +
    theme_bw()

  # ---------------------------------------------------------------------------
  # 4. Severity-modulated oscillation:
  #       f_t(t) + f_s(severity) + f_{t:s}(t, severity)
  #
  # Group is held at its typical level.
  # ---------------------------------------------------------------------------

  sev_vals <- as.numeric(
    stats::quantile(
      dat$severity,
      probs = severity_probs,
      na.rm = TRUE,
      names = FALSE
    )
  )

  nd_severity <- expand.grid(
    t = t_grid,
    severity = sev_vals,
    KEEP.OUT.ATTRS = FALSE
  )
  nd_severity$group <- rep(group_typical, nrow(nd_severity))
  nd_severity$group <- .restore_class(nd_severity$group, dat$group)

  nd_severity$process <- .predict_terms(
    fit = fit,
    newdata = nd_severity,
    include = c("t", "severity", "t:severity"),
    component = component
  )

  sev_labels <- c("Low severity", "Median severity", "High severity")
  if (length(sev_vals) != 3L) {
    sev_labels <- paste0(
      round(100 * severity_probs), "th percentile"
    )
  }

  nd_severity$severity_level <- factor(
    nd_severity$severity,
    levels = sev_vals,
    labels = sev_labels
  )

  p_severity <- ggplot(
    nd_severity,
    aes(
      x = t,
      y = process,
      color = severity_level,
      group = severity_level
    )
  ) +
    geom_line(linewidth = 1.0) +
    labs(
      title = "Severity-modulated oscillation",
      x = "Time",
      y = "Reconstructed process",
      color = "Severity"
    ) +
    theme_bw()

  # ---------------------------------------------------------------------------
  # Assemble 2 x 2 main-figure panel.
  # ---------------------------------------------------------------------------

  ((p_group | p_time) /
     (p_phase | p_severity)) +
    plot_annotation(
      title = paste0("CF", component, ": reconstructed biological processes")
    )
}

# Example:
#
# p_main <- make_single_simulation_processes_main(
#   fit = fit,
#   dat = dat,
#   component = 1
# )
#
# p_main
#
# ggsave(
#   "single_simulation_processes_main.pdf",
#   p_main,
#   width = 12,
#   height = 8
# )
