# R/preprocessing_helpers.R

is_sparse_matrix <- function(x) {
  inherits(x, "sparseMatrix")
}

as_numeric_matrix_or_sparse <- function(Y) {
  if (is_sparse_matrix(Y)) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required for sparse matrix support.")
    }
    return(Y)
  }
  
  Y <- as.matrix(Y)
  
  if (!is.numeric(Y)) {
    stop("Y_raw must be numeric.")
  }
  
  Y
}

any_negative <- function(Y) {
  if (is_sparse_matrix(Y)) {
    any(Y@x < 0)
  } else {
    any(Y < 0)
  }
}

any_nonfinite_matrix <- function(Y) {
  if (is_sparse_matrix(Y)) {
    any(!is.finite(Y@x))
  } else {
    any(!is.finite(Y))
  }
}

matrix_col_sums <- function(Y) {
  if (is_sparse_matrix(Y)) {
    Matrix::colSums(Y)
  } else {
    colSums(Y)
  }
}

matrix_row_means <- function(Y) {
  if (is_sparse_matrix(Y)) {
    Matrix::rowMeans(Y)
  } else {
    rowMeans(Y)
  }
}

matrix_row_sums_sq <- function(Y) {
  if (is_sparse_matrix(Y)) {
    Matrix::rowSums(Y^2)
  } else {
    rowSums(Y^2)
  }
}

matrix_row_vars <- function(Y) {
  n <- ncol(Y)
  
  if (n <= 1) {
    stop("Need at least two samples to compute row variances.")
  }
  
  mu <- matrix_row_means(Y)
  ss <- matrix_row_sums_sq(Y)
  
  vars <- (ss - n * mu^2) / (n - 1)
  vars[vars < 0 & vars > -1e-12] <- 0
  
  vars
}

compute_size_factors <- function(t_j) {
  if (is.null(t_j)) {
    stop("t_j cannot be NULL.")
  }
  
  if (any(!is.finite(t_j))) {
    stop("t_j contains non-finite values.")
  }
  
  if (any(t_j <= 0)) {
    stop("All entries of t_j must be positive.")
  }
  
  t_j / exp(mean(log(t_j)))
}

formula_variables <- function(formula) {
  if (is.null(formula)) {
    return(character(0))
  }
  
  all.vars(formula)
}

is_allowed_model_column <- function(x) {
  is.numeric(x) || is.integer(x) || is.factor(x) || is.ordered(x)
}

validate_formula_data <- function(data, formula, data_name, formula_name) {
  vars <- formula_variables(formula)
  
  if (length(vars) == 0) {
    return(invisible(TRUE))
  }
  
  missing_vars <- setdiff(vars, colnames(data))
  
  if (length(missing_vars) > 0) {
    stop(
      formula_name,
      " contains variables not found in ",
      data_name,
      ": ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  for (v in vars) {
    x <- data[[v]]
    
    if (!is_allowed_model_column(x)) {
      stop(
        data_name,
        "$",
        v,
        " has unsupported type: ",
        paste(class(x), collapse = "/"),
        ". Variables in ",
        formula_name,
        " must be numeric, integer, factor, or ordered factor."
      )
    }
    
    if (any(is.na(x))) {
      stop(data_name, "$", v, " contains missing values.")
    }
    
    if (is.factor(x) || is.ordered(x)) {
      if (nlevels(x) < 2) {
        stop(data_name, "$", v, " must have at least two levels.")
      }
    }
    
    if (is.numeric(x) || is.integer(x)) {
      if (!is.finite(stats::var(x)) || stats::var(x) == 0) {
        stop(data_name, "$", v, " has zero or non-finite variance.")
      }
    }
  }
  
  invisible(TRUE)
}

validate_sample_alignment <- function(Y, xdata, zdata) {
  if (nrow(xdata) != ncol(Y)) {
    stop("nrow(xdata) must equal ncol(Y).")
  }
  
  if (nrow(zdata) != ncol(Y)) {
    stop("nrow(zdata) must equal ncol(Y).")
  }
  
  if (!is.null(rownames(xdata)) && !is.null(colnames(Y))) {
    if (!all(rownames(xdata) == colnames(Y))) {
      stop("rownames(xdata) must exactly match colnames(Y).")
    }
  }
  
  if (!is.null(rownames(zdata)) && !is.null(colnames(Y))) {
    if (!all(rownames(zdata) == colnames(Y))) {
      stop("rownames(zdata) must exactly match colnames(Y).")
    }
  }
  
  invisible(TRUE)
}

sparse_safe_col_normalize <- function(Y, t_tilde_j) {
  if (is_sparse_matrix(Y)) {
    Y %*% Matrix::Diagonal(x = 1 / t_tilde_j)
  } else {
    sweep(Y, 2, t_tilde_j, "/")
  }
}

sparse_safe_log1p <- function(Y) {
  if (is_sparse_matrix(Y)) {
    Y@x <- log1p(Y@x)
    Y
  } else {
    log1p(Y)
  }
}

select_top_hvg <- function(Y, top_ftrs) {
  if (is.null(top_ftrs)) {
    return(seq_len(nrow(Y)))
  }
  
  top_ftrs <- as.integer(top_ftrs)
  
  if (length(top_ftrs) != 1 || top_ftrs <= 0) {
    stop("top_ftrs must be a single positive integer or NULL.")
  }
  
  n_select <- min(top_ftrs, nrow(Y))
  
  vars <- matrix_row_vars(Y)
  
  if (any(!is.finite(vars))) {
    stop("Some feature variances are non-finite.")
  }
  
  order(vars, decreasing = TRUE)[seq_len(n_select)]
}

row_center_scale <- function(
    Y,
    center_rows = TRUE,
    scale_rows = TRUE,
    allow_sparse_densification = FALSE
) {
  row_center <- rep(0, nrow(Y))
  row_scale <- rep(1, nrow(Y))
  
  names(row_center) <- rownames(Y)
  names(row_scale) <- rownames(Y)
  
  if (is_sparse_matrix(Y) && center_rows && !allow_sparse_densification) {
    stop(
      "center_rows = TRUE would densify sparse Y. ",
      "Set center_rows = FALSE, or explicitly set ",
      "allow_sparse_densification = TRUE."
    )
  }
  
  if (is_sparse_matrix(Y) && center_rows && allow_sparse_densification) {
    warning(
      "Converting sparse Y to dense matrix because center_rows = TRUE ",
      "and allow_sparse_densification = TRUE."
    )
    Y <- as.matrix(Y)
  }
  
  if (center_rows) {
    row_center <- rowMeans(Y)
    Y <- sweep(Y, 1, row_center, "-")
  }
  
  if (scale_rows) {
    if (is_sparse_matrix(Y)) {
      vars <- matrix_row_vars(Y)
      row_scale <- sqrt(vars)
    } else {
      row_scale <- apply(Y, 1, stats::sd)
    }
    
    bad_scale <- !is.finite(row_scale) | row_scale == 0
    
    if (any(bad_scale)) {
      warning(
        sum(bad_scale),
        " rows had zero or non-finite SD; setting scale to 1."
      )
      row_scale[bad_scale] <- 1
    }
    
    if (is_sparse_matrix(Y)) {
      Y <- Matrix::Diagonal(x = 1 / row_scale) %*% Y
    } else {
      Y <- sweep(Y, 1, row_scale, "/")
    }
  }
  
  list(
    Y = Y,
    row_center = row_center,
    row_scale = row_scale
  )
}

standardize_family <- function(family) {
  if (is.null(family)) {
    return("gaussian")
  }
  
  if (!is.character(family) || length(family) != 1) {
    stop("family must be NULL, 'gaussian', or 'poisson'.")
  }
  
  family <- match.arg(family, choices = c("gaussian", "poisson"))
  family
}

preprocess_svdx_data <- function(
    Y_raw,
    xdata,
    xformula,
    zdata = NULL,
    zformula = ~1,
    family = NULL,
    row_offsets = NULL,
    col_offsets = NULL,
    t_j = NULL,
    top_ftrs = TOP_FTRS_DEFAULT,
    normalize = NORMALIZE_DEFAULT,
    log1p_transform = LOG1P_DEFAULT,
    center_rows = CENTER_ROWS_DEFAULT,
    scale_rows = SCALE_ROWS_DEFAULT,
    allow_sparse_densification = ALLOW_SPARSE_DENSIFICATION_DEFAULT,
    feature_info = NULL,
    sample_info = NULL,
    dataset_name = NULL
) {
  if (missing(xformula)) {
    stop("xformula must be provided. There is no SVDX analysis without xformula.")
  }
  
  family_type <- standardize_family(family)
  
  Y_raw <- as_numeric_matrix_or_sparse(Y_raw)
  
  if (any_nonfinite_matrix(Y_raw)) {
    stop("Y_raw contains non-finite values.")
  }
  
  if (is.null(rownames(Y_raw))) {
    rownames(Y_raw) <- paste0("feature_", seq_len(nrow(Y_raw)))
  }
  
  if (is.null(colnames(Y_raw))) {
    colnames(Y_raw) <- paste0("sample_", seq_len(ncol(Y_raw)))
  }
  
  xdata <- as.data.frame(xdata)
  
  if (is.null(zdata)) {
    zdata <- xdata
  } else {
    zdata <- as.data.frame(zdata)
  }
  
  validate_sample_alignment(Y_raw, xdata, zdata)
  
  validate_formula_data(
    data = xdata,
    formula = xformula,
    data_name = "xdata",
    formula_name = "xformula"
  )
  
  validate_formula_data(
    data = zdata,
    formula = zformula,
    data_name = "zdata",
    formula_name = "zformula"
  )
  
  Y_has_negative <- any_negative(Y_raw)
  
  if (family_type == "poisson") {
    if (Y_has_negative) {
      stop("family = 'poisson' requires nonnegative count data.")
    }
    
    if (log1p_transform) {
      stop(
        "log1p_transform = TRUE is not allowed when family = 'poisson'. ",
        "Poisson SVDX expects count-scale data."
      )
    }
    
    if (center_rows) {
      stop(
        "center_rows = TRUE is not allowed when family = 'poisson'. ",
        "Poisson SVDX should keep Y on the count scale."
      )
    }
    
    if (scale_rows) {
      stop(
        "scale_rows = TRUE is not allowed when family = 'poisson'. ",
        "Poisson SVDX should keep Y on the count scale."
      )
    }
  }
  
  explicit_t_j <- !is.null(t_j)
  
  Y <- Y_raw
  
  normalize_requested <- normalize
  normalize_used_on_Y <- FALSE
  normalization_applied_via_col_offsets <- FALSE
  normalization_reason <- "Normalization not requested."
  
  log1p_requested <- log1p_transform
  log1p_used <- FALSE
  log1p_reason <- "log1p transform not requested."
  
  if (is.null(col_offsets)) {
    col_offsets <- rep(0, ncol(Y_raw))
    names(col_offsets) <- colnames(Y_raw)
  } else if (length(col_offsets) == ncol(Y_raw)) {
    names(col_offsets) <- colnames(Y_raw)
  } else {
    stop("col_offsets must have length ncol(Y_raw), or be NULL.")
  }
  
  if (is.null(t_j)) {
    if (normalize) {
      if (Y_has_negative) {
        stop(
          "Default column-sum normalization requires nonnegative data. ",
          "Y_raw contains negative values. Provide explicit t_j only if ",
          "normalization is scientifically appropriate, or set normalize = FALSE."
        )
      }
      
      t_j <- matrix_col_sums(Y_raw)
      names(t_j) <- colnames(Y_raw)
    } else {
      t_j <- rep(NA_real_, ncol(Y_raw))
      names(t_j) <- colnames(Y_raw)
    }
  } else {
    if (length(t_j) != ncol(Y_raw)) {
      stop("length(t_j) must equal ncol(Y_raw).")
    }
    
    names(t_j) <- colnames(Y_raw)
  }
  
  t_tilde_j <- rep(1, ncol(Y_raw))
  names(t_tilde_j) <- colnames(Y_raw)
  
  if (normalize) {
    t_tilde_j <- compute_size_factors(t_j)
    
    if (family_type == "poisson") {
      col_offsets <- col_offsets + log(t_tilde_j)
      normalization_applied_via_col_offsets <- TRUE
      normalization_reason <- paste(
        "For family='poisson', normalization was applied by adding",
        "log(t_tilde_j) to col_offsets; Y was left on the count scale."
      )
    } else {
      Y <- sparse_safe_col_normalize(Y, t_tilde_j)
      normalize_used_on_Y <- TRUE
      normalization_reason <- "Applied geometric-mean anchored column normalization to Y."
    }
  }
  
  if (family_type == "gaussian" && log1p_transform) {
    if (any_negative(Y)) {
      stop(
        "Cannot apply log1p transform because Y contains negative values. ",
        "Set log1p_transform = FALSE."
      )
    }
    
    Y <- sparse_safe_log1p(Y)
    log1p_used <- TRUE
    log1p_reason <- "Applied log1p transform."
  }
  
  selected_idx <- select_top_hvg(Y, top_ftrs = top_ftrs)
  
  Y <- Y[selected_idx, , drop = FALSE]
  selected_features <- rownames(Y)
  
  if (is.null(row_offsets)) {
    row_offsets <- rep(0, nrow(Y))
    names(row_offsets) <- rownames(Y)
  } else if (length(row_offsets) == nrow(Y_raw)) {
    row_offsets <- row_offsets[selected_idx]
    names(row_offsets) <- rownames(Y)
  } else if (length(row_offsets) == nrow(Y)) {
    names(row_offsets) <- rownames(Y)
  } else {
    stop(
      "row_offsets must have length nrow(Y_raw), ",
      "length nrow(Y after feature selection), or be NULL."
    )
  }
  
  if (family_type == "gaussian") {
    cs <- row_center_scale(
      Y,
      center_rows = center_rows,
      scale_rows = scale_rows,
      allow_sparse_densification = allow_sparse_densification
    )
    
    Y_processed <- cs$Y
    row_center <- cs$row_center
    row_scale <- cs$row_scale
  } else {
    Y_processed <- Y
    
    row_center <- rep(0, nrow(Y_processed))
    row_scale <- rep(1, nrow(Y_processed))
    
    names(row_center) <- rownames(Y_processed)
    names(row_scale) <- rownames(Y_processed)
  }
  
  if (is.null(feature_info)) {
    feature_info <- data.frame(
      feature_id = rownames(Y_raw),
      stringsAsFactors = FALSE
    )
  }
  
  if (!"feature_id" %in% colnames(feature_info)) {
    feature_info$feature_id <- rownames(Y_raw)
  }
  
  feature_info_selected <- feature_info[
    match(selected_features, feature_info$feature_id),
    ,
    drop = FALSE
  ]
  
  if (is.null(sample_info)) {
    sample_info <- data.frame(
      sample_id = colnames(Y_raw),
      stringsAsFactors = FALSE
    )
  }
  
  stopifnot(nrow(Y_processed) == length(row_offsets))
  stopifnot(ncol(Y_processed) == length(col_offsets))
  stopifnot(nrow(xdata) == ncol(Y_processed))
  stopifnot(nrow(zdata) == ncol(Y_processed))
  
  list(
    dataset_name = dataset_name,
    family = family,
    family_type = family_type,
    
    Y = Y_processed,
    
    xdata = xdata,
    zdata = zdata,
    
    xformula = xformula,
    zformula = zformula,
    
    row_offsets = row_offsets,
    col_offsets = col_offsets,
    
    feature_info = feature_info_selected,
    sample_info = sample_info,
    
    preprocessing = list(
      dataset_name = dataset_name,
      matrix_orientation = "features_by_samples",
      sparse_input = is_sparse_matrix(Y_raw),
      sparse_output = is_sparse_matrix(Y_processed),
      allow_sparse_densification = allow_sparse_densification,
      
      family = family,
      family_type = family_type,
      
      xformula_text = paste(deparse(xformula), collapse = " "),
      zformula_text = paste(deparse(zformula), collapse = " "),
      
      normalize_requested = normalize_requested,
      normalize_used_on_Y = normalize_used_on_Y,
      normalization_applied_via_col_offsets = normalization_applied_via_col_offsets,
      normalization_reason = normalization_reason,
      
      log1p_requested = log1p_requested,
      log1p_used = log1p_used,
      log1p_reason = log1p_reason,
      
      top_ftrs = top_ftrs,
      center_rows = center_rows,
      scale_rows = scale_rows,
      
      Y_raw_had_negative_values = Y_has_negative,
      
      t_j = t_j,
      t_tilde_j = t_tilde_j,
      
      selected_idx = selected_idx,
      selected_features = selected_features,
      
      row_center = row_center,
      row_scale = row_scale,
      
      offset_note = paste(
        "row_offsets and col_offsets are assumed to be supplied",
        "on the appropriate scale by the user. For family='poisson',",
        "normalization factors are incorporated by adding log(t_tilde_j)",
        "to col_offsets."
      )
    )
  )
}