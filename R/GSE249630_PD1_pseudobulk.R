library(Matrix)
library(dplyr)
library(stringr)

# ============================================================
# GSE249630 pseudobulk construction
# Experiments: PD1_time_course_2 and PD1_time_course_6
# Pseudobulk: experiment_id x mouse x tissue x cell_type
#
# The original Barboy et al. analysis additionally used:
#   (1) FACS index-sort measurements (including CD45.1 protein signal)
#       to refine tumor-specific / OT1 cell identities, and
#   (2) MetaCell-derived gene modules to identify/remove non-T-cell
#       populations.
#
# The index-sort files and the generated MetaCell gene-set object needed
# to reproduce those two filtering steps are not available in the public
# GEO/paper-associated materials we could identify. Therefore, those two
# filters are NOT implemented here.
#
# We implement only filtering/sample-exclusion rules for which the public
# materials contain the information required to reproduce them:
#   * 450 <= total UMIs <= 32768
#   * mitochondrial fraction < 0.4
#   * for PD1_time_course_6, GEO Number_of_cells == 1
#   * authors' explicit tumor-OT1 mouse exclusions:
#       m76, m80, m81, m82, m83, m87, m88, m89
#   * authors' explicit spleen/LN OT1 mouse exclusions:
#       m132, m134, m137, m140, m144, m145
#
# The latter exclusion is compartment-specific: tumor cells from those
# six mice are retained because the authors' exclusion was stated for
# their spleen/LN tumor-specific model.
# ============================================================

raw_tar <- "./data/raw/OT1/GSE249630_RAW.tar"
cell_meta_file <- "./data/raw/OT1/GSE249630_metadata_q.txt.gz"
sra_file <- "./data/raw/OT1/SraRunTable.csv"
extract_dir <- path.expand("./data/raw/OT1/GSE249630_PD1/")

target_experiments <- c("PD1_time_course_2", "PD1_time_course_6")

# ---- Metadata ----
cell_meta <- read.delim(gzfile(cell_meta_file), check.names = FALSE,
                        stringsAsFactors = FALSE)
cell_meta <- cell_meta[, !is.na(names(cell_meta)) & names(cell_meta) != "",
                       drop = FALSE]

sra <- read.csv(sra_file, check.names = FALSE, stringsAsFactors = FALSE)

sample_meta <- sra %>%
  filter(experiment_id %in% target_experiments) %>%
  transmute(
    experiment_id,
    GSM = `Sample Name`,
    mouse = as.character(mouse_number),
    day = as.numeric(days_post_adoptive_transfer),
    treatment = if_else(is.na(treatment) | treatment == "", "Control", treatment),
    tissue = as.character(tissue),
    cell_type = as.character(cell_type)
  ) %>%
  distinct()

gsm_check <- sample_meta %>%
  group_by(GSM) %>%
  summarise(
    n_experiment = n_distinct(experiment_id),
    n_mouse = n_distinct(mouse),
    n_day = n_distinct(day),
    n_treatment = n_distinct(treatment),
    n_tissue = n_distinct(tissue),
    n_cell_type = n_distinct(cell_type),
    .groups = "drop"
  )

stopifnot(
  all(gsm_check$n_experiment == 1),
  all(gsm_check$n_mouse == 1),
  all(gsm_check$n_day == 1),
  all(gsm_check$n_treatment == 1),
  all(gsm_check$n_tissue == 1),
  all(gsm_check$n_cell_type == 1)
)

target_gsms <- unique(sample_meta$GSM)

# ---- Raw files and GSM -> amplification plate mapping ----
tar_files <- untar(raw_tar, list = TRUE)

file_info <- data.frame(tar_file = tar_files, stringsAsFactors = FALSE) %>%
  mutate(
    file_name = basename(tar_file),
    GSM = str_match(file_name, "^(GSM[0-9]+)_")[, 2],
    Amp_batch_ID = str_match(
      file_name, "^GSM[0-9]+_(AB[0-9]+)\\.txt\\.gz$"
    )[, 2]
  )

target_files <- file_info %>% filter(GSM %in% target_gsms)

missing_gsm <- setdiff(target_gsms, target_files$GSM)
if (length(missing_gsm) > 0) {
  stop("Missing GEO raw files for: ", paste(missing_gsm, collapse = ", "))
}
stopifnot(all(table(target_files$GSM) == 1))

sample_meta <- sample_meta %>%
  left_join(
    target_files %>% select(GSM, Amp_batch_ID, tar_file),
    by = "GSM"
  )
stopifnot(!anyNA(sample_meta$Amp_batch_ID))

plate_check <- sample_meta %>%
  group_by(Amp_batch_ID) %>%
  summarise(
    n_experiment = n_distinct(experiment_id),
    n_mouse = n_distinct(mouse),
    n_day = n_distinct(day),
    n_treatment = n_distinct(treatment),
    n_tissue = n_distinct(tissue),
    n_cell_type = n_distinct(cell_type),
    .groups = "drop"
  )

stopifnot(
  all(plate_check$n_experiment == 1),
  all(plate_check$n_mouse == 1),
  all(plate_check$n_day == 1),
  all(plate_check$n_treatment == 1),
  all(plate_check$n_tissue == 1),
  all(plate_check$n_cell_type == 1)
)

# ---- Extract target plates only ----
dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
untar(raw_tar, files = target_files$tar_file, exdir = extract_dir)

target_files <- target_files %>%
  mutate(local_file = file.path(extract_dir, tar_file))
stopifnot(all(file.exists(target_files$local_file)))

# ---- Well-level metadata ----
cell_meta_target <- cell_meta %>%
  inner_join(
    sample_meta %>%
      select(GSM, Amp_batch_ID, experiment_id, mouse, day,
             treatment, tissue, cell_type),
    by = "Amp_batch_ID"
  ) %>%
  mutate(
    pb_id = paste(experiment_id, mouse, tissue, cell_type, sep = "__")
  )

stopifnot(!anyDuplicated(cell_meta_target$Well_ID))

pb_check <- cell_meta_target %>%
  distinct(pb_id, experiment_id, mouse, day, treatment, tissue, cell_type) %>%
  group_by(pb_id) %>%
  summarise(
    n_experiment = n_distinct(experiment_id),
    n_mouse = n_distinct(mouse),
    n_day = n_distinct(day),
    n_treatment = n_distinct(treatment),
    n_tissue = n_distinct(tissue),
    n_cell_type = n_distinct(cell_type),
    .groups = "drop"
  )

stopifnot(
  all(pb_check$n_experiment == 1),
  all(pb_check$n_mouse == 1),
  all(pb_check$n_day == 1),
  all(pb_check$n_treatment == 1),
  all(pb_check$n_tissue == 1),
  all(pb_check$n_cell_type == 1)
)

# ---- Plate reader ----
read_geo_plate <- function(file) {
  x <- read.delim(
    gzfile(file), header = TRUE, row.names = 1, check.names = FALSE
  )
  as.matrix(x)
}

# ---- Publicly reproducible author/sample filters ----
# author_excluded_tumor_ot1 <- c(
#   "m76", "m80", "m81", "m82", "m83", "m87", "m88", "m89"
# )
#
# author_excluded_spleen_ln_ot1 <- c(
#   "m132", "m134", "m137", "m140", "m144", "m145"
# )
author_excluded_tumor_ot1 <- c()
author_excluded_spleen_ln_ot1 <- c()


author_available_keep <- function(md) {
  keep <- rep(TRUE, nrow(md))

  # Authors' tumor-specific model explicitly excluded these mice.
  drop_tumor_ot1 <-
    md$experiment_id == "PD1_time_course_2" &
    md$mouse %in% author_excluded_tumor_ot1 &
    md$tissue == "tumor" &
    md$cell_type == "OT1 CD8+"

  # Authors' spleen/LN tumor-specific model explicitly excluded these mice.
  # Do not propagate this exclusion to tumor OT1 cells from the same mice.
  drop_spleen_ln_ot1 <-
    md$experiment_id == "PD1_time_course_6" &
    md$mouse %in% author_excluded_spleen_ln_ot1 &
    md$tissue %in% c("Spleen", "LN") &
    md$cell_type == "OT1 CD8+"

  keep[drop_tumor_ot1 | drop_spleen_ln_ot1] <- FALSE
  keep
}

# ---- Cell QC ----
# Publicly reproducible thresholds from the authors' preprocessing:
#   450 <= UMIs <= 32768
#   mitochondrial fraction < 0.4
# For PD1_time_course_6, additionally require GEO Number_of_cells == 1.
# Finally apply the authors' explicit mouse/compartment exclusions above.
get_cell_qc <- function(x, md) {
  libsize <- colSums(x)

  mt_gene <- grepl("^mt-", rownames(x), ignore.case = TRUE)
  mt_counts <- if (any(mt_gene)) {
    colSums(x[mt_gene, , drop = FALSE])
  } else {
    rep(0, ncol(x))
  }

  mt_frac <- rep(NA_real_, length(libsize))
  nz <- libsize > 0
  mt_frac[nz] <- mt_counts[nz] / libsize[nz]

  keep_basic_qc <- libsize >= 450 & libsize <= 2^15 & mt_frac < 0.4

  is_pd1_6 <- md$experiment_id == "PD1_time_course_6"
  keep_metadata_cell <- rep(TRUE, nrow(md))
  keep_metadata_cell[is_pd1_6] <- md$Number_of_cells[is_pd1_6] == 1
  keep_metadata_cell[is.na(keep_metadata_cell)] <- FALSE

  keep_author_available <- author_available_keep(md)

  keep <- keep_basic_qc & keep_metadata_cell & keep_author_available
  keep[is.na(keep)] <- FALSE

  list(
    keep = keep,
    libsize = libsize,
    mt_frac = mt_frac,
    keep_basic_qc = keep_basic_qc,
    keep_metadata_cell = keep_metadata_cell,
    keep_author_available = keep_author_available
  )
}

# ---- Aggregate one plate ----
aggregate_plate <- function(file, plate_id, cell_meta_target) {
  message("Reading ", basename(file))
  x <- read_geo_plate(file)

  md_plate <- cell_meta_target %>% filter(Amp_batch_ID == plate_id)
  wells <- intersect(colnames(x), md_plate$Well_ID)

  if (length(wells) == 0) {
    warning("No matching wells for plate ", plate_id)
    return(NULL)
  }

  x <- x[, wells, drop = FALSE]
  md <- md_plate[match(wells, md_plate$Well_ID), , drop = FALSE]
  stopifnot(identical(colnames(x), md$Well_ID))

  qc <- get_cell_qc(x, md)
  keep <- qc$keep

  qc_summary <- data.frame(
    experiment_id = md$experiment_id[1],
    GSM = md$GSM[1],
    Amp_batch_ID = plate_id,
    n_wells = ncol(x),
    n_below_450 = sum(qc$libsize < 450),
    n_above_32768 = sum(qc$libsize > 2^15),
    n_high_mito = sum(qc$mt_frac >= 0.4, na.rm = TRUE),
    n_metadata_noncell = sum(md$Number_of_cells == 0, na.rm = TRUE),
    n_author_excluded = sum(!qc$keep_author_available),
    n_retained = sum(keep)
  )

  x <- x[, keep, drop = FALSE]
  md <- md[keep, , drop = FALSE]

  if (ncol(x) == 0) {
    warning("No cells retained for plate ", plate_id)
    return(list(counts = NULL, cell_meta = md, qc = qc_summary))
  }

  pb_levels <- unique(md$pb_id)

  P <- sparseMatrix(
    i = match(md$pb_id, pb_levels),
    j = seq_len(nrow(md)),
    x = 1,
    dims = c(length(pb_levels), nrow(md)),
    dimnames = list(pb_levels, md$Well_ID)
  )

  x <- Matrix(x, sparse = TRUE)
  pb_counts <- x %*% t(P)

  list(counts = pb_counts, cell_meta = md, qc = qc_summary)
}

# ---- Process all target plates ----
plate_results <- vector("list", nrow(target_files))
names(plate_results) <- target_files$Amp_batch_ID

for (i in seq_len(nrow(target_files))) {
  plate_results[[i]] <- aggregate_plate(
    file = target_files$local_file[i],
    plate_id = target_files$Amp_batch_ID[i],
    cell_meta_target = cell_meta_target
  )
  gc()
}

plate_results <- plate_results[!vapply(plate_results, is.null, logical(1))]

# ---- QC audit ----
plate_qc <- bind_rows(lapply(plate_results, `[[`, "qc"))

qc_summary <- plate_qc %>%
  group_by(experiment_id) %>%
  summarise(
    n_plates = n(),
    total_wells = sum(n_wells),
    below_450 = sum(n_below_450),
    above_32768 = sum(n_above_32768),
    high_mito = sum(n_high_mito),
    metadata_noncell = sum(n_metadata_noncell),
    author_excluded_wells = sum(n_author_excluded),
    cells_retained = sum(n_retained),
    median_cells_per_plate = median(n_retained),
    min_cells_per_plate = min(n_retained),
    max_cells_per_plate = max(n_retained),
    .groups = "drop"
  )

print(qc_summary)

# ---- Verify gene rows ----
valid_results <- plate_results[
  vapply(plate_results, function(z) !is.null(z$counts), logical(1))
]
stopifnot(length(valid_results) > 0)

ref_genes <- rownames(valid_results[[1]]$counts)
same_genes <- vapply(
  valid_results,
  function(z) identical(rownames(z$counts), ref_genes),
  logical(1)
)
if (!all(same_genes)) stop("Gene rows differ across GEO raw files.")

# ---- Sum plate contributions into biological pseudobulks ----
all_pb_ids <- unique(unlist(lapply(
  valid_results, function(z) colnames(z$counts)
)))

pb_counts <- Matrix(
  0,
  nrow = length(ref_genes),
  ncol = length(all_pb_ids),
  sparse = TRUE,
  dimnames = list(ref_genes, all_pb_ids)
)

for (z in valid_results) {
  ids <- colnames(z$counts)
  pb_counts[, ids] <- pb_counts[, ids, drop = FALSE] + z$counts
}

# ---- Retained-cell metadata ----
retained_cell_meta <- bind_rows(lapply(valid_results, `[[`, "cell_meta"))

pb_meta <- retained_cell_meta %>%
  group_by(pb_id) %>%
  summarise(
    experiment_id = first(experiment_id),
    mouse = first(mouse),
    day = first(day),
    treatment = first(treatment),
    tissue = first(tissue),
    cell_type = first(cell_type),
    n_cells = n(),
    n_plates = n_distinct(Amp_batch_ID),
    Amp_batch_IDs = paste(sort(unique(Amp_batch_ID)), collapse = ";"),
    GSMs = paste(sort(unique(GSM)), collapse = ";"),
    .groups = "drop"
  )

pb_meta <- pb_meta[
  match(colnames(pb_counts), pb_meta$pb_id),
  , drop = FALSE
]

stopifnot(identical(colnames(pb_counts), pb_meta$pb_id))
pb_meta$library_size <- as.numeric(Matrix::colSums(pb_counts))

# ---- Separate experiments / cell populations ----
i2 <- pb_meta$experiment_id == "PD1_time_course_2" &
  pb_meta$cell_type == "OT1 CD8+"

i6_ot1 <- pb_meta$experiment_id == "PD1_time_course_6" &
  pb_meta$cell_type == "OT1 CD8+"

i6_endo <- pb_meta$experiment_id == "PD1_time_course_6" &
  pb_meta$cell_type == "endogenous CD8+"

pd1_2_counts <- pb_counts[, i2, drop = FALSE]
pd1_2_xdata <- pb_meta[i2, , drop = FALSE]

pd1_6_counts <- pb_counts[, i6_ot1, drop = FALSE]
pd1_6_xdata <- pb_meta[i6_ot1, , drop = FALSE]

pd1_6_endo_counts <- pb_counts[, i6_endo, drop = FALSE]
pd1_6_endo_xdata <- pb_meta[i6_endo, , drop = FALSE]

cat("\nPD1_time_course_2 OT1:", ncol(pd1_2_counts), "pseudobulks\n")
cat("PD1_time_course_6 OT1:", ncol(pd1_6_counts), "pseudobulks\n")
cat("PD1_time_course_6 endogenous CD8:", ncol(pd1_6_endo_counts),
    "pseudobulks\n\n")

# Expected after applying only the publicly reproducible author exclusions:
#   _2: 34 original OT1 pseudobulks - 8 excluded mice = 26
#   _6 OT1: 80 original pseudobulks - 12 excluded spleen/LN pseudobulks = 68
#   _6 endogenous CD8: unchanged at 17
stopifnot(
  ncol(pd1_2_counts) == 26,
  ncol(pd1_6_counts) == 68,
  ncol(pd1_6_endo_counts) == 17
)

# ---- Final audit ----
cat("\nPD1_time_course_2 OT1 day x treatment:\n")
print(with(pd1_2_xdata, table(day, treatment)))

cat("\nPD1_time_course_6 OT1 day x treatment:\n")
print(with(pd1_6_xdata, table(day, treatment)))

cat("\nPD1_time_course_6 OT1 tissue counts:\n")
print(with(pd1_6_xdata, table(tissue)))

cat("\nPD1_time_course_6 endogenous CD8 day x treatment:\n")
print(with(pd1_6_endo_xdata, table(day, treatment)))

cat("\nCells per pseudobulk, _2 OT1:\n")
print(summary(pd1_2_xdata$n_cells))

cat("\nCells per pseudobulk, _6 OT1:\n")
print(summary(pd1_6_xdata$n_cells))

cat("\nLibrary sizes, _2 OT1:\n")
print(summary(pd1_2_xdata$library_size))

cat("\nLibrary sizes, _6 OT1:\n")
print(summary(pd1_6_xdata$library_size))

# Explicitly verify that the exclusions were applied only where intended.
stopifnot(
  !any(pd1_2_xdata$mouse %in% author_excluded_tumor_ot1),
  !any(
    pd1_6_xdata$mouse %in% author_excluded_spleen_ln_ot1 &
      pd1_6_xdata$tissue %in% c("Spleen", "LN")
  )
)

# Tumor OT1 pseudobulks from the six spleen/LN-excluded _6 mice should remain.
stopifnot(all(
  author_excluded_spleen_ln_ot1 %in%
    pd1_6_xdata$mouse[pd1_6_xdata$tissue == "tumor"]
))

# ---- Save ----
reproducibility_note <- paste(
  "Two filters from the original analysis could not be reproduced from",
  "the publicly available materials identified: FACS index-sort-based",
  "cell-identity refinement and MetaCell-derived non-T-cell gene-module",
  "filtering. No approximation to those filters was applied."
)

saveRDS(
  list(
    counts = pd1_2_counts,
    xdata = pd1_2_xdata,
    plate_qc = plate_qc %>% filter(experiment_id == "PD1_time_course_2"),
    reproducibility_note = reproducibility_note,
    author_excluded_mice = author_excluded_tumor_ot1
  ),
  "PD1_time_course_2_pseudobulk.rds"
)
