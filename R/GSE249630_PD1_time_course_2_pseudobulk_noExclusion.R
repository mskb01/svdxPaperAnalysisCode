library(Matrix)
library(dplyr)
library(stringr)

# ============================================================
# GSE249630: PD1_time_course_2
# Mouse-level OT1 tumor pseudobulk, BEFORE the authors' 8-mouse exclusion.
# Publicly reproducible cell QC only:
#   450 <= total UMIs <= 32768
#   mitochondrial fraction < 0.4
# Pseudobulk: mouse (tumor, OT1 CD8+ are constant in this experiment)
# ============================================================

raw_tar        <- "./data/raw/OT1/GSE249630_RAW.tar"
cell_meta_file <- "./data/raw/OT1/GSE249630_metadata_q.txt.gz"
sra_file       <- "./data/raw/OT1/SraRunTable.csv"
extract_dir    <- "./data/raw/OT1/GSE249630_PD1_TC2/"
out_file       <- "PD1_time_course_2_pseudobulk_pre_exclusion.rds"

# ---- Metadata ----
cell_meta <- read.delim(gzfile(cell_meta_file), check.names = FALSE,
                        stringsAsFactors = FALSE)
cell_meta <- cell_meta[, !is.na(names(cell_meta)) & names(cell_meta) != "", drop = FALSE]

sra <- read.csv(sra_file, check.names = FALSE, stringsAsFactors = FALSE)

sample_meta <- sra %>%
  filter(experiment_id == "PD1_time_course_2") %>%
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

stopifnot(all(sample_meta$tissue == "tumor"),
          all(sample_meta$cell_type == "OT1 CD8+"))

# ---- Map GSM to raw plate file ----
tar_files <- untar(raw_tar, list = TRUE)
file_info <- data.frame(tar_file = tar_files, stringsAsFactors = FALSE) %>%
  mutate(
    file_name = basename(tar_file),
    GSM = str_match(file_name, "^(GSM[0-9]+)_")[, 2],
    Amp_batch_ID = str_match(file_name, "^GSM[0-9]+_(AB[0-9]+)\\.txt\\.gz$")[, 2]
  )

target_files <- file_info %>% filter(GSM %in% sample_meta$GSM)
stopifnot(setequal(sample_meta$GSM, target_files$GSM),
          all(table(target_files$GSM) == 1))

sample_meta <- sample_meta %>%
  left_join(target_files %>% select(GSM, Amp_batch_ID, tar_file), by = "GSM")

# ---- Extract target plates ----
dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
untar(raw_tar, files = target_files$tar_file, exdir = extract_dir)
target_files$local_file <- file.path(extract_dir, target_files$tar_file)
stopifnot(all(file.exists(target_files$local_file)))

# ---- Well metadata ----
cell_meta_target <- cell_meta %>%
  inner_join(sample_meta %>%
               select(GSM, Amp_batch_ID, experiment_id, mouse, day,
                      treatment, tissue, cell_type),
             by = "Amp_batch_ID") %>%
  mutate(pb_id = paste(experiment_id, mouse, tissue, cell_type, sep = "__"))

stopifnot(!anyDuplicated(cell_meta_target$Well_ID))

# ---- Read, QC, and aggregate each plate ----
read_plate <- function(file) {
  as.matrix(read.delim(gzfile(file), header = TRUE, row.names = 1,
                       check.names = FALSE))
}

aggregate_plate <- function(file, plate_id) {
  message("Reading ", basename(file))
  x <- read_plate(file)
  md_plate <- cell_meta_target %>% filter(Amp_batch_ID == plate_id)
  wells <- intersect(colnames(x), md_plate$Well_ID)
  if (!length(wells)) return(NULL)

  x <- x[, wells, drop = FALSE]
  md <- md_plate[match(wells, md_plate$Well_ID), , drop = FALSE]
  stopifnot(identical(colnames(x), md$Well_ID))

  libsize <- colSums(x)
  mt <- grepl("^mt-", rownames(x), ignore.case = TRUE)
  mt_counts <- if (any(mt)) colSums(x[mt, , drop = FALSE]) else rep(0, ncol(x))
  mt_frac <- ifelse(libsize > 0, mt_counts / libsize, NA_real_)
  keep <- libsize >= 450 & libsize <= 2^15 & mt_frac < 0.4
  keep[is.na(keep)] <- FALSE

  qc <- data.frame(
    GSM = md$GSM[1], Amp_batch_ID = plate_id, n_wells = ncol(x),
    n_below_450 = sum(libsize < 450),
    n_above_32768 = sum(libsize > 2^15),
    n_high_mito = sum(mt_frac >= 0.4, na.rm = TRUE),
    n_retained = sum(keep)
  )

  x <- x[, keep, drop = FALSE]
  md <- md[keep, , drop = FALSE]
  if (!ncol(x)) return(list(counts = NULL, cell_meta = md, qc = qc))

  pb_levels <- unique(md$pb_id)
  P <- sparseMatrix(
    i = match(md$pb_id, pb_levels),
    j = seq_len(nrow(md)), x = 1,
    dims = c(length(pb_levels), nrow(md)),
    dimnames = list(pb_levels, md$Well_ID)
  )

  list(counts = Matrix(x, sparse = TRUE) %*% t(P), cell_meta = md, qc = qc)
}

res <- lapply(seq_len(nrow(target_files)), function(i) {
  aggregate_plate(target_files$local_file[i], target_files$Amp_batch_ID[i])
})
res <- Filter(Negate(is.null), res)
plate_qc <- bind_rows(lapply(res, `[[`, "qc"))
valid <- Filter(function(z) !is.null(z$counts), res)

# ---- Sum technical plates into biological pseudobulks ----
ref_genes <- rownames(valid[[1]]$counts)
stopifnot(all(vapply(valid, function(z) identical(rownames(z$counts), ref_genes), logical(1))))

pb_ids <- unique(unlist(lapply(valid, function(z) colnames(z$counts))))
pb_counts <- Matrix(0, nrow = length(ref_genes), ncol = length(pb_ids), sparse = TRUE,
                    dimnames = list(ref_genes, pb_ids))
for (z in valid) pb_counts[, colnames(z$counts)] <- pb_counts[, colnames(z$counts), drop = FALSE] + z$counts

retained_meta <- bind_rows(lapply(valid, `[[`, "cell_meta"))
pb_meta <- retained_meta %>%
  group_by(pb_id) %>%
  summarise(
    mouse = first(mouse), day = first(day), treatment = first(treatment),
    tissue = first(tissue), cell_type = first(cell_type),
    n_cells = n(), n_plates = n_distinct(Amp_batch_ID),
    Amp_batch_IDs = paste(sort(unique(Amp_batch_ID)), collapse = ";"),
    GSMs = paste(sort(unique(GSM)), collapse = ";"),
    .groups = "drop"
  )
pb_meta <- pb_meta[match(colnames(pb_counts), pb_meta$pb_id), , drop = FALSE]
pb_meta$library_size <- as.numeric(Matrix::colSums(pb_counts))

stopifnot(identical(colnames(pb_counts), pb_meta$pb_id),
          ncol(pb_counts) == 34)

cat("\nTC2 day x treatment:\n")
print(with(pb_meta, table(day, treatment)))
cat("\nTotal pseudobulks:", ncol(pb_counts), "\n")

saveRDS(
  list(
    counts = pb_counts,
    xdata = pb_meta,
    plate_qc = plate_qc,
    note = "Pre-exclusion TC2 pseudobulk. Authors' eight-mouse exclusion intentionally NOT applied."
  ),
  out_file
)

cat("Saved:", out_file, "\n")
