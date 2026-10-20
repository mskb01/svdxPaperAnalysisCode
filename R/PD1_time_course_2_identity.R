library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# Audit whether C1 is driven by myeloid-like cells
#
# Run AFTER GSE249630_PD1_pseudobulk.R, so these objects exist:
#   target_files
#   cell_meta_target
#   read_geo_plate()
#   get_cell_qc()
#
# This script does NOT remove cells.
# It only quantifies T-cell vs myeloid marker expression among
# the retained cells and summarizes it by mouse/day/treatment.
# ============================================================

tcell_markers <- c(
  "Cd3d", "Cd3e", "Cd3g", "Trac",
  "Cd247", "Cd8a", "Cd8b1", "Lck", "Lat"
)

myeloid_markers <- c(
  "Lyz2", "C1qa", "C1qb", "C1qc",
  "Msr1", "Fcer1g", "Fcgr4", "Clec12a",
  "Ccr1", "Cd36", "Dab2", "Lgmn",
  "Alox5ap", "Cybb"
)

score_plate <- function(file, plate_id, cell_meta_target) {

  message("Scoring ", basename(file))

  x <- read_geo_plate(file)

  md_plate <- cell_meta_target %>%
    filter(Amp_batch_ID == plate_id)

  wells <- intersect(colnames(x), md_plate$Well_ID)

  if (length(wells) == 0)
    return(NULL)

  x <- x[, wells, drop = FALSE]

  md <- md_plate[
    match(wells, md_plate$Well_ID),
    ,
    drop = FALSE
  ]

  stopifnot(identical(colnames(x), md$Well_ID))

  # same cell QC used for pseudobulking
  qc <- get_cell_qc(x, md)
  keep <- qc$keep

  x <- x[, keep, drop = FALSE]
  md <- md[keep, , drop = FALSE]

  if (ncol(x) == 0)
    return(NULL)

  genes_lower <- tolower(rownames(x))

  t_idx <- match(
    tolower(tcell_markers),
    genes_lower,
    nomatch = 0
  )
  t_idx <- t_idx[t_idx > 0]

  m_idx <- match(
    tolower(myeloid_markers),
    genes_lower,
    nomatch = 0
  )
  m_idx <- m_idx[m_idx > 0]

  libsize <- colSums(x)

  t_umi <- if (length(t_idx)) {
    colSums(x[t_idx, , drop = FALSE])
  } else {
    rep(0, ncol(x))
  }

  myeloid_umi <- if (length(m_idx)) {
    colSums(x[m_idx, , drop = FALSE])
  } else {
    rep(0, ncol(x))
  }

  t_detected <- if (length(t_idx)) {
    colSums(x[t_idx, , drop = FALSE] > 0)
  } else {
    rep(0, ncol(x))
  }

  myeloid_detected <- if (length(m_idx)) {
    colSums(x[m_idx, , drop = FALSE] > 0)
  } else {
    rep(0, ncol(x))
  }

  out <- md %>%
    mutate(
      library_size = libsize,

      tcell_umi = t_umi,
      myeloid_umi = myeloid_umi,

      tcell_detected = t_detected,
      myeloid_detected = myeloid_detected,

      tcell_fraction = tcell_umi / library_size,
      myeloid_fraction = myeloid_umi / library_size,

      # diagnostic continuum; positive means more myeloid-marker
      # UMIs than T-cell-marker UMIs
      myeloid_vs_tcell =
        log2((myeloid_umi + 1) / (tcell_umi + 1)),

      # intentionally only a diagnostic flag, NOT a QC rule
      myeloid_dominant =
        myeloid_umi > tcell_umi
    )

  out
}


# ------------------------------------------------------------
# Score all retained cells from time_course_2 only
# ------------------------------------------------------------

tf2 <- target_files %>%
  filter(
    GSM %in% unique(
      cell_meta_target$GSM[
        cell_meta_target$experiment_id == "PD1_time_course_2"
      ]
    )
  )

cell_scores_2 <- lapply(seq_len(nrow(tf2)), function(i) {
  score_plate(
    file = tf2$local_file[i],
    plate_id = tf2$Amp_batch_ID[i],
    cell_meta_target = cell_meta_target
  )
}) %>%
  bind_rows()


# ------------------------------------------------------------
# Per-pseudobulk summary
# ------------------------------------------------------------

pb_scores_2 <- cell_scores_2 %>%
  group_by(
    pb_id, mouse, day, treatment, tissue, cell_type
  ) %>%
  summarise(
    n_cells = n(),

    median_tcell_fraction =
      median(tcell_fraction, na.rm = TRUE),

    median_myeloid_fraction =
      median(myeloid_fraction, na.rm = TRUE),

    median_myeloid_vs_tcell =
      median(myeloid_vs_tcell, na.rm = TRUE),

    q90_myeloid_vs_tcell =
      quantile(myeloid_vs_tcell, 0.90, na.rm = TRUE),

    frac_myeloid_dominant =
      mean(myeloid_dominant, na.rm = TRUE),

    median_tcell_detected =
      median(tcell_detected, na.rm = TRUE),

    median_myeloid_detected =
      median(myeloid_detected, na.rm = TRUE),

    .groups = "drop"
  )


# ------------------------------------------------------------
# Day x treatment summary
# ------------------------------------------------------------

day_treatment_scores_2 <- pb_scores_2 %>%
  group_by(day, treatment) %>%
  summarise(
    n_mice = n(),

    median_frac_myeloid_dominant =
      median(frac_myeloid_dominant),

    median_myeloid_vs_tcell =
      median(median_myeloid_vs_tcell),

    median_q90_myeloid_vs_tcell =
      median(q90_myeloid_vs_tcell),

    .groups = "drop"
  )

print(day_treatment_scores_2)


# ------------------------------------------------------------
# Most suspicious pseudobulks
# ------------------------------------------------------------

cat("\nPseudobulks with largest myeloid-like fraction:\n")

print(
  pb_scores_2 %>%
    arrange(desc(frac_myeloid_dominant)) %>%
    select(
      mouse, day, treatment,
      n_cells,
      frac_myeloid_dominant,
      median_myeloid_vs_tcell,
      q90_myeloid_vs_tcell,
      median_tcell_detected,
      median_myeloid_detected
    ) %>%
    head(15)
)


# ------------------------------------------------------------
# Most myeloid-like individual cells
# ------------------------------------------------------------

cat("\nMost myeloid-like retained cells:\n")

print(
  cell_scores_2 %>%
    arrange(desc(myeloid_vs_tcell)) %>%
    select(
      Well_ID, Amp_batch_ID,
      mouse, day, treatment,
      library_size,
      tcell_umi, myeloid_umi,
      tcell_detected, myeloid_detected,
      myeloid_vs_tcell
    ) %>%
    head(30)
)


# ------------------------------------------------------------
# Diagnostic plots
# ------------------------------------------------------------

p1 <- ggplot(
  cell_scores_2,
  aes(
    x = factor(day),
    y = myeloid_vs_tcell,
    fill = treatment
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  labs(
    x = "Day",
    y = "log2((myeloid marker UMIs + 1) / (T-cell marker UMIs + 1))",
    fill = "Treatment"
  ) +
  theme_bw()

ggsave(
  "~/Downloads/PD1_time_course_2_cell_myeloid_vs_tcell.pdf",
  p1,
  width = 7,
  height = 5
)


p2 <- ggplot(
  pb_scores_2,
  aes(
    x = day,
    y = frac_myeloid_dominant,
    color = treatment
  )
) +
  geom_point(size = 2) +
  geom_smooth(se = FALSE) +
  labs(
    x = "Day",
    y = "Fraction of retained cells with myeloid UMIs > T-cell UMIs",
    color = "Treatment"
  ) +
  theme_bw()

ggsave(
  "~/Downloads/PD1_time_course_2_pb_myeloid_fraction.pdf",
  p2,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# Save audit tables
# ------------------------------------------------------------

write.csv(
  cell_scores_2,
  "PD1_time_course_2_cell_identity_audit.csv",
  row.names = FALSE
)

write.csv(
  pb_scores_2,
  "PD1_time_course_2_pseudobulk_identity_audit.csv",
  row.names = FALSE
)

write.csv(
  day_treatment_scores_2,
  "PD1_time_course_2_day_treatment_identity_audit.csv",
  row.names = FALSE
)

saveRDS(
  list(
    cell_scores = cell_scores_2,
    pseudobulk_scores = pb_scores_2,
    day_treatment_scores = day_treatment_scores_2,
    tcell_markers = tcell_markers,
    myeloid_markers = myeloid_markers
  ),
  "PD1_time_course_2_identity_audit.rds"
)
