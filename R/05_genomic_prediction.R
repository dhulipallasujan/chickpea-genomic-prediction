# Genome-wide prediction on the public ICRISAT CicerSeq chickpea panel.
# Source files are downloaded locally and are excluded from Git.
options(timeout = max(600, getOption("timeout")))

required <- c("readxl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Install readxl before running this analysis.")
}

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

phenotype_file <- "data/raw/3000_Phenotyping_data_QC.xlsx"
phenotype_url <- "https://cegresources.icrisat.org/cicerseq/phenotypic_files/3000_Phenotyping_data_QC.xlsx"
if (!file.exists(phenotype_file)) {
  utils::download.file(phenotype_url, phenotype_file, mode = "wb", quiet = FALSE)
}

chromosomes <- paste0("Ca", 1:8)
marker_files <- file.path("data/raw", paste0("Cultivated_", chromosomes, ".hmp.txt.gz"))
marker_urls <- paste0("https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_",
                      chromosomes, ".hmp.txt.gz")
for (i in seq_along(marker_files)) {
  if (!file.exists(marker_files[i])) utils::download.file(marker_urls[i], marker_files[i], mode = "wb")
}

normalize_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", trimws(x)))
iupac <- c(R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC")
decode_dosage <- function(calls, ref, alt) {
  call <- toupper(trimws(calls))
  value <- rep(NA_real_, length(call))
  value[call == ref] <- 0
  value[call == alt] <- 2
  pair <- paste(sort(c(ref, alt)), collapse = "")
  heterozygous <- names(iupac)[iupac == pair]
  if (length(heterozygous)) value[call == heterozygous[1]] <- 1
  value
}

target_sheet <- "ICRISAT_2014-15"
target_trait <- "GY"
source_marker_step <- 500L
post_qc_stride <- 5L
min_maf <- 0.05
max_missing <- 0.20
folds <- 5L
repeats <- 5L
inner_folds <- 3L
source("R/gxe_utils.R")
lambda_grid <- default_lambda_grid
random_seed <- 20261020L

phenotype <- readxl::read_excel(phenotype_file, sheet = target_sheet, .name_repair = "minimal")
if (!all(c("Genotype", target_trait) %in% names(phenotype))) {
  stop("Expected Genotype and ", target_trait, " columns; found: ", paste(names(phenotype), collapse = ", "))
}
phenotype$Genotype <- trimws(as.character(phenotype$Genotype))
phenotype$y <- suppressWarnings(as.numeric(phenotype[[target_trait]]))
phenotype <- phenotype[is.finite(phenotype$y) & !is.na(phenotype$Genotype) & nzchar(phenotype$Genotype), ]
if (anyDuplicated(normalize_id(phenotype$Genotype))) {
  stop("Phenotype IDs collide after punctuation/space normalization; resolve before analysis.")
}

# Read one in every 100 HapMap rows across each chromosome, then QC the
# sampled SNPs. This phenotype-blind thinning spreads retained SNPs genome-wide.
read_thinned_chromosome <- function(path, target_ids, step) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  metadata_n <- 11L
  if (length(header) <= metadata_n) stop("No accession columns in ", path)
  sample_ids <- header[-seq_len(metadata_n)]
  canonical_samples <- normalize_id(sample_ids)
  if (anyDuplicated(canonical_samples)) stop("Duplicate marker sample IDs after normalization in ", path)
  sample_col <- match(normalize_id(target_ids), canonical_samples)
  present <- !is.na(sample_col)
  target_ids <- target_ids[present]
  sample_col <- sample_col[present] + metadata_n
  if (!length(sample_col)) stop("No phenotype accessions matched marker samples in ", path)

  rows <- list()
  marker_ids <- character()
  global_row <- 0L
  repeat {
    lines <- readLines(con, n = 400L, warn = FALSE)
    if (!length(lines)) break
    row_ix <- global_row + seq_along(lines)
    global_row <- global_row + length(lines)
    take <- which((row_ix - 1L) %% step == 0L)
    if (!length(take)) next
    parsed <- strsplit(lines[take], "\t", fixed = TRUE)
    for (one in parsed) {
      if (length(one) < max(sample_col)) next
      alleles <- strsplit(one[2], "/", fixed = TRUE)[[1]]
      if (length(alleles) != 2L || any(!alleles %in% c("A", "C", "G", "T"))) next
      dosage <- decode_dosage(one[sample_col], alleles[1], alleles[2])
      missing_rate <- mean(is.na(dosage))
      if (!is.finite(missing_rate) || missing_rate > max_missing) next
      frequency <- mean(dosage, na.rm = TRUE) / 2
      if (!is.finite(frequency) || min(frequency, 1 - frequency) < min_maf) next
      rows[[length(rows) + 1L]] <- dosage
      marker_ids <- c(marker_ids, paste(one[3], one[4], one[1], sep = ":"))
    }
  }
  if (!length(rows)) stop("No biallelic markers passed QC in ", path)
  dosage <- do.call(rbind, rows)
  colnames(dosage) <- target_ids
  rownames(dosage) <- marker_ids
  dosage
}

target_ids <- phenotype$Genotype
chromosome_matrices <- vector("list", length(marker_files))
names(chromosome_matrices) <- chromosomes
cache_file <- "data/processed/cicerseq_sampled_markers.rds"
if (file.exists(cache_file)) {
  cached <- readRDS(cache_file)
  chromosome_matrices <- cached$chromosome_matrices
  cached_step <- if (!is.null(cached$marker_step)) cached$marker_step else 500L
  if (source_marker_step != cached_step) {
    stop("Cached marker thinning step differs. Remove ", cache_file, " and rerun.")
  }
  message("Using cached QC-passed markers from ", cache_file)
} else {
  for (i in seq_along(marker_files)) {
    message("Reading markers from ", chromosomes[i], " (systematic source sampling every ", source_marker_step, " rows)")
    chromosome_matrices[[i]] <- read_thinned_chromosome(marker_files[i], target_ids, source_marker_step)
    message("  SNPs retained after QC: ", nrow(chromosome_matrices[[i]]))
  }
  dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(chromosome_matrices = chromosome_matrices, marker_step = source_marker_step), cache_file)
}
qc_passed_before_post_qc_thinning <- sum(vapply(chromosome_matrices, nrow, integer(1)))
chromosome_matrices <- lapply(chromosome_matrices, function(z) {
  z[seq.int(1L, nrow(z), by = post_qc_stride), , drop = FALSE]
})
matched_ids <- colnames(chromosome_matrices[[1]])
if (!all(vapply(chromosome_matrices, function(z) identical(colnames(z), matched_ids), logical(1)))) {
  stop("Accession matching/order differs among chromosome files.")
}
phenotype <- phenotype[match(matched_ids, phenotype$Genotype), ]
D <- do.call(rbind, chromosome_matrices)
missing_before_imputation <- mean(is.na(D))
marker_mean <- rowMeans(D, na.rm = TRUE)
for (j in seq_len(nrow(D))) D[j, is.na(D[j, ])] <- marker_mean[j]

# Standardize the QC-passed marker dosages and fit ridge regression (RR-BLUP)
# in marker space. Penalty selection is nested within each outer training set.
G <- t(D)
marker_mean <- colMeans(G)
marker_sd <- apply(G, 2L, sd)
G <- sweep(G, 2L, marker_mean, FUN = "-")
G <- sweep(G, 2L, marker_sd, FUN = "/")
G[!is.finite(G)] <- 0

n <- length(matched_ids)
if (n < folds * 4L) stop("Too few matched phenotypes for cross-validation.")
set.seed(random_seed)
fold_summary <- vector("list", repeats * folds)
prediction_list <- vector("list", repeats)
baseline_list <- vector("list", repeats)
counter <- 0L
for (r in seq_len(repeats)) {
  fold_id <- sample(rep_len(seq_len(folds), n))
  pred <- rep(NA_real_, n)
  baseline_pred <- rep(NA_real_, n)
  for (f in seq_len(folds)) {
    test <- which(fold_id == f)
    train <- which(fold_id != f)
    inner_assignment <- sample(rep_len(seq_len(inner_folds), length(train)))
    inner_rmse <- numeric(length(lambda_grid))
    for (li in seq_along(lambda_grid)) {
      inner_pred <- rep(NA_real_, length(train))
      for (ik in seq_len(inner_folds)) {
        inner_test <- which(inner_assignment == ik)
        inner_train <- which(inner_assignment != ik)
        X_inner <- G[train[inner_train], , drop = FALSE]
        inner_intercept <- mean(phenotype$y[train[inner_train]])
        inner_effect <- solve(
          crossprod(X_inner) + diag(lambda_grid[li], ncol(X_inner)),
          crossprod(X_inner, phenotype$y[train[inner_train]] - inner_intercept)
        )
        inner_pred[inner_test] <- as.numeric(
          inner_intercept + G[train[inner_test], , drop = FALSE] %*% inner_effect
        )
      }
      inner_rmse[li] <- sqrt(mean((phenotype$y[train] - inner_pred)^2))
    }
    selected_lambda <- lambda_grid[which.min(inner_rmse)]
    X_train <- G[train, , drop = FALSE]
    intercept <- mean(phenotype$y[train])
    marker_effect <- solve(crossprod(X_train) + diag(selected_lambda, ncol(X_train)),
                           crossprod(X_train, phenotype$y[train] - intercept))
    pred[test] <- as.numeric(intercept + G[test, , drop = FALSE] %*% marker_effect)
    baseline_pred[test] <- intercept
    counter <- counter + 1L
    fold_summary[[counter]] <- data.frame(
      repeat_id = r, fold = f, n_train = length(train), n_test = length(test),
      selected_lambda = selected_lambda, inner_cv_rmse = min(inner_rmse),
      pearson_r = cor(phenotype$y[test], pred[test]),
      spearman_rho = suppressWarnings(cor(phenotype$y[test], pred[test], method = "spearman")),
      rmse_rrblup = sqrt(mean((phenotype$y[test] - pred[test])^2)),
      rmse_mean = sqrt(mean((phenotype$y[test] - mean(phenotype$y[train]))^2))
    )
  }
  prediction_list[[r]] <- pred
  baseline_list[[r]] <- baseline_pred
  message("Completed outer CV repeat ", r, " of ", repeats)
}
fold_metrics <- do.call(rbind, fold_summary)
repeat_metrics <- do.call(rbind, lapply(seq_len(repeats), function(r) {
  pred <- prediction_list[[r]]
  baseline_pred <- baseline_list[[r]]
  data.frame(
    repeat_id = r,
    pearson_r = cor(phenotype$y, pred),
    spearman_rho = suppressWarnings(cor(phenotype$y, pred, method = "spearman")),
    calibration_slope = unname(coef(lm(phenotype$y ~ pred))[2]),
    rmse_rrblup = sqrt(mean((phenotype$y - pred)^2)),
    rmse_mean = sqrt(mean((phenotype$y - baseline_pred)^2))
  )
}))
repeat_metrics$rmse_reduction_pct <- 100 * (repeat_metrics$rmse_mean - repeat_metrics$rmse_rrblup) / repeat_metrics$rmse_mean
summary_table <- data.frame(
  metric = c("Pearson correlation", "Spearman rank correlation", "Calibration slope (observed ~ predicted)", "RR-BLUP RMSE", "Training-mean RMSE", "RMSE reduction (%)"),
  mean = c(mean(repeat_metrics$pearson_r), mean(repeat_metrics$spearman_rho),
           mean(repeat_metrics$calibration_slope),
           mean(repeat_metrics$rmse_rrblup), mean(repeat_metrics$rmse_mean),
           mean(repeat_metrics$rmse_reduction_pct)),
  sd_between_repeats = c(sd(repeat_metrics$pearson_r), sd(repeat_metrics$spearman_rho),
                         sd(repeat_metrics$calibration_slope),
                         sd(repeat_metrics$rmse_rrblup), sd(repeat_metrics$rmse_mean),
                         sd(repeat_metrics$rmse_reduction_pct))
)

write.csv(summary_table, "reports/genomic_prediction_summary.csv", row.names = FALSE)
write.csv(repeat_metrics, "reports/genomic_prediction_repeats.csv", row.names = FALSE)
write.csv(fold_metrics, "reports/genomic_prediction_folds.csv", row.names = FALSE)
edge_share <- boundary_share(fold_metrics$selected_lambda, lambda_grid)
if (is.finite(edge_share) && edge_share > 0.10) {
  warning(sprintf("%.0f%% of tuned penalties sit at the edge of the lambda grid; widen it.", 100 * edge_share))
}
write.csv(data.frame(
  target_environment = target_sheet, target_trait = target_trait,
  phenotype_rows_with_trait = nrow(phenotype), matched_accessions = n,
  chromosomes = length(chromosomes), qc_passed_before_post_qc_thinning = qc_passed_before_post_qc_thinning,
  markers_after_qc = nrow(D), source_marker_step = source_marker_step,
  post_qc_stride = post_qc_stride,
  min_maf = min_maf, max_marker_missing = max_missing,
  marker_missing_rate_before_imputation = missing_before_imputation,
  folds = folds, repeats = repeats, inner_folds = inner_folds,
  lambda_grid = paste(lambda_grid, collapse = ";"), seed = random_seed
), "reports/genomic_prediction_qc.csv", row.names = FALSE)
prediction_matrix <- do.call(cbind, prediction_list)
write.csv(data.frame(
  Genotype = matched_ids, Observed = phenotype$y,
  Predicted_mean = rowMeans(prediction_matrix),
  Predicted_sd_across_repeats = apply(prediction_matrix, 1L, sd)
), "reports/genomic_prediction_predictions.csv", row.names = FALSE)

# Describe marker-derived structure for QC. PCs are not used as predictors.
pca <- stats::prcomp(G, center = FALSE, scale. = FALSE)
pca_variance <- pca$sdev^2 / sum(pca$sdev^2)
write.csv(data.frame(PC = paste0("PC", seq_along(pca_variance)),
                     variance_explained = pca_variance,
                     cumulative_variance = cumsum(pca_variance)),
          "reports/marker_pca_variance.csv", row.names = FALSE)
write.csv(data.frame(Genotype = matched_ids, PC1 = pca$x[, 1], PC2 = pca$x[, 2]),
          "reports/marker_pca_scores.csv", row.names = FALSE)
png("figures/marker_population_structure_pca.png", width = 1000, height = 800, res = 140)
plot(pca$x[, 1], pca$x[, 2], pch = 16, cex = 0.55,
     col = grDevices::adjustcolor("#2D6A8A", alpha.f = 0.55),
     xlab = sprintf("PC1 (%.1f%%)", 100 * pca_variance[1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pca_variance[2]),
     main = "Cultivated chickpea panel - marker QC")
grid(col = "grey90")
dev.off()

png("figures/genomic_prediction_cv.png", width = 1500, height = 700, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1))
pred_mean <- rowMeans(prediction_matrix)
plot(phenotype$y, pred_mean, pch = 16, cex = 0.55,
     col = grDevices::adjustcolor("#2D6A8A", alpha.f = 0.45),
     xlab = "Observed grain yield (GY)", ylab = "Out-of-fold RR-BLUP prediction",
     main = "New-line prediction - ICRISAT 2014-15")
abline(lm(pred_mean ~ phenotype$y), col = "#D55E00", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "grey45")
legend("topleft", legend = c("Observed = predicted", "Fit: predicted on observed"),
       lty = c(2, 1), col = c("grey45", "#D55E00"), bty = "n", cex = 0.85)
rmse_means <- c(mean(repeat_metrics$rmse_rrblup), mean(repeat_metrics$rmse_mean))
rmse_sds <- c(sd(repeat_metrics$rmse_rrblup), sd(repeat_metrics$rmse_mean))
plot(1:2, rmse_means, pch = 19, cex = 1.5, xaxt = "n",
     xlim = c(0.6, 2.4), ylim = range(c(rmse_means - rmse_sds, rmse_means + rmse_sds)) + c(-12, 15),
     col = c("#2D6A8A", "#777777"), xlab = "", ylab = "RMSE (GY units)",
     main = "Out-of-fold error - lower is better")
axis(1, at = 1:2, labels = c("Nested RR-BLUP", "Training-mean baseline"))
arrows(1:2, rmse_means - rmse_sds, 1:2, rmse_means + rmse_sds,
       angle = 90, code = 3, length = 0.08, lwd = 1.5)
text(1:2, rmse_means, labels = sprintf("%.1f", rmse_means), pos = 3, cex = 0.85)
mtext(sprintf("RMSE reduction %.1f%% | Pearson r %.2f | Spearman rho %.2f",
              mean(repeat_metrics$rmse_reduction_pct), mean(repeat_metrics$pearson_r),
              mean(repeat_metrics$spearman_rho)), side = 3, line = 0.2, cex = 0.75)
mtext("Points: mean; whiskers: split-to-split SD. Y-axis is zoomed.", side = 1, line = 3.2, cex = 0.7)
par(op)
dev.off()

writeLines(c(
  "# Genomic prediction result notes", "",
  paste0("Target: ", target_trait, " in ", target_sheet, "."),
  paste0("Matched accessions: ", n, "; retained genome-wide SNPs: ", nrow(D), "."),
  paste0("Five-fold cross-validation repeated ", repeats, " times; split-to-split SD is not a population confidence interval."),
  "Phenotypes are accession means, not plot-level records. This is a within-environment CV1 benchmark, not independent external validation.",
  paste0("Read one marker in every ", source_marker_step, " HapMap rows, applied MAF/missingness QC, then retained every ", post_qc_stride, "th QC-passed marker for the compact model."),
  paste0("Mean RMSE reduction versus training-mean baseline: ", round(mean(repeat_metrics$rmse_reduction_pct), 1), "%."),
  paste0("Mean Pearson r: ", round(mean(repeat_metrics$pearson_r), 3),
         "; Spearman rho: ", round(mean(repeat_metrics$spearman_rho), 3),
         "; calibration slope (observed ~ predicted): ", round(mean(repeat_metrics$calibration_slope), 3), "."),
  paste0("RR-BLUP used standardized markers; ridge penalty was selected from ",
         paste(lambda_grid, collapse = ", "), " by inner three-fold CV within each training fold."),
  "CicerSeq asks users to follow the spirit of the Toronto Statement; cite the source and do not redistribute the source files."
), "reports/genomic_prediction_notes.md")

cat("Genomic prediction complete. See reports/genomic_prediction_summary.csv and figures/genomic_prediction_cv.png.\n")
