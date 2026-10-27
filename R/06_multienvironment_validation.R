# Multi-environment chickpea genomic prediction and validation.
# Uses two ICRISAT seasons with grain-yield observations; this is across-year
# validation at one location, not validation at an independent site.
options(timeout = max(600, getOption("timeout")))

source("R/gxe_utils.R")
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read the CicerSeq workbook.")

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

phenotype_file <- "data/raw/3000_Phenotyping_data_QC.xlsx"
phenotype_url <- "https://cegresources.icrisat.org/cicerseq/phenotypic_files/3000_Phenotyping_data_QC.xlsx"
if (!file.exists(phenotype_file)) utils::download.file(phenotype_url, phenotype_file, mode = "wb")

environment_levels <- c("ICRISAT_2014-15", "ICRISAT_2015-16")
trait <- "GY"
phenotype_parts <- lapply(environment_levels, function(env) {
  x <- readxl::read_excel(phenotype_file, sheet = env, .name_repair = "minimal")
  if (!all(c("Genotype", trait) %in% names(x))) {
    stop("Expected Genotype and ", trait, " in sheet ", env, ".")
  }
  data.frame(
    accession = trimws(as.character(x[["Genotype"]])),
    line_id = toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x[["Genotype"]])))),
    environment = env,
    y = suppressWarnings(as.numeric(x[[trait]])),
    stringsAsFactors = FALSE
  )
})
phenotype <- do.call(rbind, phenotype_parts)
phenotype <- phenotype[is.finite(phenotype$y) & !is.na(phenotype$line_id) & nzchar(phenotype$line_id), ]
phenotype$environment <- factor(phenotype$environment, levels = environment_levels)
cell_id <- paste(phenotype$line_id, phenotype$environment, sep = "\r")
if (anyDuplicated(cell_id)) stop("Duplicate line-by-environment grain-yield cells require review.")
if (anyDuplicated(phenotype$line_id[phenotype$environment == environment_levels[1]])) {
  stop("Duplicate accession within first environment.")
}

source_marker_step <- 500L
post_qc_stride <- 5L
min_maf <- 0.05
max_missing <- 0.20
chromosomes <- paste0("Ca", 1:8)
marker_files <- file.path("data/raw", paste0("Cultivated_", chromosomes, ".hmp.txt.gz"))
marker_urls <- paste0("https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_",
                      chromosomes, ".hmp.txt.gz")
for (i in seq_along(marker_files)) {
  if (!file.exists(marker_files[i])) utils::download.file(marker_urls[i], marker_files[i], mode = "wb")
}

normalize_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
iupac <- c(R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC")
decode_dosage <- function(calls, ref, alt) {
  call <- toupper(trimws(calls))
  value <- rep(NA_real_, length(call))
  value[call == ref] <- 0
  value[call == alt] <- 2
  pair <- paste(sort(c(ref, alt)), collapse = "")
  het <- names(iupac)[iupac == pair]
  if (length(het)) value[call == het[1]] <- 1
  value
}

read_thinned_chromosome <- function(path, target_ids) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  metadata_n <- 11L
  if (length(header) <= metadata_n) stop("No sample columns found in ", path)
  canonical_samples <- normalize_id(header[-seq_len(metadata_n)])
  if (anyDuplicated(canonical_samples)) stop("Duplicate normalized sample names in ", path)
  sample_col <- match(target_ids, canonical_samples)
  present <- !is.na(sample_col)
  ids_here <- target_ids[present]
  sample_col <- sample_col[present] + metadata_n
  if (!length(sample_col)) stop("No phenotype accessions matched markers in ", path)
  rows <- list()
  marker_ids <- character()
  global_row <- 0L
  repeat {
    lines <- readLines(con, n = 400L, warn = FALSE)
    if (!length(lines)) break
    row_ix <- global_row + seq_along(lines)
    global_row <- global_row + length(lines)
    take <- which((row_ix - 1L) %% source_marker_step == 0L)
    if (!length(take)) next
    parsed <- strsplit(lines[take], "\t", fixed = TRUE)
    for (one in parsed) {
      if (length(one) < max(sample_col)) next
      alleles <- strsplit(one[2], "/", fixed = TRUE)[[1]]
      if (length(alleles) != 2L || any(!alleles %in% c("A", "C", "G", "T"))) next
      dosage <- decode_dosage(one[sample_col], alleles[1], alleles[2])
      miss <- mean(is.na(dosage))
      maf <- mean(dosage, na.rm = TRUE) / 2
      if (!is.finite(miss) || miss > max_missing || !is.finite(maf) || min(maf, 1 - maf) < min_maf) next
      rows[[length(rows) + 1L]] <- dosage
      marker_ids <- c(marker_ids, paste(one[3], one[4], one[1], sep = ":"))
    }
  }
  if (!length(rows)) stop("No markers passed QC in ", path)
  dosage <- do.call(rbind, rows)
  colnames(dosage) <- ids_here
  rownames(dosage) <- marker_ids
  dosage
}

target_ids <- sort(unique(phenotype$line_id))
cache_file <- "data/processed/cicerseq_two_season_markers.rds"
if (file.exists(cache_file)) {
  cache <- readRDS(cache_file)
  if (!identical(cache$target_ids, target_ids) || cache$source_marker_step != source_marker_step ||
      cache$min_maf != min_maf || cache$max_missing != max_missing) {
    stop("Marker cache settings/IDs differ. Delete ", cache_file, " and rerun.")
  }
  chromosome_matrices <- cache$chromosome_matrices
  message("Using two-season marker cache: ", cache_file)
} else {
  chromosome_matrices <- vector("list", length(chromosomes))
  names(chromosome_matrices) <- chromosomes
  for (i in seq_along(chromosomes)) {
    message("Reading cultivated markers from ", chromosomes[i], " at 1-in-", source_marker_step, " row sampling")
    chromosome_matrices[[i]] <- read_thinned_chromosome(marker_files[i], target_ids)
    message("  retained after MAF/missingness QC: ", nrow(chromosome_matrices[[i]]))
  }
  if (!all(vapply(chromosome_matrices, function(z) identical(colnames(z), colnames(chromosome_matrices[[1]])), logical(1)))) {
    stop("Accession matching/order differs among chromosome files.")
  }
  saveRDS(list(chromosome_matrices = chromosome_matrices, target_ids = target_ids,
               source_marker_step = source_marker_step, min_maf = min_maf,
               max_missing = max_missing), cache_file)
}

matched_ids <- colnames(chromosome_matrices[[1]])
phenotype <- phenotype[phenotype$line_id %in% matched_ids, ]
phenotype <- phenotype[order(phenotype$line_id, phenotype$environment), ]
if (nrow(phenotype) == 0L) stop("No phenotype accessions matched the marker panel.")

marker_counts <- vapply(chromosome_matrices, nrow, integer(1))
offsets <- c(0L, head(cumsum(marker_counts), -1L))
sparse_index <- unlist(lapply(seq_along(chromosome_matrices), function(i) {
  offsets[i] + seq.int(1L, marker_counts[i], by = post_qc_stride)
}), use.names = FALSE)
D <- do.call(rbind, chromosome_matrices)
marker_missing_rate <- mean(is.na(D))
marker_means <- rowMeans(D, na.rm = TRUE)
for (j in seq_len(nrow(D))) D[j, is.na(D[j, ])] <- marker_means[j]
G_full <- t(D)
g_mean <- colMeans(G_full)
g_sd <- apply(G_full, 2L, stats::sd)
valid_marker <- is.finite(g_sd) & g_sd > 0
G_full <- G_full[, valid_marker, drop = FALSE]
G_full <- sweep(G_full, 2L, g_mean[valid_marker], "-")
G_full <- sweep(G_full, 2L, g_sd[valid_marker], "/")
G_full[!is.finite(G_full)] <- 0
sparse_index <- sparse_index[sparse_index <= ncol(G_full)]
G_small <- G_full[, sparse_index, drop = FALSE]

rownames(phenotype) <- NULL
if (!all(phenotype$line_id %in% rownames(G_full))) stop("Phenotype/marker alignment failed.")
cell_counts <- table(phenotype$line_id, phenotype$environment)
paired_ids <- rownames(cell_counts)[rowSums(cell_counts > 0) == length(environment_levels)]
if (length(paired_ids) < 100L) stop("Too few accessions observed in both seasons for CV2/GxE analysis.")

set.seed(20261021L)
pca <- stats::prcomp(G_small, center = FALSE, scale. = FALSE)
pc_n <- min(10L, ncol(pca$x))
cluster_n <- 10L
km <- stats::kmeans(pca$x[, seq_len(pc_n), drop = FALSE], centers = cluster_n,
                    nstart = 25L, iter.max = 100L)
line_cluster <- setNames(km$cluster, rownames(G_small))
saveRDS(list(phenotype = phenotype, G_small = G_small, G_full = G_full, line_cluster = line_cluster,
             environment_levels = environment_levels, paired_ids = paired_ids),
        "data/processed/analysis_objects.rds")

split_by_cluster <- function(d, k, seed) {
  assign_group_folds(line_cluster[d$line_id], k = k, seed = seed)
}
split_by_cv2 <- function(d, k, seed) {
  assign_cv2_folds(d$line_id, d$environment, k = k, seed = seed)
}

lambda_grid <- default_lambda_grid
inner_folds <- 3L
tune_lambda <- function(d, train_idx, G, model, split_function, seed, environment_fixed = TRUE) {
  train <- d[train_idx, , drop = FALSE]
  inner_id <- split_function(train, inner_folds, seed)
  rmse <- rep(NA_real_, length(lambda_grid))
  for (li in seq_along(lambda_grid)) {
    pred <- rep(NA_real_, nrow(train))
    for (f in seq_len(inner_folds)) {
      i_test <- which(inner_id == f)
      if (!length(i_test)) next
      i_train <- which(inner_id != f)
      pred[i_test] <- fit_marker_ridge(train[i_train, , drop = FALSE], train[i_test, , drop = FALSE],
                                       G, lambda_grid[li], model,
                                       env_levels = environment_levels,
                                       environment_fixed = environment_fixed)
    }
    good <- is.finite(pred)
    if (any(good)) rmse[li] <- sqrt(mean((train$y[good] - pred[good])^2))
  }
  if (all(!is.finite(rmse))) stop("Inner CV produced no valid predictions.")
  best <- which.min(rmse)
  list(lambda = lambda_grid[best], rmse = rmse[best])
}

metric_rows <- function(predictions) {
  groups <- split(predictions, interaction(predictions$repeat_id, predictions$scenario,
                                           predictions$environment, predictions$model, drop = TRUE))
  do.call(rbind, lapply(groups, function(d) {
    s <- score_prediction(d$observed, d$predicted)
    if (d$model[1] == "Environment mean") {
      s[c("pearson_r", "spearman_rho", "calibration_slope")] <- NA_real_
    }
    data.frame(scenario = d$scenario[1], repeat_id = d$repeat_id[1],
               environment = d$environment[1], model = d$model[1],
               n = s["n"], rmse = s["rmse"], pearson_r = s["pearson_r"],
               spearman_rho = s["spearman_rho"], calibration_slope = s["calibration_slope"],
               stringsAsFactors = FALSE)
  }))
}

run_nested <- function(d, G, scenario, split_function, models, repeats, seed_base,
                       cv2_baseline = FALSE, environment_fixed = TRUE) {
  prediction_parts <- list()
  fold_parts <- list()
  part <- 0L
  for (r in seq_len(repeats)) {
    outer_id <- split_function(d, 5L, seed_base + r)
    for (f in seq_len(5L)) {
      test_idx <- which(outer_id == f)
      train_idx <- which(outer_id != f)
      if (!length(test_idx)) next
      train <- d[train_idx, , drop = FALSE]
      test <- d[test_idx, , drop = FALSE]
      part <- part + 1L
      for (env in unique(as.character(test$environment))) {
        rows <- which(as.character(test$environment) == env)
        pred_env <- rep(mean(train$y[as.character(train$environment) == env]), length(rows))
        prediction_parts[[length(prediction_parts) + 1L]] <- data.frame(
          scenario = scenario, repeat_id = r, environment = env,
          line_id = test$line_id[rows], observed = test$y[rows],
          model = "Environment mean", predicted = pred_env)
      }
      if (cv2_baseline) {
        line_pred <- rep(NA_real_, nrow(test))
        for (j in seq_len(nrow(test))) {
          other <- which(train$line_id == test$line_id[j] & train$environment != test$environment[j])
          if (length(other) != 1L) stop("CV2 split failed to retain the line's other-season phenotype.")
          env <- as.character(test$environment[j])
          other_env <- as.character(train$environment[other])
          shift <- mean(train$y[as.character(train$environment) == env]) -
            mean(train$y[as.character(train$environment) == other_env])
          line_pred[j] <- train$y[other] + shift
        }
        prediction_parts[[length(prediction_parts) + 1L]] <- data.frame(
          scenario = scenario, repeat_id = r, environment = as.character(test$environment),
          line_id = test$line_id, observed = test$y,
          model = "Known-line phenotype + environment shift", predicted = line_pred)
      }
      for (model in models) {
        tuned <- tune_lambda(d, train_idx, G, model, split_function,
                             seed = seed_base + r * 100L + f,
                             environment_fixed = environment_fixed)
        predicted <- fit_marker_ridge(train, test, G, tuned$lambda, model,
                                      env_levels = environment_levels,
                                      environment_fixed = environment_fixed)
        prediction_parts[[length(prediction_parts) + 1L]] <- data.frame(
          scenario = scenario, repeat_id = r, environment = as.character(test$environment),
          line_id = test$line_id, observed = test$y, model = if (model == "main") "Marker main effects" else "Marker main + GxE",
          predicted = predicted)
        fold_parts[[length(fold_parts) + 1L]] <- data.frame(
          scenario = scenario, repeat_id = r, fold = f, model = model,
          n_train = nrow(train), n_test = nrow(test), selected_lambda = tuned$lambda,
          inner_cv_rmse = tuned$rmse)
      }
    }
    message(scenario, " repeat ", r, "/", repeats, " complete")
  }
  list(predictions = do.call(rbind, prediction_parts), folds = do.call(rbind, fold_parts))
}

core_repeats <- 10L
message("Matched lines: ", nrow(G_full), "; paired across seasons: ", length(paired_ids),
        "; markers for core models: ", ncol(G_small), "; full QC marker set: ", ncol(G_full))

cv1 <- run_nested(phenotype, G_small, "Grouped CV1 - new genetic clusters",
                  split_by_cluster, models = c("main", "gxe"), repeats = core_repeats,
                  seed_base = 20262000L, cv2_baseline = FALSE, environment_fixed = TRUE)
paired <- phenotype[phenotype$line_id %in% paired_ids, , drop = FALSE]
cv2 <- run_nested(paired, G_small, "CV2 - known line in other season",
                  split_by_cv2, models = c("main", "gxe"), repeats = core_repeats,
                  seed_base = 20263000L, cv2_baseline = TRUE, environment_fixed = TRUE)

predictions <- rbind(cv1$predictions, cv2$predictions)
fold_metrics <- rbind(cv1$folds, cv2$folds)
repeat_scores <- metric_rows(predictions)
summary_groups <- split(repeat_scores, interaction(repeat_scores$scenario,
                                                  repeat_scores$environment,
                                                  repeat_scores$model, drop = TRUE))
validation_summary <- do.call(rbind, lapply(summary_groups, function(d) {
  data.frame(scenario = d$scenario[1], environment = d$environment[1], model = d$model[1],
             n = mean(d$n), rmse_mean = mean(d$rmse), rmse_sd_between_repeats = stats::sd(d$rmse),
             pearson_r_mean = if (all(is.na(d$pearson_r))) NA_real_ else mean(d$pearson_r, na.rm = TRUE),
             pearson_r_sd = if (sum(is.finite(d$pearson_r)) < 2L) NA_real_ else stats::sd(d$pearson_r, na.rm = TRUE),
             spearman_rho_mean = if (all(is.na(d$spearman_rho))) NA_real_ else mean(d$spearman_rho, na.rm = TRUE),
             spearman_rho_sd = if (sum(is.finite(d$spearman_rho)) < 2L) NA_real_ else stats::sd(d$spearman_rho, na.rm = TRUE),
             calibration_slope_mean = if (all(is.na(d$calibration_slope))) NA_real_ else mean(d$calibration_slope, na.rm = TRUE),
             calibration_slope_sd = if (sum(is.finite(d$calibration_slope)) < 2L) NA_real_ else stats::sd(d$calibration_slope, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
write.csv(validation_summary, "reports/multienvironment_validation_summary.csv", row.names = FALSE)
write.csv(repeat_scores, "reports/multienvironment_validation_repeats.csv", row.names = FALSE)
write.csv(fold_metrics, "reports/multienvironment_validation_folds.csv", row.names = FALSE)
edge_share_cv <- boundary_share(fold_metrics$selected_lambda, lambda_grid)
if (is.finite(edge_share_cv) && edge_share_cv > 0.10) {
  warning(sprintf("%.0f%% of tuned penalties sit at the edge of the lambda grid; widen it.", 100 * edge_share_cv))
}
paired_rows <- phenotype[phenotype$line_id %in% paired_ids, , drop = FALSE]
y_first <- paired_rows[paired_rows$environment == environment_levels[1], ]
y_second <- paired_rows[paired_rows$environment == environment_levels[2], ]
season_cor <- stats::cor(y_first$y, y_second$y[match(y_first$line_id, y_second$line_id)])
# Predictions are large; store compressed and keep them out of git (see .gitignore).
write.csv(predictions, gzfile("reports/multienvironment_validation_predictions.csv.gz"), row.names = FALSE)

# Leave-one-season-out: train on one year, predict the shared accessions in the
# other year. No target-year fixed effect is estimated from withheld outcomes.
loyo_rows <- list()
for (direction in seq_along(environment_levels)) {
  train_env <- environment_levels[direction]
  test_env <- environment_levels[3L - direction]
  train_idx <- which(as.character(phenotype$environment) == train_env & phenotype$line_id %in% paired_ids)
  test_idx <- which(as.character(phenotype$environment) == test_env & phenotype$line_id %in% paired_ids)
  dtrain <- phenotype[train_idx, , drop = FALSE]
  dtest <- phenotype[test_idx, , drop = FALSE]
  split_line_groups <- function(d, k, seed) assign_group_folds(line_cluster[d$line_id], k, seed)
  tune <- tune_lambda(phenotype, train_idx, G_small, "main", split_line_groups,
                      seed = 20264000L + direction, environment_fixed = FALSE)
  pred <- fit_marker_ridge(dtrain, dtest, G_small, tune$lambda, "main",
                           env_levels = environment_levels, environment_fixed = FALSE)
  base <- rep(mean(dtrain$y), nrow(dtest))
  loyo_rows[[direction]] <- rbind(
    data.frame(direction = paste(train_env, "to", test_env), model = "Training-season mean",
               n = length(base), rmse = score_prediction(dtest$y, base)["rmse"],
               pearson_r = NA_real_, spearman_rho = NA_real_, calibration_slope = NA_real_,
               selected_lambda = NA_real_),
    data.frame(direction = paste(train_env, "to", test_env), model = "Marker main effects",
               n = nrow(dtest), rmse = score_prediction(dtest$y, pred)["rmse"],
               pearson_r = score_prediction(dtest$y, pred)["pearson_r"],
               spearman_rho = score_prediction(dtest$y, pred)["spearman_rho"],
               calibration_slope = score_prediction(dtest$y, pred)["calibration_slope"],
               selected_lambda = tune$lambda)
  )
  message("Leave-one-season-out complete: ", train_env, " -> ", test_env)
}
loyo_summary <- do.call(rbind, loyo_rows)
loyo_lambda <- loyo_summary$selected_lambda[!is.na(loyo_summary$selected_lambda)]
if (any(loyo_lambda >= max(lambda_grid) | loyo_lambda <= min(lambda_grid))) {
  warning("A leave-one-season-out penalty sits at the edge of the lambda grid; widen it.")
}
write.csv(loyo_summary, "reports/leave_one_season_out.csv", row.names = FALSE)

# Marker-density sensitivity under the same grouped CV1 split logic.
density_repeats <- 3L
density_parts <- list()
for (r in seq_len(density_repeats)) {
  outer <- split_by_cluster(phenotype, 5L, 20265000L + r)
  for (density in c("sparse_1_in_2500", "full_1_in_500")) {
    G_use <- if (density == "sparse_1_in_2500") G_small else G_full
    for (f in seq_len(5L)) {
      test_idx <- which(outer == f)
      train_idx <- which(outer != f)
      tune <- tune_lambda(phenotype, train_idx, G_use, "main", split_by_cluster,
                          seed = 20265100L + 10L * r + f, environment_fixed = TRUE)
      train <- phenotype[train_idx, , drop = FALSE]
      test <- phenotype[test_idx, , drop = FALSE]
      pred <- fit_marker_ridge(train, test, G_use, tune$lambda, "main",
                               env_levels = environment_levels, environment_fixed = TRUE)
      density_parts[[length(density_parts) + 1L]] <- data.frame(
        repeat_id = r, fold = f, density = density, selected_lambda = tune$lambda,
        environment = as.character(test$environment), line_id = test$line_id,
        observed = test$y, predicted = pred)
    }
  }
  message("Marker-density sensitivity repeat ", r, "/", density_repeats, " complete")
}
density_predictions <- do.call(rbind, density_parts)
edge_share_density <- boundary_share(vapply(density_parts, function(z) z$selected_lambda[1], numeric(1)), lambda_grid)
density_metric_parts <- lapply(split(density_predictions,
  interaction(density_predictions$repeat_id, density_predictions$density,
              density_predictions$environment, drop = TRUE)), function(d) {
    s <- score_prediction(d$observed, d$predicted)
    data.frame(repeat_id = d$repeat_id[1], density = d$density[1],
               environment = d$environment[1], n = s["n"], rmse = s["rmse"],
               pearson_r = s["pearson_r"], spearman_rho = s["spearman_rho"],
               calibration_slope = s["calibration_slope"])
  })
density_metrics <- do.call(rbind, density_metric_parts)
density_groups <- split(density_metrics, interaction(density_metrics$density,
                                                      density_metrics$environment, drop = TRUE))
density_summary <- do.call(rbind, lapply(density_groups, function(d) {
  data.frame(density = d$density[1], environment = d$environment[1],
             marker_count = if (d$density[1] == "sparse_1_in_2500") ncol(G_small) else ncol(G_full),
             rmse_mean = mean(d$rmse), rmse_sd = stats::sd(d$rmse),
             pearson_r_mean = mean(d$pearson_r), pearson_r_sd = stats::sd(d$pearson_r),
             spearman_rho_mean = mean(d$spearman_rho), spearman_rho_sd = stats::sd(d$spearman_rho))
}))
write.csv(density_summary, "reports/marker_density_sensitivity.csv", row.names = FALSE)
write.csv(density_metrics, "reports/marker_density_sensitivity_repeats.csv", row.names = FALSE)

write.csv(data.frame(
  environments = paste(environment_levels, collapse = ";"), trait = trait,
  phenotype_records = nrow(phenotype), unique_lines = length(unique(phenotype$line_id)),
  paired_lines = length(paired_ids), marker_matched_lines = nrow(G_full),
  chromosomes = length(chromosomes), source_marker_step = source_marker_step,
  post_qc_stride = post_qc_stride, markers_full_qc = ncol(G_full),
  markers_compact = ncol(G_small), marker_missing_before_imputation = marker_missing_rate,
  min_maf = min_maf, max_marker_missing = max_missing,
  groups = cluster_n, pc_count_for_clusters = pc_n, outer_folds = 5L,
  core_repeats = core_repeats, density_repeats = density_repeats,
  inner_folds = inner_folds, lambda_grid = paste(lambda_grid, collapse = ";"),
  seed = 20261021L,
  lambda_edge_share_cv = edge_share_cv, lambda_edge_share_density = edge_share_density,
  season_to_season_pearson = season_cor
), "reports/multienvironment_qc.csv", row.names = FALSE)

# Figures and the results table are built from the saved summaries by R/08_figures_and_tables.R.

writeLines(c(
  "# Multi-environment validation notes", "",
  paste0("Trait: ", trait, "; environments: ", paste(environment_levels, collapse = " and "), "."),
  paste0("", length(paired_ids), " accessions have grain-yield phenotypes in both seasons."),
  paste0("Grouped CV1 holds out all records from marker-derived PC clusters; CV2 holds out line-by-season cells while retaining the same line in the other season."),
  "The GxE model is a ridge model with marker main effects plus marker-by-season interactions; environment is fixed.",
  "Leave-one-season-out tests the same ICRISAT location across two years, using paired lines only. It is not leave-one-location-out validation.",
  "Split SDs summarize repeated random partitions, not population confidence intervals. Marker PCA clusters are defined from genotype data without phenotypes.",
  "Additional-year/site GxE, independent-location prediction, and pedigree/family external validation are not supported by these data alone.",
  "The ridge penalty grid is log-spaced from 1 to 1e5; the share of tuned penalties at a grid edge is in multienvironment_qc.csv.",
  "CicerSeq source files are downloaded locally and are not redistributed in the repository."
), "reports/multienvironment_validation_notes.md")

cat("Multi-environment genomic prediction complete. See reports/multienvironment_validation_summary.csv; run R/08_figures_and_tables.R for figures.\n")
