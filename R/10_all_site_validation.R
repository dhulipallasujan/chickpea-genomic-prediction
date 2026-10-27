# Multi-location chickpea prediction using the common yield-per-plant trait.
# Site-year sheets use YPP, except the two ICRISAT sheets, which label the
# corresponding trait YPPlnt. Locations, years, IDs, and duplicated cells are
# audited and written to reports before model comparisons.
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

normalize_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
sheets <- setdiff(readxl::excel_sheets(phenotype_file), "Description")
phenotype_parts <- lapply(sheets, function(sheet) {
  dat <- readxl::read_excel(phenotype_file, sheet = sheet, .name_repair = "minimal")
  location <- sub("_[0-9]{4}-[0-9]{2}$", "", sheet)
  year <- sub("^.*_", "", sheet)
  trait_column <- if (location == "ICRISAT") "YPPlnt" else "YPP"
  if (!all(c("Genotype", trait_column) %in% names(dat))) {
    stop("Expected Genotype and ", trait_column, " in sheet ", sheet, ".")
  }
  raw_id <- trimws(as.character(dat[["Genotype"]]))
  y <- suppressWarnings(as.numeric(dat[[trait_column]]))
  data.frame(accession = raw_id, line_id = normalize_id(raw_id),
             environment = sheet, location = location, year = year, y = y,
             source_trait = trait_column, stringsAsFactors = FALSE)
})
phenotype_raw <- do.call(rbind, phenotype_parts)
phenotype_raw <- phenotype_raw[is.finite(phenotype_raw$y) & !is.na(phenotype_raw$line_id) & nzchar(phenotype_raw$line_id), ]
cell_key <- paste(phenotype_raw$line_id, phenotype_raw$environment, sep = "\r")
duplicate_rows_collapsed <- sum(duplicated(cell_key))
# The source contains a few repeated records after punctuation normalization.
# Average them to one accession-by-site-year phenotype before CV.
phenotype <- stats::aggregate(y ~ line_id + environment + location + year,
                              data = phenotype_raw, FUN = mean)
phenotype$source_trait <- ifelse(phenotype$location == "ICRISAT", "YPPlnt", "YPP")
environment_levels <- sort(unique(phenotype$environment))
location_levels <- sort(unique(phenotype$location))
phenotype$environment <- factor(phenotype$environment, levels = environment_levels)
phenotype$location <- factor(phenotype$location, levels = location_levels)
if (anyDuplicated(paste(phenotype$line_id, phenotype$environment, sep = "\r"))) {
  stop("Duplicate accession-by-environment cells remain after aggregation.")
}

source_marker_step <- 500L
# This all-site extension uses a tractable 1-in-20 subset of QC-passed markers;
# the denser 1-in-5 panel remains available for the two-season analysis.
post_qc_stride <- 20L
min_maf <- 0.05
max_missing <- 0.20
chromosomes <- paste0("Ca", 1:8)
marker_files <- file.path("data/raw", paste0("Cultivated_", chromosomes, ".hmp.txt.gz"))
marker_urls <- paste0("https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_",
                      chromosomes, ".hmp.txt.gz")
for (i in seq_along(marker_files)) {
  if (!file.exists(marker_files[i])) utils::download.file(marker_urls[i], marker_files[i], mode = "wb")
}
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
  canonical <- normalize_id(header[-seq_len(metadata_n)])
  if (anyDuplicated(canonical)) stop("Duplicate normalized marker sample names in ", path)
  sample_col <- match(target_ids, canonical)
  present <- !is.na(sample_col)
  ids_here <- target_ids[present]
  sample_col <- sample_col[present] + metadata_n
  if (!length(sample_col)) stop("No phenotype accessions matched markers in ", path)
  rows <- list()
  marker_ids <- character()
  row_number <- 0L
  repeat {
    lines <- readLines(con, n = 400L, warn = FALSE)
    if (!length(lines)) break
    indices <- row_number + seq_along(lines)
    row_number <- row_number + length(lines)
    take <- which((indices - 1L) %% source_marker_step == 0L)
    if (!length(take)) next
    for (one in strsplit(lines[take], "\t", fixed = TRUE)) {
      if (length(one) < max(sample_col)) next
      alleles <- strsplit(one[2], "/", fixed = TRUE)[[1]]
      if (length(alleles) != 2L || any(!alleles %in% c("A", "C", "G", "T"))) next
      dosage <- decode_dosage(one[sample_col], alleles[1], alleles[2])
      maf <- mean(dosage, na.rm = TRUE) / 2
      miss <- mean(is.na(dosage))
      if (!is.finite(maf) || !is.finite(miss) || miss > max_missing || min(maf, 1 - maf) < min_maf) next
      rows[[length(rows) + 1L]] <- dosage
      marker_ids <- c(marker_ids, paste(one[3], one[4], one[1], sep = ":"))
    }
  }
  if (!length(rows)) stop("No markers passed QC in ", path)
  matrix <- do.call(rbind, rows)
  colnames(matrix) <- ids_here
  rownames(matrix) <- marker_ids
  matrix
}

target_ids <- sort(unique(phenotype$line_id))
cache_file <- "data/processed/cicerseq_all_site_ypp_markers.rds"
if (file.exists(cache_file)) {
  cache <- readRDS(cache_file)
  if (!identical(cache$target_ids, target_ids) || cache$source_marker_step != source_marker_step ||
      cache$min_maf != min_maf || cache$max_missing != max_missing) {
    stop("Marker cache settings/IDs differ. Delete ", cache_file, " and rerun.")
  }
  chromosome_matrices <- cache$chromosome_matrices
  message("Using all-site marker cache: ", cache_file)
} else {
  chromosome_matrices <- vector("list", length(chromosomes))
  names(chromosome_matrices) <- chromosomes
  for (i in seq_along(chromosomes)) {
    message("Reading cultivated markers from ", chromosomes[i], " at 1-in-", source_marker_step, " row sampling")
    chromosome_matrices[[i]] <- read_thinned_chromosome(marker_files[i], target_ids)
    message("  retained after QC: ", nrow(chromosome_matrices[[i]]))
  }
  if (!all(vapply(chromosome_matrices, function(z) identical(colnames(z), colnames(chromosome_matrices[[1]])), logical(1)))) {
    stop("Accession matching/order differs among chromosome files.")
  }
  saveRDS(list(chromosome_matrices = chromosome_matrices, target_ids = target_ids,
               source_marker_step = source_marker_step, min_maf = min_maf,
               max_missing = max_missing), cache_file)
}

matched_ids <- colnames(chromosome_matrices[[1]])
phenotype_before_match <- nrow(phenotype)
phenotype <- phenotype[phenotype$line_id %in% matched_ids, , drop = FALSE]
if (nrow(phenotype) == 0L) stop("No phenotype accessions matched the marker panel.")
marker_counts <- vapply(chromosome_matrices, nrow, integer(1))
offsets <- c(0L, head(cumsum(marker_counts), -1L))
sparse_index <- unlist(lapply(seq_along(chromosome_matrices), function(i) {
  offsets[i] + seq.int(1L, marker_counts[i], by = post_qc_stride)
}), use.names = FALSE)
D <- do.call(rbind, chromosome_matrices)
missing_before_imputation <- mean(is.na(D))
marker_means <- rowMeans(D, na.rm = TRUE)
for (j in seq_len(nrow(D))) D[j, is.na(D[j, ])] <- marker_means[j]
G_all <- t(D)
g_mean <- colMeans(G_all)
g_sd <- apply(G_all, 2L, stats::sd)
valid <- is.finite(g_sd) & g_sd > 0
G_all <- sweep(G_all[, valid, drop = FALSE], 2L, g_mean[valid], "-")
G_all <- sweep(G_all, 2L, g_sd[valid], "/")
G_all[!is.finite(G_all)] <- 0
sparse_index <- sparse_index[sparse_index <= ncol(G_all)]
G <- G_all[, sparse_index, drop = FALSE]
phenotype <- phenotype[order(phenotype$line_id, phenotype$environment), , drop = FALSE]
if (!all(phenotype$line_id %in% rownames(G))) stop("Phenotype/marker alignment failed.")

set.seed(20261021L)
pca <- stats::prcomp(G, center = FALSE, scale. = FALSE)
pc_n <- min(10L, ncol(pca$x))
cluster_n <- 10L
km <- stats::kmeans(pca$x[, seq_len(pc_n), drop = FALSE], centers = cluster_n,
                    nstart = 25L, iter.max = 100L)
line_cluster <- setNames(km$cluster, rownames(G))
split_by_cluster <- function(d, k, seed) assign_group_folds(line_cluster[d$line_id], k, seed)
split_by_cv2 <- function(d, k, seed) assign_cv2_folds(d$line_id, d$environment, k, seed)

lambda_grid <- 10^seq(0, 5, length.out = 5L)
inner_folds <- 3L
tune_lambda <- function(d, train_idx, G, model, split_function, seed,
                        environment_fixed = TRUE, env_levels = environment_levels) {
  train <- d[train_idx, , drop = FALSE]
  inner_id <- split_function(train, inner_folds, seed)
  rmse <- rep(NA_real_, length(lambda_grid))
  for (li in seq_along(lambda_grid)) {
    predicted <- rep(NA_real_, nrow(train))
    for (fold in seq_len(inner_folds)) {
      test_idx <- which(inner_id == fold)
      if (!length(test_idx)) next
      inner_train <- which(inner_id != fold)
      predicted[test_idx] <- fit_marker_ridge(train[inner_train, , drop = FALSE],
                                               train[test_idx, , drop = FALSE], G,
                                               lambda_grid[li], model, env_levels,
                                               environment_fixed)
    }
    good <- is.finite(predicted)
    if (any(good)) rmse[li] <- sqrt(mean((train$y[good] - predicted[good])^2))
  }
  if (all(!is.finite(rmse))) stop("Inner CV produced no predictions for ", model, ".")
  list(lambda = lambda_grid[which.min(rmse)], rmse = min(rmse, na.rm = TRUE))
}

score_groups <- function(predictions) {
  groups <- split(predictions, interaction(predictions$repeat_id, predictions$scenario,
                                           predictions$environment, predictions$model, drop = TRUE))
  do.call(rbind, lapply(groups, function(d) {
    s <- score_prediction(d$observed, d$predicted)
    if (d$model[1] == "Environment mean") s[c("pearson_r", "spearman_rho", "calibration_slope")] <- NA_real_
    data.frame(scenario = d$scenario[1], repeat_id = d$repeat_id[1],
               environment = d$environment[1], model = d$model[1], n = s["n"],
               rmse = s["rmse"], pearson_r = s["pearson_r"],
               spearman_rho = s["spearman_rho"], calibration_slope = s["calibration_slope"])
  }))
}
summarize_scores <- function(scores) {
  groups <- split(scores, interaction(scores$scenario, scores$environment, scores$model, drop = TRUE))
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  sd_or_na <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else stats::sd(x, na.rm = TRUE)
  do.call(rbind, lapply(groups, function(d) data.frame(
    scenario = d$scenario[1], environment = d$environment[1], model = d$model[1],
    n_mean = mean(d$n), rmse_mean = mean(d$rmse), rmse_sd = sd_or_na(d$rmse),
    pearson_r_mean = mean_or_na(d$pearson_r), pearson_r_sd = sd_or_na(d$pearson_r),
    spearman_rho_mean = mean_or_na(d$spearman_rho), spearman_rho_sd = sd_or_na(d$spearman_rho),
    calibration_slope_mean = mean_or_na(d$calibration_slope),
    calibration_slope_sd = sd_or_na(d$calibration_slope), stringsAsFactors = FALSE)))
}

run_nested <- function(d, scenario, split_function, repeats, seed_base, cv2 = FALSE) {
  predictions <- list()
  folds <- list()
  count <- 0L
  for (repeat_id in seq_len(repeats)) {
    outer <- split_function(d, 5L, seed_base + repeat_id)
    for (fold in seq_len(5L)) {
      test_idx <- which(outer == fold)
      if (!length(test_idx)) next
      train_idx <- which(outer != fold)
      train <- d[train_idx, , drop = FALSE]
      test <- d[test_idx, , drop = FALSE]
      for (env in unique(as.character(test$environment))) {
        rows <- which(as.character(test$environment) == env)
        mean_pred <- rep(mean(train$y[as.character(train$environment) == env]), length(rows))
        count <- count + 1L
        predictions[[count]] <- data.frame(scenario = scenario, repeat_id = repeat_id,
          environment = env, line_id = test$line_id[rows], observed = test$y[rows],
          model = "Environment mean", predicted = mean_pred)
      }
      if (cv2) {
        env_mean <- tapply(train$y, as.character(train$environment), mean)
        grand <- mean(train$y)
        line_pred <- vapply(seq_len(nrow(test)), function(j) {
          known <- which(train$line_id == test$line_id[j])
          if (!length(known)) stop("CV2 split did not retain a known phenotype for a held-out line.")
          adjusted <- train$y[known] - (env_mean[as.character(train$environment[known])] - grand)
          mean(adjusted) + (env_mean[as.character(test$environment[j])] - grand)
        }, numeric(1))
        count <- count + 1L
        predictions[[count]] <- data.frame(scenario = scenario, repeat_id = repeat_id,
          environment = as.character(test$environment), line_id = test$line_id,
          observed = test$y, model = "Known-line phenotype + environment shift", predicted = line_pred)
      }
      for (model in c("main", "gxe")) {
        tuned <- tune_lambda(d, train_idx, G, model, split_function,
                             seed = seed_base + repeat_id * 100L + fold,
                             environment_fixed = TRUE)
        predicted <- fit_marker_ridge(train, test, G, tuned$lambda, model,
                                      env_levels = environment_levels, environment_fixed = TRUE)
        count <- count + 1L
        predictions[[count]] <- data.frame(scenario = scenario, repeat_id = repeat_id,
          environment = as.character(test$environment), line_id = test$line_id,
          observed = test$y, model = if (model == "main") "Marker main effects" else "Marker main + GxE",
          predicted = predicted)
        folds[[length(folds) + 1L]] <- data.frame(scenario = scenario, repeat_id = repeat_id,
          fold = fold, model = model, n_train = nrow(train), n_test = nrow(test),
          selected_lambda = tuned$lambda, inner_cv_rmse = tuned$rmse)
      }
    }
    message(scenario, " repeat ", repeat_id, "/", repeats, " complete")
  }
  list(predictions = do.call(rbind, predictions), folds = do.call(rbind, folds))
}

paired_lines <- names(which(table(phenotype$line_id) >= 2L))
cv1 <- run_nested(phenotype, "Grouped CV1 - new genetic clusters", split_by_cluster,
                  repeats = 3L, seed_base = 20267000L, cv2 = FALSE)
cv2_data <- phenotype[phenotype$line_id %in% paired_lines, , drop = FALSE]
cv2 <- run_nested(cv2_data, "CV2 - known line in other environments", split_by_cv2,
                  repeats = 3L, seed_base = 20268000L, cv2 = TRUE)
predictions <- rbind(cv1$predictions, cv2$predictions)
folds <- rbind(cv1$folds, cv2$folds)
repeat_scores <- score_groups(predictions)
validation_summary <- summarize_scores(repeat_scores)
write.csv(validation_summary, "reports/all_site_validation_summary.csv", row.names = FALSE)
write.csv(repeat_scores, "reports/all_site_validation_repeats.csv", row.names = FALSE)
write.csv(folds, "reports/all_site_validation_folds.csv", row.names = FALSE)
write.csv(predictions, gzfile("reports/all_site_validation_predictions.csv.gz"), row.names = FALSE)

# Leave each location out completely, including all years recorded at that site.
# No held-out site-year fixed effect or marker-by-site coefficient is estimated.
location_holdout <- list()
location_oof <- list()
for (loc_i in seq_along(location_levels)) {
  location <- location_levels[loc_i]
  train <- phenotype[phenotype$location != location, , drop = FALSE]
  test <- phenotype[phenotype$location == location, , drop = FALSE]
  if (length(unique(train$line_id)) < inner_folds) stop("Too few training lines for location holdout.")
  split_train <- function(d, k, seed) assign_group_folds(line_cluster[d$line_id], k, seed)
  idx <- seq_len(nrow(train))
  tune <- tune_lambda(train, idx, G, "main", split_train,
                      seed = 20269000L + loc_i, environment_fixed = FALSE,
                      env_levels = environment_levels)
  predicted <- fit_marker_ridge(train, test, G, tune$lambda, "main",
                                env_levels = environment_levels, environment_fixed = FALSE)
  mean_pred <- rep(mean(train$y), nrow(test))
  model_score <- score_prediction(test$y, predicted)
  baseline_score <- score_prediction(test$y, mean_pred)
  location_oof[[loc_i]] <- rbind(
    data.frame(location = location, line_id = test$line_id, environment = test$environment,
      observed = test$y, model = "Other-location mean", predicted = mean_pred),
    data.frame(location = location, line_id = test$line_id, environment = test$environment,
      observed = test$y, model = "Marker main effects; no site effect", predicted = predicted)
  )
  location_holdout[[loc_i]] <- rbind(
    data.frame(location = location, model = "Other-location mean", n = nrow(test),
      rmse = baseline_score["rmse"], pearson_r = NA_real_, spearman_rho = NA_real_,
      calibration_slope = NA_real_, selected_lambda = NA_real_),
    data.frame(location = location, model = "Marker main effects; no site effect", n = nrow(test),
      rmse = model_score["rmse"], pearson_r = model_score["pearson_r"],
      spearman_rho = model_score["spearman_rho"], calibration_slope = model_score["calibration_slope"],
      selected_lambda = tune$lambda)
  )
  message("Leave-location-out complete: ", location, " (", nrow(test), " records)")
}
location_summary <- do.call(rbind, location_holdout)
location_predictions <- do.call(rbind, location_oof)
write.csv(location_summary, "reports/leave_one_location_out.csv", row.names = FALSE)
write.csv(location_predictions, gzfile("reports/leave_one_location_out_predictions.csv.gz"), row.names = FALSE)

env_counts <- as.data.frame(table(phenotype$environment), stringsAsFactors = FALSE)
names(env_counts) <- c("environment", "records")
env_counts$location <- as.character(sub("_[0-9]{4}-[0-9]{2}$", "", env_counts$environment))
env_counts$unique_lines <- vapply(env_counts$environment, function(e) {
  length(unique(phenotype$line_id[as.character(phenotype$environment) == e]))
}, integer(1))
env_counts$source_trait <- ifelse(env_counts$location == "ICRISAT", "YPPlnt", "YPP")
write.csv(env_counts, "reports/all_site_environment_counts.csv", row.names = FALSE)
write.csv(data.frame(
  source_workbook = basename(phenotype_file), common_trait = "yield per plant",
  environments = length(environment_levels), locations = length(location_levels),
  phenotype_records_after_aggregation = phenotype_before_match,
  phenotype_records_before_marker_match = phenotype_before_match,
  marker_matched_records = nrow(phenotype), unique_lines = length(unique(phenotype$line_id)),
  lines_in_two_or_more_environments = length(paired_lines),
  duplicate_rows_averaged = duplicate_rows_collapsed,
  chromosomes = length(chromosomes), source_marker_step = source_marker_step,
  markers_full_qc = ncol(G_all), markers_compact = ncol(G),
  marker_missing_before_imputation = missing_before_imputation,
  min_maf = min_maf, max_marker_missing = max_missing,
  pcs_for_clusters = pc_n, genetic_clusters = cluster_n, outer_folds = 5L,
  core_repeats = 3L, inner_folds = inner_folds,
  lambda_grid = paste(lambda_grid, collapse = ";"),
  lambda_boundary_share = boundary_share(folds$selected_lambda, lambda_grid),
  seed = 20261021L
), "reports/all_site_qc.csv", row.names = FALSE)

png("figures/all_site_prediction_validation.png", width = 2000, height = 1200, res = 160)
oldpar <- par(mfrow = c(2, 2), mar = c(9, 4.5, 3.5, 1))
plot_cv <- function(scenario, title, models) {
  d <- validation_summary[validation_summary$scenario == scenario & validation_summary$model %in% models, ]
  env_order <- environment_levels[environment_levels %in% d$environment]
  x <- seq_len(length(env_order))
  colors <- c("#2D6A8A", "#D55E00", "#009E73")
  y <- range(c(d$rmse_mean - d$rmse_sd, d$rmse_mean + d$rmse_sd), finite = TRUE)
  plot(NA, xlim = c(0.5, length(x) + 0.5), ylim = y,
       xaxt = "n", xlab = "", ylab = "RMSE (yield-per-plant scale)", main = title)
  for (j in seq_along(models)) {
    d_model <- d[d$model == models[j], , drop = FALSE]
    dd <- d_model[match(env_order, d_model$environment), , drop = FALSE]
    points(x + (j - (length(models) + 1) / 2) * 0.18, dd$rmse_mean,
           pch = 19, col = colors[j])
    show_error <- is.finite(dd$rmse_sd) & dd$rmse_sd > 0
    if (any(show_error)) arrows(x[show_error] + (j - (length(models) + 1) / 2) * 0.18,
           dd$rmse_mean[show_error] - dd$rmse_sd[show_error],
           x[show_error] + (j - (length(models) + 1) / 2) * 0.18,
           dd$rmse_mean[show_error] + dd$rmse_sd[show_error],
           angle = 90, code = 3, length = 0.035, col = colors[j])
  }
  axis(1, at = x, labels = env_order, las = 2, cex.axis = 0.7)
  legend("topright", legend = models, col = colors[seq_along(models)], pch = 19, cex = 0.7, bty = "n")
}
plot_cv("Grouped CV1 - new genetic clusters", "Grouped CV1 by site-year",
        c("Environment mean", "Marker main effects", "Marker main + GxE"))
plot_cv("CV2 - known line in other environments", "Known-line CV2 by site-year",
        c("Known-line phenotype + environment shift", "Marker main effects", "Marker main + GxE"))
locs <- unique(location_summary$location)
base <- location_summary[location_summary$model == "Other-location mean", ]
marker <- location_summary[location_summary$model == "Marker main effects; no site effect", ]
x <- seq_along(locs)
ylim <- range(c(base$rmse, marker$rmse), finite = TRUE)
plot(x - 0.1, base$rmse, pch = 19, col = "#999999", xaxt = "n", xlim = range(x) + c(-0.5, 0.5),
     ylim = ylim, xlab = "", ylab = "RMSE (yield-per-plant scale)", main = "Leave-one-location-out")
points(x + 0.1, marker$rmse, pch = 19, col = "#2D6A8A")
axis(1, at = x, labels = locs, las = 2, cex.axis = 0.8)
legend("topright", legend = c("Other-location mean", "Marker main effects"),
       col = c("#999999", "#2D6A8A"), pch = 19, bty = "n", cex = 0.7)
barplot(pca$sdev[1:5]^2 / sum(pca$sdev^2) * 100, names.arg = paste0("PC", 1:5),
        col = "#2D6A8A", ylab = "Variance explained (%)", main = "Marker PCA - first five PCs")
par(oldpar)
dev.off()

writeLines(c(
  "# All-site yield-per-plant validation notes", "",
  paste0("Trait harmonization: YPP in all sheets except ICRISAT, where the supplied column is YPPlnt."),
  paste0("The workbook supplied ", length(environment_levels), " site-year environments across ", length(location_levels), " locations."),
  paste0("", duplicate_rows_collapsed, " duplicate normalized accession-by-site-year rows were averaged before validation."),
  "Grouped CV1 holds all records from marker-derived genotype clusters in one fold; CV2 holds out line-by-site-year cells while retaining that line in another environment.",
  "The GxE model uses marker main effects plus marker-by-site-year interactions and fixed effects for observed site-year environments.",
  "Leave-one-location-out holds out every year from one location. Its model has no target-site effect or target-site interaction, because neither can be estimated without target-site phenotypes.",
  "Location-specific scale and management differences are therefore a substantial challenge; interpret both RMSE and ranking/calibration metrics.",
  "Three repeated split sets summarize partition variability, not population confidence intervals. Related groups are approximate PCA/k-means clusters, not recorded pedigree families.",
  paste0("The all-site marker model uses every ", post_qc_stride, "th QC-passed marker (", ncol(G),
         " markers); the full sampled QC set has ", ncol(G_all), " markers."),
  "Yield columns are named differently by site in the source workbook; site-specific trait mapping and counts are reported in all_site_environment_counts.csv.",
  paste0("The ridge penalty grid was 1 to 100,000; ", round(100 * boundary_share(folds$selected_lambda, lambda_grid), 1),
         "% of nested CV tuned folds selected a grid endpoint. Interpret tuning accordingly."),
  "CicerSeq source files are downloaded locally and are not redistributed in this repository."
), "reports/all_site_validation_notes.md")

cat("All-site prediction complete. See reports/all_site_validation_summary.csv, reports/leave_one_location_out.csv, compressed prediction files, and figures/all_site_prediction_validation.png.\n")
