# Rebuild derived validation scores from saved out-of-fold predictions.
source("R/gxe_utils.R")
obj <- readRDS("data/processed/analysis_objects.rds")
phenotype <- obj$phenotype
G <- obj$G_small
line_cluster <- obj$line_cluster
environment_levels <- obj$environment_levels
lambda_grid <- default_lambda_grid
inner_folds <- 3L

tune_lambda <- function(d, train_idx, split_function, seed) {
  train <- d[train_idx, , drop = FALSE]
  inner <- split_function(train, inner_folds, seed)
  errors <- rep(NA_real_, length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    pred <- rep(NA_real_, nrow(train))
    for (fold in seq_len(inner_folds)) {
      test <- which(inner == fold)
      if (!length(test)) next
      fit <- which(inner != fold)
      pred[test] <- fit_marker_ridge(train[fit, , drop = FALSE], train[test, , drop = FALSE],
        G, lambda_grid[i], "main", env_levels = environment_levels, environment_fixed = FALSE)
    }
    ok <- is.finite(pred)
    if (any(ok)) errors[i] <- sqrt(mean((train$y[ok] - pred[ok])^2))
  }
  best <- which.min(errors)
  list(lambda = lambda_grid[best], rmse = errors[best])
}

# Save the two leave-one-season-out prediction vectors, then rebuild their table.
paired <- names(which(table(phenotype$line_id) >= 2L))
loyo_parts <- list()
for (direction in seq_along(environment_levels)) {
  train_env <- environment_levels[direction]
  test_env <- environment_levels[3L - direction]
  train_idx <- which(as.character(phenotype$environment) == train_env & phenotype$line_id %in% paired)
  test_idx <- which(as.character(phenotype$environment) == test_env & phenotype$line_id %in% paired)
  train <- phenotype[train_idx, , drop = FALSE]
  test <- phenotype[test_idx, , drop = FALSE]
  split_groups <- function(d, k, seed) assign_group_folds(line_cluster[d$line_id], k, seed)
  tuned <- tune_lambda(phenotype, train_idx, split_groups, 20264000L + direction)
  marker <- fit_marker_ridge(train, test, G, tuned$lambda, "main",
    env_levels = environment_levels, environment_fixed = FALSE)
  loyo_parts[[direction]] <- rbind(
    data.frame(direction = paste(train_env, "to", test_env), line_id = test$line_id,
      observed = test$y, model = "Training-season mean", predicted = mean(train$y),
      selected_lambda = NA_real_),
    data.frame(direction = paste(train_env, "to", test_env), line_id = test$line_id,
      observed = test$y, model = "Marker main effects", predicted = marker,
      selected_lambda = tuned$lambda)
  )
}
loyo_predictions <- do.call(rbind, loyo_parts)
write.csv(loyo_predictions, "reports/leave_one_season_out_predictions.csv", row.names = FALSE)
loyo_groups <- split(loyo_predictions, interaction(loyo_predictions$direction, loyo_predictions$model, drop = TRUE))
loyo_summary <- do.call(rbind, lapply(loyo_groups, function(d) {
  s <- score_prediction(d$observed, d$predicted)
  data.frame(direction = d$direction[1], model = d$model[1], n = s["n"], rmse = s["rmse"],
    pearson_r = s["pearson_r"], spearman_rho = s["spearman_rho"],
    calibration_slope = s["calibration_slope"], selected_lambda = d$selected_lambda[1])
}))
write.csv(loyo_summary, "reports/leave_one_season_out.csv", row.names = FALSE)

# Verify the two-season repeat scores against the saved accession-level predictions.
core <- read.csv(gzfile("reports/multienvironment_validation_predictions.csv.gz"), stringsAsFactors = FALSE)
core_groups <- split(core, interaction(core$repeat_id, core$scenario, core$environment, core$model, drop = TRUE))
core_scores <- do.call(rbind, lapply(core_groups, function(d) {
  s <- score_prediction(d$observed, d$predicted)
  if (d$model[1] == "Environment mean") s[c("pearson_r", "spearman_rho", "calibration_slope")] <- NA_real_
  data.frame(scenario = d$scenario[1], repeat_id = d$repeat_id[1], environment = d$environment[1],
    model = d$model[1], n = s["n"], rmse = s["rmse"], pearson_r = s["pearson_r"],
    spearman_rho = s["spearman_rho"], calibration_slope = s["calibration_slope"])
}))
saved_core <- read.csv("reports/multienvironment_validation_repeats.csv", stringsAsFactors = FALSE)
core_key <- c("scenario", "repeat_id", "environment", "model")
core_check <- merge(core_scores, saved_core, by = core_key, suffixes = c("_recalc", "_saved"))
core_fields <- c("n", "rmse", "pearson_r", "spearman_rho", "calibration_slope")
for (field in core_fields) {
  delta <- abs(core_check[[paste0(field, "_recalc")]] - core_check[[paste0(field, "_saved")]])
  if (any(delta > 1e-8, na.rm = TRUE)) stop("Two-season saved-prediction check failed for ", field)
}

saved_core_summary <- read.csv("reports/multienvironment_validation_summary.csv", stringsAsFactors = FALSE)
core_summary_groups <- split(core_scores, interaction(core_scores$scenario, core_scores$environment,
  core_scores$model, drop = TRUE))
core_summary_recalc <- do.call(rbind, lapply(core_summary_groups, function(d) {
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  sd_na <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
  data.frame(scenario = d$scenario[1], environment = d$environment[1], model = d$model[1],
    n = mean(d$n), rmse_mean = mean(d$rmse), rmse_sd_between_repeats = sd_na(d$rmse),
    pearson_r_mean = mean_na(d$pearson_r), pearson_r_sd = sd_na(d$pearson_r),
    spearman_rho_mean = mean_na(d$spearman_rho), spearman_rho_sd = sd_na(d$spearman_rho),
    calibration_slope_mean = mean_na(d$calibration_slope), calibration_slope_sd = sd_na(d$calibration_slope))
}))
core_summary_check <- merge(core_summary_recalc, saved_core_summary,
  by = c("scenario", "environment", "model"), suffixes = c("_recalc", "_saved"))
for (field in setdiff(names(core_summary_recalc), c("scenario", "environment", "model"))) {
  delta <- abs(core_summary_check[[paste0(field, "_recalc")]] - core_summary_check[[paste0(field, "_saved")]])
  if (any(delta > 1e-8, na.rm = TRUE)) stop("Two-season summary check failed for ", field)
}

# Recompute all-site CV metrics and location-holdout scores from the saved OOF rows.
all_site <- read.csv(gzfile("reports/all_site_validation_predictions.csv.gz"), stringsAsFactors = FALSE)
all_groups <- split(all_site, interaction(all_site$repeat_id, all_site$scenario,
  all_site$environment, all_site$model, drop = TRUE))
all_scores <- do.call(rbind, lapply(all_groups, function(d) {
  s <- score_prediction(d$observed, d$predicted)
  if (d$model[1] == "Environment mean") s[c("pearson_r", "spearman_rho", "calibration_slope")] <- NA_real_
  data.frame(scenario = d$scenario[1], repeat_id = d$repeat_id[1], environment = d$environment[1],
    model = d$model[1], n = s["n"], rmse = s["rmse"], pearson_r = s["pearson_r"],
    spearman_rho = s["spearman_rho"], calibration_slope = s["calibration_slope"])
}))
saved_all <- read.csv("reports/all_site_validation_repeats.csv", stringsAsFactors = FALSE)
all_key <- c("scenario", "repeat_id", "environment", "model")
all_check <- merge(all_scores, saved_all, by = all_key, suffixes = c("_recalc", "_saved"))
for (field in core_fields) {
  delta <- abs(all_check[[paste0(field, "_recalc")]] - all_check[[paste0(field, "_saved")]])
  if (any(delta > 1e-8, na.rm = TRUE)) stop("All-site saved-prediction check failed for ", field)
}
all_summary <- read.csv("reports/all_site_validation_summary.csv", stringsAsFactors = FALSE)
all_summary_groups <- split(all_scores, interaction(all_scores$scenario, all_scores$environment,
  all_scores$model, drop = TRUE))
all_summary_recalc <- do.call(rbind, lapply(all_summary_groups, function(d) {
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  sd_na <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
  data.frame(scenario = d$scenario[1], environment = d$environment[1], model = d$model[1],
    n_mean = mean(d$n), rmse_mean = mean(d$rmse), rmse_sd = sd_na(d$rmse),
    pearson_r_mean = mean_na(d$pearson_r), pearson_r_sd = sd_na(d$pearson_r),
    spearman_rho_mean = mean_na(d$spearman_rho), spearman_rho_sd = sd_na(d$spearman_rho),
    calibration_slope_mean = mean_na(d$calibration_slope), calibration_slope_sd = sd_na(d$calibration_slope))
}))
all_summary_check <- merge(all_summary_recalc, all_summary,
  by = c("scenario", "environment", "model"), suffixes = c("_recalc", "_saved"))
for (field in setdiff(names(all_summary_recalc), c("scenario", "environment", "model"))) {
  delta <- abs(all_summary_check[[paste0(field, "_recalc")]] - all_summary_check[[paste0(field, "_saved")]])
  if (any(delta > 1e-8, na.rm = TRUE)) stop("All-site summary check failed for ", field)
}
location_oof <- read.csv(gzfile("reports/leave_one_location_out_predictions.csv.gz"), stringsAsFactors = FALSE)
location_groups <- split(location_oof, interaction(location_oof$location, location_oof$model, drop = TRUE))
location_recalc <- do.call(rbind, lapply(location_groups, function(d) {
  s <- score_prediction(d$observed, d$predicted)
  if (d$model[1] == "Other-location mean") s[c("pearson_r", "spearman_rho", "calibration_slope")] <- NA_real_
  data.frame(location = d$location[1], model = d$model[1], n = s["n"], rmse = s["rmse"],
    pearson_r = s["pearson_r"], spearman_rho = s["spearman_rho"],
    calibration_slope = s["calibration_slope"])
}))
saved_location <- read.csv("reports/leave_one_location_out.csv", stringsAsFactors = FALSE)
loc_check <- merge(location_recalc, saved_location, by = c("location", "model"), suffixes = c("_recalc", "_saved"))
for (field in core_fields) {
  delta <- abs(loc_check[[paste0(field, "_recalc")]] - loc_check[[paste0(field, "_saved")]])
  if (any(delta > 1e-8, na.rm = TRUE)) stop("Location saved-prediction check failed for ", field)
}

density_scores <- read.csv("reports/marker_density_sensitivity_repeats.csv", stringsAsFactors = FALSE)
density_saved <- read.csv("reports/marker_density_sensitivity.csv", stringsAsFactors = FALSE)
density_groups <- split(density_scores, interaction(density_scores$density, density_scores$environment, drop = TRUE))
density_recalc <- do.call(rbind, lapply(density_groups, function(d) data.frame(
  density = d$density[1], environment = d$environment[1],
  marker_count = if (d$density[1] == "sparse_1_in_2500") ncol(obj$G_small) else ncol(obj$G_full),
  rmse_mean = mean(d$rmse), rmse_sd = sd(d$rmse), pearson_r_mean = mean(d$pearson_r),
  pearson_r_sd = sd(d$pearson_r))))
density_check <- merge(density_recalc, density_saved,
  by = c("density", "environment", "marker_count"), suffixes = c("_recalc", "_saved"))
for (field in c("rmse_mean", "rmse_sd", "pearson_r_mean", "pearson_r_sd")) {
  delta <- abs(density_check[[paste0(field, "_recalc")]] - density_check[[paste0(field, "_saved")]])
  if (any(delta > 1e-8, na.rm = TRUE)) stop("Marker-density summary check failed for ", field)
}

writeLines(c(
  "Saved out-of-fold prediction checks passed.",
  paste("Two-season repeat groups checked:", nrow(core_check)),
  paste("Two-season summary groups checked:", nrow(core_summary_check)),
  paste("All-site repeat groups checked:", nrow(all_check)),
  paste("All-site summary groups checked:", nrow(all_summary_check)),
  paste("Leave-one-location groups checked:", nrow(loc_check)),
  paste("Marker-density summary groups checked:", nrow(density_check)),
  paste("Leave-one-season prediction rows exported:", nrow(loyo_predictions)),
  "Leave-one-season summary rebuilt from its saved prediction rows."
), "reports/prediction_verification.txt")
cat("Saved-prediction metrics verified and leave-one-season predictions exported.\n")
