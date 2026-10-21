# Cross-check the hand-written ridge against rrBLUP::mixed.solve (REML) on the
# first grouped-CV1 repeat. mixed.solve estimates the variance ratio from the
# training data, so its implied penalty also shows whether the tuning grid was wide enough.
source("R/gxe_utils.R")
if (!requireNamespace("rrBLUP", quietly = TRUE)) install.packages("rrBLUP", repos = "https://cloud.r-project.org")
if (!requireNamespace("rrBLUP", quietly = TRUE)) stop("Install rrBLUP to run the cross-check.")
objects_file <- "data/processed/analysis_objects.rds"
folds_file <- "reports/multienvironment_validation_folds.csv"
if (!file.exists(objects_file) || !file.exists(folds_file)) stop("Run R/06_multienvironment_validation.R first.")

obj <- readRDS(objects_file)
phenotype <- obj$phenotype
G <- obj$G_small
line_cluster <- obj$line_cluster
env_levels <- obj$environment_levels
folds <- read.csv(folds_file, stringsAsFactors = FALSE)
tuned <- folds[folds$scenario == "Grouped CV1 - new genetic clusters" & folds$model == "main" & folds$repeat_id == 1L, ]

# Same outer split as repeat 1 of the main grouped CV1 run.
outer <- assign_group_folds(line_cluster[phenotype$line_id], 5L, 20262000L + 1L)
rows <- list()
for (f in seq_len(5L)) {
  test_idx <- which(outer == f)
  train <- phenotype[-test_idx, , drop = FALSE]
  test <- phenotype[test_idx, , drop = FALSE]
  lambda_grid_choice <- tuned$selected_lambda[tuned$fold == f]
  own <- fit_marker_ridge(train, test, G, lambda_grid_choice, "main", env_levels = env_levels, environment_fixed = TRUE)
  x_train <- stats::model.matrix(~ factor(as.character(train$environment), levels = env_levels))
  x_test <- stats::model.matrix(~ factor(as.character(test$environment), levels = env_levels))
  z_train <- make_marker_features(G, train$line_id, train$environment, "main")
  z_test <- make_marker_features(G, test$line_id, test$environment, "main")
  fit <- rrBLUP::mixed.solve(y = train$y, Z = z_train, X = x_train)
  reference <- as.numeric(x_test %*% fit$beta + z_test %*% fit$u)
  rows[[f]] <- data.frame(
    fold = f, n_test = nrow(test),
    own_ridge_lambda = lambda_grid_choice, rrblup_reml_lambda = fit$Ve / fit$Vu,
    rmse_own_ridge = sqrt(mean((test$y - own)^2)), rmse_rrblup = sqrt(mean((test$y - reference)^2)),
    prediction_correlation = stats::cor(own, reference))
}
out <- do.call(rbind, rows)
write.csv(out, "reports/rrblup_crosscheck.csv", row.names = FALSE)
print(out)
cat("Wrote reports/rrblup_crosscheck.csv\n")
