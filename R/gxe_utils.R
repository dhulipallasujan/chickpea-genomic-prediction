# Helpers for leakage-aware group CV and ridge genomic prediction.

assign_group_folds <- function(group, k, seed = 1L) {
  group <- as.character(group)
  if (length(group) == 0L || anyNA(group) || any(!nzchar(group))) stop("Groups must be non-empty.")
  if (k < 2L || k > length(unique(group))) stop("Need at least k distinct groups.")
  counts <- table(group)
  set.seed(seed)
  tie_break <- runif(length(counts))
  group_order <- names(counts)[order(-as.integer(counts), tie_break)]
  loads <- integer(k)
  group_fold <- setNames(integer(length(group_order)), group_order)
  for (g in group_order) {
    available <- which(loads == min(loads))
    f <- sample(available, 1L)
    group_fold[[g]] <- f
    loads[f] <- loads[f] + as.integer(counts[[g]])
  }
  unname(group_fold[group])
}

assign_cv2_folds <- function(line_id, environment, k, seed = 1L) {
  line_id <- as.character(line_id)
  environment <- as.character(environment)
  if (length(line_id) != length(environment) || anyNA(line_id) || anyNA(environment)) {
    stop("Line IDs and environments must be complete and have equal length.")
  }
  cell <- paste(line_id, environment, sep = "\r")
  if (anyDuplicated(cell)) stop("Expected one record per line-by-environment cell.")
  if (k < 2L) stop("k must be at least 2.")
  set.seed(seed)
  fold <- integer(length(line_id))
  rows_by_line <- split(seq_along(line_id), line_id)
  for (rows in rows_by_line) {
    if (length(rows) == 2L && length(unique(environment[rows])) == 2L) {
      first <- sample.int(k, 1L)
      second <- sample(setdiff(seq_len(k), first), 1L)
      fold[rows] <- c(first, second)
    } else if (length(rows) > 2L) {
      fold[rows] <- sample(rep_len(seq_len(k), length(rows)))
    }
    # A singleton cell stays in training for all inner splits; it cannot
    # support a known-line CV2 prediction by itself.
  }
  fold
}

make_marker_features <- function(G, line_id, environment, model = c("main", "gxe"), env_levels = NULL) {
  model <- match.arg(model)
  idx <- match(as.character(line_id), rownames(G))
  if (anyNA(idx)) stop("A phenotype line is missing from the marker matrix.")
  x <- G[idx, , drop = FALSE]
  if (model == "main") return(x)
  if (is.null(env_levels)) env_levels <- sort(unique(as.character(environment)))
  if (length(env_levels) < 2L) stop("The marker-by-environment interaction requires at least two environments.")
  env <- factor(as.character(environment), levels = env_levels)
  env_design <- stats::model.matrix(~ env - 1L)[, -1L, drop = FALSE]
  interactions <- do.call(cbind, lapply(seq_len(ncol(env_design)), function(j) {
    sweep(x, 1L, env_design[, j], FUN = "*")
  }))
  cbind(x, interactions)
}

fit_marker_ridge <- function(train, test, G, lambda, model = c("main", "gxe"),
                             env_levels = NULL, environment_fixed = TRUE) {
  model <- match.arg(model)
  if (is.null(env_levels)) env_levels <- sort(unique(as.character(train$environment)))
  x_train <- make_marker_features(G, train$line_id, train$environment, model, env_levels)
  x_test <- make_marker_features(G, test$line_id, test$environment, model, env_levels)
  if (environment_fixed) {
    train_env <- factor(as.character(train$environment), levels = env_levels)
    test_env <- factor(as.character(test$environment), levels = env_levels)
    f_train <- stats::model.matrix(~ train_env)
    f_test <- stats::model.matrix(~ test_env)
    colnames(f_test) <- colnames(f_train)
  } else {
    f_train <- matrix(1, nrow(train), 1L, dimnames = list(NULL, "(Intercept)"))
    f_test <- matrix(1, nrow(test), 1L, dimnames = list(NULL, "(Intercept)"))
  }
  z_train <- cbind(f_train, x_train)
  z_test <- cbind(f_test, x_test)
  penalty <- rep(0, ncol(f_train) + ncol(x_train))
  penalty[seq.int(ncol(f_train) + 1L, length(penalty))] <- lambda
  system <- crossprod(z_train)
  diag(system) <- diag(system) + penalty
  coefficient <- solve(system, crossprod(z_train, train$y))
  as.numeric(z_test %*% coefficient)
}

score_prediction <- function(observed, predicted) {
  pearson <- if (stats::sd(predicted) > 0) stats::cor(observed, predicted) else NA_real_
  spearman <- if (stats::sd(predicted) > 0) suppressWarnings(stats::cor(observed, predicted, method = "spearman")) else NA_real_
  slope <- if (stats::sd(predicted) > 0) unname(stats::coef(stats::lm(observed ~ predicted))[2L]) else NA_real_
  c(n = length(observed), rmse = sqrt(mean((observed - predicted)^2)),
    pearson_r = pearson, spearman_rho = spearman, calibration_slope = slope)
}

# Ridge penalty grid: log-spaced from 1 to 1e5. An earlier grid stopped at 1000 and
# every tuned fold selected that ceiling, so the optimum was not inside the grid.
default_lambda_grid <- 10^seq(0, 5, by = 0.5)

# Share of tuned penalties that sit at either edge of the grid. A high share means
# the grid is too narrow and results depend on where it was cut off.
boundary_share <- function(selected, grid) {
  selected <- selected[is.finite(selected)]
  if (!length(selected)) return(NA_real_)
  mean(selected <= min(grid) | selected >= max(grid))
}
