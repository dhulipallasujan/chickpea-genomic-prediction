source("R/gxe_utils.R")

# Penalty grid: log-spaced, increasing, reaches 1e5.
stopifnot(length(default_lambda_grid) == 11L, !is.unsorted(default_lambda_grid),
          min(default_lambda_grid) == 1, abs(max(default_lambda_grid) - 1e5) < 1e-6)

# Edge detection.
stopifnot(boundary_share(c(1000, 1000, 10), c(1, 10, 1000)) == 2/3)
stopifnot(boundary_share(c(10, 100), c(1, 10, 100, 1000)) == 0.5)
stopifnot(is.na(boundary_share(numeric(0), c(1, 10))))

# Scoring: perfect predictions give zero RMSE and unit slope.
s <- score_prediction(c(1, 2, 3, 4), c(1, 2, 3, 4))
stopifnot(abs(s[["rmse"]]) < 1e-12, abs(s[["calibration_slope"]] - 1) < 1e-12)

# Ridge behaviour on synthetic data with two environments.
set.seed(3)
G <- matrix(rnorm(40L * 12L), nrow = 40L,
            dimnames = list(paste0("L", seq_len(40L)), paste0("S", seq_len(12L))))
dat <- expand.grid(line_id = paste0("L", seq_len(40L)), environment = c("E1", "E2"), stringsAsFactors = FALSE)
dat$y <- 2 * G[match(dat$line_id, rownames(G)), 1L] + ifelse(dat$environment == "E2", 5, 0) + rnorm(nrow(dat), sd = 0.2)
train <- dat[c(1:30, 41:70), ]
test <- dat[c(31:40, 71:80), ]

# A huge penalty shrinks marker effects to ~0, so predictions collapse to the training environment mean.
shrunk <- fit_marker_ridge(train, test, G, lambda = 1e12, model = "main", env_levels = c("E1", "E2"))
env_mean <- ifelse(test$environment == "E1", mean(train$y[train$environment == "E1"]), mean(train$y[train$environment == "E2"]))
stopifnot(max(abs(shrunk - env_mean)) < 1e-3)

# A small penalty recovers the real marker signal.
fitted <- fit_marker_ridge(train, test, G, lambda = 0.01, model = "main", env_levels = c("E1", "E2"))
stopifnot(stats::cor(fitted, test$y) > 0.9)

# Environment fixed effect is not penalized: a constant shift between seasons is recovered.
stopifnot(abs(mean(fitted[test$environment == "E2"]) - mean(fitted[test$environment == "E1"]) - 5) < 1)

cat("Ridge behaviour checks passed.\n")
