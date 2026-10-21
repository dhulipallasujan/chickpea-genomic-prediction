source("R/gxe_utils.R")

# Group CV keeps each genetic cluster entirely in one fold.
line <- paste0("L", seq_len(60L))
cluster <- rep(paste0("C", seq_len(10L)), each = 6L)
fold <- assign_group_folds(cluster, k = 5L, seed = 42L)
stopifnot(all(fold %in% 1:5))
stopifnot(all(vapply(split(fold, cluster), function(z) length(unique(z)) == 1L, logical(1))))

# CV2 cell holds are valid only if that line's other environment remains in train.
cv2_data <- expand.grid(line_id = paste0("L", seq_len(40L)),
                        environment = c("E1", "E2"), stringsAsFactors = FALSE)
cv2_fold <- assign_cv2_folds(cv2_data$line_id, cv2_data$environment, k = 5L, seed = 7L)
stopifnot(all(vapply(split(cv2_fold, cv2_data$line_id), function(z) length(unique(z)) == 2L, logical(1))))
for (f in 1:5) {
  held <- cv2_data[cv2_fold == f, , drop = FALSE]
  training <- cv2_data[cv2_fold != f, , drop = FALSE]
  stopifnot(all(held$line_id %in% training$line_id))
  for (i in seq_len(nrow(held))) {
    stopifnot(any(training$line_id == held$line_id[i] & training$environment != held$environment[i]))
  }
}

# Marker main-effect and marker-by-environment ridge fits return finite predictions.
set.seed(10)
G <- matrix(rnorm(40L * 12L), nrow = 40L,
            dimnames = list(paste0("L", seq_len(40L)), paste0("S", seq_len(12L))))
dat <- cv2_data
dat$y <- as.numeric(G[match(dat$line_id, rownames(G)), 1L]) +
  ifelse(dat$environment == "E2", 1, 0) + rnorm(nrow(dat), sd = 0.2)
pred_main <- fit_marker_ridge(dat[1:60, ], dat[61:80, ], G, lambda = 10,
                              model = "main", env_levels = c("E1", "E2"))
pred_gxe <- fit_marker_ridge(dat[1:60, ], dat[61:80, ], G, lambda = 10,
                             model = "gxe", env_levels = c("E1", "E2"))
stopifnot(length(pred_main) == 20L, all(is.finite(pred_main)))
stopifnot(length(pred_gxe) == 20L, all(is.finite(pred_gxe)))

cat("Validation split and ridge helper checks passed.\n")
