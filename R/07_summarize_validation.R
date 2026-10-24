# Rebuild repeat-level validation metrics from saved out-of-fold predictions.
# This keeps reporting independent from the relatively slow model-fitting stage.
source("R/gxe_utils.R")
prediction_file <- "reports/multienvironment_validation_predictions.csv.gz"
if (!file.exists(prediction_file)) stop("Run R/06_multienvironment_validation.R first.")
predictions <- read.csv(prediction_file, stringsAsFactors = FALSE)
groups <- split(predictions, interaction(predictions$repeat_id, predictions$scenario,
                                         predictions$environment, predictions$model, drop = TRUE))
repeat_scores <- do.call(rbind, lapply(groups, function(d) {
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
summary_groups <- split(repeat_scores, interaction(repeat_scores$scenario,
                                                   repeat_scores$environment,
                                                   repeat_scores$model, drop = TRUE))
validation_summary <- do.call(rbind, lapply(summary_groups, function(d) {
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  sd_or_na <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else stats::sd(x, na.rm = TRUE)
  data.frame(scenario = d$scenario[1], environment = d$environment[1], model = d$model[1],
             n = mean(d$n), rmse_mean = mean(d$rmse), rmse_sd_between_repeats = stats::sd(d$rmse),
             pearson_r_mean = mean_or_na(d$pearson_r), pearson_r_sd = sd_or_na(d$pearson_r),
             spearman_rho_mean = mean_or_na(d$spearman_rho), spearman_rho_sd = sd_or_na(d$spearman_rho),
             calibration_slope_mean = mean_or_na(d$calibration_slope),
             calibration_slope_sd = sd_or_na(d$calibration_slope), stringsAsFactors = FALSE)
}))
write.csv(repeat_scores, "reports/multienvironment_validation_repeats.csv", row.names = FALSE)
write.csv(validation_summary, "reports/multienvironment_validation_summary.csv", row.names = FALSE)
cat("Validation summaries rebuilt from saved out-of-fold predictions.\n")
