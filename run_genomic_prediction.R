# Run the full analysis from the repository root.
scripts <- c("R/05_genomic_prediction.R", "R/06_multienvironment_validation.R",
             "R/07_summarize_validation.R", "R/08_figures_and_tables.R",
             "R/09_rrblup_crosscheck.R")
if (!all(file.exists(scripts))) stop("Run this file from the repository root; required analysis scripts are missing.")
executable <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
for (script in scripts) {
  cat("\n==> ", script, "\n", sep = "")
  status <- system2(file.path(R.home("bin"), executable), script)
  if (status != 0) stop("Genomic prediction workflow stopped in ", script, " (exit ", status, ").")
}
