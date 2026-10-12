## Runs the DualPulse pipeline from the repository root, each stage in its own R process.
##   Rscript run.R              all nine stages, in order
##   Rscript run.R 5 6          stages 5 and 6
##   Rscript run.R 6:46,47      sections 46 and 47 of stage 6 (sections are named in each stage's header)
## Paths are set in config_local.R (see config_template.R).

args   <- commandArgs(trailingOnly = TRUE)
stages <- sort(list.files("scripts", pattern = "^[0-9]{2}_.*\\.R$"))
if (!length(stages)) stop("Run from the repository root (the folder holding run.R).")
if (!length(args)) args <- as.character(seq_along(stages))
rscript <- file.path(R.home("bin"), "Rscript")
t0 <- Sys.time()
for (a in args) {
  k    <- as.integer(sub(":.*", "", a))
  only <- if (grepl(":", a)) sub(".*:", "", a) else ""
  if (is.na(k) || k < 1 || k > length(stages)) stop("no stage ", a)
  message(sprintf("\n#### Stage %d: %s%s", k, stages[k], if (nzchar(only)) paste0(" (sections ", only, ")") else ""))
  Sys.setenv(DP_ONLY = only)
  status <- system2(rscript, shQuote(file.path("scripts", stages[k])))
  if (status != 0) stop(sprintf("stage %d stopped with exit status %d", k, status))
}
message(sprintf("\nDone in %.0f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
