## DualPulse setup: paths and the section runner. Sourced by run.R and at the top of every stage script.
## Run everything from the repository root.
##
## Paths. ROOT holds data/ (inputs), output/ (results) and manuscript/ (figures and tables the paper reads);
## by default it is the repository itself. The passenger counts and the timetable live with the PEKA
## project (PEKA_DIR). Set your own paths in config_local.R (see config_template.R); it is not tracked by git.

REPO <- normalizePath(getwd(), winslash = "/")
if (!file.exists(file.path(REPO, "R", "setup.R"))) stop("Run DualPulse from the repository root (the folder holding run.R).")

ROOT          <- REPO                                     # data/, output/, manuscript/
PEKA_DIR      <- file.path(REPO, "peka")                  # PEKA project: output/tables/peka_2025_MM.rds, gtfs/
PEKA_RAW      <- file.path(REPO, "peka_raw")              # monthly 7z archives of the PEKA visit records (ingest only)
SEVENZIP      <- "7z"                                     # 7-Zip executable (ingest only)
CITY_DIR      <- file.path(REPO, "data", "external")     # poz.gpkg: city boundary for the 200 m grid
RCN_GPKG      <- file.path(REPO, "data", "external", "rcn_poznan_metro.gpkg")   # property transactions (section 19)
OSM_EXTRACT   <- file.path(REPO, "data", "external", "poznan_agglomeration.osm") # OSM extract (section 14, first run only)
PEKA_MONTHS   <- 1:10                                     # January-October 2025
WORKERS       <- 3L                                       # parallel R processes in sections ingest and 40 (~3 GB RAM each)
if (file.exists(file.path(REPO, "config_local.R"))) source(file.path(REPO, "config_local.R"))
if (!exists("PEKA_TABLES"))   PEKA_TABLES   <- file.path(PEKA_DIR, "output", "tables")
if (!exists("GTFS_DIR"))      GTFS_DIR      <- file.path(PEKA_DIR, "gtfs", "extracted_202510")   # October 2025 feed
if (!exists("GTFS_FALLBACK")) GTFS_FALLBACK <- file.path(PEKA_DIR, "gtfs", "extracted")          # later feed, ingest fallback
invisible(lapply(file.path(ROOT, c("output/stage1", "output/stage2", "output/figures", "manuscript/figures", "manuscript/tables")),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

## Each former script runs as one section, in a fresh function frame and with only the base packages
## attached, as if it were run on its own with Rscript. DP_ONLY (comma-separated section ids) restricts a run.
detach_packages <- function() {
  base <- paste0("package:", c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base"))
  for (p in setdiff(grep("^package:", search(), value = TRUE), base)) try(detach(p, character.only = TRUE), silent = TRUE)
}
section <- function(id, title, expr) {
  only <- strsplit(Sys.getenv("DP_ONLY", ""), ",")[[1]]
  if (length(only) && !id %in% only) return(invisible(FALSE))
  detach_packages()
  message(sprintf("\n== [%s] %s", id, title)); t0 <- Sys.time()
  fn <- function() NULL; body(fn) <- substitute(expr); environment(fn) <- globalenv()
  fn()
  message(sprintf("== [%s] done in %.1f min", id, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(TRUE)
}
