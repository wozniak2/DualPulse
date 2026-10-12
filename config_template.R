## Copy this file to config_local.R (not tracked by git) and set your own paths.
## Variables left out keep the defaults of R/setup.R.

ROOT        <- "path/to/project"                   # holds data/ (inputs); output/ and manuscript/ are written here
PEKA_DIR    <- "path/to/PEKA"                      # output/tables/peka_2025_MM.rds and gtfs/extracted_202510/
PEKA_RAW    <- "path/to/PEKA_raw"                  # MM_csv_2025.7z monthly archives (ingest section only)
SEVENZIP    <- "C:/Program Files/7-Zip/7z.exe"     # 7-Zip (ingest section only)
CITY_DIR    <- "path/to/folder"                    # holds poz.gpkg, the city boundary
RCN_GPKG    <- "path/to/rcn_poznan_metro.gpkg"     # property transactions (section 19)
OSM_EXTRACT <- "path/to/poznan_agglomeration.osm"  # OpenStreetMap extract (section 14; cached after the first run)
WORKERS     <- 3L                                  # parallel R processes for the ingest and section 40
