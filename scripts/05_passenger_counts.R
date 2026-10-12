## DualPulse | stage 5: Passenger counts: ingest, cleaning, flows and same-trip indices
##
## Ten months of PEKA vehicle location and door-count records: ingest, trip-level cleaning, stop-to-stop
## flows by three estimators, expansion to all operated visits, same-trip indices and descriptive checks
## (main text 2.4, supplementary S3).
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   ingest  Ingest of the PEKA monthly archives
##   40      Trip-level cleaning, stop-to-stop flows and interchange visits
##   41      Expansion and stop-level measures
##   42      Same-trip indices
##   22      Counter coverage (October 2025)
##   23      Does usage mirror the timetable? (October 2025)
##   58      Descriptive statistics of the counts
##
## Run from the repository root: Rscript run.R 5   (one section: Rscript run.R 5:58)

if (!exists("section")) source("R/setup.R")

## ── ingest · Ingest of the PEKA monthly archives ───────────────────────────────────────────────────────────
section("ingest", "Ingest of the PEKA monthly archives", {
  ## Fast multi-month ingest of PEKA PUNKTUALNOSC AVL/APC data (optimised
  ## version of 01_ingest_month.R; same output schema, so downstream scripts
  ## that read peka_2025_MM.rds work unchanged).
  ##
  ## Usage: Rscript 09_ingest_months_fast.R [months] [workers]
  ##   months : comma-separated, default 1,2,3,4,5,6,7,8,9
  ##   workers: parallel R processes, default 3 (about 3 GB RAM each)
  ##
  ## Per month:
  ##   1. extracts <MM>_csv_2025.7z with 7-Zip into Downloads/PEKA/<MM>_csv_2025
  ##      (skipped when the daily CSVs are already there);
  ##   2. reads the daily files with fread, NR_LINI forced to character at read
  ##      time (the October run failed on mixed integer/character line numbers);
  ##   3. keeps executed visits, fixes Windows-1250 text by converting unique
  ##      values only, computes delay (midnight rollover) and scheduled hour,
  ##      English weekday names independent of the system locale, APC fields;
  ##   4. adds stop-cluster coordinates and route mode by in-place joins, from
  ##      the October 2025 feed with the 2026 feed as fallback (stop matching by
  ##      SYM_PRZYST = stop_code without trailing digits, as in 01); lines in
  ##      neither feed get the mode of the vehicles that run them;
  ##   5. saves output/tables/peka_2025_MM.rds (full visit records) and
  ##      output/tables/peka_2025_MM_cells.rds: one row per date x stop cluster
  ##      x line x scheduled hour with visits, counter-equipped visits, counted
  ##      boardings and alightings, delay sums and on-time counts, dwell sum --
  ##      enough for volume, timing, day-type, seasonal and reliability measures
  ##      without reloading the visit records.
  ## Month 10 (already ingested by 01) gets only its cells table.
  ## Log: output/tables/ingest_log_MM.txt per month; summary: ingest_summary_2025.csv.
  ## Combining months: visits after midnight stay in the file of their duty day
  ## (no visit appears in two files), so cells tables can be stacked; one stray
  ## 2024 date occurs in September, so filter dates to the study period.

  suppressPackageStartupMessages({library(data.table); library(parallel)})
  ## months and parallel workers come from the configuration (R/setup.R); months already ingested are skipped
  months  <- PEKA_MONTHS[!file.exists(file.path(PEKA_TABLES, sprintf("peka_2025_%02d_cells.rds", PEKA_MONTHS)))]
  workers <- WORKERS
  if (!length(months)) return(invisible(message("   PEKA tables present for all months; ingest skipped")))
  cfg <- list(raw_root = PEKA_RAW, sevenzip = SEVENZIP,
              tables = PEKA_TABLES,
              gtfs = c(GTFS_DIR,
                       GTFS_FALLBACK))

  ingest_month <- function(m, cfg) {
    suppressPackageStartupMessages(library(data.table))
    setDTthreads(4)
    mm <- sprintf("%02d", m)
    tab_dir <- cfg$tables
    out_rds <- file.path(tab_dir, sprintf("peka_2025_%s.rds", mm))
    log_f   <- file.path(tab_dir, sprintf("ingest_log_%s.txt", mm))
    logm <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", file = log_f, append = TRUE, sep = "")
    cat("", file = log_f)
    t0 <- Sys.time()

    ## lines absent from both feeds (summer temporary tram lines, replacement buses,
    ## discontinued suburban lines): mode of the vehicles that run them, learnt from
    ## the same month's lines of known mode (fleet numbers do not overlap: trams
    ## 142-955, buses 1001-9903 in August 2025)
    fill_mode <- function(dt, logm) {
      dt[mode == "unknown", mode := NA_character_]
      vm <- dt[mode %in% c("tram", "bus"), .N, by = .(NB_WOZU, mode)][order(-N)][!duplicated(NB_WOZU), .(NB_WOZU, vmode = mode)]
      unk <- dt[is.na(mode) & !is.na(NR_LINI) & NR_LINI != "", .(NR_LINI, NB_WOZU)][vm, vmode := i.vmode, on = "NB_WOZU"]
      lm_ <- unk[!is.na(vmode), .N, by = .(NR_LINI, vmode)][order(-N)][!duplicated(NR_LINI)]
      dt[lm_, mode := fifelse(is.na(mode), i.vmode, mode), on = "NR_LINI"]
      dt[is.na(mode), mode := "unknown"]
      logm(sprintf("mode filled from vehicles for %d lines; mode known for %.2f%% of visits", nrow(lm_), 100 * mean(dt$mode != "unknown")))
    }

    if (m == 10 && file.exists(out_rds)) {
      ## visit records written by 01 are left unchanged; the vehicle-based mode fill
      ## is applied only to the cells table built below
      dt <- readRDS(out_rds); setDT(dt); logm("month 10: existing visit records read (", nrow(dt), " rows)")
      fill_mode(dt, logm)
    } else {
      ## 1. extract
      raw_dir <- file.path(cfg$raw_root, sprintf("%s_csv_2025", mm))
      ndays <- as.integer(format(as.Date(sprintf("2025-%s-01", mm)) + 31 - as.integer(format(as.Date(sprintf("2025-%s-01", mm)) + 31, "%d")), "%d"))
      have <- if (dir.exists(raw_dir)) length(list.files(raw_dir, pattern = "^PUNKTUALNOSC_.*\\.CSV$")) else 0L
      if (have < ndays) {
        arc <- file.path(cfg$raw_root, sprintf("%s_csv_2025.7z", mm))
        st <- system2(cfg$sevenzip, c("x", shQuote(arc, type = "cmd"), paste0("-o", shQuote(raw_dir, type = "cmd")), "-y", "-bso0", "-bsp0"))
        if (st != 0) stop("7-Zip failed for month ", mm)
        logm("extracted ", arc, " in ", round(as.numeric(Sys.time() - t0, units = "secs")), " s")
      }
      files <- sort(list.files(raw_dir, pattern = "^PUNKTUALNOSC_.*\\.CSV$", full.names = TRUE))
      logm(length(files), " daily files (calendar days: ", ndays, ")")

      ## 2. read
      col_names <- c("DT_KARTY","NB_WOZU","NR_WYJKWI","NR_WYJPOJ","NR_KURSOWK","LP_KURSU","ID_KURSU","NR_LINI","LP_PRZYST","ID_PRZYST",
                     "NUM_SLUPKA","NAZ_PRZYST","RJ_KIEDY","DATA","KIEDY","POSTOJ","SYM_PRZYST","SYM_SLUPKA","WESZLO_NIS","WESZLO_WYS",
                     "WYSZLO_NIS","WYSZLO_WYS","OPIS_TABL","NR_IBIS","NE_KIER","ID_KURSU_BUSMAN","WAR_TRASY")
      dt <- rbindlist(lapply(files, function(f) fread(f, header = FALSE, sep = ";", quote = "\"", col.names = col_names, fill = TRUE,
                                                       na.strings = "", encoding = "Latin-1", showProgress = FALSE,
                                                       colClasses = list(character = c("NR_LINI", "NR_KURSOWK", "SYM_PRZYST", "SYM_SLUPKA", "RJ_KIEDY", "KIEDY", "WAR_TRASY")))),
                      use.names = TRUE, fill = TRUE)
      logm("rows loaded: ", nrow(dt), " (", round(as.numeric(Sys.time() - t0, units = "secs")), " s)")

      ## 3. clean and derive
      dt <- dt[!is.na(DATA)]
      logm("executed visits: ", nrow(dt))
      fixenc <- function(v) { u <- unique(v); iconv(u, from = "WINDOWS-1250", to = "UTF-8")[match(v, u)] }
      dt[, `:=`(NAZ_PRZYST = fixenc(NAZ_PRZYST), OPIS_TABL = fixenc(OPIS_TABL))]
      ss <- as.integer(substr(dt$RJ_KIEDY, 1, 2)) * 3600L + as.integer(substr(dt$RJ_KIEDY, 4, 5)) * 60L
      as_ <- as.integer(substr(dt$KIEDY, 1, 2)) * 3600L + as.integer(substr(dt$KIEDY, 4, 5)) * 60L + as.integer(substr(dt$KIEDY, 7, 8))
      dl <- as_ - ss
      dl[!is.na(dl) & dl >  43200L] <- dl[!is.na(dl) & dl >  43200L] - 86400L
      dl[!is.na(dl) & dl < -43200L] <- dl[!is.na(dl) & dl < -43200L] + 86400L
      dt[, `:=`(delay_sec = dl, sched_hour = as.integer(substr(RJ_KIEDY, 1, 2)))]
      rm(ss, as_, dl)
      wd <- data.table(DATA = unique(dt$DATA))
      wd[, weekday := c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")[as.integer(format(as.Date(DATA), "%u"))]]
      dt[wd, weekday := i.weekday, on = "DATA"]
      dt[, `:=`(is_weekend = weekday %in% c("Saturday", "Sunday"), has_apc = !is.na(WESZLO_NIS),
                boardings = WESZLO_NIS + WESZLO_WYS, alightings = WYSZLO_NIS + WYSZLO_WYS)]

      ## 4. GTFS enrichment (October 2025 feed first, 2026 feed as fallback)
      sc <- rbindlist(lapply(seq_along(cfg$gtfs), function(k) {
        s <- fread(file.path(cfg$gtfs[k], "stops.txt"), encoding = "UTF-8")[, prefix := sub("[0-9]+$", "", stop_code)]
        s[, .(stop_lat = mean(stop_lat), stop_lon = mean(stop_lon), stop_name_gtfs = stop_name[1], n_posts = .N, pri = k), by = prefix]
      }))[order(pri)][!duplicated(prefix)]
      dt[sc, `:=`(stop_lat = i.stop_lat, stop_lon = i.stop_lon, stop_name_gtfs = i.stop_name_gtfs, n_posts = i.n_posts), on = c(SYM_PRZYST = "prefix")]
      rl <- rbindlist(lapply(seq_along(cfg$gtfs), function(k) {
        r <- fread(file.path(cfg$gtfs[k], "routes.txt"), encoding = "UTF-8", colClasses = list(character = "route_short_name"))
        r[, .(route_type = route_type[1], pri = k), by = route_short_name]
      }))[order(pri)][!duplicated(route_short_name)]
      rl[, mode := fifelse(route_type == 0, "tram", fifelse(route_type == 3, "bus", "other"))]
      dt[rl, mode := i.mode, on = c(NR_LINI = "route_short_name")]
      logm(sprintf("GTFS stop match %.1f%%, route-mode match %.1f%%", 100 * mean(!is.na(dt$stop_lat)), 100 * mean(!is.na(dt$mode))))
      fill_mode(dt, logm)

      ## same column order as the October file written by 01
      setcolorder(dt, c("NR_LINI", "SYM_PRZYST", setdiff(names(dt), c("NR_LINI", "SYM_PRZYST"))))
      saveRDS(dt, out_rds)
      logm("saved ", out_rds, " (", round(file.size(out_rds) / 2^20), " MB)")
    }

    ## 5. cells table
    cl <- dt[!is.na(SYM_PRZYST) & NR_LINI != "",
             .(visits = .N, counted = sum(has_apc), b = sum(boardings, na.rm = TRUE), a = sum(alightings, na.rm = TRUE),
               delay_sum = sum(delay_sec[!is.na(delay_sec) & abs(delay_sec) < 3600]), delay_n = sum(!is.na(delay_sec) & abs(delay_sec) < 3600),
               ontime_n = sum(!is.na(delay_sec) & delay_sec >= -60 & delay_sec <= 180), dwell_sum = sum(POSTOJ, na.rm = TRUE), mode = mode[1]),
             by = .(DATA, SYM_PRZYST, NR_LINI, sched_hour)]
    cl[, DATA := as.IDate(DATA)]
    saveRDS(cl, file.path(tab_dir, sprintf("peka_2025_%s_cells.rds", mm)))
    logm("cells: ", nrow(cl), " rows; total ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
    data.table(month = mm, visits = nrow(dt), days = uniqueN(dt$DATA), apc_share = mean(dt$has_apc), stop_match = mean(!is.na(dt$stop_lat)),
               mode_match = mean(dt$mode != "unknown"), cells = nrow(cl), minutes = round(as.numeric(Sys.time() - t0, units = "mins"), 1))
  }

  t0 <- Sys.time()
  cl <- makeCluster(min(workers, length(months)))
  res <- parLapplyLB(cl, months, ingest_month, cfg = cfg)
  stopCluster(cl)
  res <- rbindlist(res)
  print(res, digits = 3)
  fwrite(res, file.path(cfg$tables, "ingest_summary_2025.csv"))
  cat("All months done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 40 · Trip-level cleaning, stop-to-stop flows and interchange visits ────────────────────────────────────
section("40", "Trip-level cleaning, stop-to-stop flows and interchange visits", {
  ## Stage 2, step 1: one pass over the visit records of each month, January-
  ## October 2025 (manuscript §2.6). Replaces the former scripts 40 (trip
  ## cleaning), 44 (stop-to-stop flows) and 47 (15-minute timing), which repeated
  ## the same cleaning.
  ##
  ## Input: PEKA/output/tables/peka_2025_MM.rds (PEKA scripts 01 and 09).
  ## Cleaning (Furth et al. 2006): one record per trip and stop index (trip =
  ## duty date x course id); direction of a trip = its terminal stop cluster;
  ## vehicle-days whose counted boardings and alightings differ by more than 10%
  ## are screened out; within each counted trip, boardings and alightings are
  ## rescaled to their mean total; trips with partial counter data, a residual
  ## imbalance above max(5, 30%) before rescaling, or an implied on-board load
  ## below -3 after rescaling are screened out. Lines missing a mode get the mode
  ## of the vehicles that run them.
  ## Outputs per month (output/stage2/...):
  ##   base/peka_base_2025_MM.rds  date x stop cluster x pole x line x direction x
  ##                               scheduled hour: operated visits, valid counted
  ##                               visits, balanced boardings and alightings (all lines)
  ##   od/od_2025_MM.rds           stop-to-stop flows by line and direction on each valid counted trip
  ##                               (complex pairs, "outside" beyond the city), by month, day type and
  ##                               period of the boarding. flow/fkm: Tsygalnitsky recursion = null-seed
  ##                               biproportional fit, which reproduces boardings, alightings and loads
  ##                               (Furth & Navick 1992); flow_dist/fkm_dist: biproportional fit with a
  ##                               seed proportional to ride length (Navick & Furth 1994);
  ##                               flow_prop/fkm_prop: boardings assigned to the stops further along in
  ##                               proportion to their alightings (the earlier rule, which does not
  ##                               reproduce alightings or loads). fkm = flow x distance along the route
  ##                               (cumulative steps between stop-cluster centroids); hour_d =
  ##                               scheduled hour at the alighting stop
  ##   hubv/hubvisits_2025_MM.rds  every visit at the 64 interchanges (all lines, counted or not):
  ##                               time to the second, line, direction, trip first/last stop,
  ##                               balanced counts on valid counted trips (for transfer estimation)
  ##   fine/fine15_2025_MM.rds     day lines, regular working days: counted boardings
  ##                               and alightings in 15-minute bins of actual time by
  ##                               line x direction x complex, with a term-time
  ##                               Monday-Thursday flag
  ## Summaries: base/cleaning_summary.csv, od/od_summary.csv
  ## Usage: Rscript 40_trip_level_outputs.R [months] [workers]

  suppressPackageStartupMessages({library(data.table); library(parallel)})
  ## months and parallel workers come from the configuration (R/setup.R)
  months  <- PEKA_MONTHS
  workers <- WORKERS
  cfg <- list(peka = PEKA_TABLES,
              gtfs = GTFS_DIR,
              s2 = file.path(ROOT, "output/stage2"))
  for (k in c("base", "od", "fine", "hubv")) dir.create(file.path(cfg$s2, k), showWarnings = FALSE, recursive = TRUE)

  month_pass <- function(m, cfg) {
    suppressPackageStartupMessages({library(data.table); library(sf)}); setDTthreads(4)
    t0 <- Sys.time(); mm <- sprintf("%02d", m)
    C <- st_read(file.path(cfg$s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
    gt <- unique(fread(file.path(cfg$gtfs, "stops.txt"), encoding = "UTF-8")[, .(code = sub("[0-9]+$", "", stop_code), stop_name)])
    gs <- fread(file.path(cfg$gtfs, "stops.txt"), encoding = "UTF-8")[, code := sub("[0-9]+$", "", stop_code)]
    gt <- gt[stop_name %in% C$stop_name][!duplicated(code)]
    hubs <- fread(file.path(cfg$s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[interchange %in% TRUE, stop_name]
    ## stop-cluster centroids (all clusters, in and beyond the city), metres in the complexes' CRS
    cc <- gs[, .(lon = mean(stop_lon), lat = mean(stop_lat)), by = code]
    xy <- st_coordinates(st_transform(st_as_sf(cc, coords = c("lon", "lat"), crs = 4326), st_crs(C)))
    cc[, `:=`(X = xy[, 1], Y = xy[, 2])]
    x <- readRDS(file.path(cfg$peka, sprintf("peka_2025_%s.rds", mm))); setDT(x)
    x <- x[!is.na(DATA) & !is.na(SYM_PRZYST) & !is.na(NR_LINI) & NR_LINI != "",
           .(DT_KARTY, ID_KURSU, LP_PRZYST, NB_WOZU, DATA, KIEDY, NR_LINI, SYM_PRZYST, SYM_SLUPKA, sched_hour, has_apc, boardings, alightings, mode)]
    n_raw <- nrow(x)
    x <- unique(x, by = c("DT_KARTY", "ID_KURSU", "LP_PRZYST"))

    ## ── cleaning ───────────────────────────────────────────────────────────
    x[mode == "unknown", mode := NA_character_]
    vm <- x[mode %in% c("tram", "bus"), .N, by = .(NB_WOZU, mode)][order(-N)][!duplicated(NB_WOZU), .(NB_WOZU, vmode = mode)]
    lm_ <- x[is.na(mode)][vm, on = "NB_WOZU", nomatch = 0][, .N, by = .(NR_LINI, vmode)][order(-N)][!duplicated(NR_LINI)]
    x[lm_, mode := fifelse(is.na(mode), i.vmode, mode), on = "NR_LINI"][is.na(mode), mode := "unknown"]
    x[, tid := .GRP, by = .(DT_KARTY, ID_KURSU)]
    setorder(x, tid, LP_PRZYST)
    x[, dir := SYM_PRZYST[.N], by = tid]
    x[, `:=`(b = fcoalesce(as.numeric(boardings), 0), a = fcoalesce(as.numeric(alightings), 0))]
    vd <- x[has_apc == TRUE, .(B = sum(b), A = sum(a)), by = .(NB_WOZU, DT_KARTY)]
    vd[, vd_ok := pmax(B, A) == 0 | abs(B - A) <= 0.10 * pmax(B, A)]
    x[vd, vd_ok := i.vd_ok, on = .(NB_WOZU, DT_KARTY)]
    tr <- x[, .(all_apc = all(has_apc), B = sum(b), A = sum(a)), by = tid]
    tr[, imb_ok := abs(B - A) <= pmax(5, 0.30 * pmax(B, A))]
    tr[, `:=`(fb = fifelse(B > 0, (B + A) / 2 / B, 1), fa = fifelse(A > 0, (B + A) / 2 / A, 1))]
    x[tr, `:=`(all_apc = i.all_apc, imb_ok = i.imb_ok, fb = i.fb, fa = i.fa), on = "tid"]
    x[, `:=`(bb = b * fb, ab = a * fa)]
    x[, load := cumsum(bb - ab), by = tid]
    x[x[, .(load_ok = min(load) >= -3), by = tid], load_ok := i.load_ok, on = "tid"]
    x[, counted := has_apc & all_apc & (vd_ok %in% TRUE) & imb_ok & load_ok]
    ap <- x[has_apc == TRUE]
    diag <- data.table(month = mm, visits = nrow(x), duplicates_removed = n_raw - nrow(x), apc_visits = nrow(ap),
                       apc_share = nrow(ap) / nrow(x), screened_vehicle_day = mean(!(ap$vd_ok %in% TRUE)),
                       screened_partial_trip = mean(!ap$all_apc), screened_trip_imbalance = mean(!ap$imb_ok), screened_negative_load = mean(!ap$load_ok),
                       valid_counted_share_of_apc = mean(ap$counted), valid_counted_share_of_all = mean(x$counted),
                       trips = uniqueN(x$tid), directions_per_line = x[, uniqueN(dir), by = NR_LINI][, median(V1)])

    ## ── (1) base table ─────────────────────────────────────────────────────
    base <- x[, .(visits = .N, counted = sum(counted), B = sum(bb[counted]), A = sum(ab[counted])),
              by = .(DATA, SYM_PRZYST, SYM_SLUPKA, NR_LINI, dir, sched_hour, mode)]
    saveRDS(base, file.path(cfg$s2, "base", sprintf("peka_base_2025_%s.rds", mm)))
    diag[, base_rows := nrow(base)]; rm(base)

    ## ── (4) every visit at an interchange (for transfer estimation) ───────
    hub_codes <- gt[stop_name %in% hubs, code]
    ends <- x[, .(lp_min = min(LP_PRZYST), lp_max = max(LP_PRZYST)), by = tid]
    hv <- x[SYM_PRZYST %in% hub_codes, .(DATA, code = SYM_PRZYST, NR_LINI, dir, mode, tid, LP_PRZYST, sched_hour, KIEDY, has_apc, counted, bb, ab)]
    hv[ends, `:=`(first = LP_PRZYST == i.lp_min, last = LP_PRZYST == i.lp_max), on = "tid"]
    hv[gt, cx := i.stop_name, on = "code"]
    hv[, sec := as.integer(substr(KIEDY, 1, 2)) * 3600L + as.integer(substr(KIEDY, 4, 5)) * 60L + as.integer(substr(KIEDY, 7, 8))]
    hv[counted == FALSE, `:=`(bb = NA_real_, ab = NA_real_)]
    saveRDS(hv[, .(DATA, cx, NR_LINI, dir, mode, tid, sched_hour, sec, first, last, has_apc, counted, B = bb, A = ab)],
            file.path(cfg$s2, "hubv", sprintf("hubvisits_2025_%s.rds", mm)))
    diag[, hub_visits := nrow(hv)]; rm(hv, ends)

    ## valid counted trips only for flows and fine timing
    ok <- x[, .(ok = all(counted)), by = tid][ok == TRUE, tid]
    v <- x[tid %in% ok]; rm(x, ap); gc()
    v[gt, cx := i.stop_name, on = c(SYM_PRZYST = "code")]
    hol <- as.IDate(c("2025-01-01", "2025-01-06", "2025-04-20", "2025-04-21", "2025-05-01", "2025-05-03", "2025-06-08", "2025-06-19", "2025-08-15"))
    bridge <- as.IDate(c("2025-01-02", "2025-01-03", "2025-05-02", "2025-06-20")); dst <- as.IDate(c("2025-03-30", "2025-10-26"))
    dow <- as.integer(format(as.Date(v$DATA), "%u"))
    v[, daytype := fifelse(DATA %in% hol | dow == 7, "sunhol", fifelse(DATA %in% bridge, "bridge", fifelse(dow == 6, "sat", fifelse(dow == 5, "fri", "monthu"))))]

    ## ── (2) stop-to-stop flows ─────────────────────────────────────────────
    w <- copy(v)[is.na(cx), cx := "outside"]
    w[, period := as.character(cut(sched_hour, c(-1, 3, 5, 8, 14, 18, 99), labels = c("night", "early", "am_peak", "midday", "pm_peak", "evening")))]
    w[, S := rev(cumsum(rev(ab))) - ab, by = tid]
    ## distance along the route: cumulative straight-line steps between stop-cluster centroids
    w[cc, `:=`(X = i.X, Y = i.Y), on = c(SYM_PRZYST = "code")]
    w[, `:=`(X = nafill(nafill(X, "locf"), "nocb"), Y = nafill(nafill(Y, "locf"), "nocb")), by = tid]
    w[, km := cumsum(fcoalesce(sqrt((X - shift(X))^2 + (Y - shift(Y))^2), 0)) / 1000, by = tid]
    ## Tsygalnitsky recursion, identical to the null-seed biproportional fit (Furth & Navick 1992): at each
    ## stop the alightings are drawn from everyone on board in proportion to their numbers, so the flows
    ## reproduce boardings, alightings and on-board loads. p = share of the arriving load that alights;
    ## where everyone alights (z), riders from earlier stops cannot travel further.
    w[, L := shift(cumsum(bb - ab), fill = 0), by = tid]
    w[, p := fifelse(ab <= 0, 0, fifelse(L > ab, ab / L, 1))]
    w[, z := p >= 1 - 1e-9]
    w[, `:=`(C = cumsum(fifelse(z, 0, log1p(-p))), g = cumsum(z)), by = tid]
    w[, `:=`(Cp = shift(C, fill = 0), gp = shift(g, fill = 0L)), by = tid]
    O <- w[bb > 0 & S > 0, .(tid, lp_o = LP_PRZYST, o = cx, NR_LINI, dir, bb, S, km_o = km, C_o = C, g_o = g, daytype, period, DATA)]
    Dd <- w[ab > 0, .(tid, lp_d = LP_PRZYST, d = cx, ab, km_d = km, hd = sched_hour, p, Cp, gp)]
    P <- Dd[O, on = .(tid, lp_d > lp_o), .(o = i.o, d = x.d, NR_LINI = i.NR_LINI, dir = i.dir, daytype = i.daytype, period = i.period, DATA = i.DATA,
                                           f_prop = i.bb * x.ab / i.S, f_tsy = i.bb * x.p * exp(x.Cp - i.C_o) * (x.gp == i.g_o),
                                           km = x.km_d - i.km_o, hour_d = x.hd, tid = i.tid, lp_o = i.lp_o, lp_d = x.lp_d, bo = i.bb, ad = x.ab),
            allow.cartesian = TRUE, nomatch = 0]
    ## distance-seeded biproportional fit (seed proportional to ride length; Navick & Furth 1994), as a sensitivity variant
    P[, go := frank(.SD, ties.method = "dense"), .SDcols = c("tid", "lp_o")]
    P[, gd := frank(.SD, ties.method = "dense"), .SDcols = c("tid", "lp_d")]
    f <- pmax(P$km, 0.2); bo <- P$bo; ad <- P$ad; go <- P$go; gd <- P$gd
    for (it in 1:40) {
      f <- f * bo / as.vector(rowsum(f, go))[go]
      f <- f * ad / as.vector(rowsum(f, gd))[gd]
    }
    P[, f_dist := f]; rm(f, bo, ad, go, gd)
    ## reproduction of the counts by each estimator (share of alightings misallocated; row error of the distance fit)
    ae <- function(v) { s <- P[, .(s = sum(get(v)), ad = ad[1]), by = gd]; sum(abs(s$s - s$ad)) / sum(s$ad) }
    re <- P[, .(s = sum(f_dist), bo = bo[1]), by = go][, sum(abs(s - bo)) / sum(bo)]
    diag[, `:=`(alight_err_prop = ae("f_prop"), alight_err_tsy = ae("f_tsy"), alight_err_dist = ae("f_dist"), board_err_dist = re)]
    P <- P[DATA >= as.IDate("2025-01-01") & DATA <= as.IDate("2025-10-31")][, month := as.integer(format(as.Date(DATA), "%m"))]
    od <- P[, .(flow = sum(f_tsy), fkm = sum(f_tsy * km), flow_prop = sum(f_prop), fkm_prop = sum(f_prop * km),
                flow_dist = sum(f_dist), fkm_dist = sum(f_dist * km), trips = uniqueN(tid)), by = .(o, d, NR_LINI, dir, month, daytype, period, hour_d)]
    saveRDS(od, file.path(cfg$s2, "od", sprintf("od_2025_%s.rds", mm)))
    diag[, `:=`(valid_trips = length(ok), boardings_assigned = round(sum(P$f_tsy)), od_rows = nrow(od),
                coords_share = w[, mean(SYM_PRZYST %in% cc$code)], mean_ride_km = sum(P$f_tsy * P$km) / sum(P$f_tsy),
                mean_ride_km_prop = sum(P$f_prop * P$km) / sum(P$f_prop), mean_ride_km_dist = sum(P$f_dist * P$km) / sum(P$f_dist))]; rm(w, O, Dd, P, od); gc()

    ## ── (3) 15-minute timing on regular working days (day lines) ──────────
    f <- v[!is.na(cx) & !is.na(KIEDY) & !grepl("^2[0-9][0-9]$", NR_LINI) & DATA >= as.IDate("2025-01-01") & DATA <= as.IDate("2025-10-31") &
           !DATA %in% c(hol, bridge, dst) & daytype %in% c("monthu", "fri")]
    f[, month := as.integer(format(as.Date(DATA), "%m"))]
    f[, term_monthu := daytype == "monthu" & ((month %in% 3:5) | (month == 6 & DATA <= as.IDate("2025-06-27")) | month == 10) &
                       !(DATA >= as.IDate("2025-04-17") & DATA <= as.IDate("2025-04-22"))]
    f[, bin := (as.integer(substr(KIEDY, 1, 2)) * 60L + as.integer(substr(KIEDY, 4, 5))) %/% 15L]
    saveRDS(f[, .(B = sum(bb), A = sum(ab), visits = .N), by = .(NR_LINI, dir, cx, month, term_monthu, bin)],
            file.path(cfg$s2, "fine", sprintf("fine15_2025_%s.rds", mm)))

    diag[, minutes := round(as.numeric(Sys.time() - t0, units = "mins"), 1)]
    diag
  }

  t0 <- Sys.time()
  cl <- makeCluster(min(workers, length(months)))
  res <- rbindlist(parLapplyLB(cl, months, month_pass, cfg = cfg))
  stopCluster(cl)
  setorder(res, month)
  ODC <- c("valid_trips", "boardings_assigned", "od_rows", "coords_share", "mean_ride_km", "mean_ride_km_prop", "mean_ride_km_dist",
           "alight_err_prop", "alight_err_tsy", "alight_err_dist", "board_err_dist")
  fwrite(res[, !ODC, with = FALSE], file.path(cfg$s2, "base", "cleaning_summary.csv"))
  fwrite(res[, c("month", ODC), with = FALSE], file.path(cfg$s2, "od", "od_summary.csv"))
  print(res, digits = 3)
  cat("All months done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 41 · Expansion and stop-level measures ─────────────────────────────────────────────────────────────────
section("41", "Expansion and stop-level measures", {
  ## Stage 2, step 2: expansion of the cleaned counts to all operated visits,
  ## line x stop cells for the models, and stop-level usage measures
  ## (manuscript §2.6), January-October 2025.
  ##
  ## Input: output/stage2/base/peka_base_2025_MM.rds (script 40), stop clusters
  ## mapped to the 607 complexes through the October 2025 GTFS.
  ## Calendar (as script 33): public holidays count as Sundays; bridge days kept
  ## apart; daylight-saving changeover days excluded; term time = March to
  ## 27 June and October without the Easter break (17-22 April); summer =
  ## July-August; Wielkopolska winter school break 20 January-2 February;
  ## September = schools without universities.
  ## Diversion days: stop-days with fewer than half of the complex's median
  ## operated visits for their day type are excluded throughout.
  ## Expansion: per-visit balanced boardings and alightings from valid counted
  ## visits, taken from the first stratum with at least five counted visits:
  ##   L1 pole x line x direction x hour x day type x month
  ##   L2 complex x line x hour x day type x month (poles and directions merged)
  ##   L3 complex x line x period of the day x day type x month
  ##   L4 complex x line x period x day type, across months
  ## visits still uncovered are left out; per-visit means are never pooled across
  ## lines. Complexes with fewer than 10% of their visits counted are flagged.
  ## Outputs (output/stage2/):
  ##   stage2_line_stop_cells.rds  complex x line x direction x month x day type x
  ##                               hour: departures, expanded boardings/alightings,
  ##                               counted visits, expansion level shares
  ##   stage2_complex_day.rds      complex x date totals (all lines; day lines apart)
  ##   stage2_complex_measures.csv primary and secondary measures per complex
  ##   stage2_expansion_summary.csv, stage2_coverage_check.csv

  suppressPackageStartupMessages({library(data.table); library(sf)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  gtfs <- GTFS_DIR
  t0 <- Sys.time()

  ## ── base table, complexes, calendar ──────────────────────────────────────
  x <- rbindlist(lapply(sprintf("%02d", 1:10), function(m) readRDS(file.path(s2, "base", sprintf("peka_base_2025_%s.rds", m)))))
  x <- x[DATA >= as.IDate("2025-01-01") & DATA <= as.IDate("2025-10-31")]
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  gt <- unique(fread(file.path(gtfs, "stops.txt"), encoding = "UTF-8")[, .(code = sub("[0-9]+$", "", stop_code), stop_name)])
  gt <- gt[stop_name %in% C$stop_name][!duplicated(code)]
  x[gt, stop_name := i.stop_name, on = c(SYM_PRZYST = "code")]
  cat(sprintf("Operated visits: %.1f M; at city complexes %.1f%%\n", sum(x$visits) / 1e6, 100 * x[!is.na(stop_name), sum(visits)] / sum(x$visits)))
  x <- x[!is.na(stop_name)]
  D <- data.table(DATA = seq(as.IDate("2025-01-01"), as.IDate("2025-10-31"), by = 1))
  hol    <- as.IDate(c("2025-01-01", "2025-01-06", "2025-04-20", "2025-04-21", "2025-05-01", "2025-05-03", "2025-06-08", "2025-06-19", "2025-08-15"))
  bridge <- as.IDate(c("2025-01-02", "2025-01-03", "2025-05-02", "2025-06-20"))
  dst    <- as.IDate(c("2025-03-30", "2025-10-26"))
  D[, `:=`(dow = as.integer(format(as.Date(DATA), "%u")), month = as.integer(format(as.Date(DATA), "%m")))]
  D[, daytype := fifelse(DATA %in% hol | dow == 7, "sunhol", fifelse(DATA %in% bridge, "bridge", fifelse(dow == 6, "sat", fifelse(dow == 5, "fri", "monthu"))))]
  D[, season := fcase(DATA >= as.IDate("2025-01-20") & DATA <= as.IDate("2025-02-02"), "ferie",
                      month %in% 1:2, "winter_school",
                      DATA >= as.IDate("2025-04-17") & DATA <= as.IDate("2025-04-22"), "easter_break",
                      (month %in% 3:5) | (month == 6 & DATA <= as.IDate("2025-06-27")) | month == 10, "term",
                      month %in% 7:8, "summer", month == 6, "summer_break_june", month == 9, "september")]
  D <- D[!DATA %in% dst]
  x <- merge(x, D[, .(DATA, month, daytype, season)], by = "DATA")
  x[, night := grepl("^2[0-9][0-9]$", NR_LINI)]
  x[, period := cut(sched_hour, c(-1, 3, 5, 8, 14, 18, 99), labels = c("night", "early", "am_peak", "midday", "pm_peak", "evening"))]

  ## ── diversion days ───────────────────────────────────────────────────────
  sd <- x[, .(visits = sum(visits)), by = .(stop_name, DATA, daytype)]
  sd[, diversion := visits < 0.5 * median(visits), by = .(stop_name, daytype)]
  x[sd, diversion := i.diversion, on = .(stop_name, DATA)]
  cat(sprintf("Stop-days flagged as diversions: %.1f%% (%.1f%% of visits)\n", 100 * mean(sd$diversion), 100 * x[diversion == TRUE, sum(visits)] / sum(x$visits)))
  x <- x[diversion == FALSE]

  ## ── expansion ────────────────────────────────────────────────────────────
  lev <- list(L1 = c("SYM_SLUPKA", "NR_LINI", "dir", "sched_hour", "daytype", "month"),
              L2 = c("stop_name", "NR_LINI", "sched_hour", "daytype", "month"),
              L3 = c("stop_name", "NR_LINI", "period", "daytype", "month"),
              L4 = c("stop_name", "NR_LINI", "period", "daytype"))
  x[, `:=`(level = NA_integer_, mb = NA_real_, ma = NA_real_)]
  for (k in seq_along(lev)) {
    s <- x[, .(c = sum(counted), b = sum(B), a = sum(A)), by = c(lev[[k]])][c >= 5]
    x[s, `:=`(level = fifelse(is.na(level), k, level), mb = fifelse(is.na(mb), i.b / i.c, mb), ma = fifelse(is.na(ma), i.a / i.c, ma)), on = lev[[k]]]
  }
  x[, `:=`(Bx = visits * mb, Ax = visits * ma)]
  ex <- x[, .(visits = sum(visits)), by = .(month, level)][, share := visits / sum(visits), by = month]
  exw <- dcast(ex, month ~ level, value.var = "share", fill = 0); setnames(exw, as.character(c(1:4, NA))[as.character(c(1:4, NA)) %in% names(exw)], paste0("L", c(1:4, "_uncovered"))[as.character(c(1:4, NA)) %in% names(exw)])
  print(exw, digits = 3)
  fwrite(exw, file.path(s2, "stage2_expansion_summary.csv"))

  ## ── line x stop cells (model data) ───────────────────────────────────────
  cells <- x[, .(departures = sum(visits), counted = sum(counted), B = sum(Bx, na.rm = TRUE), A = sum(Ax, na.rm = TRUE),
                 uncovered = sum(visits[is.na(level)]), l1 = sum(visits[level %in% 1L]), days = uniqueN(DATA)),
             by = .(stop_name, NR_LINI, dir, mode, night, month, daytype, season, sched_hour, period)]
  saveRDS(cells, file.path(s2, "stage2_line_stop_cells.rds"))
  ## expansion factors for the stop-to-stop flows of script 40 (counted trips only): operated over validly
  ## counted visits of the boarding stratum, with two fallbacks (used by od_expand.R)
  ec <- copy(cells)[, period := as.character(period)]
  odx <- list(e1 = ec[, .(dep = sum(departures), cnt = sum(counted)), by = .(o = stop_name, NR_LINI, dir, month, daytype, period)][cnt > 0, .(o, NR_LINI, dir, month, daytype, period, e1 = dep / cnt)],
              e2 = ec[, .(dep = sum(departures), cnt = sum(counted)), by = .(o = stop_name, NR_LINI, dir, daytype, period)][cnt > 0, .(o, NR_LINI, dir, daytype, period, e2 = dep / cnt)],
              e3 = ec[, .(dep = sum(departures), cnt = sum(counted)), by = .(NR_LINI, dir, month, daytype, period)][cnt > 0, .(NR_LINI, dir, month, daytype, period, e3 = dep / cnt)])
  saveRDS(odx, file.path(s2, "od_expansion_factors.rds")); rm(ec, odx)
  cat(sprintf("Line x stop cells: %d rows; lines %d; directions %d; complexes %d\n", nrow(cells), uniqueN(cells$NR_LINI), uniqueN(cells[, paste(NR_LINI, dir)]), uniqueN(cells$stop_name)))

  ## ── complex-day totals ───────────────────────────────────────────────────
  cd <- x[, .(visits = sum(visits), counted = sum(counted), V = sum(Bx + Ax, na.rm = TRUE), Vday = sum((Bx + Ax)[!night], na.rm = TRUE),
              uncovered = sum(visits[is.na(level)]), tram_visits = sum(visits[mode == "tram"])), by = .(stop_name, DATA, month, daytype, season)]
  saveRDS(cd, file.path(s2, "stage2_complex_day.rds"))

  ## ── stop-level measures ──────────────────────────────────────────────────
  W <- cd[daytype == "monthu"]
  mv <- function(s, dt = "monthu") cd[season %in% s & daytype == dt, .(v = mean(V), visits = mean(visits), vpv = sum(V) / sum(visits), days = .N), by = stop_name]
  tm <- mv("term"); sm <- mv("summer"); se <- mv("september"); fe <- mv("ferie"); ws <- mv("winter_school")
  fr <- mv("term", "fri"); sa <- mv("term", "sat"); su <- mv("term", "sunhol")
  cv_m <- W[season == "term", .(v = mean(V)), by = .(stop_name, month)][, .(term_cv = sd(v) / mean(v)), by = stop_name]
  hr <- x[night == FALSE & season == "term", .(B = sum(Bx, na.rm = TRUE), A = sum(Ax, na.rm = TRUE)), by = .(stop_name, daytype, sched_hour)]
  tw <- hr[daytype == "monthu"]
  tim <- tw[, .(tidal = sum(B[sched_hour %in% 6:8]) / (sum(B[sched_hour %in% 6:8]) + sum(A[sched_hour %in% 6:8])) -
                        sum(B[sched_hour %in% 15:17]) / (sum(B[sched_hour %in% 15:17]) + sum(A[sched_hour %in% 15:17])),
                onset_b = sum(B[sched_hour >= 3 & sched_hour < 6]) / sum(B), onset_a = sum(A[sched_hour >= 3 & sched_hour < 6]) / sum(A),
                conc = (sum((B + A)[sched_hour %in% 7:8]) / 2) / (sum((B + A)[sched_hour %in% 10:13]) / 4),
                late_b = sum(B[sched_hour >= 21 | sched_hour < 3]) / sum(B), Bd = sum(B), Ad = sum(A)), by = stop_name]
  dd <- cd[season == "term" & daytype %in% c("monthu", "sat"), .(days = .N), by = .(stop_name, daytype)]
  sat <- merge(hr[daytype == "sat", .(Bs = sum(B), As = sum(A)), by = stop_name], dd[daytype == "sat", .(stop_name, ds = days)], by = "stop_name")
  sat <- merge(sat, merge(tim[, .(stop_name, Bd, Ad)], dd[daytype == "monthu", .(stop_name, dw = days)], by = "stop_name"), by = "stop_name")
  sat[, `:=`(sat = ((Bs + As) / ds) / ((Bd + Ad) / dw), sat_b = (Bs / ds) / (Bd / dw), sat_a = (As / ds) / (Ad / dw))]
  cov <- cd[season == "term" & daytype == "monthu", .(counted_share = sum(counted) / sum(visits), uncovered_share = sum(uncovered) / sum(visits)), by = stop_name]
  M <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(
    tm[, .(stop_name, volume = v, visits_term = visits, vpv_term = vpv, days_term = days)],
    tim[, .(stop_name, tidal, onset_b, onset_a, conc, late_b)], sat[, .(stop_name, sat, sat_b, sat_a)],
    sm[, .(stop_name, vol_summer = v, visits_summer = visits, vpv_summer = vpv)], se[, .(stop_name, vol_sept = v)],
    fe[, .(stop_name, vol_ferie = v)], ws[, .(stop_name, vol_winter_school = v)],
    fr[, .(stop_name, vol_fri = v)], su[, .(stop_name, vol_sun = v)], cv_m, cov))
  M[, `:=`(summer_ratio = vol_summer / volume, summer_ratio_pv = vpv_summer / vpv_term, sept_ratio = vol_sept / volume,
           ferie_ratio = vol_ferie / vol_winter_school, fri_ratio = vol_fri / volume, sun_ratio = vol_sun / volume,
           svc_summer = visits_summer / visits_term)]
  M[, `:=`(service_changed = svc_summer < 0.75 | svc_summer > 1.33, low_coverage = counted_share < 0.10)]
  fwrite(M, file.path(s2, "stage2_complex_measures.csv"))
  cat(sprintf("\nComplexes with term-time measures: %d | low counter coverage (<10%% counted): %d | summer service changed: %d\n",
              M[!is.na(volume), .N], M[low_coverage == TRUE, .N], M[service_changed == TRUE, .N]))
  print(round(M[volume >= 100, sapply(.SD, quantile, c(.1, .5, .9), na.rm = TRUE),
                .SDcols = c("volume", "tidal", "onset_b", "onset_a", "conc", "sat", "summer_ratio", "sept_ratio", "ferie_ratio", "term_cv", "counted_share")], 3))

  ## ── does counter coverage depend on land use or deprivation? ─────────────
  cvn <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  mun <- fread(file.path(s2, "municipal_catchments_network.csv"), encoding = "UTF-8")[, .(stop_name, muni_share)]
  k <- Reduce(function(a, b) merge(a, b, by = "stop_name"), list(M[, .(stop_name, counted_share)], cvn, mun))[residents > 0]
  sp <- function(v) cor(k$counted_share, k[[v]], method = "spearman", use = "complete.obs")
  cc <- data.table(variable = c("dep", "muni_share", "ilr_res_vs_nonres", "landuse_mix", "wc_tree", "dist_km", "residents"),
                   rho_with_counted_share = sapply(c("dep", "muni_share", "ilr_res_vs_nonres", "landuse_mix", "wc_tree", "dist_km", "residents"), sp))
  cc[, rho_net_of_distance := sapply(variable, function(v) if (v == "dist_km") NA_real_ else {
    z <- k[is.finite(counted_share) & is.finite(get(v))]; cor(resid(lm(rank(counted_share) ~ splines::ns(dist_km, 3), z)), resid(lm(rank(get(v)) ~ splines::ns(dist_km, 3), z)))
  })]
  cat("\nCounter coverage (share of term-time Monday-Thursday visits counted) vs catchment characteristics, walking catchments:\n"); print(cc, digits = 2)
  fwrite(cc, file.path(s2, "stage2_coverage_check.csv"))
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 42 · Same-trip indices ─────────────────────────────────────────────────────────────────────────────────
section("42", "Same-trip indices", {
  ## Stage 2, step 3: same-trip stop indices (manuscript §2.7, Identification).
  ##
  ## Every comparison is between stops served by the same trips. For each
  ## measure, a stop's value is computed within each trip group (line x
  ## direction, and period where the measure is defined within a period) as a
  ## rate per operated departure or a share within the trip group, expressed as
  ## the stop's deviation from the weighted mean of the other stops in that group,
  ## and averaged over the stop's groups. This is the fixed-effects contrast of a
  ## model with group fixed effects, computed directly; the timetable of each
  ## group is held fixed by construction. Only groups with at least two stops
  ## enter. Rates use covered departures (operated visits minus the 2-3% left
  ## unexpanded); small counts get +0.5.
  ##
  ## Indices (term-time Monday-Thursday unless stated; day lines for timing):
  ##   volume    log((B + A + 0.5) / departures), groups line x direction x period (all lines)
  ##   tidal     logit boarding share 06-09 minus logit boarding share 15-18, groups line x direction
  ##   onset_b/a logit share of boardings (alightings) before 06:00, groups line x direction that run before 06:00
  ##   conc(_b/_a) log of the per-departure rate 07-09 over 10-14, groups line x direction
  ##   sat(_b/_a)  log of the per-departure rate on Saturdays over Monday-Thursday, groups line x direction
  ##   summer    log of the per-departure rate in July-August over term time (Monday-Thursday, all lines)
  ##   sept, ferie  the same for September and for the winter school break (vs the school weeks of Jan-Feb)
  ## Monthly versions (persistence): volume, tidal, onset_b, conc_b, sat_b for each month (Monday-Thursday
  ## and Saturdays of that month).
  ## Outputs: output/stage2/stage2_stop_indices.csv, stage2_stop_indices_monthly.csv

  suppressPackageStartupMessages(library(data.table))
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  t0 <- Sys.time()
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))
  cells[, dep := departures - uncovered]
  cells <- cells[dep > 0]

  ## deviation of each stop from its trip-group peers, averaged over the stop's groups
  index <- function(dt, gcols) {
    dt <- dt[is.finite(q) & w > 0]
    dt[, n_in := .N, by = gcols]
    dt <- dt[n_in >= 2]
    dt[, dev := q - weighted.mean(q, w), by = gcols]
    dt[, .(index = weighted.mean(dev, w), groups = .N, weight = sum(w)), by = stop_name]
  }
  lg <- function(k, n) qlogis((k + 0.5) / (n + 1))

  build <- function(C, label) {
    Tm <- C[daytype == "monthu"]; Dm <- Tm[night == FALSE]
    out <- list()
    ## volume
    v <- Tm[, .(y = sum(B + A), dep = sum(dep)), by = .(NR_LINI, dir, period, stop_name)][, `:=`(q = log((y + 0.5) / dep), w = dep)]
    out$volume <- index(v, c("NR_LINI", "dir", "period"))
    ## tidal
    am <- Dm[sched_hour %in% 6:8, .(B = sum(B), A = sum(A)), by = .(NR_LINI, dir, stop_name)][, `:=`(q = lg(B, B + A), w = B + A)]
    pm <- Dm[sched_hour %in% 15:17, .(B = sum(B), A = sum(A)), by = .(NR_LINI, dir, stop_name)][, `:=`(q = lg(B, B + A), w = B + A)]
    ti <- merge(index(am, c("NR_LINI", "dir")), index(pm, c("NR_LINI", "dir")), by = "stop_name", suffixes = c("_am", "_pm"))
    out$tidal <- ti[, .(stop_name, index = index_am - index_pm, groups = pmin(groups_am, groups_pm), weight = weight_am + weight_pm)]
    ## onset (lines that run before 06:00)
    early <- Dm[sched_hour < 6, .(dep_early = sum(dep)), by = .(NR_LINI, dir)]
    on <- merge(Dm[, .(Be = sum(B[sched_hour < 6]), Ball = sum(B), Ae = sum(A[sched_hour < 6]), Aall = sum(A)), by = .(NR_LINI, dir, stop_name)],
                early[dep_early > 0], by = c("NR_LINI", "dir"))
    out$onset_b <- index(copy(on)[, `:=`(q = lg(Be, Ball), w = Ball)], c("NR_LINI", "dir"))
    out$onset_a <- index(copy(on)[, `:=`(q = lg(Ae, Aall), w = Aall)], c("NR_LINI", "dir"))
    ## concentration
    cc <- Dm[, .(Bp = sum(B[sched_hour %in% 7:8]), Ap = sum(A[sched_hour %in% 7:8]), dp = sum(dep[sched_hour %in% 7:8]),
                 Bm = sum(B[sched_hour %in% 10:13]), Am = sum(A[sched_hour %in% 10:13]), dm = sum(dep[sched_hour %in% 10:13])), by = .(NR_LINI, dir, stop_name)][dp > 0 & dm > 0]
    rr <- function(yp, ym) log((yp + 0.5) / cc$dp) - log((ym + 0.5) / cc$dm)
    out$conc   <- index(copy(cc)[, `:=`(q = rr(Bp + Ap, Bm + Am), w = dp + dm)], c("NR_LINI", "dir"))
    out$conc_b <- index(copy(cc)[, `:=`(q = rr(Bp, Bm), w = dp + dm)], c("NR_LINI", "dir"))
    out$conc_a <- index(copy(cc)[, `:=`(q = rr(Ap, Am), w = dp + dm)], c("NR_LINI", "dir"))
    ## Saturday
    sa <- merge(C[daytype == "sat" & night == FALSE, .(Bs = sum(B), As = sum(A), ds = sum(dep)), by = .(NR_LINI, dir, stop_name)],
                Dm[, .(Bw = sum(B), Aw = sum(A), dw = sum(dep)), by = .(NR_LINI, dir, stop_name)], by = c("NR_LINI", "dir", "stop_name"))
    rs <- function(ys, yw) log((ys + 0.5) / sa$ds) - log((yw + 0.5) / sa$dw)
    out$sat   <- index(copy(sa)[, `:=`(q = rs(Bs + As, Bw + Aw), w = ds + dw)], c("NR_LINI", "dir"))
    out$sat_b <- index(copy(sa)[, `:=`(q = rs(Bs, Bw), w = ds + dw)], c("NR_LINI", "dir"))
    out$sat_a <- index(copy(sa)[, `:=`(q = rs(As, Aw), w = ds + dw)], c("NR_LINI", "dir"))
    out
  }
  season_ratio <- function(num, den) {
    a <- cells[season %in% num & daytype == "monthu", .(y1 = sum(B + A), d1 = sum(dep)), by = .(NR_LINI, dir, stop_name)]
    b <- cells[season %in% den & daytype == "monthu", .(y0 = sum(B + A), d0 = sum(dep)), by = .(NR_LINI, dir, stop_name)]
    m <- merge(a, b, by = c("NR_LINI", "dir", "stop_name"))[, `:=`(q = log((y1 + 0.5) / d1) - log((y0 + 0.5) / d0), w = d1 + d0)]
    index(m, c("NR_LINI", "dir"))
  }

  ## ── term-time indices ────────────────────────────────────────────────────
  res <- build(cells[season == "term"], "term")
  res$summer <- season_ratio("summer", "term")
  res$sept   <- season_ratio("september", "term")
  res$ferie  <- season_ratio("ferie", "winter_school")
  W <- Reduce(function(a, b) merge(a, b, by = "stop_name", all = TRUE),
              lapply(names(res), function(k) setnames(copy(res[[k]]), c("index", "groups", "weight"), c(k, paste0(k, "_groups"), paste0(k, "_weight")))))
  fwrite(W, file.path(s2, "stage2_stop_indices.csv"))
  cat("Stops with each index:\n"); print(sapply(names(res), function(k) sum(is.finite(W[[k]]))))
  cat("\nIndex distributions (deviation from same-trip peers):\n")
  print(round(sapply(W[, names(res), with = FALSE], quantile, c(.1, .5, .9), na.rm = TRUE), 3))

  ## ── monthly indices (persistence) ────────────────────────────────────────
  mon <- rbindlist(lapply(1:10, function(m) {
    r <- build(cells[month == m & season != "easter_break"], m)
    rbindlist(lapply(c("volume", "tidal", "onset_b", "conc_b", "sat_b"), function(k) r[[k]][, .(stop_name, measure = k, month = m, index, weight)]))
  }))
  fwrite(mon, file.path(s2, "stage2_stop_indices_monthly.csv"))
  cat("\nMonthly indices:", nrow(mon), "rows. Correlation of each month with term time (Spearman, median over months):\n")
  cm <- merge(mon, melt(W[, .(stop_name, volume, tidal, onset_b, conc_b, sat_b)], id.vars = "stop_name", variable.name = "measure", value.name = "term"), by = c("stop_name", "measure"))
  print(cm[, .(rho = cor(index, term, method = "spearman", use = "complete.obs")), by = .(measure, month)][, .(median_rho = round(median(rho), 2), min_rho = round(min(rho), 2)), by = measure])
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 22 · Counter coverage (October 2025) ───────────────────────────────────────────────────────────────────
section("22", "Counter coverage (October 2025)", {
  ## Door-counter (APC) coverage diagnostics for the usage measures (manuscript
  ## §2.6), on the processed October 2025 PEKA month (PEKA project, script 01).
  ##
  ## Counters are fitted to part of the fleet, so counted visits are expanded to
  ## all operated visits within strata. This script reports how much of the
  ## operated service falls in strata with no counted visit at the finest
  ## stratification and after collapsing to complex x line x period x day type,
  ## the coverage by line and by stop complex, and how much early-morning
  ## (before 05:00) boarding comes from the night network (lines 2xx).
  ##
  ## Inputs: PEKA/output/tables/peka_2025_10.rds, October 2025 GTFS stops,
  ## output/stage2/stop_complex_attributes.csv (21).
  ## Output: output/stage2/apc_coverage_summary.csv

  suppressPackageStartupMessages(library(data.table))
  peka <- PEKA_DIR
  out_dir <- ROOT
  x <- readRDS(file.path(PEKA_TABLES, "peka_2025_10.rds")); setDT(x)
  cx <- fread(file.path(out_dir, "output/stage2/stop_complex_attributes.csv"), encoding = "UTF-8")
  gt <- fread(file.path(GTFS_DIR, "stops.txt"), encoding = "UTF-8")[, .(code = sub("[0-9]+$", "", stop_code), stop_name)]
  city_codes <- unique(gt[stop_name %in% cx$stop_name, code])

  res <- list()
  res$overall_coverage <- mean(x$has_apc)
  ln <- x[, .(visits = .N, cov = mean(has_apc)), by = NR_LINI]
  res$lines <- nrow(ln); res$line_cov_median <- median(ln$cov); res$lines_below_5pct <- sum(ln$cov < 0.05)
  cat(sprintf("Network: %.1f%% of operated visits counted | %d lines, median line coverage %.0f%%, %d lines below 5%%\n",
              100 * res$overall_coverage, res$lines, 100 * res$line_cov_median, res$lines_below_5pct))

  xc <- x[SYM_PRZYST %in% city_codes]
  xc[, hr := floor(as.numeric(sched_hour))]
  xc[, period := cut(hr, c(-1, 5, 9, 14, 18, 21, 24), labels = c("early", "am", "midday", "pm", "evening", "late"))]
  zero <- function(keys) { st <- xc[, .(v = .N, a = sum(has_apc)), by = keys]; c(strata_zero = mean(st$a == 0), visits_zero = sum(st[a == 0]$v) / sum(st$v)) }
  fine <- zero(c("SYM_SLUPKA", "NR_LINI", "OPIS_TABL", "hr", "is_weekend"))
  coarse <- zero(c("SYM_PRZYST", "NR_LINI", "period", "is_weekend"))
  res$fine_strata_zero <- fine[1]; res$fine_visits_zero <- fine[2]; res$coarse_strata_zero <- coarse[1]; res$coarse_visits_zero <- coarse[2]
  cat(sprintf("City complexes, one month: pole x line x destination x hour x day type -> %.0f%% of strata (%.0f%% of operated visits) have no counted visit;\n  complex x line x period x day type -> %.0f%% of strata (%.0f%% of visits)\n",
              100 * fine[1], 100 * fine[2], 100 * coarse[1], 100 * coarse[2]))

  cc <- xc[, .(cov = mean(has_apc)), by = SYM_PRZYST]
  res$complex_cov_median <- median(cc$cov); res$complexes_below_10pct <- sum(cc$cov < 0.10)
  cat(sprintf("Complexes: median coverage %.0f%%, %d of %d below 10%%\n", 100 * res$complex_cov_median, res$complexes_below_10pct, nrow(cc)))

  xc[, night := grepl("^2[0-9][0-9]$", as.character(NR_LINI))]
  res$night_line_coverage <- xc[night == TRUE, mean(has_apc)]
  res$early_before5_night_visits <- xc[hr < 5, mean(night)]
  res$early_before5_night_boardings <- xc[hr < 5 & has_apc == TRUE, sum(boardings[night], na.rm = TRUE) / sum(boardings, na.rm = TRUE)]
  res$h5_6_night_visits <- xc[hr == 5, mean(night)]
  cat(sprintf("Night lines (2xx): counter coverage %.0f%% | operated visits before 05:00 on night lines %.0f%% | counted boardings before 05:00 on night lines %.0f%% | visits 05:00-06:00 on night lines %.0f%%\n",
              100 * res$night_line_coverage, 100 * res$early_before5_night_visits, 100 * res$early_before5_night_boardings, 100 * res$h5_6_night_visits))
  fwrite(as.data.table(lapply(res, unname)), file.path(out_dir, "output/stage2/apc_coverage_summary.csv"))
})

## ── 23 · Does usage mirror the timetable? (October 2025) ───────────────────────────────────────────────────
section("23", "Does usage mirror the timetable? (October 2025)", {
  ## How much of a stop's usage rhythm mirrors its timetable rhythm?
  ## (October 2025 PEKA month; city stop complexes; day lines only.)
  ##
  ## For each complex and each timing measure of §2.6, compares the usage
  ## version (expanded boardings) with the same measure computed on operated
  ## departures (the timetable actually run): onset share (before 06:00),
  ## concentration ratio (07-09 vs 10-14, per hour) and Saturday ratio
  ## (Saturday vs Monday-Thursday, per day). Also checks within-line variation:
  ## the share of variation in per-visit boardings that lies between stops of
  ## the same line, direction and hour, where the timetable is identical.
  ##
  ## Expansion here is a simple within complex x line x hour x day-type
  ## scaling of counted visits (the full §2.6 cleaning is not applied).
  ## Output: output/stage2/supply_vs_usage_rhythm.csv

  suppressPackageStartupMessages(library(data.table))
  peka <- PEKA_DIR
  out_dir <- ROOT
  x <- readRDS(file.path(PEKA_TABLES, "peka_2025_10.rds")); setDT(x)
  cx <- fread(file.path(out_dir, "output/stage2/stop_complex_attributes.csv"), encoding = "UTF-8")
  gt <- fread(file.path(GTFS_DIR, "stops.txt"), encoding = "UTF-8")[, .(code = sub("[0-9]+$", "", stop_code), stop_name)]
  x <- x[SYM_PRZYST %in% unique(gt[stop_name %in% cx$stop_name, code]) & !grepl("^2[0-9][0-9]$", as.character(NR_LINI))]
  x[, `:=`(hr = floor(as.numeric(sched_hour)), dow = as.integer(format(as.Date(DATA), "%u")))]
  x[, daytype := fifelse(dow <= 4, "monthu", fifelse(dow == 6, "sat", "other"))]
  x <- x[daytype != "other"]
  ndays <- x[, .(nd = uniqueN(as.Date(DATA))), by = daytype]

  ## crude expansion: mean boardings per counted visit x operated visits, within complex x line x hour x day type
  cell <- x[, .(visits = .N, counted = sum(has_apc), b_counted = sum(boardings[has_apc], na.rm = TRUE)), by = .(SYM_PRZYST, NR_LINI, hr, daytype)]
  cell[, b_exp := fifelse(counted > 0, b_counted / counted * visits, NA_real_)]
  ch <- cell[, .(visits = sum(visits), b = sum(b_exp, na.rm = TRUE), cov = sum(counted) / sum(visits)), by = .(SYM_PRZYST, hr, daytype)]
  ch <- merge(ch, ndays, by = "daytype")
  m <- function(dt, d, h) dt[daytype == d & hr %in% h]
  per <- ch[, .(
    dep_onset = sum(visits[daytype == "monthu" & hr < 6]) / sum(visits[daytype == "monthu"]),
    use_onset = sum(b[daytype == "monthu" & hr < 6]) / sum(b[daytype == "monthu"]),
    dep_conc  = (sum(visits[daytype == "monthu" & hr %in% 7:8]) / 2) / (sum(visits[daytype == "monthu" & hr %in% 10:13]) / 4),
    use_conc  = (sum(b[daytype == "monthu" & hr %in% 7:8]) / 2) / (sum(b[daytype == "monthu" & hr %in% 10:13]) / 4),
    dep_sat   = (sum(visits[daytype == "sat"]) / nd[daytype == "sat"][1]) / (sum(visits[daytype == "monthu"]) / nd[daytype == "monthu"][1]),
    use_sat   = (sum(b[daytype == "sat"]) / nd[daytype == "sat"][1]) / (sum(b[daytype == "monthu"]) / nd[daytype == "monthu"][1]),
    cov = sum(cov * visits) / sum(visits), b_day = sum(b[daytype == "monthu"]) / nd[daytype == "monthu"][1]), by = SYM_PRZYST]
  per <- per[cov >= 0.10 & b_day >= 50 & is.finite(dep_conc) & is.finite(use_conc)]
  sp <- function(a, b) round(cor(a, b, method = "spearman", use = "complete.obs"), 2)
  res <- data.table(measure = c("onset share", "concentration ratio", "Saturday ratio"),
                    rho_usage_vs_timetable = c(sp(per$use_onset, per$dep_onset), sp(per$use_conc, per$dep_conc), sp(per$use_sat, per$dep_sat)),
                    timetable_range_p10_p90 = c(paste(round(quantile(per$dep_onset, c(.1, .9), na.rm = TRUE), 3), collapse = "-"),
                                                paste(round(quantile(per$dep_conc, c(.1, .9), na.rm = TRUE), 2), collapse = "-"),
                                                paste(round(quantile(per$dep_sat, c(.1, .9), na.rm = TRUE), 2), collapse = "-")))
  cat("Complexes analysed (>=10% counted, >=50 boardings per weekday):", nrow(per), "\n"); print(res)

  ## within-line comparison: per-visit boardings, share of variance between stops within line x direction x hour x day type
  v <- x[has_apc == TRUE, .(bpv = mean(boardings, na.rm = TRUE), n = .N), by = .(SYM_PRZYST, NR_LINI, OPIS_TABL, hr, daytype)][n >= 5]
  v[, lb := log1p(bpv)]
  tot <- var(v$lb); within <- v[, .(w = sum((lb - mean(lb))^2)), by = .(NR_LINI, OPIS_TABL, hr, daytype)][, sum(w)] / (nrow(v) - 1)
  cat(sprintf("Per-visit boardings (log): %.0f%% of variance lies between stops served by the same trips (line x direction x hour x day type)\n", 100 * within / tot))
  fwrite(cbind(res, n_complexes = nrow(per), within_trip_share = within / tot), file.path(out_dir, "output/stage2/supply_vs_usage_rhythm.csv"))
})

## ── 58 · Descriptive statistics of the counts ──────────────────────────────────────────────────────────────
section("58", "Descriptive statistics of the counts", {
  ## Descriptive statistics of the passenger-count (AVL/APC) data, January-October
  ## 2025, for supplementary Table S3: (a) by month: operated stop visits, share
  ## counted, share of counted visits passing the trip-level checks, valid counted
  ## trips, counted boardings on them, expansion levels and expanded boardings;
  ## (b) by stop group (term-time Monday-Thursday): departures, counter coverage,
  ## expanded boardings per stop and per departure, for all complexes and for the
  ## thirds of deprivation among stops with residents that are not flagged.
  ## Inputs: output/stage2/base/cleaning_summary.csv, od/od_summary.csv (script 40);
  ##         stage2_line_stop_cells.rds, stage2_complex_day.rds,
  ##         stage2_expansion_summary.csv (script 41); catchment_covariates_network.csv.
  ## Outputs: output/stage2/apc_descriptives_monthly.csv, apc_descriptives_groups.csv,
  ##          manuscript/tables/tab_apc_monthly.tex, tab_apc_groups.tex (table rows)

  suppressPackageStartupMessages(library(data.table))
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  TERM <- c(3:6, 10)

  ## ── (a) by month ─────────────────────────────────────────────────────────
  cl <- fread(file.path(s2, "base/cleaning_summary.csv"))[, .(month = as.integer(month), visits, apc_share, valid_counted_share_of_apc)]
  od <- fread(file.path(s2, "od/od_summary.csv"))[, .(month = as.integer(month), valid_trips, boardings_assigned)]
  ex <- fread(file.path(s2, "stage2_expansion_summary.csv"))
  setnames(ex, 2, "uncovered")   # second column: share of visits left uncovered (header "NA")
  ex <- ex[, .(month = as.integer(month), l1 = L1, uncovered)]
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))
  eb <- cells[, .(expanded_B = sum(B), departures = sum(departures), counted = sum(counted)), by = month]
  A <- Reduce(function(a, b) merge(a, b, by = "month"), list(cl, od, ex, eb))[order(month)]
  tot <- A[, .(month = NA_integer_, visits = sum(visits), apc_share = weighted.mean(apc_share, visits), valid_counted_share_of_apc = weighted.mean(valid_counted_share_of_apc, visits * apc_share),
               valid_trips = sum(valid_trips), boardings_assigned = sum(boardings_assigned), l1 = weighted.mean(l1, departures), uncovered = weighted.mean(uncovered, departures),
               expanded_B = sum(expanded_B), departures = sum(departures), counted = sum(counted))]
  A <- rbind(A, tot)
  fwrite(A, file.path(s2, "apc_descriptives_monthly.csv"))
  print(A[, .(month, visits_M = round(visits / 1e6, 2), counted = round(100 * apc_share), pass = round(100 * valid_counted_share_of_apc), trips_k = round(valid_trips / 1e3),
              counted_B_M = round(boardings_assigned / 1e6, 1), l1 = round(100 * l1), uncovered = round(100 * uncovered, 1), expanded_B_M = round(expanded_B / 1e6, 1))])

  ## ── (b) by stop group, term-time Monday-Thursday ─────────────────────────
  cd <- readRDS(file.path(s2, "stage2_complex_day.rds"))[daytype == "monthu" & month %in% TERM]
  ndays <- uniqueN(cd$DATA)
  cat("Term-time Monday-Thursday dates:", ndays, "\n")
  city_day <- cells[daytype == "monthu" & month %in% TERM, sum(B)] / ndays
  cat(sprintf("Expanded boardings on an average term-time Monday-Thursday: %.0f\n", city_day))
  mode_b <- cells[month %in% TERM & daytype == "monthu", .(B = sum(B), departures = sum(departures), counted = sum(counted)), by = mode][, `:=`(B_share = B / sum(B), cov = counted / departures)]
  print(mode_b)
  st <- merge(cd[, .(dep_day = sum(visits) / uniqueN(DATA), cov = sum(counted) / sum(visits)), by = stop_name],
              cells[daytype == "monthu" & month %in% TERM, .(B_day = sum(B) / ndays), by = stop_name], by = "stop_name")
  st[, bpd := B_day / dep_day]
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[, .(stop_name, dep, residents, flagged)]
  st <- merge(st, cv, by = "stop_name", all.x = TRUE)
  q <- st[residents > 0 & flagged == FALSE, quantile(dep, c(1/3, 2/3))]
  st[residents > 0 & flagged == FALSE, third := fcase(dep <= q[1], "Least deprived third", dep >= q[2], "Most deprived third", default = "Middle third")]
  summ <- function(D, g) D[, .(group = g, stops = .N, dep_day = median(dep_day), cov = median(cov), B_day = median(B_day), B_q1 = quantile(B_day, 0.25), B_q3 = quantile(B_day, 0.75), bpd = median(bpd))]
  G <- rbind(summ(st, "All stop complexes"), summ(st[!is.na(third)], "Stops with residents, not flagged"),
             summ(st[third == "Least deprived third"], "Least deprived third"), summ(st[third == "Middle third"], "Middle third"), summ(st[third == "Most deprived third"], "Most deprived third"))
  fwrite(G, file.path(s2, "apc_descriptives_groups.csv"))
  print(G, digits = 3)

  ## ── LaTeX rows ───────────────────────────────────────────────────────────
  mn <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October")
  f1 <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = "{,}")
  La <- A[, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", ifelse(is.na(month), "\\midrule January--October", mn[month]),
                    f1(visits / 1e6, 2), f1(100 * apc_share, 0), f1(100 * valid_counted_share_of_apc, 0), f1(valid_trips / 1e3, 0),
                    f1(boardings_assigned / 1e6, 1), f1(100 * l1, 0), f1(100 * uncovered, 1), f1(expanded_B / 1e6, 1))]
  Lb <- G[, sprintf("%s & %d & %s & %s & %s (%s--%s) & %s \\\\", group, stops, f1(dep_day, 0), f1(100 * cov, 0), f1(B_day, 0), f1(B_q1, 0), f1(B_q3, 0), f1(bpd, 1))]
  dir.create(file.path(out_dir, "manuscript/tables"), showWarnings = FALSE)
  writeLines(La, file.path(out_dir, "manuscript/tables/tab_apc_monthly.tex"), useBytes = TRUE)
  writeLines(Lb, file.path(out_dir, "manuscript/tables/tab_apc_groups.tex"), useBytes = TRUE)
  cat("Saved apc_descriptives_monthly.csv, apc_descriptives_groups.csv and the table rows\n")
})

