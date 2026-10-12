## DualPulse | stage 8: The 22 gaps: decomposition, inference and checks
##
## Slopes of every measure on deprivation raw, net of position and net of land use and population, their
## decomposition and small-cluster inference, and the checks built on them (main text 3.5, supplementary
## S7, S8, S10, S11).
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   56      Gap ladder, decomposition and inference
##   55      The extremes of deprivation
##   59      Why deprived areas sit where they do
##   63      Ride length: line mix and estimators
##   65      Concentration indices and access poverty
##   66      Spatial-noise placebo
##   67      Counterfactual geographies
##   68      Synthetic test of the destination estimators
##   73      Destinations offered and reached
##   74      Offer and use: what the network offers and how riders use it
##   75      Peripheral deprivation: what the network offers deprived areas far from the centre
##   76      Early boardings: the early offer, its use, and what sector stands for
##   77      What location stands for: jobs nearby and the structure of service
##
## Run from the repository root: Rscript run.R 8   (one section: Rscript run.R 8:73)

if (!exists("section")) source("R/setup.R")

## ── 56 · Gap ladder, decomposition and inference ───────────────────────────────────────────────────────────
section("56", "Gap ladder, decomposition and inference", {
  ## How public transit usage differs between more and less deprived areas, and
  ## how much of each difference position on the network, land use and
  ## deprivation itself account for: the central figure and table of the paper.
  ##
  ## For every stop-level measure of the study, the slope on catchment
  ## deprivation (both standardised over stops with residents, not flagged) in
  ## three steps:
  ##   raw                     no adjustment
  ##   net of network position location (spline of distance to the centre and
  ##                           sector, the first two harmonics of the bearing from
  ##                           the Old Market Square) and service at the stop (log
  ##                           lines, log departures, tram share, share of
  ##                           departures on lines that pass the centre); variants
  ##                           without service (network_noserv, landuse_noserv)
  ##   net of land use too     + residential, multi/single, services/production
  ##                           and retail/workplace balances, land-use mix,
  ##                           floor-area density, tree cover, age shares, log
  ##                           residents, log catchment area
  ## Primary specification ("position", "full"): the six same-trip usage indices
  ## already hold each trip's timetable fixed, and service partly responds to
  ## demand, so for them position = location only (the no-service variants); for
  ## the destination, route, access and response measures position includes
  ## service at the stop.
  ## A fourth estimate uses the finer indicator of deprivation, the share of
  ## catchment residents in municipal buildings (script 37), within osiedla:
  ## osiedle fixed effects, net of position on the network and land use.
  ## Inference (42 osiedla, unbalanced): CR2 standard errors clustered by osiedle
  ## with Satterthwaite (Bell-McCaffrey) degrees of freedom (Pustejovsky & Tipton
  ## 2018; Imbens & Kolesar 2016); checks: CV3 with t(G-1), restricted wild
  ## cluster bootstrap (WCR-C, Webb weights, 9,999 draws; MacKinnon, Nielsen &
  ## Webb 2023), partial leverage and the effective number of clusters G* (Carter,
  ## Schnepel & Steigerwald 2017), leave-one-osiedle-out. Benjamini-Hochberg q
  ## across the 22 measures in each step; Romano-Wolf step-down p-values for the
  ## remainders from the osiedle bootstrap.
  ## Measurement error: the score averages two rounds that correlate at r = 0.60
  ## across osiedla, reliability R = 2r/(1+r) = 0.75; with error-free controls a
  ## net slope is attenuated by lambda = (R - R2)/(1 - R2), R2 = share of the
  ## score's variance explained by the controls. Disattenuated slopes for R = 0.75
  ## and R = 0.60; re-estimates with the three-round mean, each single round and
  ## the catchments that draw >= 90% of residents from one osiedle (script 62).
  ## Equivalence: the full-step slope is "equivalent to zero" when its 90%
  ## interval lies within +-0.20 SD per SD (two one-sided tests).
  ## Decomposition (Gelbach 2016) of the raw slope into the parts accounted for by
  ## location, service, land use and population, and the remainder; 1,999
  ## bootstrap resamples of osiedla, shared across measures (joint draws give
  ## intervals for the headline ratios).
  ## Measures: usage (same-trip indices, term-time Monday-Thursday, script 42);
  ## destinations (46) and interchange use and connections (from flows and the
  ## transfer rates of 47); routes (51); access in the morning peak (54);
  ## responses to leave days, service cuts and weather (52).
  ## Outputs: output/stage2/gap_ladder.csv, gap_inference.csv, gap_loo.csv,
  ##          gap_measurement.csv, gap_decomposition.csv, gap_headline.csv;
  ##          output/figures/Fig_gap_ladder.png, Fig_gap_loo.png;
  ##          manuscript/tables/tab_usage_differences.tex (table rows)

  suppressPackageStartupMessages({library(data.table); library(splines); library(ggplot2)})
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R"))   # flows expanded from counted to all operated trips
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  MON <- c(3:6, 10)
  PER <- function(h) as.character(cut(h, c(-1, 3, 5, 8, 14, 18, 99), labels = c("night", "early", "am_peak", "midday", "pm_peak", "evening")))

  ## ── stops and covariates ────────────────────────────────────────────────
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  HUBS <- cv[interchange %in% TRUE, stop_name]
  td <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")[, .(stop_name, b1 = ilr_res_vs_nonres, b2 = ilr_multi_vs_single,
          b3 = ilr_serv_vs_prod, b4 = ilr_retail_vs_work, landuse_mix, fsi = log(fsi + 0.01), wc_tree)]
  d <- merge(cv[residents > 0 & flagged == FALSE, .(stop_name, dep, dist_km, osiedle, share_0_18, share_61plus, lres = log1p(as.numeric(residents)), larea = log(catchment_area_m2))],
             td, by = "stop_name")
  ## service at the stop, and the share of departures on lines that pass the centre
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))[daytype == "monthu" & month %in% MON & night == FALSE]
  sv <- cells[, .(l_lines = log(uniqueN(NR_LINI)), l_dep = log(sum(departures)), tram_share = sum(departures[mode == "tram"]) / sum(departures)), by = .(stop_name)]
  segdep <- cells[, .(dep = sum(departures)), by = .(o = stop_name, NR_LINI, dir)]; rm(cells); gc()
  od <- rbindlist(lapply(sprintf("%02d", MON), function(m) expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))))[daytype == "monthu" & o != "outside",
          .(flow = sum(flow), flow_c = sum(flow_c)), by = .(o, d, NR_LINI, dir, hour_d)]))
  od <- merge(od, cv[, .(d = stop_name, centre_km = dist_km)], by = "d", all.x = TRUE)
  reach <- merge(od[, .(reach = any(centre_km <= 1.5, na.rm = TRUE)), by = .(o, NR_LINI, dir)], segdep, by = c("o", "NR_LINI", "dir"))[, .(direct_centre = sum(dep * reach) / sum(dep)), by = .(stop_name = o)]
  d <- merge(merge(d, sv, by = "stop_name"), reach, by = "stop_name")
  d <- d[complete.cases(d)]
  ## sector: direction of the complex from the Old Market Square, as the first and second harmonics of its bearing
  xy <- fread(file.path(s2, "stop_complex_attributes.csv"), encoding = "UTF-8")[, .(stop_name, x, y)]
  ctr <- sf::st_coordinates(sf::st_transform(sf::st_sfc(sf::st_point(c(16.9335, 52.4083)), crs = 4326), 2180))
  xy[, th := atan2(y - ctr[2], x - ctr[1])][, `:=`(sec_s1 = sin(th), sec_c1 = cos(th), sec_s2 = sin(2 * th), sec_c2 = cos(2 * th), d_check = sqrt((x - ctr[1])^2 + (y - ctr[2])^2) / 1000)]
  d <- merge(d, xy[, .(stop_name, sec_s1, sec_c1, sec_s2, sec_c2, d_check)], by = "stop_name")
  cat(sprintf("Sector check: correlation of recomputed distance with dist_km = %.3f (stops: %d)\n", cor(d$d_check, d$dist_km), nrow(d)))
  ## municipal housing: share of catchment residents in municipal buildings (script 37)
  d <- merge(d, fread(file.path(s2, "municipal_catchments_network.csv"), encoding = "UTF-8")[, .(stop_name, muni_share)], by = "stop_name", all.x = TRUE)
  ## deprivation by round and the osiedle mix of each catchment (script 62)
  d <- merge(d, fread(file.path(s2, "catchment_dep_rounds.csv"), encoding = "UTF-8")[, .(stop_name, dep3, dep2019, dep2021, top_share)], by = "stop_name", all.x = TRUE)
  ## interchange use and estimated connections per 100 rides (transfer rates of script 47)
  RT <- fread(file.path(s2, "hub_transfer_rates.csv"))[w == 10 & !is.na(T_eb), .(d = hub, p = period, T_eb)]
  od[, p := PER(hour_d)]
  od <- merge(od, RT, by = c("d", "p"), all.x = TRUE)[is.na(T_eb), T_eb := 0]
  HC <- od[, .(boardings = sum(flow_c), hub_share = 100 * sum(flow[d %in% HUBS]) / sum(flow), conn_100 = 100 * sum(flow * T_eb) / sum(flow)), by = .(stop_name = o)][boardings >= 2000]

  ## ── measures ─────────────────────────────────────────────────────────────
  I  <- fread(file.path(s2, "stage2_stop_indices.csv"), encoding = "UTF-8")[, .(stop_name, volume, tidal, onset_b, conc_b, sat_b, summer)]
  DP <- fread(file.path(s2, vtag("dest_mech_profiles.csv")), encoding = "UTF-8")[transfers == "none" & comp == "observed" & window == "weekday all day",
          .(stop_name, ride_km, dest_dep, to_centre = 100 * to_centre, same_osiedle = 100 * same_osiedle)]
  RO <- fread(file.path(s2, vtag("routes_origin.csv")), encoding = "UTF-8")[comp == "observed" & window == "Mon-Thu", .(stop_name = o, eff_dest, dir_conc)]
  AC <- fread(file.path(s2, "access_origin.csv"), encoding = "UTF-8")[window == "Weekday 07-09 (morning peak)", .(stop_name, jobs30 = 100 * jobs30, transfers_jobs, t_centre)]
  BH <- dcast(fread(file.path(s2, "behaviour_stop_measures.csv"), encoding = "UTF-8")[measure %in% c("bridge", "irregularity", "dependence", "weather_rain", "night_vs_day")],
              stop_name ~ measure, value.var = "value")
  X <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(d, I, DP, HC[, .(stop_name, hub_share, conn_100)], RO, AC, BH))
  M <- data.table(
    msr = c("volume", "tidal", "onset_b", "conc_b", "sat_b", "summer",
            "ride_km", "dest_dep", "to_centre", "same_osiedle",
            "eff_dest", "dir_conc", "hub_share", "conn_100",
            "jobs30", "transfers_jobs", "t_centre",
            "bridge", "irregularity", "dependence", "weather_rain", "night_vs_day"),
    dim = rep(c("How much and when", "Where to", "Routes and transfers", "Access by the timetable", "Responses"), c(6, 4, 4, 3, 5)),
    label = c("Volume per departure", "Tidal index (morning boarding share)", "Boardings before 06:00", "Morning-peak concentration", "Saturday ratio", "Summer ratio",
              "Ride length (km)", "Deprivation of destinations (z)", "Alighting near the centre (%)", "Alighting in own osiedle (%)",
              "Effective number of destinations", "Directional concentration (0-1)", "Alighting at interchanges (%)", "Connections per 100 rides",
              "City jobs within 30 min (%)", "Transfers to jobs", "Time to the centre (min)",
              "Bridge days (log ratio)", "Day-to-day irregularity", "Riders kept when service is cut", "Rain (log riders per departure)", "Night relative to day volume"))
  M[, serv_in := dim != "How much and when"]   # primary specification: service at the stop in position, except for the same-trip usage indices

  ## ── inference helpers ────────────────────────────────────────────────────
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  source(file.path(REPO, "R/cr_inference.R"))   # crinf(), wcr()
  ## position = location (distance spline and sector) + service at the stop; variants without service
  LOC  <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2")
  SERV <- c("l_lines", "l_dep", "tram_share", "direct_centre")
  LUSE <- c("b1", "b2", "b3", "b4", "landuse_mix", "fsi", "wc_tree")
  POPV <- c("share_0_18", "share_61plus", "lres", "larea")
  NET <- c(LOC, SERV)
  LUP <- c(NET, LUSE, POPV)
  LUP_NS <- c(LOC, LUSE, POPV)
  STEPS <- list(raw = character(0), network = NET, landuse = LUP, network_noserv = LOC, landuse_noserv = LUP_NS)
  q3 <- quantile(d$dep, c(1/3, 2/3))
  X[, third := fcase(dep <= q3[1], "least", dep >= q3[2], "most", default = "middle")]
  fwrite(X[, .(stop_name, dep, third)], file.path(s2, vtag("gap_thirds_stops.csv")))   # the comparison groups, for the map of script 60
  saveRDS(X, file.path(s2, vtag("gap_stops.rds")))                                     # stops, covariates and measures, for script 63
  prep <- function(y) {   # standardised measure, exposure and covariates on the stops where the measure exists
    Z <- X[is.finite(get(y))]
    Z[, `:=`(yy = zs(get(y)), dz = zs(dep), mz = zs(muni_share), d3 = zs(dep3), d19 = zs(dep2019), d21 = zs(dep2021))]
    for (v in c(SERV, LUSE, POPV)) Z[, (v) := zs(get(v))]
    Z[]
  }
  RELS <- c(0.75, 0.60)   # reliability of the two-round score (Spearman-Brown on r = 0.60), and a pessimistic value

  ## ── the ladder ───────────────────────────────────────────────────────────
  set.seed(20261009)
  RES <- rbindlist(lapply(seq_len(nrow(M)), function(i) {
    y <- M$msr[i]; Z <- prep(y)
    steps <- c(lapply(STEPS, function(ctrl) list("dz", ctrl)), list(municipal = list("mz", c(LUP, "factor(osiedle)"))))
    rbindlist(lapply(names(steps), function(st) { m <- steps[[st]]
      ZZ <- Z[is.finite(get(m[[1]]))]
      f <- lm(reformulate(c(m[[1]], m[[2]]), "yy"), ZZ); ci <- crinf(f, ZZ$osiedle, m[[1]])
      r2 <- if (length(m[[2]]) && m[[1]] == "dz") summary(lm(reformulate(m[[2]], "dz"), ZZ))$r.squared else 0
      data.table(step = st, est = ci$b, se = ci$se2, df = ci$df, p = 2 * pt(-abs(ci$b / ci$se2), ci$df),
                 lo = ci$b - qt(0.975, ci$df) * ci$se2, hi = ci$b + qt(0.975, ci$df) * ci$se2,
                 lo90 = ci$b - qt(0.95, ci$df) * ci$se2, hi90 = ci$b + qt(0.95, ci$df) * ci$se2,
                 p_g2 = 2 * pt(-abs(ci$b / ci$se2), ci$G - 2), p_cv3 = 2 * pt(-abs(ci$b / ci$se3), ci$G - 1), p_cv1 = 2 * pt(-abs(ci$b / ci$se1), ci$G - 1),
                 p_wcr = if (st == "municipal") NA_real_ else wcr(f, ZZ$osiedle, m[[1]]),
                 G = ci$G, Gstar = ci$Gstar, lev_max = ci$lev_max, lev_os = ci$lev_os, r2_ctrl = r2)
    }))[, `:=`(msr = y, n = nrow(Z), least = Z[third == "least", median(get(y))], most = Z[third == "most", median(get(y))])]
  }))
  RES <- merge(RES, M, by = "msr")
  ## primary specification: no service for the same-trip usage indices
  PRIM <- rbind(RES[(serv_in & step == "network") | (!serv_in & step == "network_noserv")][, step := "position"],
                RES[(serv_in & step == "landuse") | (!serv_in & step == "landuse_noserv")][, step := "full"])
  RES <- rbind(RES, PRIM)
  RES[, q := p.adjust(p, "BH"), by = step]
  RES[, q_g2 := p.adjust(p_g2, "BH"), by = step]
  ## measurement error: attenuation factor of each slope under reliability R, and the disattenuated slope and interval
  for (R in RELS) { tag <- sprintf("%02d", round(100 * R))
    RES[, (paste0("lambda_", tag)) := fifelse(step == "raw", R, (R - r2_ctrl) / (1 - r2_ctrl))]
    RES[step == "municipal", (paste0("lambda_", tag)) := NA_real_]
    RES[, (paste0("est_dis_", tag)) := est / get(paste0("lambda_", tag))]
    RES[, (paste0("lo90_dis_", tag)) := lo90 / get(paste0("lambda_", tag))][, (paste0("hi90_dis_", tag)) := hi90 / get(paste0("lambda_", tag))] }
  ## equivalence (two one-sided tests at 5%): 90% interval within +-DELTA; "differs" when the 95% interval excludes zero
  DELTA <- 0.20
  eqc <- function(lo, hi, lo90, hi90) fifelse(lo > 0 | hi < 0, fifelse(lo90 > -DELTA & hi90 < DELTA, "differs, smaller than 0.2", "differs"),
                                             fifelse(lo90 > -DELTA & hi90 < DELTA, "equivalent", "inconclusive"))
  RES[, equiv := eqc(lo, hi, lo90, hi90)]
  RES[, equiv_dis75 := eqc(lo / lambda_75, hi / lambda_75, lo90_dis_75, hi90_dis_75)]
  fwrite(RES, file.path(s2, vtag("gap_ladder.csv")))

  W <- dcast(RES, dim + label + msr + n + least + most ~ step, value.var = c("est", "p", "q"))
  W[, msr := factor(msr, levels = M$msr)]; setorder(W, msr)
  cat("Slope on deprivation (SD per SD): raw / position / full (primary specification); CR2-Satterthwaite p; medians least / most deprived third\n")
  print(W[, .(label, n, least = signif(least, 3), most = signif(most, 3), raw = sprintf("%.2f (p=%.3f)", est_raw, p_raw), position = sprintf("%.2f (p=%.3f)", est_position, p_position),
              full = sprintf("%.2f (p=%.3f, q=%.2f)", est_full, p_full, q_full), municipal = sprintf("%.2f (p=%.3f)", est_municipal, p_municipal))], row.names = FALSE)
  cat("\nInference checks at the full step: Satterthwaite df, p (CR2-Satterthwaite / old t(G-2) / CV3 / WCR bootstrap), G*, largest partial leverage\n")
  IN <- RES[step == "full", .(msr, est = round(est, 2), df = round(df, 1), p = signif(p, 2), p_g2 = signif(p_g2, 2), p_cv3 = signif(p_cv3, 2), p_wcr = signif(p_wcr, 2),
                              q = signif(q, 2), G, Gstar = round(Gstar, 1), lev_max = round(lev_max, 2), lev_os)]
  print(IN, row.names = FALSE)
  fwrite(RES[, .(msr, step, est, se, df, p, p_g2, p_cv1, p_cv3, p_wcr, q, q_g2, G, Gstar, lev_max, lev_os, r2_ctrl)], file.path(s2, vtag("gap_inference.csv")))
  cat("\nMeasures with q < 0.05 at each step (CR2-Satterthwaite):\n"); print(RES[q < 0.05, .(step, label, est = round(est, 2), p = signif(p, 2), q = signif(q, 2))][order(step)], row.names = FALSE)
  cat("\nMeasurement error at the full step: share of the score's variance explained by the controls, attenuation factor and disattenuated slope (R = 0.75; R = 0.60)\n")
  print(RES[step %in% c("raw", "position", "full"), .(msr, step, est = round(est, 2), r2_ctrl = round(r2_ctrl, 2), lambda_75 = round(lambda_75, 2), dis_75 = round(est_dis_75, 2),
                                                      lambda_60 = round(lambda_60, 2), dis_60 = round(est_dis_60, 2))][step == "full"], row.names = FALSE)
  cat("\nEquivalence at the full step (margin +-0.20 SD per SD; 90% interval):\n")
  print(RES[step == "full", .(msr, est = round(est, 2), lo90 = round(lo90, 2), hi90 = round(hi90, 2), equiv, equiv_dis75)], row.names = FALSE)
  print(RES[step == "full", .N, by = equiv]); print(RES[step == "full", .N, by = equiv_dis75])

  ## ── robustness of the exposure: three-round mean, single rounds, pure catchments ──
  ROB <- rbindlist(lapply(seq_len(nrow(M)), function(i) {
    y <- M$msr[i]; Z <- prep(y); ctrl <- if (M$serv_in[i]) LUP else LUP_NS
    one <- function(ZZ, ex, lab) { ZZ <- ZZ[is.finite(get(ex))]; f <- lm(reformulate(c(ex, ctrl), "yy"), ZZ); ci <- crinf(f, ZZ$osiedle, ex)
      data.table(msr = y, variant = lab, est = ci$b, p = 2 * pt(-abs(ci$b / ci$se2), ci$df), n = nrow(ZZ)) }
    rbind(one(Z, "dz", "two-round mean (main)"), one(Z, "d3", "three-round mean"), one(Z, "d19", "2019 round"), one(Z, "d21", "2021 round"),
          one(Z[top_share >= 0.9], "dz", "pure catchments (>= 90% one osiedle)"))
  }))
  ## connections per 100 rides with transfer rates estimated for 5-, 10- and 15-minute windows (script 47): unshrunk rates
  ## clipped to 0-1 (script 47 shrinks only the 10-minute rates), so the three windows are compared like with like
  for (wn in c(5, 10, 15)) {
    RTw <- fread(file.path(s2, "hub_transfer_rates.csv"), encoding = "UTF-8")[w == wn & !is.na(T), .(d = hub, p = period, Tw = pmin(pmax(T, 0), 1))]
    HCw <- merge(od, RTw, by = c("d", "p"), all.x = TRUE)[is.na(Tw), Tw := 0][, .(boardings = sum(flow_c), cw = 100 * sum(flow * Tw) / sum(flow)), by = .(stop_name = o)][boardings >= 2000]
    Z <- prep("conn_100"); Z <- merge(Z, HCw[, .(stop_name, cw)], by = "stop_name"); Z[, yy := zs(cw)]
    for (st in c("raw", "full")) { f <- lm(reformulate(c("dz", if (st == "raw") character(0) else LUP), "yy"), Z); ci <- crinf(f, Z$osiedle, "dz")
      ROB <- rbind(ROB, data.table(msr = "conn_100", variant = sprintf("transfer window %d min, unshrunk (%s)", wn, st), est = ci$b, p = 2 * pt(-abs(ci$b / ci$se2), ci$df), n = nrow(Z))) }
    cat(sprintf("Connections per 100 rides, %d-minute window: medians least / most %.2f / %.2f\n", wn, Z[third == "least", median(cw)], Z[third == "most", median(cw)]))
  }
  fwrite(ROB, file.path(s2, vtag("gap_measurement.csv")))
  cat("\nFull-step slope with other versions of the exposure:\n")
  print(dcast(ROB, msr ~ variant, value.var = "est")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)][order(match(msr, M$msr))], row.names = FALSE)

  ## ── leave one osiedle out (full step) ────────────────────────────────────
  LOO <- rbindlist(lapply(seq_len(nrow(M)), function(i) {
    y <- M$msr[i]; Z <- prep(y); ctrl <- if (M$serv_in[i]) LUP else LUP_NS
    rbindlist(lapply(unique(Z$osiedle), function(o) { ZZ <- Z[osiedle != o]; f <- lm(reformulate(c("dz", ctrl), "yy"), ZZ); ci <- crinf(f, ZZ$osiedle, "dz")
      data.table(msr = y, dropped = o, est = ci$b, p = 2 * pt(-abs(ci$b / ci$se2), ci$df)) }))
  }))
  fwrite(LOO, file.path(s2, vtag("gap_loo.csv")))
  cat("\nLeave-one-osiedle-out, full step: range of slopes and the osiedle whose omission moves the slope most\n")
  LS <- merge(LOO, RES[step == "full", .(msr, est_all = est, p_all = p)], by = "msr")[, .(est_all = round(est_all[1], 2), min = round(min(est), 2), max = round(max(est), 2),
         most_influential = dropped[which.max(abs(est - est_all))], p_max = signif(max(p), 2), p_min = signif(min(p), 2)), by = msr][order(match(msr, M$msr))]
  print(LS, row.names = FALSE)
  FOC <- c("volume", "dest_dep", "dir_conc", "night_vs_day")
  gl <- ggplot(merge(LOO[msr %in% FOC], M[, .(msr, label)], by = "msr"), aes(est, reorder(dropped, est))) + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_point(size = 1.1, colour = "#2a78d6") + geom_vline(data = merge(RES[step == "full" & msr %in% FOC, .(msr, est)], M[, .(msr, label)], by = "msr"), aes(xintercept = est), colour = "#e34948", linetype = "dashed") +
    facet_wrap(~label, nrow = 1, scales = "free_x") + labs(x = "Slope net of position, land use and population, with one osiedle left out (SD per SD)", y = "Osiedle left out") +
    theme_minimal(base_size = 8.5) + theme(panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, vtag("Fig_gap_loo.png")), gl, width = 12, height = 7.5, dpi = 220, bg = "#ffffff")

  ## ── position on the network of the least and most deprived thirds ───────
  PM <- X[third != "middle", .(stops = .N, dist_km = median(dist_km), lines = median(exp(l_lines)), departures = median(exp(l_dep)),
                               tram_share = median(tram_share), direct_centre = median(direct_centre)), by = third]
  cat("\nPosition of the thirds (medians; departures = Mon-Thu day departures, March-June and October):\n"); print(PM)
  fwrite(PM, file.path(s2, vtag("gap_position_thirds.csv")))

  ## ── decomposition (Gelbach 2016), joint osiedle bootstrap ────────────────
  ## raw slope = remainder (slope net of all covariates of the primary specification) + the parts accounted for by
  ## location (distance spline and sector), service at the stop, land use and population; the bootstrap also
  ## returns the position-step slope, so the headline ratios get intervals from the same draws
  SECT <- c("sec_s1", "sec_c1", "sec_s2", "sec_c2")
  gel <- function(Z, ctrl, pos) {
    f <- lm(reformulate(c("dz", ctrl), "yy"), Z); b <- coef(f); b[is.na(b)] <- 0
    Wm <- model.matrix(f); Wm <- Wm[, setdiff(colnames(Wm), c("(Intercept)", "dz")), drop = FALSE]
    part <- apply(Wm, 2, function(w) cov(w, Z$dz) / var(Z$dz)) * b[colnames(Wm)]
    grp <- ifelse(grepl("dist_km", colnames(Wm), fixed = TRUE) | colnames(Wm) %in% SECT, "location",
                  ifelse(colnames(Wm) %in% SERV, "service", ifelse(colnames(Wm) %in% LUSE, "landuse", "population")))
    pt <- tapply(part, grp, sum); g0 <- function(k) if (k %in% names(pt)) pt[[k]] else 0
    c(raw = coef(lm(yy ~ dz, Z))[["dz"]], location = g0("location"), service = g0("service"), landuse = g0("landuse"), population = g0("population"),
      remainder = b[["dz"]], position = coef(lm(reformulate(c("dz", pos), "yy"), Z))[["dz"]])
  }
  OS <- sort(unique(X$osiedle)); NB <- 1999
  RS <- replicate(NB, sample(OS, replace = TRUE), simplify = FALSE)
  BOOT <- list()
  GB <- rbindlist(lapply(seq_len(nrow(M)), function(i) {
    y <- M$msr[i]; Z <- prep(y); ctrl <- if (M$serv_in[i]) LUP else LUP_NS; pos <- if (M$serv_in[i]) NET else LOC
    e <- gel(Z, ctrl, pos); stopifnot(abs(e[["raw"]] - sum(e[c("location", "service", "landuse", "population", "remainder")])) < 1e-8)
    ZL <- split(Z, by = "osiedle")
    bt <- t(vapply(RS, function(k) gel(rbindlist(ZL[k[k %in% names(ZL)]]), ctrl, pos), numeric(7)))
    BOOT[[y]] <<- bt
    data.table(msr = y, part = names(e), est = e, lo = apply(bt, 2, quantile, 0.025), hi = apply(bt, 2, quantile, 0.975), sd = apply(bt, 2, sd))
  }))
  fwrite(GB, file.path(s2, vtag("gap_decomposition.csv")))
  GW <- dcast(GB[part != "position"], msr ~ part, value.var = "est")
  cat("\nDecomposition of the raw slope (SD per SD; primary specification):\n"); print(merge(M[, .(msr, label)], GW, by = "msr")[order(match(msr, M$msr)), .(label, raw = round(raw, 2), location = round(location, 2),
    service = round(service, 2), landuse = round(landuse, 2), population = round(population, 2), remainder = round(remainder, 2))], row.names = FALSE)
  cat("Parts whose 95% bootstrap interval excludes zero:\n"); print(GB[!part %in% c("raw", "remainder", "position") & (lo > 0 | hi < 0), .(msr, part, est = round(est, 2), lo = round(lo, 2), hi = round(hi, 2))], row.names = FALSE)

  ## Romano-Wolf step-down p-values for the 22 remainders (bootstrap distribution of the centred, studentised remainders)
  rem_hat <- GB[part == "remainder", setNames(est, msr)]; rem_sd <- GB[part == "remainder", setNames(sd, msr)]
  TB <- sapply(names(rem_hat), function(k) abs(BOOT[[k]][, "remainder"] - rem_hat[[k]]) / rem_sd[[k]])
  t_obs <- abs(rem_hat) / rem_sd; ord <- order(-t_obs); prw <- numeric(length(t_obs))
  for (s in seq_along(ord)) { idx <- ord[s:length(ord)]; mx <- apply(TB[, idx, drop = FALSE], 1, max); prw[ord[s]] <- mean(mx >= t_obs[ord[s]]) }
  prw <- cummax(prw[ord])[order(ord)]
  RW <- data.table(msr = names(rem_hat), remainder = rem_hat, p_boot = 2 * pmin(sapply(names(rem_hat), function(k) mean(BOOT[[k]][, "remainder"] <= 0)),
                                                                                   sapply(names(rem_hat), function(k) mean(BOOT[[k]][, "remainder"] >= 0))), p_rw = prw)
  RES <- merge(RES, RW[, .(msr, step = "full", p_rw)], by = c("msr", "step"), all.x = TRUE)
  fwrite(RES, file.path(s2, vtag("gap_ladder.csv")))
  cat("\nRemainders: Romano-Wolf step-down p (familywise across 22 measures):\n"); print(RW[order(p_rw), .(msr, remainder = round(remainder, 2), p_boot = signif(p_boot, 2), p_rw = signif(p_rw, 2))], row.names = FALSE)

  ## ── headline: how much of the gaps position removes, usage and access separately, with joint-bootstrap intervals ──
  ## access is computed from the timetable and the stop's location, so position explains it largely by construction
  H <- dcast(RES[step %in% c("raw", "position", "full", "network", "network_noserv")], msr + dim ~ step, value.var = "est")
  H[, group := fifelse(dim == "Access by the timetable", "access", "usage")]
  hl <- function(raw, pos, grp, sel15) { u <- grp == "usage"; a <- !u
    c(usage_agg_removed = 1 - sum(abs(pos[u])) / sum(abs(raw[u])), usage_sum_raw = sum(abs(raw[u])), usage_sum_pos = sum(abs(pos[u])),
      usage_med_removed_015 = median((1 - pos / raw)[u & sel15]), access_agg_removed = 1 - sum(abs(pos[a])) / sum(abs(raw[a]))) }
  setorder(H, msr); sel15 <- abs(H$raw) >= 0.15
  hl0 <- hl(H$raw, H$position, H$group, sel15)
  hb <- t(sapply(seq_len(NB), function(b) { r <- sapply(H$msr, function(k) BOOT[[k]][b, "raw"]); p <- sapply(H$msr, function(k) BOOT[[k]][b, "position"]); hl(r, p, H$group, sel15) }))
  HL <- data.table(quantity = names(hl0), est = hl0, lo = apply(hb, 2, quantile, 0.025), hi = apply(hb, 2, quantile, 0.975))
  fwrite(HL, file.path(s2, vtag("gap_headline.csv")))
  cat("\nHeadline (position = primary specification; 95% joint osiedle-bootstrap intervals):\n"); print(HL, digits = 2)
  cat("Usage measures with |raw| >= 0.15:", paste(H[sel15 & group == "usage", msr], collapse = ", "), "\n")

  ## ── figure ───────────────────────────────────────────────────────────────
  F <- RES[step %in% c("raw", "position", "full", "municipal")]
  F[, `:=`(label = factor(label, levels = rev(M$label)), dim = factor(dim, levels = unique(M$dim)),
           step = factor(step, levels = c("raw", "position", "full", "municipal"),
                         labels = c("Deprivation: raw", "Net of position on the network", "Net of land use and population too", "Municipal-housing share, within osiedla")))]
  pal <- c("Deprivation: raw" = "#2a78d6", "Net of position on the network" = "#eb6834", "Net of land use and population too" = "#1baf7a", "Municipal-housing share, within osiedla" = "#4a3aa7")
  g <- ggplot(F, aes(est, label, colour = step, group = factor(step, levels = rev(levels(step))))) +
    annotate("rect", xmin = -DELTA, xmax = DELTA, ymin = -Inf, ymax = Inf, fill = "#f0efec") + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.75), size = 0.25, linewidth = 0.45) +
    facet_grid(dim ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_colour_manual(values = pal, name = NULL) +
    labs(x = "Slope: SD of the measure per SD of the indicator (95% interval, CR2-Satterthwaite by osiedle; band: +-0.2)", y = NULL) + guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
                                          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1),
                                          panel.spacing.y = unit(0.6, "lines"), plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, vtag("Fig_gap_ladder.png")), g, width = 11, height = 11, dpi = 220, bg = "#ffffff")

  ## ── LaTeX rows for the combined table ───────────────────────────────────
  W <- dcast(RES, dim + label + msr + n + least + most ~ step, value.var = c("est", "p", "q"))
  W[, msr := factor(msr, levels = M$msr)]; setorder(W, msr)
  nz <- function(s) sub("^-(0\\.0+)$", "\\1", s)   # no negative zeros
  fnum <- function(k, v) { if (!is.finite(v)) return("--"); if (k %in% c("to_centre", "same_osiedle", "hub_share", "jobs30")) nz(sprintf("%.0f", v))
    else if (k %in% c("ride_km", "eff_dest", "t_centre")) nz(sprintf("%.1f", v)) else nz(sprintf("%.2f", v)) }
  fs <- function(e, p) paste0(gsub("-", "$-$", nz(sprintf("%.2f", e))), ifelse(p < 0.01, "$^{**}$", ifelse(p < 0.05, "$^{*}$", "")))
  dir.create(file.path(out_dir, "manuscript/tables"), showWarnings = FALSE)
  L <- c()
  for (dm in unique(M$dim)) {
    L <- c(L, sprintf("\\multicolumn{9}{l}{\\emph{%s}} \\\\", dm))
    for (k in M[dim == dm, msr]) { r <- W[msr == k]
      lab <- gsub("%", "\\\\%", r$label)
      ## parts: * when the 95% bootstrap interval excludes zero; raw and remainder: CR2-Satterthwaite p, as in the figure
      fb <- function(pt) { if (pt == "service" && !M[msr == k, serv_in]) return("--")
        g <- GB[msr == k & part == pt]; paste0(gsub("-", "$-$", nz(sprintf("%.2f", g$est))), ifelse(g$lo > 0 | g$hi < 0, "$^{*}$", "")) }
      L <- c(L, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", lab, gsub("-", "$-$", fnum(k, r$least)), gsub("-", "$-$", fnum(k, r$most)),
                        fs(r$est_raw, r$p_raw), fb("location"), fb("service"), fb("landuse"), fb("population"), fs(r$est_full, r$p_full))) }
  }
  writeLines(L, file.path(out_dir, "manuscript/tables", vtag("tab_usage_differences.tex")), useBytes = TRUE)
  cat("Saved gap_ladder.csv, gap_inference.csv, gap_measurement.csv, gap_loo.csv, gap_decomposition.csv, gap_headline.csv, Fig_gap_ladder.png, Fig_gap_loo.png, tab_usage_differences.tex\n")
})

## ── 55 · The extremes of deprivation ───────────────────────────────────────────────────────────────────────
section("55", "The extremes of deprivation", {
  ## Do the extremes tell a different story? The 5% most and 5% least deprived
  ## stop catchments (10% as a check) compared with the thirds used elsewhere,
  ## on every stop-level measure of the study. (Exploratory, descriptive.)
  ##
  ## Deprivation is measured for 42 osiedla, so the extremes come from a handful
  ## of osiedla (listed); differences between them are comparisons of a few
  ## places, reported descriptively. Check: the 5% extremes without the two
  ## osiedla whose rates are driven by single institutions (Główna, Krzesiny-
  ## Pokrzywno-Garaszewo; as in Stage 1).
  ## places, reported descriptively. Every measure is standardised over all stops
  ## with residents (not flagged). For each group: median of the measure (raw),
  ## median residual from a spline of distance to the centre (net of location),
  ## and median residual from distance plus land use and population (net of
  ## location and land use: residential balance, multi/single, services/
  ## production, retail/workplaces, mix, floor-area density, tree cover, age
  ## shares, log residents). Differences are most minus least deprived group.
  ## Measures: same-trip usage indices (script 42), destination profile and ride
  ## length (46), route heterogeneity (51), accessibility (54), responses to leave
  ## days, service cuts and weather (52).
  ## Outputs: output/stage2/extremes_groups.csv, extremes_differences.csv;
  ##          output/figures/Fig_extremes.png

  suppressPackageStartupMessages({library(data.table); library(splines); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[residents > 0 & flagged == FALSE,
          .(stop_name, dep, dist_km, osiedle, share_0_18, share_61plus, lres = log1p(as.numeric(residents)))]
  td <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")[, .(stop_name, b1 = ilr_res_vs_nonres, b2 = ilr_multi_vs_single,
          b3 = ilr_serv_vs_prod, b4 = ilr_retail_vs_work, landuse_mix, fsi = log(fsi + 0.01), wc_tree)]
  d <- merge(cv, td, by = "stop_name")
  d <- d[complete.cases(d)]
  q <- quantile(d$dep, c(0.05, 0.10, 1/3, 2/3, 0.90, 0.95))
  d[, `:=`(g5 = fcase(dep <= q[1], "least", dep >= q[6], "most", default = NA_character_),
           g10 = fcase(dep <= q[2], "least", dep >= q[5], "most", default = NA_character_),
           g33 = fcase(dep <= q[3], "least", dep >= q[4], "most", default = NA_character_))]
  ## check: extremes without the two osiedla whose rates are driven by single institutions (as in Stage 1)
  INST <- c("Główna", "Krzesiny-Pokrzywno-Garaszewo")
  qx <- quantile(d[!osiedle %in% INST, dep], c(0.05, 0.95))
  d[, g5x := fifelse(osiedle %in% INST, NA_character_, fcase(dep <= qx[1], "least", dep >= qx[2], "most", default = NA_character_))]
  cat("Stops with residents (not flagged):", nrow(d), "\n")
  for (g in c("g5", "g10", "g5x")) {
    cat("\n", g, "groups: osiedla and distance to the centre\n")
    print(d[!is.na(get(g)), .(stops = .N, dist_km = round(median(dist_km), 1)), by = .(group = get(g), osiedle)][order(group, -stops)], row.names = FALSE)
  }

  ## ── measures ─────────────────────────────────────────────────────────────
  I  <- fread(file.path(s2, "stage2_stop_indices.csv"), encoding = "UTF-8")[, .(stop_name, volume, tidal, onset_b, conc_b, sat_b, summer)]
  DP <- fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & comp == "observed" & window == "weekday all day",
          .(stop_name, ride_km, dest_dep, to_centre, same_osiedle)]
  RT <- fread(file.path(s2, "routes_origin.csv"), encoding = "UTF-8")[comp == "observed" & window == "Mon-Thu", .(stop_name = o, eff_dest, dir_conc)]
  AC <- fread(file.path(s2, "access_origin.csv"), encoding = "UTF-8")[window == "Weekday 07-09 (morning peak)", .(stop_name, jobs30, transfers_jobs, t_centre)]
  BH <- dcast(fread(file.path(s2, "behaviour_stop_measures.csv"), encoding = "UTF-8")[measure %in% c("bridge", "irregularity", "dependence", "weather_rain", "night_vs_day")],
              stop_name ~ measure, value.var = "value")
  X <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(d, I, DP, RT, AC, BH))
  LAB <- c(volume = "Usage: volume per departure", tidal = "Usage: tidal index", onset_b = "Usage: early boardings (before 06:00)",
           conc_b = "Usage: morning-peak concentration", sat_b = "Usage: Saturday ratio", summer = "Usage: summer ratio",
           ride_km = "Destinations: ride length", dest_dep = "Destinations: deprivation", to_centre = "Destinations: share to the centre",
           same_osiedle = "Destinations: share in own osiedle", eff_dest = "Routes: effective destinations", dir_conc = "Routes: directional concentration",
           jobs30 = "Access: jobs within 30 min", transfers_jobs = "Access: transfers to jobs", t_centre = "Access: time to the centre",
           bridge = "Responses: bridge days", irregularity = "Responses: day-to-day irregularity", dependence = "Responses: riders kept when service is cut",
           weather_rain = "Responses: rain", night_vs_day = "Responses: night vs day")
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  LU <- c("b1", "b2", "b3", "b4", "landuse_mix", "fsi", "wc_tree", "share_0_18", "share_61plus", "lres")
  RES <- rbindlist(lapply(names(LAB), function(y) {
    Z <- X[is.finite(get(y))]; Z[, z := zs(get(y))]
    Z[, r_loc := resid(lm(z ~ ns(dist_km, 3), Z))]
    Z[, r_lu := resid(lm(reformulate(c("ns(dist_km, 3)", LU), "z"), Z))]
    rbindlist(lapply(c("g5", "g5x", "g10", "g33"), function(g) {
      m <- Z[!is.na(get(g)), .(n = .N, osiedla = uniqueN(osiedle), raw = median(z), loc = median(r_loc), lu = median(r_lu)), by = .(group = get(g))]
      w <- dcast(melt(m, id.vars = "group"), variable ~ group)
      data.table(measure = y, cut = g, n_least = w[variable == "n", least], n_most = w[variable == "n", most],
                 os_least = w[variable == "osiedla", least], os_most = w[variable == "osiedla", most],
                 raw = w[variable == "raw", most - least], net_location = w[variable == "loc", most - least], net_location_landuse = w[variable == "lu", most - least])
    }))
  }))
  fwrite(RES, file.path(s2, "extremes_differences.csv"))
  fwrite(d[, .(stop_name, osiedle, dep, dist_km, g5, g5x, g10, g33)], file.path(s2, "extremes_groups.csv"))
  cat("\nMost minus least deprived group, in SD of the measure (medians): raw / net of location / net of location and land use\n")
  P <- dcast(RES, measure ~ cut, value.var = c("raw", "net_location", "net_location_landuse"))
  P[, measure := factor(LAB[measure], levels = LAB)]; setorder(P, measure)
  print(P[, .(measure, thirds_raw = round(raw_g33, 2), thirds_loc = round(net_location_g33, 2), thirds_lu = round(net_location_landuse_g33, 2), p5x_raw = round(raw_g5x, 2), p5x_loc = round(net_location_g5x, 2), p5x_lu = round(net_location_landuse_g5x, 2),
              p10_raw = round(raw_g10, 2), p10_loc = round(net_location_g10, 2), p10_lu = round(net_location_landuse_g10, 2),
              p5_raw = round(raw_g5, 2), p5_loc = round(net_location_g5, 2), p5_lu = round(net_location_landuse_g5, 2))], row.names = FALSE)
  cat("\nStops per group (5% cut): ", RES[cut == "g5", paste(unique(paste(n_least, n_most)), collapse = ", ")], "| osiedla:", RES[cut == "g5", paste(unique(paste(os_least, os_most)), collapse = ", ")], "\n")
  ## location share of the raw difference, by cut (median over measures with |raw| >= 0.3)
  cat("\nHow much of the raw gap location removes (median of 1 - |net of location| / |raw|, measures with |raw| >= 0.3):\n")
  print(RES[abs(raw) >= 0.3, .(measures = .N, location_removes = round(median(1 - abs(net_location) / abs(raw)), 2),
                               location_landuse_removes = round(median(1 - abs(net_location_landuse) / abs(raw)), 2)), by = cut])

  ## ── figure ───────────────────────────────────────────────────────────────
  F <- melt(RES, id.vars = c("measure", "cut"), measure.vars = c("raw", "net_location", "net_location_landuse"))
  F[, `:=`(lab = factor(LAB[measure], levels = rev(LAB)),
           cut = factor(cut, levels = c("g33", "g10", "g5", "g5x"), labels = c("Thirds", "10% extremes", "5% extremes", "5% extremes without the two institution osiedla")),
           variable = factor(variable, levels = c("raw", "net_location", "net_location_landuse"), labels = c("Raw", "Net of location", "Net of location and land use")))]
  g <- ggplot(F, aes(value, lab, colour = cut, shape = cut)) + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_point(size = 2, position = position_dodge(width = 0.6)) + facet_wrap(~ variable, nrow = 1) +
    scale_colour_manual(values = c("Thirds" = "#9a9893", "10% extremes" = "#eb6834", "5% extremes" = "#e34948", "5% extremes without the two institution osiedla" = "#2a78d6"), name = NULL) +
    scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
    labs(x = "Most minus least deprived group (SD of the measure; medians)", y = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
                                          plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, "Fig_extremes.png"), g, width = 13, height = 7, dpi = 220, bg = "#ffffff")
  cat("Saved extremes_groups.csv, extremes_differences.csv, Fig_extremes.png\n")
})

## ── 59 · Why deprived areas sit where they do ──────────────────────────────────────────────────────────────
section("59", "Why deprived areas sit where they do", {
  ## Why deprived areas sit where they do on the network (manuscript §3.1).
  ## A. Service: is the better service of deprived stops more than their location
  ##    predicts? Slopes of log lines, log departures, tram share and the share of
  ##    departures on lines that pass the centre on catchment deprivation, raw,
  ##    net of a distance spline, and net of distance and sector (first two
  ##    harmonics of the bearing from the Old Market Square); CR2 by osiedle;
  ##    428 stops with residents, not flagged; term-time Monday-Thursday.
  ## B. Housing: across the 42 osiedla, how much of deprivation distance and
  ##    sector account for (R2), how central the compact mixed-use (pre-war)
  ##    blocks of the land-use typology are, how much of the link between
  ##    deprivation and distance their share accounts for (partial correlation),
  ##    and mean deprivation by octant around the centre.
  ## Outputs: output/stage2/why_position_service.csv, why_position_osiedla.csv
  suppressPackageStartupMessages({library(data.table); library(splines); library(sf)})
  out_dir <- ROOT
  s1 <- file.path(out_dir, "output/stage1"); s2 <- file.path(out_dir, "output/stage2")
  MON <- c(3:6, 10)
  cr2 <- function(fit, cl) {
    X <- model.matrix(fit); e <- resid(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]
    B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) { i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); Am <- ev$vectors %*% diag(ifelse(ev$values > 1e-8, 1 / sqrt(pmax(ev$values, 1e-8)), 0), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (Am %*% e[i]); meat <- meat + u %*% t(u) }
    setNames(sqrt(diag(B %*% meat %*% B)), colnames(X))
  }
  zs <- function(v) (v - mean(v)) / sd(v)
  ctr <- st_coordinates(st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180))
  sect <- function(x, y) { th <- atan2(y - ctr[2], x - ctr[1]); data.table(s1 = sin(th), c1 = cos(th), s2 = sin(2 * th), c2 = cos(2 * th), th = th) }

  ## ── A. Is the service advantage of deprived stops more than geometry? ───
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[residents > 0 & flagged == FALSE, .(stop_name, dep, dist_km, osiedle)]
  xy <- fread(file.path(s2, "stop_complex_attributes.csv"), encoding = "UTF-8")[, .(stop_name, x, y)]
  cv <- merge(cv, cbind(xy[, .(stop_name)], sect(xy$x, xy$y)), by = "stop_name")
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))[daytype == "monthu" & month %in% MON & night == FALSE]
  sv <- cells[, .(l_lines = log(uniqueN(NR_LINI)), l_dep = log(sum(departures)), tram_share = sum(departures[mode == "tram"]) / sum(departures)), by = stop_name]
  segdep <- cells[, .(dep = sum(departures)), by = .(o = stop_name, NR_LINI, dir)]; rm(cells); gc()
  od <- rbindlist(lapply(sprintf("%02d", MON), function(m) readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m)))[daytype == "monthu" & o != "outside", .(flow = sum(flow)), by = .(o, d, NR_LINI, dir)]))
  od <- merge(od, fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[, .(d = stop_name, centre_km = dist_km)], by = "d", all.x = TRUE)
  reach <- merge(od[, .(reach = any(centre_km <= 1.5, na.rm = TRUE)), by = .(o, NR_LINI, dir)], segdep, by = c("o", "NR_LINI", "dir"))[, .(direct_centre = sum(dep * reach) / sum(dep)), by = .(stop_name = o)]
  A <- merge(merge(cv, sv, by = "stop_name"), reach, by = "stop_name")
  A[, dz := zs(dep)]
  G <- uniqueN(A$osiedle)
  cat("A. Slope of service at the stop on catchment deprivation (SD per SD), CR2 by osiedle; stops:", nrow(A), "\n")
  SV <- rbindlist(lapply(c("l_lines", "l_dep", "tram_share", "direct_centre"), function(v) {
    A[, yy := zs(get(v))]
    rbindlist(lapply(list(raw = "dz", distance = c("dz", "ns(dist_km, 3)"), distance_sector = c("dz", "ns(dist_km, 3)", "s1", "c1", "s2", "c2")), function(rhs) {
      f <- lm(reformulate(rhs, "yy"), A); b <- coef(f)[["dz"]]; se <- cr2(f, A$osiedle)[["dz"]]
      data.table(est = b, se = se, p = 2 * pt(-abs(b / se), G - 2)) }), idcol = "adjustment")[, service := v]
  }))
  fwrite(SV, file.path(s2, "why_position_service.csv"))
  print(dcast(SV[, .(service, adjustment, v = sprintf("%5.2f (p=%.3f)", est, p))], service ~ adjustment, value.var = "v"))
  cat("  Share of the variance in service explained by distance and sector alone (R2):\n")
  for (v in c("l_lines", "l_dep", "tram_share", "direct_centre")) cat(sprintf("    %-14s %.2f\n", v, summary(lm(reformulate(c("ns(dist_km, 3)", "s1", "c1", "s2", "c2"), v), A))$r.squared))

  ## ── B. Why is deprivation central, and in which sectors? (42 osiedla) ───
  P <- fread(file.path(s1, "s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, dep, dist_km)]
  TY <- fread(file.path(s1, "s1ms_types.csv"))[, .(type, name)]
  TS <- fread(file.path(s1, "osiedle_type_shares.csv"), encoding = "UTF-8")
  cat("\nB. Types in osiedle_type_shares:", paste(sort(unique(TS$type)), collapse = ","), "| s1ms types:", paste(TY$type, TY$name, sep = "=", collapse = "; "), "\n")
  TW <- dcast(TS, osiedle ~ paste0("t", type), value.var = "share", fill = 0)
  P <- merge(P, TW, by = "osiedle")
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  pp <- st_transform(st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE), 2180)
  pp$pop <- pp[[grep("ZAMELDOWANYCH", names(pp), value = TRUE)[1]]]; pp$osiedle <- osi$name[st_nearest_feature(pp, osi)]
  pc <- data.table(st_coordinates(pp), pop = pp$pop, osiedle = pp$osiedle)[pop > 0, .(x = weighted.mean(X, pop), y = weighted.mean(Y, pop)), by = osiedle]
  P <- merge(P, cbind(pc[, .(osiedle)], sect(pc$x, pc$y)), by = "osiedle")
  cat("Check against s1ms_types share_rho (compact blocks):", round(cor(P$t1, P$dep, method = "spearman"), 3), "\n")
  r2 <- function(rhs) summary(lm(reformulate(rhs, "dep"), P))$r.squared
  cat(sprintf("R2 of deprivation across 42 osiedla: distance %.2f | sector %.2f | distance + sector %.2f | compact-block share %.2f | distance + sector + compact %.2f\n",
              r2("ns(dist_km, 3)"), r2(c("s1", "c1", "s2", "c2")), r2(c("ns(dist_km, 3)", "s1", "c1", "s2", "c2")), r2("t1"), r2(c("ns(dist_km, 3)", "s1", "c1", "s2", "c2", "t1"))))
  cat(sprintf("Spearman: compact share ~ distance %.2f; deprivation ~ distance %.2f; deprivation ~ distance net of compact share (partial, linear) %.2f\n",
              cor(P$t1, P$dist_km, method = "spearman"), cor(P$dep, P$dist_km, method = "spearman"),
              cor(resid(lm(dep ~ t1, P)), resid(lm(dist_km ~ t1, P)), method = "spearman")))
  P[, octant := cut((th * 180 / pi + 360 + 22.5) %% 360, seq(0, 360, 45), labels = c("E", "NE", "N", "NW", "W", "SW", "S", "SE"), include.lowest = TRUE)]
  print(P[, .(osiedla = .N, mean_dep = round(mean(dep), 2), compact_share = round(mean(t1), 2), dist_km = round(mean(dist_km), 1)), by = octant][order(-mean_dep)])
  print(P[order(-dep)][1:10, .(osiedle, dep = round(dep, 2), dist_km = round(dist_km, 1), compact = round(t1, 2), octant)])
  fwrite(P[, .(osiedle, dep, dist_km, compact_share = t1, estates_share = t2, single_family_share = t3, loose_share = t4, bearing_deg = (th * 180 / pi + 360) %% 360, octant)],
         file.path(s2, "why_position_osiedla.csv"))
})

## ── 63 · Ride length: line mix and estimators ──────────────────────────────────────────────────────────────
section("63", "Ride length: line mix and estimators", {
  ## Ride length under three OD estimators, split into a part fixed by the counts and a part that depends
  ## on the estimator (reviewer check of the destination estimator).
  ## On a line, the counts fix the mean ride: boardings and alightings fix the on-board load of every
  ## segment, and sum(load x segment length) / boardings is the mean ride, whatever OD matrix reproduces
  ## both margins (Furth & Navick 1992). Each origin's ride length is therefore split into
  ##   line mix     the boarding-weighted mean of the mean rides of the line x direction x period
  ##                combinations its riders board (fixed by the counts under any margin-consistent
  ##                estimator: it says which lines a stop's riders use, not where they alight), and
  ##   within line  the rest, which depends on where along each line the stop lies and on the estimator.
  ## Estimators (script 40): tsy (Tsygalnitsky = null-seed biproportional, main), dist (distance-seeded
  ## biproportional fit), prop (proportional to downstream alightings, the earlier rule).
  ## Sample and covariates: the stops of script 56 (gap_stops.rds), term-time Monday-Thursday, origins
  ## with at least 2,000 counted boardings. Slopes on deprivation: raw and net of position (location +
  ## service at the stop), CR2-Satterthwaite by osiedle.
  ## Output: output/stage2/ride_length_linemix.csv

  suppressPackageStartupMessages({library(data.table); library(splines)})
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R")); source(file.path(REPO, "R/cr_inference.R"))
  s2 <- file.path(out_dir, "output/stage2"); MON <- c(3:6, 10)
  X <- readRDS(file.path(s2, "gap_stops.rds"))
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  NET <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2", "l_lines", "l_dep", "tram_share", "direct_centre")

  OUT <- rbindlist(lapply(c("tsy", "dist", "prop"), function(v) {
    od <- rbindlist(lapply(sprintf("%02d", MON), function(m) expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))), variant = v)[daytype == "monthu" & o != "outside",
            .(flow = sum(flow), fkm = sum(fkm), flow_c = sum(flow_c)), by = .(o, NR_LINI, dir, period)]))
    lm_ <- od[, .(rbar = sum(fkm) / sum(flow)), by = .(NR_LINI, dir, period)]                 # mean ride of each line x direction x period
    od <- merge(od, lm_, by = c("NR_LINI", "dir", "period"))
    R <- od[, .(boardings_c = sum(flow_c), ride_km = sum(fkm) / sum(flow), linemix_km = sum(flow * rbar) / sum(flow)), by = .(stop_name = o)][boardings_c >= 2000]
    R[, within_km := ride_km - linemix_km]
    Z <- merge(X[, c("stop_name", "dep", "third", "osiedle", "dist_km", "sec_s1", "sec_c1", "sec_s2", "sec_c2", "l_lines", "l_dep", "tram_share", "direct_centre"), with = FALSE], R, by = "stop_name")
    Z[, dz := zs(dep)]; for (k in c("l_lines", "l_dep", "tram_share", "direct_centre")) Z[, (k) := zs(get(k))]
    sd_ride <- sd(Z$ride_km)
    rbindlist(lapply(c("ride_km", "linemix_km", "within_km"), function(y) {
      Z[, yy := get(y) / sd_ride]            # all three components on the scale of the SD of ride length, so the slopes add up
      rbindlist(lapply(list(raw = character(0), position = NET), function(ctrl) {
        f <- lm(reformulate(c("dz", ctrl), "yy"), Z); ci <- crinf(f, Z$osiedle, "dz")
        data.table(estimator = v, component = y, step = if (length(ctrl)) "position" else "raw", est = ci$b, p = 2 * pt(-abs(ci$b / ci$se2), ci$df), df = ci$df)
      }))[, `:=`(least = Z[third == "least", median(get(y))], most = Z[third == "most", median(get(y))], mean_all = mean(Z[[y]]), n = nrow(Z))]
    }))
  }))
  fwrite(OUT, file.path(s2, "ride_length_linemix.csv"))
  cat("Ride length by OD estimator: medians (km) of the least and most deprived thirds, and slopes on deprivation (in SD of ride length per SD), raw and net of position\n")
  print(dcast(OUT, estimator + component + least + most + mean_all ~ step, value.var = c("est", "p"))[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)][order(estimator, component)], row.names = FALSE)
  cat("Saved ride_length_linemix.csv\n")
})

## ── 65 · Concentration indices and access poverty ──────────────────────────────────────────────────────────
section("65", "Concentration indices and access poverty", {
  ## Equity block: how service, access and rides are distributed over residents ranked by deprivation.
  ## Concentration index (Wagstaff, Paci & van Doorslaer 1991): residents ranked from the most to the least
  ## deprived catchment (fractional rank, weighted by residents); CI = 2 cov(y, rank) / mean(y). Negative values
  ## mean the quantity is concentrated among residents of more deprived catchments. For a bounded share (jobs
  ## within 30 minutes, 0-100%) the Erreygers (2009) correction E = 4 mean(y) / (max - min) x CI is reported.
  ## Conditional on position (indirect standardisation): CI(y) - CI(y expected from position alone), where the
  ## expectation comes from a resident-weighted regression of y on deprivation and position (distance spline,
  ## sector, service at the stop) evaluated at mean deprivation. A value near zero means that position accounts
  ## for the whole concentration.
  ## Access poverty: share of residents whose stop reaches fewer than half of the resident-weighted median
  ## share of the city's jobs within 30 minutes (FGT headcount and gap), overall and by thirds of deprivation.
  ## Sample: the stops of script 56 (residents, not flagged). Intervals: 999 bootstrap resamples of osiedla.
  ## Outputs: output/stage2/equity_block.csv, equity_access_poverty.csv; output/figures/Fig_concentration.png

  suppressPackageStartupMessages({library(data.table); library(splines); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  X <- readRDS(file.path(s2, "gap_stops.rds"))
  X[, `:=`(residents = round(expm1(lres)), departures = exp(l_dep))]   # day departures, term-time Mon-Thu (indices are scale-free)
  X[, dep_per_1000 := 1000 * departures / residents]
  NETV <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2", "l_lines", "l_dep", "tram_share", "direct_centre")
  LOCV <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2")
  Y <- data.table(var = c("departures", "dep_per_1000", "jobs30", "ride_km", "conn_100"),
                  label = c("Departures at the stop", "Departures per 1,000 residents", "City jobs within 30 min (%)", "Ride length (km)", "Connections per 100 rides"),
                  bounded = c(FALSE, FALSE, TRUE, FALSE, FALSE), pos = list(LOCV, LOCV, NETV, NETV, NETV))
  wcov <- function(a, b, w) sum(w * (a - sum(w * a) / sum(w)) * (b - sum(w * b) / sum(w))) / sum(w)
  frank_w <- function(dep, w) { o <- order(-dep); cw <- cumsum(w[o]) / sum(w); r <- numeric(length(w)); r[o] <- cw - w[o] / sum(w) / 2; r }   # most deprived first
  cindex <- function(y, r, w) 2 * wcov(y, r, w) / (sum(w * y) / sum(w))
  one <- function(Z, v, bounded, pos) {
    Z <- Z[is.finite(get(v)) & residents > 0]; w <- Z$residents; r <- frank_w(Z$dep, w); y <- Z[[v]]
    ci <- cindex(y, r, w)
    f <- lm(reformulate(c("dep", pos), v), Z, weights = residents); Z0 <- copy(Z)[, dep := sum(w * dep) / sum(w)]
    yhat <- predict(f, Z0); ci_pos <- cindex(yhat, r, w)
    k <- if (bounded) 4 * (sum(w * y) / sum(w)) / 100 else 1
    c(ci = k * ci, ci_expected_from_position = k * ci_pos, ci_conditional = k * (ci - ci_pos), mean = sum(w * y) / sum(w), n = nrow(Z))
  }
  set.seed(20261009); OS <- unique(X$osiedle); ZL <- split(X, by = "osiedle")
  EQ <- rbindlist(lapply(seq_len(nrow(Y)), function(i) {
    e <- one(X, Y$var[i], Y$bounded[i], Y$pos[[i]])
    bt <- t(replicate(999, one(rbindlist(ZL[sample(OS, replace = TRUE)]), Y$var[i], Y$bounded[i], Y$pos[[i]])))
    data.table(var = Y$var[i], label = Y$label[i], index = if (Y$bounded[i]) "Erreygers" else "Wagstaff", quantity = names(e), est = e,
               lo = apply(bt, 2, quantile, 0.025), hi = apply(bt, 2, quantile, 0.975))
  }))
  fwrite(EQ, file.path(s2, "equity_block.csv"))
  cat("Concentration indices over residents ranked from the most deprived (negative = concentrated in more deprived catchments); 95% osiedle-bootstrap intervals\n")
  print(dcast(EQ[quantity %in% c("ci", "ci_expected_from_position", "ci_conditional")], label + index ~ quantity, value.var = "est")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
  print(EQ[quantity %in% c("ci", "ci_conditional"), .(label, quantity, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3))])

  ## access poverty: relative line = half the resident-weighted median of jobs within 30 minutes
  wmed <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]] }
  Z <- X[is.finite(jobs30) & residents > 0]; line <- 0.5 * wmed(Z$jobs30, Z$residents)
  Z[, `:=`(poor = jobs30 < line, gap = pmax(0, line - jobs30) / line)]
  AP <- rbind(Z[, .(group = "all", residents = sum(residents), headcount = sum(residents * poor) / sum(residents), gap = sum(residents * gap) / sum(residents))],
              Z[, .(residents = sum(residents), headcount = sum(residents * poor) / sum(residents), gap = sum(residents * gap) / sum(residents)), by = .(group = third)])
  AP[, line_pct := line]
  far <- Z[third == "most" & dist_km > 6, .(group = "most deprived third beyond 6 km", residents = sum(residents), headcount = sum(residents * poor) / sum(residents), gap = sum(residents * gap) / sum(residents), line_pct = line)]
  AP <- rbind(AP, far)
  fwrite(AP, file.path(s2, "equity_access_poverty.csv"))
  cat(sprintf("\nAccess poverty line: %.1f%% of the city's jobs within 30 minutes (half the resident-weighted median)\n", line)); print(AP, digits = 3)

  ## concentration curves
  CC <- rbindlist(lapply(c("departures", "dep_per_1000", "jobs30", "ride_km"), function(v) {
    Z <- X[is.finite(get(v)) & residents > 0][order(-dep)]
    data.table(var = v, cum_res = c(0, cumsum(Z$residents) / sum(Z$residents)), cum_y = c(0, cumsum(Z$residents * Z[[v]]) / sum(Z$residents * Z[[v]])))
  }))
  CC <- merge(CC, Y[, .(var, label)], by = "var")
  gc <- ggplot(CC, aes(cum_res, cum_y, colour = label)) + geom_abline(slope = 1, intercept = 0, colour = "#9a9893", linetype = "dashed") + geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c("Departures at the stop" = "#eb6834", "Departures per 1,000 residents" = "#4a3aa7", "City jobs within 30 min (%)" = "#2a78d6", "Ride length (km)" = "#1baf7a"), name = NULL) +
    coord_equal() + labs(x = "Cumulative share of residents, from the most to the least deprived catchment", y = "Cumulative share of the quantity") +
    theme_minimal(base_size = 10) + theme(legend.position = "right", panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, "Fig_concentration.png"), gc, width = 8.5, height = 5.2, dpi = 220, bg = "#ffffff")
  cat("Saved equity_block.csv, equity_access_poverty.csv, Fig_concentration.png\n")
})

## ── 66 · Spatial-noise placebo ─────────────────────────────────────────────────────────────────────────────
section("66", "Spatial-noise placebo", {
  ## Spatial dependence checks for the slopes of script 56 (full step, primary specification).
  ## (1) Spatial-noise placebo (Conley & Kelly 2025): simulated deprivation fields with the spatial structure of
  ##     the real score (Gaussian, exponential covariance between osiedle centroids, range chosen so that the
  ##     fields' mean Moran's I matches the observed score's) are mapped to catchments with the residents'
  ##     weights (script 62) and pushed through the same model. With calibrated inference about 5% of the
  ##     placebo slopes reach p < 0.05; a higher share means the real p-values are too small.
  ## (2) Moran's I of the 42 osiedle means of the full-step residuals, with permutation p-values.
  ## Weights: the five nearest osiedle centroids (row-standardised).
  ## Outputs: output/stage2/spatial_placebo.csv, spatial_moran.csv

  suppressPackageStartupMessages({library(data.table); library(splines); library(sf)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  source(file.path(REPO, "R/cr_inference.R"))
  X <- readRDS(file.path(s2, "gap_stops.rds")); RES <- fread(file.path(s2, "gap_ladder.csv"), encoding = "UTF-8")
  Wc <- fread(file.path(s2, "catchment_osiedle_weights.csv"), encoding = "UTF-8")
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  cen <- st_coordinates(st_centroid(st_geometry(osi))) / 1000; ON <- osi$name
  D <- as.matrix(dist(cen)); Wd <- t(apply(D, 1, function(r) { w <- numeric(length(r)); w[order(r)[2:6]] <- 1; w })); Wd <- Wd / rowSums(Wd)   # 5 nearest neighbours
  moran <- function(z) { z <- z - mean(z); as.numeric(length(z) * t(z) %*% Wd %*% z / sum(Wd) / sum(z^2)) }
  P <- fread(file.path(out_dir, "output/stage1/s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, dep)]
  I_obs <- moran(P$dep[match(ON, P$osiedle)])
  set.seed(20261009)
  sim_fields <- function(r, k) { L <- t(chol(exp(-D / r) + diag(1e-8, nrow(D)))); L %*% matrix(rnorm(nrow(D) * k), nrow(D), k) }
  grid <- c(0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16)
  mI <- sapply(grid, function(r) mean(apply(sim_fields(r, 400), 2, moran)))
  r_best <- grid[which.min(abs(mI - I_obs))]
  cat(sprintf("Observed Moran's I of the osiedle score: %.3f; simulated mean I by range (km): %s; chosen range %.1f km\n", I_obs, paste(sprintf("%.1f:%.2f", grid, mI), collapse = " "), r_best))

  ## catchment weights matrix (stops x osiedla)
  Wm <- dcast(Wc[stop_name %in% X$stop_name], stop_name ~ osiedle, value.var = "w", fill = 0)
  for (o in setdiff(ON, names(Wm))) Wm[, (o) := 0]
  Wmat <- as.matrix(Wm[, ON, with = FALSE]); rownames(Wmat) <- Wm$stop_name
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  LOC <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2"); SERV <- c("l_lines", "l_dep", "tram_share", "direct_centre")
  LUSE <- c("b1", "b2", "b3", "b4", "landuse_mix", "fsi", "wc_tree"); POPV <- c("share_0_18", "share_61plus", "lres", "larea")
  M <- unique(RES[step == "full", .(msr, dim, label, p_real = p)]); M[, serv_in := dim != "How much and when"]
  NS <- 1000; FS <- sim_fields(r_best, NS)
  PL <- rbindlist(lapply(seq_len(nrow(M)), function(i) {
    y <- M$msr[i]; Z <- X[is.finite(get(y))]; Z[, yy := zs(get(y))]; for (v in c(SERV, LUSE, POPV)) Z[, (v) := zs(get(v))]
    ctrl <- if (M$serv_in[i]) c(LOC, SERV, LUSE, POPV) else c(LOC, LUSE, POPV)
    Wz <- Wmat[Z$stop_name, , drop = FALSE]
    ps <- vapply(seq_len(NS), function(s) { Z[, dz := zs(as.vector(Wz %*% FS[, s]))]
      f <- lm(reformulate(c("dz", ctrl), "yy"), Z); ci <- crinf(f, Z$osiedle, "dz"); 2 * pt(-abs(ci$b / ci$se2), ci$df) }, numeric(1))
    ## Moran's I of osiedle-mean residuals of the real full-step model
    Z[, dz := zs(dep)]; f <- lm(reformulate(c("dz", ctrl), "yy"), Z); Z[, e := resid(f)]
    om <- Z[, .(e = mean(e)), by = osiedle]; om <- om[match(ON, osiedle)]; ok <- !is.na(om$e)
    Wsub <- Wd[ok, ok]; Wsub <- Wsub / rowSums(Wsub); zz <- om$e[ok] - mean(om$e[ok])
    mI_res <- as.numeric(sum(ok) * t(zz) %*% Wsub %*% zz / sum(Wsub) / sum(zz^2))
    perm <- replicate(999, { zp <- sample(zz); as.numeric(sum(ok) * t(zp) %*% Wsub %*% zp / sum(Wsub) / sum(zp^2)) })
    cat(sprintf("%-15s placebo rejection at 5%%: %.3f | Moran's I of osiedle residuals %.2f (p %.3f)\n", y, mean(ps < 0.05), mI_res, mean(perm >= mI_res)))
    data.table(msr = y, rej05 = mean(ps < 0.05), rej01 = mean(ps < 0.01), p_real = M$p_real[i], p_placebo_rank = mean(ps <= M$p_real[i]), moran_res = mI_res, moran_p = mean(perm >= mI_res))
  }))
  PL[, range_km := r_best][, moran_score := I_obs]
  fwrite(PL, file.path(s2, "spatial_placebo.csv"))
  cat(sprintf("\nMean placebo rejection rate at 5%% across the 22 measures: %.3f (range %.3f-%.3f)\n", mean(PL$rej05), min(PL$rej05), max(PL$rej05)))
  cat("p_placebo_rank = share of placebo p-values at or below the real p-value (a spatially calibrated p-value)\n")
  print(PL[order(p_real), .(msr, p_real = signif(p_real, 2), p_placebo_rank = signif(p_placebo_rank, 2), rej05 = round(rej05, 3), moran_res = round(moran_res, 2), moran_p = round(moran_p, 3))], row.names = FALSE)
  cat("Saved spatial_placebo.csv\n")
})

## ── 67 · Counterfactual geographies ────────────────────────────────────────────────────────────────────────
section("67", "Counterfactual geographies", {
  ## Counterfactual geographies on the existing network: how much of the raw gaps the radial network produces
  ## given where deprivation sits, and what displacement would do to deprived residents' access.
  ## (1) Allocations of the same 42 osiedle scores: observed; "mirror city" (scores re-assigned so that
  ##     deprivation rises with the residents' distance from the centre, as where poverty is peripheral);
  ##     999 random allocations. Catchments take their residents' mean score (weights of script 62). For each,
  ##     the raw slope of the spatial measures (SD per SD) and the medians of the least and most deprived thirds.
  ## (2) Displacement scenario: a share (10%, 25%, 50%) of the residents of the most deprived inner catchments
  ##     (most deprived third, within 4 km of the centre) moves to (a) the outer large estates (least deprived
  ##     third, beyond 4 km) or (b) the peripheral deprived enclaves (most deprived third, beyond 6 km), spread
  ##     over the destination catchments in proportion to their residents; service and access are held fixed.
  ##     Reported: the resident-weighted mean share of jobs within 30 minutes, time to the centre and the
  ##     access-poverty headcount (below half the median, as in script 65) of the deprived population.
  ## Inputs: gap_stops.rds (script 56), catchment_osiedle_weights.csv (62), osiedle profiles (stage 1)
  ## Outputs: output/stage2/counterfactual_allocations.csv, counterfactual_displacement.csv

  suppressPackageStartupMessages({library(data.table); library(sf)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  X <- readRDS(file.path(s2, "gap_stops.rds")); X[, residents := round(expm1(lres))]
  Wc <- fread(file.path(s2, "catchment_osiedle_weights.csv"), encoding = "UTF-8")
  P <- fread(file.path(out_dir, "output/stage1/s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, dep, dist_km)]
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  MS <- c(ride_km = "Ride length", dest_dep = "Deprivation of destinations", conn_100 = "Connections per 100 rides", hub_share = "Alighting at interchanges",
          jobs30 = "City jobs within 30 min", t_centre = "Time to the centre", transfers_jobs = "Transfers to jobs", l_dep = "Departures (log)", dist_km = "Distance to the centre")
  Wm <- dcast(Wc[stop_name %in% X$stop_name], stop_name ~ osiedle, value.var = "w", fill = 0)
  ON <- P$osiedle; for (o in setdiff(ON, names(Wm))) Wm[, (o) := 0]
  Wmat <- as.matrix(Wm[, ON, with = FALSE]); rownames(Wmat) <- Wm$stop_name; Wmat <- Wmat[X$stop_name, ]
  alloc <- function(score) {   # catchment deprivation from an osiedle-level allocation of the scores
    dz <- as.vector(Wmat %*% score); q <- quantile(dz, c(1/3, 2/3))
    rbindlist(lapply(names(MS), function(m) { ok <- is.finite(X[[m]])
      data.table(msr = m, slope = coef(lm(zs(X[[m]][ok]) ~ zs(dz[ok])))[[2]],
                 least = median(X[[m]][ok & dz <= q[1]]), most = median(X[[m]][ok & dz >= q[2]])) }))
  }
  obs <- alloc(P$dep)
  mirror <- alloc(sort(P$dep)[rank(P$dist_km, ties.method = "first")])   # most deprived score to the most distant osiedle
  set.seed(20261009)
  rnd <- rbindlist(lapply(1:999, function(b) alloc(sample(P$dep))[, b := b]))
  CA <- merge(merge(obs[, .(msr, observed = slope, obs_least = least, obs_most = most)], mirror[, .(msr, mirror = slope, mir_least = least, mir_most = most)], by = "msr"),
              rnd[, .(random_mean = mean(slope), random_lo = quantile(slope, 0.025), random_hi = quantile(slope, 0.975), share_random_beyond_observed = mean(abs(slope) >= abs(obs$slope[match(.BY$msr, obs$msr)]))), by = msr], by = "msr")
  CA[, label := MS[msr]]
  fwrite(CA, file.path(s2, "counterfactual_allocations.csv"))
  cat(sprintf("Osiedle score vs distance: observed Spearman %.2f; mirror %.2f\n", cor(P$dep, P$dist_km, method = "spearman"), cor(sort(P$dep)[rank(P$dist_km, ties.method = "first")], P$dist_km, method = "spearman")))
  cat("Raw slopes on catchment deprivation (SD per SD) under the observed, mirror and random allocations of the same scores; medians least / most deprived third\n")
  print(CA[, .(label, observed = round(observed, 2), mirror = round(mirror, 2), random = sprintf("%.2f [%.2f, %.2f]", random_mean, random_lo, random_hi),
               obs = sprintf("%.2f / %.2f", obs_least, obs_most), mirror_thirds = sprintf("%.2f / %.2f", mir_least, mir_most))], row.names = FALSE)

  ## displacement scenario
  wmed <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]] }
  A <- X[is.finite(jobs30) & residents > 0]
  line <- 0.5 * wmed(A$jobs30, A$residents)
  src <- A[third == "most" & dist_km <= 4]; destA <- A[third == "least" & dist_km > 4]; destB <- A[third == "most" & dist_km > 6]
  deprived <- A[third == "most"]
  stat <- function(D) D[, .(residents = sum(residents), jobs30 = sum(residents * jobs30) / sum(residents), t_centre = sum(residents * t_centre) / sum(residents),
                            headcount = sum(residents * (jobs30 < line)) / sum(residents))]
  DS <- rbind(stat(deprived)[, `:=`(scenario = "observed", moved = 0)],
    rbindlist(lapply(c(0.10, 0.25, 0.50), function(s) rbindlist(lapply(list(`to the outer estates` = destA, `to the peripheral enclaves` = destB), function(dd) {
      moved <- s * sum(src$residents)
      stay <- copy(deprived)[stop_name %in% src$stop_name, residents := residents * (1 - s)]
      arrive <- copy(dd)[, residents := moved * residents / sum(residents)]
      stat(rbind(stay, arrive, fill = TRUE)) }), idcol = "scenario")[, moved := s])))
  DS[, line_pct := line]
  fwrite(DS, file.path(s2, "counterfactual_displacement.csv"))
  cat(sprintf("\nDisplacement of residents of the most deprived third within 4 km (%s residents in %d catchments); access of the deprived population (most deprived third, after moves)\n",
              format(sum(src$residents), big.mark = ","), nrow(src)))
  print(DS[, .(scenario, moved, jobs30 = round(jobs30, 1), t_centre = round(t_centre, 1), headcount = round(headcount, 3))], row.names = FALSE)
  cat(sprintf("Destinations: outer estates %d catchments (median %.1f%% jobs), peripheral enclaves %d catchments (median %.1f%% jobs)\n",
              nrow(destA), median(destA$jobs30), nrow(destB), median(destB$jobs30)))
  cat("Saved counterfactual_allocations.csv, counterfactual_displacement.csv\n")
})

## ── 68 · Synthetic test of the destination estimators ──────────────────────────────────────────────────────
section("68", "Synthetic test of the destination estimators", {
  ## Synthetic recovery test of the OD estimators (no ground-truth OD exists for Poznań).
  ## On the real valid counted trips of October 2025 (stop sequences, boardings and distances along the route),
  ## known "true" flows are generated: riders boarding at stop i choose a stop j further along the trip with
  ## probability proportional to a_j x exp(-(log d_ij - log m)^2 / (2 x 0.6^2)), a log-normal decay around a modal
  ## ride m = 1.5, 3 and 5 km, where a_j is the stop's counted alightings on that trip (+0.5) and d_ij the distance along the
  ## route (at least 0.2 km). The decay cannot be split into origin and destination factors, so no estimator is exact by
  ## construction (an exponential decay would make the distance-seeded fit exact). The true flows imply
  ## alightings at every stop; the three estimators of script 40 (Tsygalnitsky = null-seed biproportional,
  ## distance-seeded biproportional, proportional to downstream alightings) are applied to these boardings and
  ## alightings, and each origin complex's estimated mean ride is compared with its true mean ride, as is the
  ## raw slope of ride length on catchment deprivation (stops of script 56).
  ## Output: output/stage2/synthetic_od_test.csv

  suppressPackageStartupMessages({library(data.table); library(sf)})
  out_dir <- ROOT
  cfg <- list(peka = PEKA_TABLES, gtfs = GTFS_DIR,
              s2 = file.path(out_dir, "output/stage2"))
  setDTthreads(8)
  ## ── valid counted trips of October, cleaned as in script 40 ─────────────
  C <- st_read(file.path(cfg$s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  gs <- fread(file.path(cfg$gtfs, "stops.txt"), encoding = "UTF-8")[, code := sub("[0-9]+$", "", stop_code)]
  gt <- unique(gs[, .(code, stop_name)])[stop_name %in% C$stop_name][!duplicated(code)]
  cc <- gs[, .(lon = mean(stop_lon), lat = mean(stop_lat)), by = code]
  xy <- st_coordinates(st_transform(st_as_sf(cc, coords = c("lon", "lat"), crs = 4326), st_crs(C))); cc[, `:=`(X = xy[, 1], Y = xy[, 2])]
  x <- readRDS(file.path(cfg$peka, "peka_2025_10.rds")); setDT(x)
  x <- x[!is.na(DATA) & !is.na(SYM_PRZYST) & !is.na(NR_LINI) & NR_LINI != "", .(DT_KARTY, ID_KURSU, LP_PRZYST, NB_WOZU, DATA, NR_LINI, SYM_PRZYST, has_apc, boardings, alightings)]
  x <- unique(x, by = c("DT_KARTY", "ID_KURSU", "LP_PRZYST"))
  x[, tid := .GRP, by = .(DT_KARTY, ID_KURSU)]; setorder(x, tid, LP_PRZYST)
  x[, `:=`(b = fcoalesce(as.numeric(boardings), 0), a = fcoalesce(as.numeric(alightings), 0))]
  vd <- x[has_apc == TRUE, .(B = sum(b), A = sum(a)), by = .(NB_WOZU, DT_KARTY)][, vd_ok := pmax(B, A) == 0 | abs(B - A) <= 0.10 * pmax(B, A)]
  x[vd, vd_ok := i.vd_ok, on = .(NB_WOZU, DT_KARTY)]
  tr <- x[, .(all_apc = all(has_apc), B = sum(b), A = sum(a)), by = tid][, imb_ok := abs(B - A) <= pmax(5, 0.30 * pmax(B, A))]
  tr[, `:=`(fb = fifelse(B > 0, (B + A) / 2 / B, 1), fa = fifelse(A > 0, (B + A) / 2 / A, 1))]
  x[tr, `:=`(all_apc = i.all_apc, imb_ok = i.imb_ok, fb = i.fb, fa = i.fa), on = "tid"]
  x[, `:=`(bb = b * fb, ab = a * fa)][, load := cumsum(bb - ab), by = tid]
  x[x[, .(load_ok = min(load) >= -3), by = tid], load_ok := i.load_ok, on = "tid"]
  x[, counted := has_apc & all_apc & (vd_ok %in% TRUE) & imb_ok & load_ok]
  ok <- x[, .(ok = all(counted)), by = tid][ok == TRUE, tid]
  w <- x[tid %in% ok, .(tid, LP_PRZYST, SYM_PRZYST, bb, ab)]; rm(x); gc()
  w[gt, cx := i.stop_name, on = c(SYM_PRZYST = "code")][is.na(cx), cx := "outside"]
  w[cc, `:=`(X = i.X, Y = i.Y), on = c(SYM_PRZYST = "code")]
  w[, `:=`(X = nafill(nafill(X, "locf"), "nocb"), Y = nafill(nafill(Y, "locf"), "nocb")), by = tid]
  w[, km := cumsum(fcoalesce(sqrt((X - shift(X))^2 + (Y - shift(Y))^2), 0)) / 1000, by = tid]
  cat("Valid counted trips (October):", uniqueN(w$tid), "\n")

  ## ── pairs on each trip ───────────────────────────────────────────────────
  O <- w[bb > 0, .(tid, lp_o = LP_PRZYST, o = cx, bo = bb, km_o = km)]
  Dd <- w[, .(tid, lp_d = LP_PRZYST, att = ab + 0.5, km_d = km)]
  P <- Dd[O, on = .(tid, lp_d > lp_o), .(tid = i.tid, lp_o = i.lp_o, lp_d = x.lp_d, o = i.o, bo = i.bo, att = x.att, km = x.km_d - i.km_o), allow.cartesian = TRUE, nomatch = 0]
  P[, go := frank(.SD, ties.method = "dense"), .SDcols = c("tid", "lp_o")]; P[, gd := frank(.SD, ties.method = "dense"), .SDcols = c("tid", "lp_d")]
  setorder(P, tid, lp_o, lp_d)
  cat("Origin-destination pairs:", nrow(P), "\n")

  estimate <- function(P) {   # P has bo (boardings at the origin) and a_syn (alightings at the destination) on each pair
    ## proportional to downstream alightings
    S <- P[, .(S = sum(a_syn[!duplicated(lp_d)])), by = .(tid, lp_o)]   # downstream alightings of each origin (each destination once)
    P[S, S := i.S, on = .(tid, lp_o)]
    P[, f_prop := bo * a_syn / S]
    ## Tsygalnitsky recursion on the stop sequence of each trip
    st <- unique(rbind(P[, .(tid, lp = lp_o, b = bo, a = NA_real_)], P[, .(tid, lp = lp_d, b = NA_real_, a = a_syn)]), by = NULL)
    st <- st[, .(b = max(c(0, b), na.rm = TRUE), a = max(c(0, a), na.rm = TRUE)), by = .(tid, lp)][order(tid, lp)]
    st[, L := shift(cumsum(b - a), fill = 0), by = tid][, p := fifelse(a <= 0, 0, fifelse(L > a, a / L, 1))][, z := p >= 1 - 1e-9]
    st[, `:=`(Cc = cumsum(fifelse(z, 0, log1p(-p))), g = cumsum(z)), by = tid][, `:=`(Cp = shift(Cc, fill = 0), gp = shift(g, fill = 0L)), by = tid]
    P[st, `:=`(C_o = i.Cc, g_o = i.g), on = .(tid, lp_o = lp)]; P[st, `:=`(p_d = i.p, Cp_d = i.Cp, gp_d = i.gp), on = .(tid, lp_d = lp)]
    P[, f_tsy := bo * p_d * exp(Cp_d - C_o) * (gp_d == g_o)]
    ## distance-seeded biproportional fit
    f <- pmax(P$km, 0.2); for (it in 1:40) { f <- f * P$bo / as.vector(rowsum(f, P$go))[P$go]; f <- f * P$a_syn / as.vector(rowsum(f, P$gd))[P$gd] }
    P[, f_dist := f]
    P[]
  }
  X56 <- readRDS(file.path(cfg$s2, "gap_stops.rds"))[, .(stop_name, dep, osiedle)]
  zs <- function(v) (v - mean(v)) / sd(v)
  RES <- rbindlist(lapply(c(1.5, 3, 5), function(lam) {
    P[, u := att * exp(-(log(pmax(km, 0.2)) - log(lam))^2 / (2 * 0.6^2))][, t_true := bo * u / sum(u), by = go]
    A <- P[, .(a_syn = sum(t_true)), by = gd]; P[A, a_syn := i.a_syn, on = "gd"]
    ## destinations with no true flow keep zero alightings
    P[is.na(a_syn), a_syn := 0]
    estimate(P)
    OR <- P[o != "outside", .(true = sum(t_true * km) / sum(t_true), tsy = sum(f_tsy * km) / sum(f_tsy), dist = sum(f_dist * km) / sum(f_dist), prop = sum(f_prop * km) / sum(f_prop),
                              boardings = sum(t_true)), by = .(stop_name = o)][boardings >= 400]
    Z <- merge(OR, X56, by = "stop_name")
    sl <- function(v) coef(lm(zs(Z[[v]]) ~ zs(Z$dep)))[[2]]
    sl_km <- function(v) coef(lm(Z[[v]] ~ zs(Z$dep)))[[2]]
    rbindlist(lapply(c("tsy", "dist", "prop"), function(e) data.table(mode_km = lam, estimator = e, origins = nrow(Z),
      mean_ride_true = P[, sum(t_true * km) / sum(t_true)], mean_ride_est = P[, sum(get(paste0("f_", e)) * km) / sum(get(paste0("f_", e)))],
      origin_bias_km = mean(Z[[e]] - Z$true), origin_mae_km = mean(abs(Z[[e]] - Z$true)), origin_cor = cor(Z[[e]], Z$true),
      slope_true_sd = sl("true"), slope_est_sd = sl(e), slope_true_km = sl_km("true"), slope_est_km = sl_km(e))))
  }))
  fwrite(RES, file.path(cfg$s2, "synthetic_od_test.csv"))
  cat("Recovery of known flows on October's counted trips: mean ride (km), origin-level error and correlation, and the raw slope of ride length on deprivation\n")
  print(RES[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)], row.names = FALSE)
  cat("Saved synthetic_od_test.csv\n")
})

## ── 73 · Destinations offered and reached ──────────────────────────────────────────────────────────────────
section("73", "Destinations offered and reached", {
  ## Is the remaining gap in the deprivation of destinations made by the network or by riders' use of it?
  ## The ladder of script 56 (raw, net of position, net of land use and population; CR2-Satterthwaite by
  ## osiedle) is applied to the deprivation of the destinations (a) reached by riders, (b) offered by the
  ## network (riders assigned to every line and direction at their stop in proportion to its departures) and
  ## (c) their difference (the part made by the lines and directions riders use); weekday, all day.
  ## Inputs: gap_stops.rds (56), dest_mech_profiles.csv (46). Output: output/stage2/destination_network_check.csv

  suppressPackageStartupMessages({library(data.table); library(splines)})
  out_dir <- ROOT
  source(file.path(REPO, "R/cr_inference.R"))
  s2 <- file.path(out_dir, "output/stage2")
  X <- readRDS(file.path(s2, "gap_stops.rds"))
  DP <- dcast(fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & window == "weekday all day" & comp %in% c("observed", "network")],
              stop_name ~ comp, value.var = "dest_dep")
  Z <- merge(X[, !c("dest_dep"), with = FALSE], DP, by = "stop_name")
  Z[, use := observed - network]
  zs <- function(v) (v - mean(v)) / sd(v)
  LOC <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2"); SERV <- c("l_lines", "l_dep", "tram_share", "direct_centre")
  LUSE <- c("b1", "b2", "b3", "b4", "landuse_mix", "fsi", "wc_tree"); POPV <- c("share_0_18", "share_61plus", "lres", "larea")
  for (v in c(SERV, LUSE, POPV)) Z[, (v) := zs(get(v))]
  Z[, dz := zs(dep)]
  s_obs <- sd(Z$observed)   # all three on the scale of the SD of the reached destinations, so (b) + (c) = (a)
  OUT <- rbindlist(lapply(c("observed", "network", "use"), function(y) { Z[, yy := get(y) / s_obs]
    rbindlist(lapply(list(raw = character(0), position = c(LOC, SERV), full = c(LOC, SERV, LUSE, POPV)), function(ctrl) {
      f <- lm(reformulate(c("dz", ctrl), "yy"), Z); ci <- crinf(f, Z$osiedle, "dz")
      data.table(component = y, est = ci$b, lo = ci$b - qt(0.975, ci$df) * ci$se2, hi = ci$b + qt(0.975, ci$df) * ci$se2, p = 2 * pt(-abs(ci$b / ci$se2), ci$df)) }), idcol = "step") }))
  fwrite(OUT, file.path(s2, "destination_network_check.csv"))
  cat(sprintf("Origins: %d. Slopes of destination deprivation on origin deprivation (in SD of the reached destinations' deprivation):\n", nrow(Z)))
  print(dcast(OUT, component ~ step, value.var = c("est", "p"))[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
  cat("Saved destination_network_check.csv\n")
})


## ── 74 · Offer and use: what the network offers and how riders use it ───────────────────────────────────────
section("74", "Offer and use: what the network offers and how riders use it", {
  ## Position on the network, as distance and sector, is a proxy for what the network offers at a stop. This
  ## section measures the offer directly and splits every origin measure into what the network offers and how
  ## riders use it. The offer is the network baseline of sections 46 and 51: each origin's boardings spread
  ## over every line and direction at the stop in proportion to its departures, carried by the same
  ## destination estimate. Use is the observed measure minus the offer; it reflects the lines and directions
  ## riders board, which the counts observe directly.
  ## (a) Split of the gap: slopes of observed, offered and use on origin deprivation, all in SD of the observed
  ##     measure, so offered + use = observed; raw, net of position (location and service) and net of land
  ##     use and population too.
  ## (b) Use given the offer: the observed measure on deprivation with the offer held fixed (regression
  ##     adjustment), alone, with position and with position, land use and population.
  ## (c) How much of the variation between origins the offer and location capture (R2, deprivation left out).
  ## CR2 standard errors by osiedle with Satterthwaite df; Benjamini-Hochberg q across the nine measures within
  ## each window, part and step. Windows: Monday-Thursday as a whole (main) and its morning peak, when most
  ## boardings start from home (check).
  ## Inputs: gap_stops.rds (56), dest_mech_profiles.csv (46), routes_origin.csv (51).
  ## Outputs: output/stage2/offer_use.csv, offer_use_r2.csv; output/figures/Fig_offer_use.png

  suppressPackageStartupMessages({library(data.table); library(splines); library(ggplot2)})
  out_dir <- ROOT
  source(file.path(REPO, "R/cr_inference.R"))
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  X <- readRDS(file.path(s2, "gap_stops.rds"))
  K <- c(ride_km = "Ride length", dest_dep = "Deprivation of destinations", to_centre = "Share carried to the centre",
         same_osiedle = "Share within the origin's osiedle", hub = "Share alighting at interchanges", office = "Office floor at destinations",
         eff_dest = "Effective destinations", dir_conc = "Directional concentration", spread_km = "Spread of destinations")
  DM <- fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & comp %in% c("observed", "network")]
  DM[, win := fcase(window == "weekday all day", "weekday", window == "weekday morning peak (06-09)", "morning peak")]
  RO <- fread(file.path(s2, "routes_origin.csv"), encoding = "UTF-8")[comp %in% c("observed", "network")]
  RO[, win := fcase(window == "Mon-Thu", "weekday", window == "Mon-Thu morning peak", "morning peak")]
  setnames(RO, "o", "stop_name")
  long <- rbind(melt(DM[!is.na(win)], id.vars = c("stop_name", "comp", "win"), measure.vars = intersect(names(K), names(DM)), variable.name = "msr"),
                melt(RO[!is.na(win)], id.vars = c("stop_name", "comp", "win"), measure.vars = intersect(names(K), names(RO)), variable.name = "msr"))
  W <- dcast(long, stop_name + win + msr ~ comp, value.var = "value")[!is.na(observed) & !is.na(network)]
  zs <- function(v) (v - mean(v)) / sd(v)
  LOC <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2"); SERV <- c("l_lines", "l_dep", "tram_share", "direct_centre")
  LUSE <- c("b1", "b2", "b3", "b4", "landuse_mix", "fsi", "wc_tree"); POPV <- c("share_0_18", "share_61plus", "lres", "larea")
  CV <- c("stop_name", "dep", "dist_km", "osiedle", "sec_s1", "sec_c1", "sec_s2", "sec_c2", SERV, LUSE, POPV)
  STEPS <- list(raw = character(0), position = c(LOC, SERV), full = c(LOC, SERV, LUSE, POPV))
  slope <- function(Z, y, ctrl) {
    f <- lm(reformulate(c("dz", ctrl), y), Z); ci <- crinf(f, Z$osiedle, "dz")
    data.table(est = ci$b, se = ci$se2, df = ci$df, lo = ci$b - qt(0.975, ci$df) * ci$se2, hi = ci$b + qt(0.975, ci$df) * ci$se2,
               p = 2 * pt(-abs(ci$b / ci$se2), ci$df), n = nobs(f))
  }
  OUT <- list(); R2 <- list()
  for (w in c("weekday", "morning peak")) for (m in names(K)) {
    Z <- merge(X[, ..CV], W[win == w & msr == m, .(stop_name, observed, network)], by = "stop_name")
    if (nrow(Z) < 50) next
    for (v in c(SERV, LUSE, POPV)) Z[, (v) := zs(get(v))]
    s_obs <- sd(Z$observed)
    Z[, `:=`(dz = zs(dep), y_obs = observed / s_obs, y_net = network / s_obs, y_use = (observed - network) / s_obs, offer = zs(network))]
    for (st in names(STEPS)) {
      for (pt_ in c("observed", "offered", "use"))
        OUT[[length(OUT) + 1]] <- cbind(data.table(win = w, msr = m, analysis = "split", part = pt_, step = st),
                                        slope(Z, c(observed = "y_obs", offered = "y_net", use = "y_use")[[pt_]], STEPS[[st]]))
      OUT[[length(OUT) + 1]] <- cbind(data.table(win = w, msr = m, analysis = "given offer", part = "use", step = st),
                                      slope(Z, "y_obs", c("offer", STEPS[[st]])))
    }
    r2 <- function(rhs) summary(lm(reformulate(rhs, "y_obs"), Z))$r.squared
    R2[[length(R2) + 1]] <- data.table(win = w, msr = m, n = nrow(Z), offer = r2("offer"), location = r2(LOC), position = r2(c(LOC, SERV)),
                                       offer_location = r2(c("offer", LOC)))
  }
  OUT <- rbindlist(OUT); R2 <- rbindlist(R2)
  OUT[, q := p.adjust(p, "BH"), by = .(win, analysis, part, step)]
  OUT[, label := K[as.character(msr)]]
  fwrite(OUT, file.path(s2, "offer_use.csv")); fwrite(R2, file.path(s2, "offer_use_r2.csv"))

  f2 <- function(x) sprintf("%.2f", x)
  for (w in c("weekday", "morning peak")) {
    cat(sprintf("\n== %s: slopes on origin deprivation (SD of the observed measure per SD)\n", w))
    T1 <- dcast(OUT[win == w & ((analysis == "split" & step == "raw") | analysis == "given offer")],
                label ~ paste(analysis, part, step), value.var = "est")
    print(T1[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
    cat("use given the offer, q (raw / position / full):\n")
    print(dcast(OUT[win == w & analysis == "given offer"], label ~ step, value.var = "q")[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 2) else x)])
    cat("R2 of the offer and of location:\n"); print(R2[win == w, .(msr, n, offer = round(offer, 2), location = round(location, 2), position = round(position, 2), offer_location = round(offer_location, 2))])
  }

  ## figure: (a) the raw gap split into offer and use; (b) use given the offer, three steps
  P <- OUT[win == "weekday"]
  P[, label := factor(label, levels = rev(K))]
  pa <- P[analysis == "split" & step == "raw"][, part := factor(part, levels = c("observed", "offered", "use"),
                                                                  labels = c("Observed gap", "Offered by the network", "Use of the network"))]
  pb <- P[analysis == "given offer"][, step := factor(step, levels = c("raw", "position", "full"),
                                                        labels = c("Offer held fixed", "+ position", "+ land use and population"))]
  th <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = "bottom",
                                              plot.title = element_text(face = "bold", size = 10), legend.title = element_blank())
  dodge <- position_dodge(width = 0.65)
  g1 <- ggplot(pa, aes(est, label, colour = part)) + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_linerange(aes(xmin = lo, xmax = hi), position = dodge, linewidth = 0.6) + geom_point(position = dodge, size = 2) +
    scale_colour_manual(values = c("#2b2b28", "#8f8d86", "#eb6834")) + labs(title = "(a) The gap: offered + use = observed", x = "Slope on origin deprivation (SD per SD)", y = NULL) + th
  g2 <- ggplot(pb, aes(est, label, colour = step)) + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_linerange(aes(xmin = lo, xmax = hi), position = dodge, linewidth = 0.6) + geom_point(position = dodge, size = 2) +
    scale_colour_manual(values = c("#eb6834", "#2a78d6", "#6b6a64")) + labs(title = "(b) Use with the offer held fixed", x = "Slope on origin deprivation (SD per SD)", y = NULL) + th +
    theme(axis.text.y = element_blank())
  if (requireNamespace("patchwork", quietly = TRUE)) {
    g <- patchwork::wrap_plots(g1, g2, widths = c(1, 1))
    ggsave(file.path(fig, "Fig_offer_use.png"), g, width = 11, height = 5.2, dpi = 200, bg = "white")
  }
  cat("\nSaved offer_use.csv, offer_use_r2.csv, Fig_offer_use.png\n")
})

## ── 75 · Peripheral deprivation: what the network offers deprived areas far from the centre ─────────────────
section("75", "Peripheral deprivation: what the network offers deprived areas far from the centre", {
  ## The advantage of deprived areas is central. Where deprivation lies far from the centre, does the network
  ## offer it less than it offers other areas equally far out, as the spatial-mismatch hypothesis expects?
  ## Peripheral: complexes at least 5 km from the Old Market Square (4 and 6 km as checks). Deprived: the most
  ## deprived third of section 56. Comparisons at equal distance, and at equal distance south of the centre
  ## (where all peripheral deprived osiedla lie; sine of the bearing below zero). Offer: jobs within 30 min,
  ## time to the centre and transfers to jobs by the timetable (section 54), and the ride length, effective
  ## destinations and destination deprivation offered by the network baseline (sections 46, 51). Use: the
  ## observed ride length and effective destinations, and connections per 100 rides (section 47).
  ## Few osiedla hold peripheral deprivation, so cluster-robust slopes are unreliable. Inference is exact
  ## randomisation at the osiedle level: each osiedle with peripheral complexes is one unit (median over its
  ## peripheral complexes); the statistic is the mean of the deprived osiedla minus the mean of the others,
  ## and its p-value the share of the ways of choosing as many osiedla as deprived (all of them up to 100,000, otherwise 100,000 random draws) that give a difference at
  ## least as large in absolute value. An osiedle counts as deprived when most of its peripheral complexes
  ## are in the most deprived third. Check: without Krzesiny-Pokrzywno-Garaszewo (single-institution score).
  ## Inputs: gap_stops.rds (56), dest_mech_profiles.csv (46), routes_origin.csv (51).
  ## Outputs: output/stage2/periphery_stops.csv, periphery_osiedla.csv, periphery_tests.csv;
  ##          output/figures/Fig_periphery.png

  suppressPackageStartupMessages({library(data.table); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  X <- readRDS(file.path(s2, "gap_stops.rds"))
  DM <- dcast(fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & window == "weekday all day" & comp == "network"],
              stop_name ~ comp, value.var = c("ride_km", "dest_dep"))
  RO <- fread(file.path(s2, "routes_origin.csv"), encoding = "UTF-8")[window == "Mon-Thu" & comp == "network", .(stop_name = o, eff_dest_network = eff_dest)]
  X <- merge(merge(X, DM, by = "stop_name", all.x = TRUE), RO, by = "stop_name", all.x = TRUE)
  X[, south := sec_s1 < 0]
  M <- c(jobs30 = "Jobs within 30 min (%)", t_centre = "Time to the centre (min)", transfers_jobs = "Transfers to jobs",
         ride_km_network = "Ride length offered (km)", eff_dest_network = "Effective destinations offered", dest_dep_network = "Deprivation of destinations offered",
         ride_km = "Ride length (km)", eff_dest = "Effective destinations", conn_100 = "Connections per 100 rides")
  KIND <- c(rep("offer", 6), rep("use", 3)); names(KIND) <- names(M)
  ex <- "Krzesiny-Pokrzywno-Garaszewo"

  ## stop-level medians: deprived and other complexes beyond the threshold, all and south of the centre
  ST <- rbindlist(lapply(c(4, 5, 6), function(km) rbindlist(lapply(c("all", "south"), function(ar) {
    P <- X[dist_km >= km & (ar == "all" | south)]
    P[, grp := fifelse(third == "most", "deprived", "other")]
    rbindlist(lapply(names(M), function(m) P[!is.na(get(m)), .(km = km, area = ar, msr = m, kind = KIND[[m]], stops = .N, osiedla = uniqueN(osiedle),
                                                            median = median(get(m))), by = grp]))
  }))))
  fwrite(ST, file.path(s2, "periphery_stops.csv"))

  ## osiedle level: median over each osiedle's peripheral complexes
  osi_level <- function(km, ar, drop = character(0)) {
    P <- X[dist_km >= km & (ar == "all" | south) & !osiedle %in% drop]
    O <- P[, c(list(complexes = .N, deprived = mean(third == "most") > 0.5, km_median = median(dist_km)),
               lapply(.SD, function(v) median(v, na.rm = TRUE))), by = osiedle, .SDcols = names(M)]
    O[]
  }
  ## randomisation test of mean(deprived) - mean(other) over the choices of the deprived osiedla:
  ## all of them up to 100,000, otherwise 100,000 random draws
  ri <- function(v, d) {
    ok <- !is.na(v); v <- v[ok]; d <- d[ok]; k <- sum(d); n <- length(v)
    if (k < 1 || k >= n) return(list(diff = NA_real_, p = NA_real_, k = k, n = n, perms = 0))
    obs <- mean(v[d]) - mean(v[!d])
    cmb <- if (choose(n, k) <= 1e5) combn(n, k) else { set.seed(75); replicate(1e5, sample.int(n, k)) }
    st <- colSums(matrix(v[cmb], nrow = k)); s <- sum(v)
    dd <- st / k - (s - st) / (n - k)
    list(diff = obs, p = mean(abs(dd) >= abs(obs) - 1e-12), k = k, n = n, perms = ncol(cmb))
  }
  TS <- list(); OS <- list()
  for (km in c(4, 5, 6)) for (ar in c("all", "south")) for (chk in c("main", "without Krzesiny")) {
    O <- osi_level(km, ar, if (chk == "main") character(0) else ex)
    if (km == 5 && chk == "main") OS[[length(OS) + 1]] <- cbind(data.table(km = km, area = ar), O)
    for (m in names(M)) { r <- ri(O[[m]], O$deprived)
      TS[[length(TS) + 1]] <- data.table(km = km, area = ar, check = chk, msr = m, kind = KIND[[m]], deprived_osiedla = r$k, osiedla = r$n,
                                         deprived = mean(O[[m]][O$deprived], na.rm = TRUE), other = mean(O[[m]][!O$deprived], na.rm = TRUE),
                                         diff = r$diff, p_exact = r$p, permutations = r$perms) }
  }
  TS <- rbindlist(TS); OS <- rbindlist(OS)
  fwrite(TS, file.path(s2, "periphery_tests.csv")); fwrite(OS, file.path(s2, "periphery_osiedla.csv"))

  cat("Deprived osiedla beyond 5 km (osiedle level):\n")
  print(OS[deprived == TRUE, .(area, osiedle, complexes, km = round(km_median, 1), jobs30 = round(jobs30, 1), t_centre = round(t_centre), conn_100 = round(conn_100, 2))])
  cat("\nStop-level medians, 5 km:\n")
  print(dcast(ST[km == 5], msr + kind ~ area + grp, value.var = "median")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
  cat("\nOsiedle-level exact randomisation tests, 5 km:\n")
  print(TS[km == 5, .(area, check, msr, k = deprived_osiedla, n = osiedla, deprived = round(deprived, 2), other = round(other, 2), diff = round(diff, 2), p = round(p_exact, 3))])
  cat("\nThreshold checks (south of the centre, main): jobs30, t_centre, transfers_jobs\n")
  print(TS[area == "south" & check == "main" & msr %in% c("jobs30", "t_centre", "transfers_jobs"), .(km, msr, k = deprived_osiedla, n = osiedla, diff = round(diff, 2), p = round(p_exact, 3))])

  ## figure: osiedle-level offer against distance, peripheral osiedla, deprived highlighted
  F <- OS[area == "all"]
  FL <- melt(F, id.vars = c("osiedle", "deprived", "km_median", "complexes"), measure.vars = c("jobs30", "t_centre", "eff_dest_network", "conn_100"))
  FL[, variable := factor(M[as.character(variable)], levels = M[c("jobs30", "t_centre", "eff_dest_network", "conn_100")])]
  g <- ggplot(FL, aes(km_median, value)) + geom_point(aes(colour = deprived, size = complexes), alpha = 0.85) +
    ggrepel::geom_text_repel(data = FL[deprived == TRUE], aes(label = osiedle), size = 2.6, colour = "#e34948", min.segment.length = 0, seed = 1) +
    facet_wrap(~variable, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = c(`TRUE` = "#e34948", `FALSE` = "#8f8d86"), labels = c(`TRUE` = "Most deprived third", `FALSE` = "Other"), name = NULL) +
    scale_size_area(max_size = 4, guide = "none") + labs(x = "Distance from the Old Market Square (km, median of the osiedle's complexes)", y = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(file.path(fig, "Fig_periphery.png"), g, width = 12, height = 4.2, dpi = 200, bg = "white")
  cat("\nSaved periphery_stops.csv, periphery_osiedla.csv, periphery_tests.csv, Fig_periphery.png\n")
})

## ── 76 · Early boardings: the early offer, its use, and what "sector" stands for ──────────────────────────────
section("76", "Early boardings: the early offer, its use, and what sector stands for", {
  ## Early boardings follow the sector of the city rather than deprivation (section 56). The onset index of
  ## section 42 compares a stop's share of boardings before 06:00 with that of other stops on the same lines,
  ## but early trips need not serve a whole line: pull-outs from depots and short early workings serve only
  ## part of it. The index is therefore split, with the same groups and weights, into
  ##   offer  the stop's share of its departures that run before 06:00 (logit), against the other stops of
  ##          the line and direction: which stops early trips serve;
  ##   use    early boardings per early departure against all-day boardings per departure: how much the early
  ##          trips are used where they run (stops with early departures only).
  ## Then: how much of the onset index each captures, whether sector still adds to the offer, the slopes of
  ## the three on deprivation (CR2-Satterthwaite by osiedle, raw / net of location / net of location, land use
  ## and population, as for the usage indices in section 56), means by octant, and the land use of the
  ## destinations of early riders (04-06) against midday riders, from section 46.
  ## Inputs: stage2_line_stop_cells.rds (41), gap_stops.rds (56), dest_mech_profiles.csv (46).
  ## Outputs: output/stage2/early_offer_use.csv, early_offer_use_tests.csv, early_sector.csv;
  ##          output/figures/Fig_early_offer.png

  suppressPackageStartupMessages({library(data.table); library(splines); library(ggplot2)})
  out_dir <- ROOT
  source(file.path(REPO, "R/cr_inference.R"))
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))
  cells[, dep := departures - uncovered]
  Dm <- cells[season == "term" & daytype == "monthu" & night == FALSE & dep > 0]
  rm(cells); invisible(gc())
  index <- function(dt, gcols) {
    dt <- dt[is.finite(q) & w > 0]
    dt[, n_in := .N, by = gcols]; dt <- dt[n_in >= 2]
    dt[, dev := q - weighted.mean(q, w), by = gcols]
    dt[, .(index = weighted.mean(dev, w), groups = .N, weight = sum(w)), by = stop_name]
  }
  lg <- function(k, n) qlogis((k + 0.5) / (n + 1))
  early <- Dm[sched_hour < 6, .(dep_early = sum(dep)), by = .(NR_LINI, dir)]
  on <- merge(Dm[, .(Be = sum(B[sched_hour < 6]), Ball = sum(B), De = sum(dep[sched_hour < 6]), Dall = sum(dep)), by = .(NR_LINI, dir, stop_name)],
              early[dep_early > 0], by = c("NR_LINI", "dir"))
  I <- Reduce(function(a, b) merge(a, b, by = "stop_name", all = TRUE), list(
    setnames(index(copy(on)[, `:=`(q = lg(Be, Ball), w = Ball)], c("NR_LINI", "dir"))[, .(stop_name, index)], "index", "onset"),
    setnames(index(copy(on)[, `:=`(q = lg(De, Dall), w = Ball)], c("NR_LINI", "dir"))[, .(stop_name, index)], "index", "offer"),
    setnames(index(copy(on)[De > 0][, `:=`(q = log((Be + 0.5) / De) - log((Ball + 0.5) / Dall), w = Ball)], c("NR_LINI", "dir"))[, .(stop_name, index)], "index", "use")))
  G <- readRDS(file.path(s2, "gap_stops.rds"))
  X <- merge(G, I, by = "stop_name")
  cat(sprintf("Stops: %d. Recomputed onset against section 42: r = %.3f\n", nrow(X), cor(X$onset, X$onset_b, use = "complete.obs")))
  fwrite(X[, .(stop_name, osiedle, dep, dist_km, onset_b, onset, offer, use)], file.path(s2, "early_offer_use.csv"))

  ## how much of the onset index the offer captures, and whether sector still adds to it
  Z <- X[!is.na(onset_b) & !is.na(offer)]
  r2 <- function(rhs) summary(lm(reformulate(rhs, "onset_b"), Z))$r.squared
  LOC <- "ns(dist_km, 3)"; SEC <- c("sec_s1", "sec_c1", "sec_s2", "sec_c2")
  R <- data.table(model = c("distance", "distance + sector", "offer", "distance + offer", "distance + offer + sector"),
                  r2 = c(r2(LOC), r2(c(LOC, SEC)), r2("offer"), r2(c(LOC, "offer")), r2(c(LOC, "offer", SEC))))
  cat(sprintf("\nCorrelation of the onset index with the offer: %.2f, with use: %.2f\n", cor(Z$onset_b, Z$offer), cor(Z$onset_b, Z$use, use = "complete.obs")))
  cat("R2 of the onset index:\n"); print(R[, r2 := round(r2, 3)])

  ## slopes on deprivation, as for the usage indices of section 56 (position = location alone)
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  LUSE <- c("b1", "b2", "b3", "b4", "landuse_mix", "fsi", "wc_tree"); POPV <- c("share_0_18", "share_61plus", "lres", "larea")
  STEPS <- list(raw = character(0), location = c(LOC, SEC), full = c(LOC, SEC, LUSE, POPV))
  TS <- rbindlist(lapply(c("onset_b", "offer", "use"), function(y) {
    W <- X[!is.na(get(y))]; W[, `:=`(yy = zs(get(y)), dz = zs(dep))]
    rbindlist(lapply(names(STEPS), function(st) { f <- lm(reformulate(c("dz", STEPS[[st]]), "yy"), W); ci <- crinf(f, W$osiedle, "dz")
      data.table(measure = y, step = st, est = ci$b, lo = ci$b - qt(0.975, ci$df) * ci$se2, hi = ci$b + qt(0.975, ci$df) * ci$se2,
                 p = 2 * pt(-abs(ci$b / ci$se2), ci$df), n = nobs(f)) }))
  }))
  fwrite(TS, file.path(s2, "early_offer_use_tests.csv"))
  cat("\nSlopes on deprivation (SD per SD; CR2-Satterthwaite):\n")
  print(dcast(TS, measure ~ step, value.var = "est")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
  print(dcast(TS, measure ~ step, value.var = "p")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

  ## octants, and where early riders are carried
  X[, ang := (atan2(sec_s1, sec_c1) * 180 / pi) %% 360]
  X[, octant := factor(c("E", "NE", "N", "NW", "W", "SW", "S", "SE")[floor(((ang + 22.5) %% 360) / 45) + 1], levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))]
  DM <- fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & comp == "observed" &
         window %in% c("weekday early (04-06)", "weekday midday (09-15)")]
  DM[, win := fifelse(grepl("early", window), "early", "midday")]
  DW <- dcast(DM, stop_name ~ win, value.var = c("industry", "z_production", "jobs_res"))
  X <- merge(X, DW, by = "stop_name", all.x = TRUE)
  O <- X[, .(stops = .N, deprivation = mean(dep), km = median(dist_km), onset = mean(onset_b, na.rm = TRUE), offer = mean(offer, na.rm = TRUE),
             use = mean(use, na.rm = TRUE), industry_early = mean(industry_early, na.rm = TRUE), industry_midday = mean(industry_midday, na.rm = TRUE),
             production_early = mean(z_production_early, na.rm = TRUE), jobs_per_resident_early = mean(jobs_res_early, na.rm = TRUE)), by = octant][order(octant)]
  fwrite(O, file.path(s2, "early_sector.csv"))
  cat("\nBy octant (direction from the Old Market Square):\n"); print(O[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
  cat(sprintf("\nDestinations of early (04-06) and midday riders, mean over origins: industrial floor %.3f vs %.3f; production zoning %.3f vs %.3f; jobs per resident %.2f vs %.2f\n",
              mean(X$industry_early, na.rm = TRUE), mean(X$industry_midday, na.rm = TRUE), mean(X$z_production_early, na.rm = TRUE),
              mean(X$z_production_midday, na.rm = TRUE), mean(X$jobs_res_early, na.rm = TRUE), mean(X$jobs_res_midday, na.rm = TRUE)))
  Y <- X[!is.na(onset_b) & !is.na(offer) & !is.na(industry_early) & !is.na(z_production_early) & !is.na(jobs_res_early)]
  r2y <- function(rhs) summary(lm(reformulate(rhs, "onset_b"), Y))$r.squared
  cat(sprintf("Same %d origins, R2 of onset: distance %.3f; + sector %.3f; + early destinations' land use %.3f; + offer %.3f; + offer + sector %.3f\n", nrow(Y),
              r2y(LOC), r2y(c(LOC, SEC)), r2y(c(LOC, "industry_early", "z_production_early", "jobs_res_early")), r2y(c(LOC, "offer")), r2y(c(LOC, "offer", SEC))))

  ## figure: onset, offer and use by octant
  P <- melt(O, id.vars = "octant", measure.vars = c("onset", "offer", "use"))
  P[, variable := factor(variable, levels = c("onset", "offer", "use"), labels = c("Early boardings (onset index)", "Early offer: share of departures before 06:00", "Use: early boardings per early departure"))]
  g <- ggplot(P, aes(octant, value, fill = variable)) + geom_hline(yintercept = 0, colour = "#9a9893") + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = c("#2b2b28", "#8f8d86", "#eb6834"), name = NULL) +
    labs(x = "Direction from the Old Market Square", y = "Mean same-trip index (deviation from stops on the same lines)") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  ggsave(file.path(fig, "Fig_early_offer.png"), g, width = 8.5, height = 4.5, dpi = 200, bg = "white")
  cat("\nSaved early_offer_use.csv, early_offer_use_tests.csv, early_sector.csv, Fig_early_offer.png\n")
})

## ── 77 · What location stands for: jobs nearby and the structure of service ─────────────────────────────────
section("77", "What location stands for: jobs nearby and the structure of service", {
  ## Location (distance spline and sector) removes most gaps in section 56. Here it is replaced by named
  ## mechanisms: jobs within 1 and 2 km of the stop complex (modelled 200 m employment grid, cells whose centroid
  ## lies within the radius; log), the structure of service at the stop (log lines, log departures, tram share,
  ## share of departures on lines through the centre) and, for origins, the effective number of line-directions
  ## offered (section 51). For each measure: the slope on deprivation raw, net of location, net of the
  ## mechanisms and net of both (CR2-Satterthwaite by osiedle), and the share of variation each set captures.
  ## Jobs nearby and distance are closely correlated, so this gives the location part a meaning; it does not
  ## separate distance from the centre and proximity to jobs.
  ## Inputs: gap_stops.rds (56), stop_complexes.gpkg (21), data/workplace_grid.gpkg, routes_origin.csv (51).
  ## Output: output/stage2/location_mechanisms.csv

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines)})
  out_dir <- ROOT
  source(file.path(REPO, "R/cr_inference.R"))
  s2 <- file.path(out_dir, "output/stage2")
  G <- readRDS(file.path(s2, "gap_stops.rds"))
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  W <- st_transform(st_read(file.path(out_dir, "data/workplace_grid.gpkg"), quiet = TRUE), st_crs(C))
  wc <- st_coordinates(st_centroid(st_geometry(W))); cc <- st_coordinates(st_geometry(C))
  d <- sqrt(outer(cc[, 1], wc[, 1], "-")^2 + outer(cc[, 2], wc[, 2], "-")^2)
  J <- data.table(stop_name = C$stop_name, jobs1 = as.vector((d <= 1000) %*% W$workplaces), jobs2 = as.vector((d <= 2000) %*% W$workplaces))
  R <- fread(file.path(s2, "routes_origin.csv"), encoding = "UTF-8")[window == "Mon-Thu" & comp == "network", .(stop_name = o, eff_lines_net = eff_lines)]
  X <- merge(merge(G, J, by = "stop_name"), R, by = "stop_name", all.x = TRUE)
  X[, `:=`(lj1 = log1p(jobs1), lj2 = log1p(jobs2))]
  cat(sprintf("Jobs within 2 km: median %.0f (most deprived third) against %.0f (least); correlation with distance %.2f, with deprivation %.2f\n",
              X[third == "most", median(jobs2)], X[third == "least", median(jobs2)], cor(X$lj2, X$dist_km), cor(X$lj2, X$dep)))
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  LOC <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2")
  MECH <- c("lj1", "lj2", "l_lines", "l_dep", "tram_share", "direct_centre")
  ORIG <- c("ride_km", "dest_dep", "to_centre", "same_osiedle", "hub_share", "conn_100", "eff_dest", "dir_conc")
  ms <- c("volume", "tidal", "onset_b", "conc_b", "sat_b", "summer", ORIG, "jobs30", "transfers_jobs", "t_centre",
          "bridge", "irregularity", "dependence", "weather_rain", "night_vs_day")
  ms <- intersect(ms, names(X))
  OUT <- rbindlist(lapply(ms, function(m) {
    Z <- X[!is.na(get(m))]; mech <- MECH
    if (m %in% ORIG) { Z <- Z[!is.na(eff_lines_net)]; mech <- c(MECH, "eff_lines_net") }
    Z[, `:=`(y = zs(get(m)), dz = zs(dep))]
    sl <- function(ctrl) { f <- lm(reformulate(c("dz", ctrl), "y"), Z); ci <- crinf(f, Z$osiedle, "dz")
      c(b = ci$b, p = 2 * pt(-abs(ci$b / ci$se2), ci$df)) }
    r2 <- function(ctrl) summary(lm(reformulate(ctrl, "y"), Z))$r.squared
    a <- sl(character(0)); b <- sl(LOC); e <- sl(mech); f <- sl(c(LOC, mech))
    data.table(measure = m, n = nrow(Z), raw = a[["b"]], p_raw = a[["p"]], net_location = b[["b"]], p_location = b[["p"]],
               net_mechanisms = e[["b"]], p_mechanisms = e[["p"]], net_both = f[["b"]], p_both = f[["p"]],
               r2_location = r2(LOC), r2_mechanisms = r2(mech), r2_both = r2(c(LOC, mech)))
  }))
  fwrite(OUT, file.path(s2, "location_mechanisms.csv"))
  print(OUT[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
  cat("\nSaved location_mechanisms.csv\n")
})
