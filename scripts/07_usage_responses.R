## DualPulse | stage 7: Usage models, early boardings, responses and modes
##
## The tests stated before modelling, the sharpened early-boarding estimate, responses to leave days,
## service cuts and weather, and tram against bus (supplementary S3, S5, S6, S9).
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   43      Tests stated before modelling
##   48      Early boardings, sharpened
##   52      Responses to leave days, service cuts and weather
##   61      Tram and bus
##
## Run from the repository root: Rscript run.R 7   (one section: Rscript run.R 7:61)

if (!exists("section")) source("R/setup.R")

## ── 43 · Tests stated before modelling ─────────────────────────────────────────────────────────────────────
section("43", "Tests stated before modelling", {
  ## Stage 2 models (manuscript §2.7), January-October 2025, walking catchments.
  ##
  ## Outcomes: the same-trip stop indices of script 42 (each stop relative to the
  ## other stops served by the same trips, so the timetable is held fixed).
  ## Exposures and covariates: land-use descriptors measured in each walking
  ## catchment (script 15 run on exclusive_400m_network), catchment deprivation,
  ## age shares, residents and paid-parking share (script 34), municipal-housing
  ## share (script 37). Skewed descriptors are transformed (log floor-area density,
  ## log1p gap and setback, log block area, square-root intersection density and
  ## university/dormitory shares); all variables are standardised, so coefficients
  ## are SD of the index per SD of the characteristic.
  ## Common covariates: natural spline of distance to the centre (3 df), shares of
  ## residents aged 0-18 and 61+, log residents, log catchment area (complexes
  ## without residents: age shares set to the mean, with an indicator).
  ## (1) Activity path: index ~ within-osiedle deviations of every land-use term +
  ##     between-osiedle means of the residential vs non-residential balance,
  ##     land-use mix and tree cover + covariates (Mundlak); Conley standard errors
  ##     (Bartlett kernel, 1 km); osiedle fixed effects as a check. All complexes.
  ## (2) Deprivation path: index ~ catchment deprivation + residential balance,
  ##     mix and tree cover + the other five floor-area balances + principal
  ##     components of the remaining core descriptors (intensity and configuration,
  ##     >= 80% of their variance) + university and dormitory shares + paid-parking
  ##     share + covariates. Complexes with residents. CR2 standard errors clustered
  ##     by osiedle with t(G - 2) reference; 42-point second stage (osiedle effects
  ##     on osiedle deprivation, HC2); permutation of osiedle deprivation within
  ##     tertiles of distance to the centre (999).
  ## (3) Primary tests (Table 1, 23 cells) with Benjamini-Hochberg control.
  ## (4) Gelbach decomposition of the deprivation gradient (base: deprivation +
  ##     covariates; parts: composition, mix, intensity, configuration, green,
  ##     parking), bootstrap over osiedla (999).
  ## (5) Persistence: models refitted on monthly indices.
  ## (6) Robustness: flagged complexes excluded; low counter coverage excluded;
  ##     summer index without changed summer service; municipal-housing share in
  ##     place of deprivation with osiedle fixed effects; stop totals (script 41)
  ##     instead of same-trip indices.
  ## Outputs: output/stage2/stage2_primary_tests.csv, stage2_activity_models.csv,
  ##          stage2_deprivation_models.csv, stage2_gelbach.csv, stage2_persistence.csv,
  ##          stage2_robustness.csv

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines)})
  set.seed(20261007)
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  t0 <- Sys.time()

  ## ── data ─────────────────────────────────────────────────────────────────
  I   <- fread(file.path(s2, "stage2_stop_indices.csv"), encoding = "UTF-8")
  cvn <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[, .(stop_name, dist_km, flagged, residents, dep, osiedle_dom = osiedle,
                                                                                          share_0_18, share_61plus, spp_share)]
  td  <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")
  mun <- fread(file.path(s2, "municipal_catchments_network.csv"), encoding = "UTF-8")[, .(stop_name, muni_share)]
  m41 <- fread(file.path(s2, "stage2_complex_measures.csv"), encoding = "UTF-8")
  P   <- fread(file.path(out_dir, "output/stage1/s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, os_dep = dep, os_dist = dist_km)]
  C   <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  loc <- data.table(stop_name = C$stop_name, os_loc = osi$name[st_nearest_feature(C, osi)], X = st_coordinates(C)[, 1], Y = st_coordinates(C)[, 2])
  d <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(loc, I, cvn, td[, .(stop_name, fsi, med_gap, med_setback, med_block_area_ha, intersection_dens,
               deadend_share, wc_tree, b1 = ilr_res_vs_nonres, b2 = ilr_multi_vs_single, b3 = ilr_serv_vs_prod, b4 = ilr_retail_vs_work, b5 = ilr_office_vs_eduhealth,
               b6 = ilr_industry_vs_other, landuse_mix, share_university, share_dormitory, area_km2)], mun,
               m41[, .(stop_name, low_coverage, service_changed, volume_total = log(volume), tidal_total = tidal, onset_b_total = onset_b, conc_total = conc,
                       sat_total = sat, summer_total = summer_ratio)]))
  d[, osiedle := fifelse(!is.na(osiedle_dom) & osiedle_dom != "", osiedle_dom, os_loc)]
  d[, `:=`(no_res = residents %in% c(0, NA), lres = log1p(fcoalesce(as.numeric(residents), 0)), larea = log(area_km2))]
  d[no_res == TRUE, `:=`(share_0_18 = NA_real_, share_61plus = NA_real_)]
  for (v in c("share_0_18", "share_61plus")) d[is.na(get(v)), (v) := mean(d[[v]], na.rm = TRUE)]
  d[, `:=`(fsi = log(fsi + 0.01), med_gap = log1p(med_gap), med_setback = log1p(med_setback), med_block_area_ha = log(med_block_area_ha),
           intersection_dens = sqrt(intersection_dens), share_university = sqrt(share_university), share_dormitory = sqrt(share_dormitory),
           muni = sqrt(muni_share))]
  LU <- c("b1", "b2", "b3", "b4", "b5", "b6", "landuse_mix", "fsi", "med_gap", "med_setback", "med_block_area_ha", "intersection_dens", "deadend_share",
          "wc_tree", "share_university", "share_dormitory")
  BETWEEN <- c("b1", "landuse_mix", "wc_tree")
  COV <- c("share_0_18", "share_61plus", "lres", "larea", "no_res")
  d <- d[complete.cases(d[, c(LU, "dist_km", "larea"), with = FALSE])]
  cat("Complexes with complete land-use descriptors:", nrow(d), "| with residents:", d[no_res == FALSE, .N], "| osiedla:", uniqueN(d$osiedle), "\n")
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  for (v in c(LU, "dep", "spp_share", "muni", "share_0_18", "share_61plus", "lres", "larea")) d[, (v) := zs(get(v))]
  for (v in LU) d[, (paste0("w_", v)) := get(v) - ave(get(v), osiedle)]
  for (v in BETWEEN) d[, (paste0("m_", v)) := ave(get(v), osiedle)]
  ## principal components of the remaining core descriptors (intensity and configuration)
  core_rest <- c("fsi", "med_gap", "med_setback", "med_block_area_ha", "intersection_dens", "deadend_share")
  pca <- prcomp(d[, core_rest, with = FALSE])
  npc <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.8)[1]
  PCS <- paste0("PC", seq_len(npc)); d[, (PCS) := as.data.table(pca$x[, seq_len(npc), drop = FALSE])]
  cat("Principal components kept for intensity and configuration:", npc, "\n")
  DEPADJ <- c("b1", "landuse_mix", "wc_tree", "b2", "b3", "b4", "b5", "b6", PCS, "share_university", "share_dormitory", "spp_share")
  D <- as.matrix(dist(d[, .(X, Y)]))
  OUT <- c("volume", "tidal", "onset_b", "onset_a", "conc", "conc_b", "conc_a", "sat", "sat_b", "sat_a", "summer", "sept", "ferie")

  ## ── estimators ───────────────────────────────────────────────────────────
  rhs_cov <- function(dd) { cv <- COV; if (uniqueN(dd$no_res) < 2) cv <- setdiff(cv, "no_res"); c(cv, "ns(dist_km, 3)") }
  conley <- function(fit, idx, cutoff = 1000) {
    X <- model.matrix(fit); e <- resid(fit); B <- solve(crossprod(X)); K <- 1 - D[idx, idx] / cutoff; K[K < 0] <- 0
    Xe <- X * e; V <- B %*% (t(Xe) %*% K %*% Xe) %*% B; sqrt(diag(V))
  }
  cr2 <- function(fit, cl) {
    X <- model.matrix(fit); e <- resid(fit); B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) {
      i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); A <- ev$vectors %*% diag(1 / sqrt(pmax(ev$values, 1e-10)), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (A %*% e[i]); meat <- meat + u %*% t(u)
    }
    sqrt(diag(B %*% meat %*% B))
  }
  activity <- function(dd, y) {
    dd <- dd[is.finite(get(y))]; dd[, yy := zs(get(y))]
    f <- lm(reformulate(c(paste0("w_", LU), paste0("m_", BETWEEN), rhs_cov(dd)), "yy"), dd)
    used <- as.integer(rownames(model.matrix(f))); if (length(used) < nrow(dd)) cat("  [", y, "] rows dropped by lm:", nrow(dd) - length(used), "\n")
    idx <- match(dd$stop_name[used], d$stop_name); se <- conley(f, idx)
    k <- paste0("w_", LU)
    fe <- lm(reformulate(c(LU, rhs_cov(dd), "factor(osiedle)"), "yy"), dd)
    data.table(outcome = y, term = LU, est = coef(f)[k], se = se[k], est_osiedle_fe = coef(fe)[LU], n = nrow(dd))
  }
  deprivation <- function(dd, y, expo = "dep", adj = DEPADJ, perm = TRUE) {
    dd <- dd[no_res == FALSE & is.finite(get(y)) & is.finite(get(expo))]; dd[, yy := zs(get(y))]
    f <- lm(reformulate(c(expo, adj, rhs_cov(dd)), "yy"), dd)
    used <- as.integer(rownames(model.matrix(f))); se <- cr2(f, dd$osiedle[used])[expo]; G <- uniqueN(dd$osiedle[used]); tval <- coef(f)[[expo]] / se
    out <- data.table(outcome = y, exposure = expo, est = coef(f)[[expo]], se_cr2 = se, t = tval, df = G - 2, n = nrow(dd), osiedla = G)
    if (perm) {
      ## 42-point second stage: osiedle effects net of stop-level characteristics, on osiedle deprivation
      f2 <- lm(reformulate(c(adj, rhs_cov(dd), "0 + factor(osiedle)"), "yy"), dd)
      oe <- coef(f2)[grepl("^factor\\(osiedle\\)", names(coef(f2)))]
      so <- data.table(osiedle = sub("^factor\\(osiedle\\)", "", names(oe)), eff = unname(oe))
      so <- merge(merge(so, dd[, .(n = .N), by = osiedle], by = "osiedle"), P, by = "osiedle")[is.finite(eff)]
      s2f <- lm(eff ~ os_dep, so, weights = n)
      se2 <- sqrt(diag(sandwich::vcovHC(s2f, type = "HC2")))[["os_dep"]]
      out[, `:=`(second_stage_slope = coef(s2f)[["os_dep"]], second_stage_p = 2 * pt(-abs(coef(s2f)[["os_dep"]] / se2), nrow(so) - 2))]
      ## permutation of osiedle deprivation within distance tertiles
      pt_ <- copy(P)[, terc := cut(os_dist, quantile(os_dist, c(0, 1/3, 2/3, 1)), include.lowest = TRUE)]
      raw <- dd$dep; os_map <- match(dd$osiedle, pt_$osiedle)
      null <- replicate(999, {
        pp <- copy(pt_)[, os_dep := sample(os_dep), by = terc]
        dd[, dep_p := zs(pp$os_dep[os_map])]
        coef(lm(reformulate(c("dep_p", adj, rhs_cov(dd)), "yy"), dd))[["dep_p"]]
      })
      out[, perm_p := mean(abs(null) >= abs(est))]
    }
    out
  }

  ## ── (1)-(2) models for every index ───────────────────────────────────────
  AM <- rbindlist(lapply(OUT, function(y) activity(d, y)))
  AM[, `:=`(z = est / se, p_two = 2 * pnorm(-abs(est / se)))]
  fwrite(AM, file.path(s2, "stage2_activity_models.csv"))
  DM <- rbindlist(lapply(OUT, function(y) deprivation(d, y)))
  DM[, p_two := 2 * pt(-abs(t), df)]
  fwrite(DM, file.path(s2, "stage2_deprivation_models.csv"))
  cat("\nDeprivation path (SD of index per SD of catchment deprivation; CR2 t with G-2 df; 42-point second stage; spatially constrained permutation):\n")
  print(DM[, .(outcome, n, est = round(est, 3), se = round(se_cr2, 3), p = round(p_two, 3), stage2_slope = round(second_stage_slope, 3), stage2_p = round(second_stage_p, 3), perm_p = round(perm_p, 3))])

  ## ── (3) primary tests ────────────────────────────────────────────────────
  tab1 <- data.table(
    id = c(paste0("A", 1:17), paste0("D", 1:6)),
    outcome = c(rep("volume", 7), "tidal", "tidal", "onset_a", "conc_a", "conc_a", "sat_a", "sat_a", "summer", "summer", "summer",
                "volume", "tidal", "onset_b", "conc_b", "sat_b", "summer"),
    term = c("fsi", "landuse_mix", "intersection_dens", "med_gap", "med_setback", "med_block_area_ha", "wc_tree", "b1", "b2", "b3", "landuse_mix", "b1",
             "b4", "b3", "share_university", "share_dormitory", "b4", rep("dep", 6)),
    sign = c(1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 0, 1, 1, 0, 1, 0))
  pt1 <- rbind(merge(tab1[term != "dep"], AM[, .(outcome, term, est, se, stat = z)], by = c("outcome", "term"))[, ref := "normal"],
               merge(tab1[term == "dep"], DM[, .(outcome, est, se = se_cr2, stat = t, df)], by = "outcome")[, ref := "t"], fill = TRUE)
  pt1[, p := fifelse(ref == "normal",
                     fifelse(sign == 0, 2 * pnorm(-abs(stat)), pnorm(-sign * stat)),
                     fifelse(sign == 0, 2 * pt(-abs(stat), df), pt(-sign * stat, df)))]
  pt1[, q := p.adjust(p, "BH")][, supported := q < 0.05 & (sign == 0 | sign(est) == sign)]
  setorder(pt1, id)
  fwrite(pt1, file.path(s2, "stage2_primary_tests.csv"))
  cat("\nPrimary tests (Table 1): one-sided where a sign is expected; BH across the 23 tests\n")
  print(pt1[, .(id, outcome, term, expected = fifelse(sign == 0, "two-sided", fifelse(sign > 0, "+", "-")), est = round(est, 3), se = round(se, 3), p = signif(p, 2), q = signif(q, 2), supported)])

  ## ── (4) Gelbach decomposition ────────────────────────────────────────────
  GR <- list(composition = c("b1", "b2", "b3", "b4", "b5", "b6", "share_university", "share_dormitory"), mix = "landuse_mix", intensity = "fsi",
             configuration = c("med_gap", "med_setback", "med_block_area_ha", "intersection_dens", "deadend_share"), green = "wc_tree", parking = "spp_share")
  gel <- function(dd, y) {
    dd <- dd[no_res == FALSE & is.finite(get(y))]; dd[, yy := zs(get(y))]
    base <- coef(lm(reformulate(c("dep", rhs_cov(dd)), "yy"), dd))[["dep"]]
    full <- coef(lm(reformulate(c("dep", unname(unlist(GR)), rhs_cov(dd)), "yy"), dd))
    gam <- vapply(unname(unlist(GR)), function(k) coef(lm(reformulate(c("dep", rhs_cov(dd)), k), dd))[["dep"]], numeric(1))
    parts <- vapply(GR, function(g) sum(gam[g] * full[g]), numeric(1))
    c(base = base, remainder = full[["dep"]], parts)
  }
  PRIM <- c("volume", "tidal", "onset_b", "conc_b", "sat_b", "summer")
  ixo <- split(seq_len(nrow(d)), d$osiedle); on <- names(ixo)
  GL <- rbindlist(lapply(PRIM, function(y) {
    e <- gel(d, y)
    bt <- t(replicate(999, { k <- sample(on, replace = TRUE); gel(d[unlist(ixo[k])], y) }))
    data.table(outcome = y, part = names(e), est = e, lo = apply(bt, 2, quantile, .025, na.rm = TRUE), hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
  }))
  fwrite(GL, file.path(s2, "stage2_gelbach.csv"))
  GL[, s := sprintf("%5.2f [%5.2f,%5.2f]", est, lo, hi)]
  cat("\nGelbach decomposition of the deprivation gradient (SD units):\n")
  print(dcast(GL, outcome ~ factor(part, levels = c("base", names(GR), "remainder")), value.var = "s"), width = 250)

  ## ── (5) persistence ──────────────────────────────────────────────────────
  mon <- fread(file.path(s2, "stage2_stop_indices_monthly.csv"), encoding = "UTF-8")
  PS <- rbindlist(lapply(1:10, function(m) {
    mm <- dcast(mon[month == m], stop_name ~ measure, value.var = "index")
    dm <- merge(d[, !c("volume", "tidal", "onset_b", "conc_b", "sat_b"), with = FALSE], mm, by = "stop_name")
    rbind(rbindlist(lapply(c("volume", "tidal"), function(y) activity(dm, y)))[, .(month = m, outcome, term, est)],
          rbindlist(lapply(c("volume", "tidal", "onset_b", "conc_b", "sat_b"), function(y) deprivation(dm, y, perm = FALSE)))[, .(month = m, outcome, term = "dep", est)])
  }))
  PSs <- merge(PS, pt1[, .(id, outcome, term, pooled = est)], by = c("outcome", "term"))[, .(months = .N, same_sign = sum(sign(est) == sign(pooled)),
                                                                                               min = min(est), max = max(est)), by = .(id, outcome, term)]
  setorder(PSs, id); fwrite(PSs, file.path(s2, "stage2_persistence.csv"))
  cat("\nPersistence: months in which the monthly estimate has the sign of the pooled estimate:\n"); print(PSs, digits = 2)

  ## ── (6) robustness ───────────────────────────────────────────────────────
  rb <- function(dd, label) {
    rbind(merge(tab1[term != "dep"], rbindlist(lapply(unique(tab1[term != "dep", outcome]), function(y) activity(dd, y)))[, .(outcome, term, est)], by = c("outcome", "term")),
          merge(tab1[term == "dep"], rbindlist(lapply(unique(tab1[term == "dep", outcome]), function(y) deprivation(dd, y, perm = FALSE)))[, .(outcome, est)], by = "outcome"))[, check := label]
  }
  RB <- rbind(rb(d[flagged == FALSE], "without flagged complexes"), rb(d[low_coverage == FALSE], "without low counter coverage"))
  RB <- rbind(RB, deprivation(d[service_changed == FALSE], "summer", perm = FALSE)[, .(outcome, est, id = "D6", term = "dep", check = "summer without changed service")], fill = TRUE)
  MU <- rbindlist(lapply(PRIM, function(y) {
    dd <- d[no_res == FALSE & is.finite(get(y))]; dd[, yy := zs(get(y))]
    f <- lm(reformulate(c("muni", DEPADJ, rhs_cov(dd), "factor(osiedle)"), "yy"), dd)
    used <- as.integer(rownames(model.matrix(f)))
    data.table(outcome = y, est = coef(f)[["muni"]], se = cr2(f, dd$osiedle[used])[["muni"]])
  }))[, `:=`(p_two = 2 * pt(-abs(est / se), uniqueN(d$osiedle) - 2), check = "municipal-housing share within osiedla")]
  TOT <- rbindlist(lapply(c("volume_total", "tidal_total", "onset_b_total", "conc_total", "sat_total", "summer_total"), function(y) deprivation(d, y, perm = FALSE)))[, check := "stop totals"]
  fwrite(rbind(RB, MU, TOT[, .(outcome, est, se = se_cr2, check)], fill = TRUE), file.path(s2, "stage2_robustness.csv"))
  cat("\nRobustness of the primary estimates:\n")
  print(dcast(merge(pt1[, .(id, pooled = round(est, 3))], RB[, .(id, check, est = round(est, 3))], by = "id"), id + pooled ~ check, value.var = "est"))
  cat("\nMunicipal-housing share within osiedla (osiedle fixed effects, CR2):\n"); print(MU[, .(outcome, est = round(est, 3), se = round(se, 3), p = round(p_two, 3))])
  cat("\nDeprivation path on stop totals (not same-trip):\n"); print(TOT[, .(outcome, est = round(est, 3), se = round(se_cr2, 3), p = round(2 * pt(-abs(t), df), 3))])
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 48 · Early boardings, sharpened ────────────────────────────────────────────────────────────────────────
section("48", "Early boardings, sharpened", {
  ## Sharpening the deprivation-path estimate for early boardings (Table 1, D3)
  ## with four pre-stated refinements, applied to every measure where they make
  ## sense and reported in full (exploratory; nothing is selected after the fact):
  ##   1. more data: onset (and volume, tidal, concentration) over all regular
  ##      working days, January-October (Monday-Friday without holidays, bridge
  ##      and daylight-saving days), next to the primary term-time Monday-Thursday;
  ##   2. precision weights: each stop index carries a sampling variance from its
  ##      counts (effective counts = expanded counts x counted / operated visits);
  ##      the deprivation model is also fitted by weighted least squares with
  ##      weights 1 / (sampling variance + between-stop variance, method of moments);
  ##   3. finer timing from actual departure times in 15-minute bins (script 47):
  ##      share of boardings before 06:00, share in 05:00-05:59, and the onset
  ##      quantiles (minute by which 5% / 10% of a stop's daily boardings on a line
  ##      and direction have occurred), each relative to the same trips;
  ##   4. residential stops only (at least 80% residential floor in the catchment).
  ## Deprivation model, covariates and inference as in script 43 (CR2 by osiedle,
  ## t with G - 2 df). For the onset variants also: the unemployment and poverty
  ## components separately (catchment means of the pooled 2019/2021 z-scores) and
  ## the 42-point second stage.
  ## Expected signs (Table 1): onset shares +, onset quantiles - (earlier); tidal +,
  ## Saturday +; volume, concentration and summer two-sided.
  ## Outputs: output/stage2/sharpen_indices.csv, sharpen_results.csv,
  ##          output/figures/Fig_sharpen_onset.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  t0 <- Sys.time()

  ## ── step 1: indices with sampling variances ──────────────────────────────
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))
  cells[, dp := departures - uncovered]; cells <- cells[dp > 0]
  cells[, eff := pmin(1, counted / dp)]                                   # share of operated visits actually counted
  index_v <- function(dt, gcols) {
    dt <- dt[is.finite(q) & w > 0 & is.finite(v)]
    dt[, n_in := .N, by = gcols]; dt <- dt[n_in >= 2]
    dt[, dev := q - weighted.mean(q, w), by = gcols]
    dt[, .(index = sum(w * dev) / sum(w), var = sum(w^2 * v) / sum(w)^2), by = stop_name]
  }
  lgt <- function(k, n) qlogis((k + 0.5) / (n + 1))
  vlg <- function(k, n) 1 / (k + 0.5) + 1 / (n - k + 0.5)
  day_measures <- function(C) {
    D <- C[night == FALSE]; out <- list()
    v <- C[, .(y = sum(B + A), dp = sum(dp), e = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, period, stop_name)]
    out$volume <- index_v(v[, `:=`(q = log((y + 0.5) / dp), w = dp, v = 1 / (y * e + 0.5))], c("NR_LINI", "dir", "period"))
    am <- D[sched_hour %in% 6:8, .(B = sum(B), A = sum(A), e = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)]
    pm <- D[sched_hour %in% 15:17, .(B = sum(B), A = sum(A), e = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)]
    ia <- index_v(am[, `:=`(q = lgt(B, B + A), w = B + A, v = vlg(B * e, (B + A) * e))], c("NR_LINI", "dir"))
    ip <- index_v(pm[, `:=`(q = lgt(B, B + A), w = B + A, v = vlg(B * e, (B + A) * e))], c("NR_LINI", "dir"))
    out$tidal <- merge(ia, ip, by = "stop_name")[, .(stop_name, index = index.x - index.y, var = var.x + var.y)]
    early <- D[sched_hour < 6, .(de = sum(dp)), by = .(NR_LINI, dir)][de > 0]
    on <- merge(D[, .(Be = sum(B[sched_hour < 6]), Ball = sum(B), e = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)], early, by = c("NR_LINI", "dir"))
    out$onset_b <- index_v(on[, `:=`(q = lgt(Be, Ball), w = Ball, v = vlg(Be * e, Ball * e))], c("NR_LINI", "dir"))
    cc <- D[, .(Bp = sum(B[sched_hour %in% 7:8]), dpp = sum(dp[sched_hour %in% 7:8]), Bm = sum(B[sched_hour %in% 10:13]), dpm = sum(dp[sched_hour %in% 10:13]),
                e = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)][dpp > 0 & dpm > 0]
    out$conc_b <- index_v(cc[, `:=`(q = log((Bp + 0.5) / dpp) - log((Bm + 0.5) / dpm), w = dpp + dpm, v = 1 / (Bp * e + 0.5) + 1 / (Bm * e + 0.5))], c("NR_LINI", "dir"))
    out
  }
  term <- cells[season == "term" & daytype == "monthu"]
  work <- cells[daytype %in% c("monthu", "fri")]
  R <- list()
  for (k in c("volume", "tidal", "onset_b", "conc_b")) { R[[paste0(k, "|term")]] <- day_measures(term)[[k]] }
  dw <- day_measures(work); for (k in names(dw)) R[[paste0(k, "|workdays")]] <- dw[[k]]
  ## Saturday and summer (single definitions)
  sa <- merge(cells[season == "term" & daytype == "sat" & night == FALSE, .(Bs = sum(B), ds = sum(dp), es = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)],
              term[night == FALSE, .(Bw = sum(B), dw = sum(dp), ew = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)], by = c("NR_LINI", "dir", "stop_name"))
  R[["sat_b|term"]] <- index_v(sa[, `:=`(q = log((Bs + 0.5) / ds) - log((Bw + 0.5) / dw), w = ds + dw, v = 1 / (Bs * es + 0.5) + 1 / (Bw * ew + 0.5))], c("NR_LINI", "dir"))
  su <- merge(cells[season == "summer" & daytype == "monthu", .(y1 = sum(B + A), d1 = sum(dp), e1 = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)],
              term[, .(y0 = sum(B + A), d0 = sum(dp), e0 = sum(eff * dp) / sum(dp)), by = .(NR_LINI, dir, stop_name)], by = c("NR_LINI", "dir", "stop_name"))
  R[["summer|term"]] <- index_v(su[, `:=`(q = log((y1 + 0.5) / d1) - log((y0 + 0.5) / d0), w = d1 + d0, v = 1 / (y1 * e1 + 0.5) + 1 / (y0 * e0 + 0.5))], c("NR_LINI", "dir"))

  ## fine timing (counted trips, actual time)
  fine <- rbindlist(lapply(sprintf("%02d", 1:10), function(m) readRDS(file.path(s2, "fine", sprintf("fine15_2025_%s.rds", m)))))
  fine_measures <- function(F) {
    g <- F[, .(B = sum(B), visits = sum(visits)), by = .(NR_LINI, dir, cx, bin)]
    early <- g[bin < 24, .(ve = sum(visits)), by = .(NR_LINI, dir)][ve > 0]
    tot <- merge(g[, .(Ball = sum(B), Be = sum(B[bin < 24]), B5 = sum(B[bin %in% 20:23])), by = .(NR_LINI, dir, cx)], early, by = c("NR_LINI", "dir"))
    setnames(tot, "cx", "stop_name")
    out <- list()
    out$onset_actual <- index_v(copy(tot)[, `:=`(q = lgt(Be, Ball), w = Ball, v = vlg(Be, Ball))], c("NR_LINI", "dir"))
    out$window05 <- index_v(copy(tot)[, `:=`(q = lgt(B5, Ball), w = Ball, v = vlg(B5, Ball))], c("NR_LINI", "dir"))
    ## onset quantiles: minute by which 5% / 10% of the day's boardings have occurred (linear within the 15-minute bin)
    g <- g[order(NR_LINI, dir, cx, bin)][, `:=`(cum = cumsum(B), tot = sum(B)), by = .(NR_LINI, dir, cx)]
    qq <- function(p) g[tot >= 20, { k <- which(cum >= p * tot[1])[1]; prev <- if (k > 1) cum[k - 1] else 0
                                   .(q = 15 * bin[k] + 15 * (p * tot[1] - prev) / max(B[k], 1e-9), w = tot[1]) }, by = .(NR_LINI, dir, stop_name = cx)]
    out$q05 <- index_v(qq(0.05)[, v := 1 / w], c("NR_LINI", "dir"))
    out$q10 <- index_v(qq(0.10)[, v := 1 / w], c("NR_LINI", "dir"))
    out
  }
  ft <- fine_measures(fine[term_monthu == TRUE]); for (k in names(ft)) R[[paste0(k, "|term")]] <- ft[[k]]
  fw <- fine_measures(fine); for (k in names(fw)) R[[paste0(k, "|workdays")]] <- fw[[k]]
  IX <- rbindlist(lapply(names(R), function(k) R[[k]][, `:=`(measure = sub("\\|.*", "", k), window = sub(".*\\|", "", k))]))
  fwrite(IX, file.path(s2, "sharpen_indices.csv"))
  cat("Indices built:", paste(names(R), collapse = ", "), "\n")

  ## ── covariates (as script 43) and deprivation components ─────────────────
  cvn <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  td  <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")
  d <- merge(cvn[residents > 0, .(stop_name, dist_km, residents, dep, osiedle, share_0_18, share_61plus, spp_share)],
             td[, .(stop_name, fsi, med_gap, med_setback, med_block_area_ha, intersection_dens, deadend_share, wc_tree, b1 = ilr_res_vs_nonres,
                    b2 = ilr_multi_vs_single, b3 = ilr_serv_vs_prod, b4 = ilr_retail_vs_work, b5 = ilr_office_vs_eduhealth, b6 = ilr_industry_vs_other,
                    landuse_mix, share_university, share_dormitory, area_km2)], by = "stop_name")
  d <- d[complete.cases(d)]
  d[, residential := b1 >= sqrt(1/2) * log(4)]                              # at least 80% residential floor
  ## unemployment and poverty components: population-weighted catchment means of the pooled 2019/2021 z-scores (signed: + = more deprived)
  AS <- fread(file.path(out_dir, "output/stage1/area_status_osiedla.csv"), encoding = "UTF-8")[, .(osiedle, unemp = -A_z, pov = -B_z)]
  vc  <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  pp  <- st_transform(st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE), 2180)
  pp$pop <- pp[[grep("ZAMELDOWANYCH", names(pp), value = TRUE)[1]]]; pp$os <- osi$name[st_nearest_feature(pp, osi)]
  jn <- as.data.table(st_drop_geometry(st_join(pp[, c("pop", "os")], st_transform(vc[, "stop_name"], 2180), join = st_within)))[!is.na(stop_name) & pop > 0]
  jn <- merge(jn, AS, by.x = "os", by.y = "osiedle")
  d <- merge(d, jn[, .(unemp = sum(pop * unemp) / sum(pop), pov = sum(pop * pov) / sum(pop)), by = stop_name], by = "stop_name")
  cat(sprintf("Stops with residents and complete covariates: %d (residential: %d) | rho(unemployment, poverty components) = %.2f\n",
              nrow(d), d[residential == TRUE, .N], cor(d$unemp, d$pov, method = "spearman")))
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  d[, `:=`(fsi = log(fsi + 0.01), med_gap = log1p(med_gap), med_setback = log1p(med_setback), med_block_area_ha = log(med_block_area_ha),
           intersection_dens = sqrt(intersection_dens), share_university = sqrt(share_university), share_dormitory = sqrt(share_dormitory),
           lres = log(residents), larea = log(area_km2))]
  num <- c("fsi", "med_gap", "med_setback", "med_block_area_ha", "intersection_dens", "deadend_share", "wc_tree", paste0("b", 1:6), "landuse_mix",
           "share_university", "share_dormitory", "spp_share", "share_0_18", "share_61plus", "lres", "larea", "dep", "unemp", "pov")
  d[, (paste0(num, "_z")) := lapply(.SD, zs), .SDcols = num]
  pca <- prcomp(d[, paste0(c("fsi", "med_gap", "med_setback", "med_block_area_ha", "intersection_dens", "deadend_share"), "_z"), with = FALSE])
  npc <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.8)[1]; PCS <- paste0("PC", seq_len(npc)); d[, (PCS) := as.data.table(pca$x[, seq_len(npc), drop = FALSE])]
  ADJ <- c(paste0(c("b1", "landuse_mix", "wc_tree", "b2", "b3", "b4", "b5", "b6", "share_university", "share_dormitory", "spp_share",
                    "share_0_18", "share_61plus", "lres", "larea"), "_z"), PCS, "ns(dist_km, 3)")
  P <- fread(file.path(out_dir, "output/stage1/s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, os_dep = dep)]

  ## ── step 2: deprivation model, unweighted and precision-weighted ─────────
  cr2 <- function(fit, cl, w = NULL) {
    X <- model.matrix(fit); e <- resid(fit); if (!is.null(w)) { X <- X * sqrt(w); e <- e * sqrt(w) }
    B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) {
      i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); A <- ev$vectors %*% diag(1 / sqrt(pmax(ev$values, 1e-10)), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (A %*% e[i]); meat <- meat + u %*% t(u)
    }
    sqrt(diag(B %*% meat %*% B))
  }
  fitdep <- function(dd, expo, weighted) {
    dd <- dd[is.finite(index) & is.finite(var)]
    fml <- reformulate(c(expo, ADJ), "index")
    f0 <- lm(fml, dd)
    if (!weighted) { f <- f0; w <- NULL } else {
      tau2 <- max(var(resid(f0)) - mean(dd$var), 0.1 * var(resid(f0)))
      dd[, wt := 1 / (var + tau2)]; f <- lm(fml, dd, weights = wt); w <- dd$wt
    }
    used <- as.integer(rownames(model.matrix(f)))
    se <- cr2(f, dd$osiedle[used], if (is.null(w)) NULL else w[used])[[expo]]
    sdy <- sd(dd$index)
    data.table(est = coef(f)[[expo]] / sdy, se = se / sdy, t = coef(f)[[expo]] / se, df = uniqueN(dd$osiedle[used]) - 2, n = length(used))
  }
  sign_of <- c(volume = 0, tidal = 1, onset_b = 1, conc_b = 0, sat_b = 1, summer = 0, onset_actual = 1, window05 = 1, q05 = -1, q10 = -1)
  grid <- CJ(key = names(R), sample = c("all stops with residents", "residential stops"), weighting = c("unweighted", "precision-weighted"), sorted = FALSE)
  res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    k <- grid$key[i]; dd <- merge(d, R[[k]], by = "stop_name")
    if (grid$sample[i] == "residential stops") dd <- dd[residential == TRUE]
    if (grepl("^q", k) && grid$weighting[i] == "precision-weighted") return(NULL)  # no analytic sampling variance for quantiles
    ex <- c("dep_z", if (grepl("^(onset|window|q[0-9])", k)) c("unemp_z", "pov_z"))
    rbindlist(lapply(ex, function(e) cbind(data.table(measure = sub("\\|.*", "", k), window = sub(".*\\|", "", k), sample = grid$sample[i],
                                                      weighting = grid$weighting[i], exposure = sub("_z$", "", e)), fitdep(dd, e, grid$weighting[i] == "precision-weighted"))))
  }))
  res[, expected := sign_of[measure]]
  res[, p := fifelse(expected == 0, 2 * pt(-abs(t), df), pt(-expected * t, df))]
  fwrite(res, file.path(s2, "sharpen_results.csv"))

  ## 42-point second stage for the onset variants (all stops, unweighted adjustment, osiedle effects on osiedle deprivation)
  s2s <- rbindlist(lapply(grep("^(onset|window|q[0-9])", names(R), value = TRUE), function(k) {
    dd <- merge(d, R[[k]], by = "stop_name")[is.finite(index)]
    f2 <- lm(reformulate(c(ADJ, "0 + factor(osiedle)"), "index"), dd)
    oe <- coef(f2)[grepl("^factor\\(osiedle\\)", names(coef(f2)))]
    so <- merge(merge(data.table(osiedle = sub("^factor\\(osiedle\\)", "", names(oe)), eff = unname(oe)), dd[, .(n = .N), by = osiedle], by = "osiedle"), P, by = "osiedle")[is.finite(eff)]
    s <- lm(eff ~ os_dep, so, weights = n); se <- sqrt(diag(sandwich::vcovHC(s, type = "HC2")))[["os_dep"]]
    data.table(measure = sub("\\|.*", "", k), window = sub(".*\\|", "", k), slope = coef(s)[["os_dep"]] / sd(dd$index), p_two = 2 * pt(-abs(coef(s)[["os_dep"]] / se), nrow(so) - 2), osiedla = nrow(so))
  }))
  fwrite(s2s, file.path(s2, "sharpen_second_stage.csv"))

  ## ── report ───────────────────────────────────────────────────────────────
  cat("\nDeprivation path, SD of index per SD of exposure (CR2, t with G-2 df; p one-sided where Table 1 states a sign):\n")
  show <- res[exposure == "dep", .(measure, window, sample, weighting, est = round(est, 3), se = round(se, 3), p = signif(p, 2), n)]
  print(dcast(show, measure + window ~ sample + weighting, value.var = c("est", "p")), width = 250)
  cat("\nOnset variants, by exposure (all stops with residents):\n")
  print(dcast(res[grepl("^(onset|window|q[0-9])", measure) & sample == "all stops with residents",
                  .(measure, window, weighting, exposure, s = sprintf("%5.3f (p %.3f)", est, p))], measure + window + weighting ~ exposure, value.var = "s"), width = 250)
  cat("\n42-point second stage for the onset variants (two-sided):\n"); print(s2s, digits = 3)

  ## figure: estimates for the onset family with 95% intervals
  fo <- res[grepl("^(onset|window|q[0-9])", measure) & exposure == "dep"]
  fo[, lab := factor(paste0(c(onset_b = "Share before 06:00 (scheduled hour)", onset_actual = "Share before 06:00 (actual time)", window05 = "Share 05:00-05:59",
                              q05 = "5% onset quantile (minutes; - = earlier)", q10 = "10% onset quantile (minutes; - = earlier)")[measure], ", ", window))]
  fo[, `:=`(lo = est - qt(0.975, df) * se, hi = est + qt(0.975, df) * se)]
  g <- ggplot(fo, aes(est, lab, colour = weighting, shape = sample)) + geom_vline(xintercept = 0, colour = "#898781", linewidth = 0.3) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), size = 0.35) +
    scale_colour_manual(values = c("unweighted" = "#898781", "precision-weighted" = "#e34948"), name = NULL) +
    scale_shape_manual(values = c("all stops with residents" = 16, "residential stops" = 17), name = NULL) +
    labs(x = "SD of the same-trip index per SD of catchment deprivation (95% interval, CR2)", y = NULL,
         title = "Early boardings and deprivation: refined measures, weights and samples",
         caption = "Positive = more early boardings in more deprived catchments for shares; negative = earlier onset for quantiles. Stops served by the same trips compared; January-October 2025.") +
    theme_minimal(base_size = 10) + theme(plot.background = element_rect(fill = "#ffffff", colour = NA), legend.position = "top",
                                          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), plot.caption = element_text(colour = "#898781", size = 7.5))
  ggsave(file.path(out_dir, "output/figures/Fig_sharpen_onset.png"), g, width = 10, height = 6.5, dpi = 200, bg = "#ffffff")
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 52 · Responses to leave days, service cuts and weather ─────────────────────────────────────────────────
section("52", "Responses to leave days, service cuts and weather", {
  ## Behavioural differences between riders from more and less deprived
  ## catchments (exploratory, secondary to the main analysis). January-October
  ## 2025, day lines unless stated; counted visits from script 40 (base tables),
  ## night lines from script 41 cells, daily weather from IMGW (Poznań-Ławica,
  ## station 330; data/raw/weather/s_d_330_2025.csv).
  ##
  ## (1) Leave days (job flexibility). Same-trip ratio of counted boardings per
  ##     counted visit on days many workers take off, against regular days of the
  ##     same weekday: bridge days (2-3 January, 2 May, 20 June) against the other
  ##     Thursdays/Fridays of their month (outside school breaks), and the
  ##     long-weekend days 28-30 April against the Monday-Wednesdays of 7-16 April
  ##     and 5-14 May. Whole day and 06-09 (outbound commute).
  ## (2) Regularity. Term-time Monday-Thursday: daily log boardings per counted
  ##     visit for each stop, line and direction, net of line-direction means and
  ##     city-wide daily shocks; day-to-day variance of the stop's daily mean minus
  ##     its expected sampling variance (Poisson) = excess day-to-day variability.
  ## (3) Night travel. Night lines (2xx), term months: same-trip volume index
  ##     (boardings + alightings per departure within night line x direction) on
  ##     weeknights (Monday-Thursday) and weekend nights (Saturday and Sunday
  ##     00-04, i.e. Friday and Saturday nights); night intensity relative to the
  ##     day (log per-departure volume, night minus day lines).
  ## (4) Service dependence. Daily riders per departure at each stop against the
  ##     stop's daily departures (day lines), net of stop x day-type and date
  ##     effects: elasticity b. b = -1: riders stay when service is cut (they
  ##     concentrate on the remaining departures); b = 0: riders fall in
  ##     proportion. Variation comes from the 2025 track works, diversions and line
  ##     changes. Stops with enough variation (SD of net log departures >= 0.05).
  ## (5) Weather. Daily log riders per departure against rain (>= 1 mm), heavy
  ##     rain (>= 10 mm), frost (mean temperature < 0 C), snow cover and heat
  ##     (maximum >= 30 C), with stop x day type x month and school-season effects;
  ##     per-stop sensitivities.
  ## Each stop-level measure is related to catchment deprivation with the
  ## deprivation-path model of script 43 (land-use adjustment, covariates, CR2 by
  ## osiedle, t with G - 2 df), and with the covariates only.
  ## Outputs: output/stage2/behaviour_stop_measures.csv, behaviour_tests.csv;
  ##          output/figures/Fig_behaviour.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  t0 <- Sys.time()
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  gt <- unique(fread(file.path(GTFS_DIR, "stops.txt"), encoding = "UTF-8")[, .(code = sub("[0-9]+$", "", stop_code), stop_name)])
  gt <- gt[stop_name %in% C$stop_name][!duplicated(code)]
  hol <- as.IDate(c("2025-01-01", "2025-01-06", "2025-04-20", "2025-04-21", "2025-05-01", "2025-05-03", "2025-06-08", "2025-06-19", "2025-08-15"))
  bridge <- as.IDate(c("2025-01-02", "2025-01-03", "2025-05-02", "2025-06-20")); dst <- as.IDate(c("2025-03-30", "2025-10-26"))
  season <- function(D) fcase(D >= as.IDate("2025-01-20") & D <= as.IDate("2025-02-02"), "ferie", D >= as.IDate("2025-04-17") & D <= as.IDate("2025-04-22"), "easter",
                              D >= as.IDate("2025-06-28") & D <= as.IDate("2025-08-31"), "summer", D >= as.IDate("2025-09-01") & D <= as.IDate("2025-09-30"), "september",
                              D >= as.IDate("2025-03-01"), "term", default = "winter")
  dtype <- function(D) { w <- as.integer(format(as.Date(D), "%u")); fcase(D %in% hol | w == 7, "sunhol", D %in% bridge, "bridge", w == 6, "sat", w == 5, "fri", default = "monthu") }

  ## ── counted visits, per month ────────────────────────────────────────────
  SD <- list(); LN <- list()
  for (m in sprintf("%02d", 1:10)) {
    b <- readRDS(file.path(s2, "base", sprintf("peka_base_2025_%s.rds", m)))
    b <- b[!grepl("^2[0-9][0-9]$", NR_LINI) & sched_hour %in% 5:21 & !DATA %in% dst]
    b[gt, cx := i.stop_name, on = c(SYM_PRZYST = "code")]; b <- b[!is.na(cx)]
    SD[[m]] <- b[, .(N = sum(visits), Cn = sum(counted), B = sum(B)), by = .(cx, DATA)]                       # stop-day totals
    LN[[m]] <- b[counted > 0, .(Cn = sum(counted), B = sum(B), B_am = sum(B[sched_hour %in% 6:8]), C_am = sum(counted[sched_hour %in% 6:8])),
                 by = .(cx, NR_LINI, dir, DATA)]                                                              # stop x line x direction x day
    rm(b); gc()
  }
  SD <- rbindlist(SD); LN <- rbindlist(LN)
  SD[, `:=`(daytype = dtype(DATA), season = season(DATA), month = as.integer(format(as.Date(DATA), "%m")))]
  LN[, `:=`(daytype = dtype(DATA), season = season(DATA), month = as.integer(format(as.Date(DATA), "%m")), wd = as.integer(format(as.Date(DATA), "%u")))]
  cat("Stop-days:", nrow(SD), "| stop-line-direction-days:", nrow(LN), "\n")

  ## same-trip index helper: deviation from the weighted mean of the group, averaged per stop
  st_index <- function(dt, gcols) {
    dt <- dt[is.finite(q) & w > 0]; dt[, n_in := .N, by = gcols]; dt <- dt[n_in >= 2]
    dt[, dev := q - weighted.mean(q, w), by = gcols]; dt[, .(index = weighted.mean(dev, w), groups = .N), by = .(stop_name = cx)]
  }
  lr <- function(b, c) log((b + 0.5) / c)

  ## ── (1) leave days ───────────────────────────────────────────────────────
  leave <- function(test, ref, am = FALSE) {
    a <- LN[DATA %in% test, .(B = sum(if (am) B_am else B), Cn = sum(if (am) C_am else Cn)), by = .(cx, NR_LINI, dir)]
    r <- LN[DATA %in% ref, .(Br = sum(if (am) B_am else B), Cr = sum(if (am) C_am else Cn)), by = .(cx, NR_LINI, dir)]
    x <- merge(a, r, by = c("cx", "NR_LINI", "dir"))[Cn >= 3 & Cr >= 6]
    x[, `:=`(q = lr(B, Cn) - lr(Br, Cr), w = 1 / (1 / Cn + 1 / Cr))]
    st_index(x, c("NR_LINI", "dir"))
  }
  alld <- sort(unique(LN$DATA))
  ref_bridge <- alld[(format(alld, "%m") %in% c("01", "05", "06")) & as.integer(format(as.Date(alld), "%u")) %in% 4:5 &
                     !alld %in% c(hol, bridge) & !season(alld) %in% c("ferie", "summer") & !(alld >= as.IDate("2025-01-01") & alld <= as.IDate("2025-01-06"))]
  maj <- as.IDate(c("2025-04-28", "2025-04-29", "2025-04-30"))
  ref_maj <- alld[as.integer(format(as.Date(alld), "%u")) %in% 1:3 & ((alld >= as.IDate("2025-04-07") & alld <= as.IDate("2025-04-16")) | (alld >= as.IDate("2025-05-05") & alld <= as.IDate("2025-05-14")))]
  cat("Reference days: bridge", length(ref_bridge), "| long weekend", length(ref_maj), "\n")
  M <- list(bridge = leave(bridge, ref_bridge), bridge_am = leave(bridge, ref_bridge, TRUE), longweekend = leave(maj, ref_maj), longweekend_am = leave(maj, ref_maj, TRUE))

  ## ── (2) regularity ───────────────────────────────────────────────────────
  R <- LN[daytype == "monthu" & season == "term" & Cn >= 2]
  R[, y := lr(B, Cn)]
  R[, y := y - weighted.mean(y, Cn), by = .(cx, NR_LINI, dir)]                                             # line-direction means
  R[, y := y - weighted.mean(y, Cn), by = DATA]                                                            # city-wide daily shocks
  R[, sv := 1 / (B + 0.5)]
  RD <- R[, .(v = weighted.mean(y, Cn), sv = sum(Cn^2 * sv) / sum(Cn)^2, Cn = sum(Cn), B = sum(B)), by = .(cx, DATA)]
  REG <- RD[, .(days = .N, total_var = var(v), samp_var = mean(sv), B_day = mean(B)), by = .(stop_name = cx)][days >= 40]
  REG[, `:=`(excess_var = total_var - samp_var, excess_sd = sign(total_var - samp_var) * sqrt(abs(total_var - samp_var)))]
  cat("Regularity: stops", nrow(REG), "| median total SD", round(median(sqrt(REG$total_var)), 3), "| median sampling SD", round(median(sqrt(REG$samp_var)), 3),
      "| median excess SD", round(median(REG$excess_sd), 3), "\n")
  M$irregularity <- REG[, .(stop_name, index = excess_sd)]

  ## ── (3) night travel (cells: night lines) ───────────────────────────────
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))[month %in% c(3:6, 10)]
  cells[, dp := departures - uncovered]; cells <- cells[dp > 0]
  ng <- function(dt) { x <- dt[, .(y = sum(B + A), dep = sum(dp)), by = .(cx = stop_name, NR_LINI, dir)][, `:=`(q = log((y + 0.5) / dep), w = dep)]; st_index(x, c("NR_LINI", "dir")) }
  M$night_weekday <- ng(cells[night == TRUE & daytype == "monthu"])
  M$night_weekend <- ng(cells[night == TRUE & daytype %in% c("sat", "sunhol") & sched_hour <= 4])
  nd <- merge(cells[night == TRUE, .(yn = sum(B + A), dn = sum(dp)), by = .(stop_name)], cells[night == FALSE & daytype == "monthu", .(yd = sum(B + A), dd = sum(dp)), by = .(stop_name)], by = "stop_name")
  M$night_vs_day <- nd[dn >= 50, .(stop_name, index = log((yn + 0.5) / dn) - log((yd + 0.5) / dd))]
  rm(cells); gc()

  ## ── (4) service dependence ──────────────────────────────────────────────
  dm2 <- function(D, cols, f1, f2, it = 10) { for (k in seq_len(it)) for (v in cols) { D[, (v) := get(v) - mean(get(v)), by = f1]; D[, (v) := get(v) - mean(get(v)), by = f2] }; D }
  S <- SD[Cn >= 10 & B > 0 & N >= 20 & daytype != "bridge"]
  S[, `:=`(y = log(B / Cn), x = log(N), sd_ = paste(cx, daytype))]
  S <- dm2(S, c("y", "x"), "sd_", "DATA")
  DEPN <- S[, .(days = .N, sdx = sd(x), b = sum(x * y) / sum(x^2), prec = sum(x^2)), by = .(stop_name = cx)][days >= 60]
  cat("Service dependence: stops with SD of net log departures >= 0.05:", DEPN[sdx >= 0.05, .N], "of", nrow(DEPN), "| median elasticity", round(DEPN[sdx >= 0.05, median(b)], 3), "\n")
  M$dependence <- DEPN[sdx >= 0.05, .(stop_name, index = -b)]                                             # higher = more riders kept when service is cut

  ## ── (5) weather ──────────────────────────────────────────────────────────
  W <- fread(file.path(out_dir, "data/raw/weather/s_d_330_2025.csv"), header = FALSE, encoding = "Latin-1")
  ## status "9" = phenomenon absent (value empty -> 0), "8" = not measured (missing)
  W <- W[, .(DATA = as.IDate(sprintf("%s-%02d-%02d", V3, as.integer(V4), as.integer(V5))), tmax = V6, tmin = V8, tmean = V10,
             prcp = fifelse(is.na(V14) & as.character(V15) == "9", 0, as.numeric(V14)), snow = fifelse(is.na(V17) & as.character(V18) == "9", 0, as.numeric(V17)))]
  cat("Weather: missing precipitation days", W[is.na(prcp), .N], "| missing snow days", W[is.na(snow), .N], "\n")
  W[is.na(prcp), prcp := 0][is.na(snow), snow := 0]
  W[, `:=`(rain = as.integer(prcp >= 1), heavy = as.integer(prcp >= 10), frost = as.integer(tmean < 0), snowc = as.integer(snow > 0), heat = as.integer(tmax >= 30))]
  cat("Weather days Jan-Oct: rain", W[DATA <= as.IDate("2025-10-31"), sum(rain)], "| heavy", W[DATA <= as.IDate("2025-10-31"), sum(heavy)], "| frost",
      W[DATA <= as.IDate("2025-10-31"), sum(frost)], "| snow cover", W[DATA <= as.IDate("2025-10-31"), sum(snowc)], "| heat", W[DATA <= as.IDate("2025-10-31"), sum(heat)], "\n")
  S2 <- merge(SD[Cn >= 10 & B > 0 & N >= 20 & daytype != "bridge"], W, by = "DATA")
  S2[, `:=`(y = log(B / Cn), fe = paste(cx, daytype, month, season))]
  WX <- c("rain", "heavy", "frost", "snowc", "heat")
  S2[, (WX) := lapply(.SD, as.numeric), .SDcols = WX]                                                    # integer columns would truncate the demeaned values
  for (v in c("y", WX)) S2[, (v) := as.numeric(get(v)) - mean(as.numeric(get(v))), by = fe]
  WS <- S2[, { X <- as.matrix(.SD[, WX, with = FALSE]); keep <- colSums(X^2) > 0
               keep[is.na(keep)] <- FALSE
               cf <- rep(NA_real_, length(WX)); if (sum(keep) > 0) cf[keep] <- tryCatch(qr.coef(qr(X[, keep, drop = FALSE]), y), error = function(e) rep(NA_real_, sum(keep)))
               as.list(setNames(cf, WX)) }, by = .(stop_name = cx)]
  cat("City-wide mean weather sensitivities (log riders per departure):\n"); print(WS[, lapply(.SD, median, na.rm = TRUE), .SDcols = WX], digits = 3)
  for (v in WX) M[[paste0("weather_", v)]] <- WS[is.finite(get(v)), .(stop_name, index = get(v))]

  ## ── stop-level data and the deprivation model (as script 43) ────────────
  cvn <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[, .(stop_name, dist_km, flagged, residents, dep, osiedle_dom = osiedle, share_0_18, share_61plus, spp_share)]
  td  <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  loc <- data.table(stop_name = C$stop_name, os_loc = osi$name[st_nearest_feature(C, osi)])
  d <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(loc, cvn, td[, .(stop_name, fsi, med_gap, med_setback, med_block_area_ha, intersection_dens,
               deadend_share, wc_tree, b1 = ilr_res_vs_nonres, b2 = ilr_multi_vs_single, b3 = ilr_serv_vs_prod, b4 = ilr_retail_vs_work, b5 = ilr_office_vs_eduhealth,
               b6 = ilr_industry_vs_other, landuse_mix, share_university, share_dormitory, area_km2)]))
  d[, osiedle := fifelse(!is.na(osiedle_dom) & osiedle_dom != "", osiedle_dom, os_loc)]
  d <- d[residents > 0]
  d[, `:=`(lres = log1p(as.numeric(residents)), larea = log(area_km2), fsi = log(fsi + 0.01), med_gap = log1p(med_gap), med_setback = log1p(med_setback),
           med_block_area_ha = log(med_block_area_ha), intersection_dens = sqrt(intersection_dens), share_university = sqrt(share_university), share_dormitory = sqrt(share_dormitory))]
  LU <- c("b1", "b2", "b3", "b4", "b5", "b6", "landuse_mix", "fsi", "med_gap", "med_setback", "med_block_area_ha", "intersection_dens", "deadend_share", "wc_tree", "share_university", "share_dormitory")
  d <- d[complete.cases(d[, c(LU, "dist_km", "larea", "dep", "share_0_18", "share_61plus"), with = FALSE])]
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  for (v in c(LU, "dep", "spp_share", "share_0_18", "share_61plus", "lres", "larea")) d[, (v) := zs(get(v))]
  pca <- prcomp(d[, .(fsi, med_gap, med_setback, med_block_area_ha, intersection_dens, deadend_share)])
  npc <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= 0.8)[1]; PCS <- paste0("PC", seq_len(npc)); d[, (PCS) := as.data.table(pca$x[, seq_len(npc), drop = FALSE])]
  DEPADJ <- c("b1", "landuse_mix", "wc_tree", "b2", "b3", "b4", "b5", "b6", PCS, "share_university", "share_dormitory", "spp_share")
  COV <- c("share_0_18", "share_61plus", "lres", "larea", "ns(dist_km, 3)")
  cr2 <- function(fit, cl) {
    X <- model.matrix(fit); e <- resid(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]
    B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) { i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); Am <- ev$vectors %*% diag(1 / sqrt(pmax(ev$values, 1e-10)), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (Am %*% e[i]); meat <- meat + u %*% t(u) }
    setNames(sqrt(diag(B %*% meat %*% B)), colnames(X))
  }
  fitdep <- function(dd, adj) {
    f <- lm(reformulate(c("dep", adj, COV), "yy"), dd); used <- as.integer(rownames(model.matrix(f)))
    se <- cr2(f, dd$osiedle[used])[["dep"]]; G <- uniqueN(dd$osiedle[used])
    data.table(est = coef(f)[["dep"]], se = se, p = 2 * pt(-abs(coef(f)[["dep"]] / se), G - 2), n = length(used), osiedla = G)
  }
  d[, third := cut(frank(dep, ties.method = "first"), quantile(seq_len(.N), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("least", "middle", "most"))]
  SM <- rbindlist(lapply(names(M), function(k) M[[k]][, .(stop_name, measure = k, value = index)]))
  fwrite(SM, file.path(s2, "behaviour_stop_measures.csv"))
  TS <- rbindlist(lapply(names(M), function(k) {
    dd <- merge(d, M[[k]], by = "stop_name")[is.finite(index) & flagged == FALSE]
    if (nrow(dd) < 50) return(NULL)
    dd[, yy := zs(index)]
    rbind(fitdep(dd, DEPADJ)[, model := "deprivation path (land use adjusted)"], fitdep(dd, character(0))[, model := "covariates only"])[,
          `:=`(measure = k, least = dd[third == "least", median(index)], most = dd[third == "most", median(index)])]
  }))
  fwrite(TS, file.path(s2, "behaviour_tests.csv"))
  cat("\nBehavioural measures: slope on catchment deprivation (SD per SD; CR2 by osiedle, t(G-2)); medians of the raw measure by deprivation third\n")
  print(dcast(TS[, .(measure, model, s = sprintf("%5.2f (p=%.3f)", est, p), n, least = signif(least, 3), most = signif(most, 3))], measure + n + least + most ~ model, value.var = "s"),
        width = 200, row.names = FALSE)

  ## ── figure ───────────────────────────────────────────────────────────────
  LAB <- c(bridge = "Bridge days: boardings vs regular days", bridge_am = "Bridge days: boardings 06-09", longweekend = "28-30 April: boardings vs regular days",
           longweekend_am = "28-30 April: boardings 06-09", irregularity = "Day-to-day irregularity (excess SD)", night_weekday = "Night lines, weeknights (same-trip)",
           night_weekend = "Night lines, Friday and Saturday nights", night_vs_day = "Night relative to day volume", dependence = "Riders kept when service is cut",
           weather_rain = "Rain (>= 1 mm)", weather_heavy = "Heavy rain (>= 10 mm)", weather_frost = "Frost (mean < 0 C)", weather_snowc = "Snow cover",
           weather_heat = "Heat (max >= 30 C)")
  FG <- TS[model == "deprivation path (land use adjusted)"][, lab := factor(LAB[measure], levels = rev(LAB))]
  g <- ggplot(FG, aes(est, lab)) + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), colour = "#e34948", size = 0.3, linewidth = 0.5) +
    labs(x = "Slope on catchment deprivation (SD of the stop measure per SD; 95% interval, CR2)", y = NULL) +
    theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, "Fig_behaviour.png"), g, width = 8.5, height = 5.2, dpi = 220, bg = "#ffffff")
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 61 · Tram and bus ──────────────────────────────────────────────────────────────────────────────────────
section("61", "Tram and bus", {
  ## Do the differences between more and less deprived areas differ between tram
  ## and bus lines? (Exploratory.)
  ## (1) Mode reliance: at stops served by both modes, the logit of the tram
  ##     share of boardings minus the logit of the tram share of departures
  ##     (riders' use of trams beyond what the timetable offers); term-time
  ##     Monday-Thursday, day lines, expanded counts (script 41).
  ## (2) Mode-specific usage indices: the six same-trip indices of script 42
  ##     (volume, tidal, onset, concentration, Saturday and summer, boardings
  ##     for the timing measures) computed separately on tram and on bus lines.
  ## For each: medians in the least and most deprived thirds, slopes on catchment
  ## deprivation raw, net of position (distance spline, sector, log lines, log
  ## departures, share of departures on lines that pass the centre) and net of
  ## land use and population too; CR2 by osiedle, t with G - 2 df. A stacked
  ## model tests whether the slope differs between tram and bus lines (mode x
  ## deprivation interaction, covariates by mode).
  ## Outputs: output/stage2/bus_tram_slopes.csv, bus_tram_reliance.csv

  suppressPackageStartupMessages({library(data.table); library(splines); library(sf)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  MON <- c(3:6, 10)
  cr2 <- function(fit, cl) {
    X <- model.matrix(fit); e <- resid(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]
    B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) { i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); Am <- ev$vectors %*% diag(ifelse(ev$values > 1e-8, 1 / sqrt(pmax(ev$values, 1e-8)), 0), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (Am %*% e[i]); meat <- meat + u %*% t(u) }
    setNames(sqrt(diag(B %*% meat %*% B)), colnames(X))
  }
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)

  ## ── stops and covariates (as in scripts 56 and 59) ──────────────────────
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  td <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")[, .(stop_name, b1 = ilr_res_vs_nonres, b2 = ilr_multi_vs_single,
          b3 = ilr_serv_vs_prod, b4 = ilr_retail_vs_work, landuse_mix, fsi = log(fsi + 0.01), wc_tree)]
  d <- merge(cv[residents > 0 & flagged == FALSE, .(stop_name, dep, dist_km, osiedle, share_0_18, share_61plus, lres = log1p(as.numeric(residents)), larea = log(catchment_area_m2))], td, by = "stop_name")
  xy <- fread(file.path(s2, "stop_complex_attributes.csv"), encoding = "UTF-8")[, .(stop_name, x, y)]
  ctr <- st_coordinates(st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180))
  xy[, th := atan2(y - ctr[2], x - ctr[1])][, `:=`(sec_s1 = sin(th), sec_c1 = cos(th), sec_s2 = sin(2 * th), sec_c2 = cos(2 * th))]
  d <- merge(d, xy[, .(stop_name, sec_s1, sec_c1, sec_s2, sec_c2)], by = "stop_name")
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))
  cells[, dp := departures - uncovered]; cells <- cells[dp > 0]
  T0 <- cells[daytype == "monthu" & month %in% MON & night == FALSE]
  sv <- T0[, .(l_lines = log(uniqueN(NR_LINI)), l_dep = log(sum(departures))), by = .(stop_name)]
  segdep <- T0[, .(dep = sum(departures)), by = .(o = stop_name, NR_LINI, dir)]
  od <- rbindlist(lapply(sprintf("%02d", MON), function(m) readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m)))[daytype == "monthu" & o != "outside", .(n = .N), by = .(o, d, NR_LINI, dir)]))
  od <- merge(od, cv[, .(d = stop_name, centre_km = dist_km)], by = "d", all.x = TRUE)
  reach <- merge(od[, .(reach = any(centre_km <= 1.5, na.rm = TRUE)), by = .(o, NR_LINI, dir)], segdep, by = c("o", "NR_LINI", "dir"))[, .(direct_centre = sum(dep * reach) / sum(dep)), by = .(stop_name = o)]
  d <- merge(merge(d, sv, by = "stop_name"), reach, by = "stop_name")
  d <- d[complete.cases(d)]
  q3 <- quantile(d$dep, c(1/3, 2/3)); d[, third := fcase(dep <= q3[1], "least", dep >= q3[2], "most", default = "middle")]
  LOC <- c("ns(dist_km, 3)", "sec_s1", "sec_c1", "sec_s2", "sec_c2"); SERV <- c("l_lines", "l_dep", "direct_centre")
  LUSE <- c("b1", "b2", "b3", "b4", "landuse_mix", "fsi", "wc_tree"); POPV <- c("share_0_18", "share_61plus", "lres", "larea")
  STEPS <- list(raw = character(0), position = c(LOC, SERV), full = c(LOC, SERV, LUSE, POPV))
  slope <- function(Z, y, rhs) { Z <- Z[is.finite(get(y))]; Z[, `:=`(yy = zs(get(y)), dz = zs(dep))]; G <- uniqueN(Z$osiedle)
    f <- lm(reformulate(c("dz", rhs), "yy"), Z); b <- coef(f)[["dz"]]; se <- cr2(f, Z$osiedle)[["dz"]]
    data.table(est = b, se = se, p = 2 * pt(-abs(b / se), G - 2), n = nrow(Z), osiedla = G) }

  ## ── (1) mode reliance ────────────────────────────────────────────────────
  mr <- T0[, .(Bt = sum(B[mode == "tram"]), Bb = sum(B[mode == "bus"]), Dt = sum(dp[mode == "tram"]), Db = sum(dp[mode == "bus"])), by = stop_name][Dt > 0 & Db > 0]
  mr[, `:=`(tram_board_share = Bt / (Bt + Bb), tram_dep_share = Dt / (Dt + Db))]
  mr[, reliance := qlogis((Bt + 0.5) / (Bt + Bb + 1)) - qlogis((Dt + 0.5) / (Dt + Db + 1))]
  R1 <- merge(d, mr, by = "stop_name")
  cat(sprintf("(1) Stops served by both modes (with residents, not flagged): %d in %d osiedla\n", nrow(R1), uniqueN(R1$osiedle)))
  print(R1[third != "middle", .(stops = .N, tram_dep_share = round(median(tram_dep_share), 2), tram_board_share = round(median(tram_board_share), 2), reliance = round(median(reliance), 2)), by = third])
  REL <- rbindlist(lapply(names(STEPS), function(s) slope(R1, "reliance", STEPS[[s]])[, step := s]))
  cat("Slope of tram reliance on deprivation:\n"); print(REL[, .(step, est = round(est, 2), p = round(p, 3), n)])
  fwrite(REL, file.path(s2, "bus_tram_reliance.csv"))

  ## ── (2) mode-specific same-trip indices (as script 42) ──────────────────
  index <- function(dt, gcols) { dt <- dt[is.finite(q) & w > 0]; dt[, n_in := .N, by = gcols]; dt <- dt[n_in >= 2]
    dt[, dev := q - weighted.mean(q, w), by = gcols]; dt[, .(index = weighted.mean(dev, w)), by = stop_name] }
  lg <- function(k, n) qlogis((k + 0.5) / (n + 1))
  build <- function(C) {
    C <- C[season == "term"]; Tm <- C[daytype == "monthu"]; Dm <- Tm[night == FALSE]; out <- list()
    v <- Tm[, .(y = sum(B + A), dep = sum(dp)), by = .(NR_LINI, dir, period, stop_name)][, `:=`(q = log((y + 0.5) / dep), w = dep)]
    out$volume <- index(v, c("NR_LINI", "dir", "period"))
    am <- Dm[sched_hour %in% 6:8, .(B = sum(B), A = sum(A)), by = .(NR_LINI, dir, stop_name)][, `:=`(q = lg(B, B + A), w = B + A)]
    pm <- Dm[sched_hour %in% 15:17, .(B = sum(B), A = sum(A)), by = .(NR_LINI, dir, stop_name)][, `:=`(q = lg(B, B + A), w = B + A)]
    ti <- merge(index(am, c("NR_LINI", "dir")), index(pm, c("NR_LINI", "dir")), by = "stop_name", suffixes = c("_am", "_pm"))
    out$tidal <- ti[, .(stop_name, index = index_am - index_pm)]
    early <- Dm[sched_hour >= 3 & sched_hour < 6, .(dep_early = sum(dp)), by = .(NR_LINI, dir)]
    on <- merge(Dm[, .(Be = sum(B[sched_hour >= 3 & sched_hour < 6]), Ball = sum(B)), by = .(NR_LINI, dir, stop_name)], early[dep_early > 0], by = c("NR_LINI", "dir"))
    out$onset_b <- index(on[, `:=`(q = lg(Be, Ball), w = Ball)], c("NR_LINI", "dir"))
    cc <- Dm[, .(Bp = sum(B[sched_hour %in% 7:8]), dpp = sum(dp[sched_hour %in% 7:8]), Bm = sum(B[sched_hour %in% 10:13]), dm = sum(dp[sched_hour %in% 10:13])), by = .(NR_LINI, dir, stop_name)][dpp > 0 & dm > 0]
    out$conc_b <- index(cc[, `:=`(q = log((Bp + 0.5) / dpp) - log((Bm + 0.5) / dm), w = dpp + dm)], c("NR_LINI", "dir"))
    sa <- merge(C[daytype == "sat" & night == FALSE, .(Bs = sum(B), ds = sum(dp)), by = .(NR_LINI, dir, stop_name)], Dm[, .(Bw = sum(B), dw = sum(dp)), by = .(NR_LINI, dir, stop_name)], by = c("NR_LINI", "dir", "stop_name"))
    out$sat_b <- index(sa[, `:=`(q = log((Bs + 0.5) / ds) - log((Bw + 0.5) / dw), w = ds + dw)], c("NR_LINI", "dir"))
    out
  }
  summer <- function(C) {
    a <- C[season == "summer" & daytype == "monthu", .(y1 = sum(B + A), d1 = sum(dp)), by = .(NR_LINI, dir, stop_name)]
    b <- C[season == "term" & daytype == "monthu", .(y0 = sum(B + A), d0 = sum(dp)), by = .(NR_LINI, dir, stop_name)]
    index(merge(a, b, by = c("NR_LINI", "dir", "stop_name"))[, `:=`(q = log((y1 + 0.5) / d1) - log((y0 + 0.5) / d0), w = d1 + d0)], c("NR_LINI", "dir"))
  }
  IX <- rbindlist(lapply(c("tram", "bus"), function(m) { C <- cells[mode == m]; r <- build(C); r$summer <- summer(C)
    rbindlist(lapply(names(r), function(k) r[[k]][, .(stop_name, measure = k, index)]))[, mode := m] }))
  IXw <- dcast(IX, stop_name + mode ~ measure, value.var = "index")
  D2 <- merge(d, IXw, by = "stop_name")
  MEAS <- c("volume", "tidal", "onset_b", "conc_b", "sat_b", "summer")
  S <- rbindlist(lapply(c("tram", "bus"), function(m) rbindlist(lapply(MEAS, function(y) rbindlist(lapply(names(STEPS), function(s)
    slope(D2[mode == m], y, STEPS[[s]])[, `:=`(step = s, measure = y)]))))[, mode := m]))
  ## difference between modes: stacked, measure standardised within mode, covariates by mode, interaction tested
  DIFF <- rbindlist(lapply(MEAS, function(y) rbindlist(lapply(names(STEPS), function(s) {
    Z <- D2[is.finite(get(y))]; Z[, `:=`(yy = zs(get(y)), dz = zs(dep)), by = mode]; Z[, tram := as.numeric(mode == "tram")]
    rhs <- if (length(STEPS[[s]])) paste0("mode:(", paste(STEPS[[s]], collapse = " + "), ")") else character(0)
    f <- lm(reformulate(c("mode", "dz", "dz:tram", rhs), "yy"), Z); b <- coef(f)[["dz:tram"]]; se <- cr2(f, Z$osiedle)[["dz:tram"]]; G <- uniqueN(Z$osiedle)
    data.table(measure = y, step = s, diff_tram_minus_bus = b, se_diff = se, p_diff = 2 * pt(-abs(b / se), G - 2)) }))))
  S <- merge(S, DIFF, by = c("measure", "step"), all.x = TRUE)
  med <- D2[third != "middle", lapply(.SD, median, na.rm = TRUE), by = .(mode, third), .SDcols = MEAS]
  fwrite(S, file.path(s2, "bus_tram_slopes.csv"))
  cat("\n(2) Stops with mode-specific indices: tram", D2[mode == "tram", uniqueN(stop_name)], "| bus", D2[mode == "bus", uniqueN(stop_name)], "\n")
  print(dcast(S[, .(measure, mode, step, v = sprintf("%5.2f%s", est, ifelse(p < 0.01, "**", ifelse(p < 0.05, "* ", "  "))))], measure + mode ~ step, value.var = "v")[order(measure, mode)])
  cat("\nTram minus bus slope (stacked; covariates by mode):\n")
  print(dcast(DIFF[, .(measure, step, v = sprintf("%5.2f (p=%.3f)", diff_tram_minus_bus, p_diff))], measure ~ step, value.var = "v"))
  cat("\nMedians by mode and third:\n"); print(med, digits = 2)
})

