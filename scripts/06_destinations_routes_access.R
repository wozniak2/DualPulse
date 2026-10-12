## DualPulse | stage 6: Where riders are carried, along which routes, and what the network offers
##
## Destinations against a network baseline, interchanges and transfers, route heterogeneity, route
## networks and timetable-based access (main text 3.3-3.4, supplementary S7-S8).
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   45      Destination profiles and alighting density
##   46      Destinations: network, timing and the lines riders use
##   47      Interchanges and transfers
##   50      Destination and interchange maps, and their reasons
##   51      Route heterogeneity
##   53      Route networks traced along the lines
##   54      Access by the timetable
##
## Run from the repository root: Rscript run.R 6   (one section: Rscript run.R 6:54)

if (!exists("section")) source("R/setup.R")

## ── 45 · Destination profiles and alighting density ────────────────────────────────────────────────────────
section("45", "Destination profiles and alighting density", {
  ## Where are riders from more and less deprived catchments carried?
  ## (Exploratory; flows from script 40, January-October 2025.) Replaces the
  ## former scripts 45 (destination profiles) and 46 (density maps).
  ##
  ## For every origin complex, the boardings of term-time Monday-Thursday are
  ## followed to their estimated alighting stops on the same trip (proportional
  ## downstream assignment, script 40). Flow-weighted destination profile per
  ## origin, all day and for the morning peak (06-09):
  ##   trip_km        straight-line distance from origin to destination complex
  ##   dest_centre_km destination distance to the centre (Stary Rynek)
  ##   to_centre      share alighting within 1.5 km of the centre
  ##   to_outside     share alighting beyond the city
  ##   same_osiedle   share alighting in the origin's osiedle
  ##   dest_b1        residential vs non-residential floor at the destination (low = workplaces)
  ##   dest_jobs      jobs per resident at the destination (modelled workplace grid)
  ##   dest_industry  industry and storage floor share at the destination
  ##   dest_dep       catchment deprivation at the destination
  ##   dest_muni      municipal-housing share at the destination
  ## Origins: complexes with residents, not terminals, interchanges or boundary
  ## complexes, at least 2,000 assigned term-time boardings (500 in the morning
  ## peak; walking catchments).
  ## (1) Comparisons: medians by thirds of origin deprivation; slope of each
  ## profile measure on origin deprivation (SD units) with base controls (origin
  ## distance spline, age shares, log residents, log catchment area), + origin
  ## land use, and within osiedla; cluster bootstrap over osiedla (999).
  ## (2) Density maps: each destination complex carries the morning-peak
  ## alightings it receives from the most and the least deprived third of
  ## origins, smoothed with a Gaussian kernel (bandwidth 500 m; 300 and 800 m as
  ## checks) onto a 100 m grid clipped to the city; each surface sums to 100% of
  ## its group's alightings.
  ## Outputs: output/stage2/dest_profiles.csv, dest_tests.csv;
  ##          output/figures/Fig_destination_density.png, Fig_destination_density_bandwidths.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines); library(terra); library(ggplot2); library(patchwork)})
  set.seed(20261007)
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R"))   # flows expanded from counted to all operated trips
  s2 <- file.path(out_dir, "output/stage2")
  od <- rbindlist(lapply(sprintf("%02d", 1:10), function(m) expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))))))
  ## term time = March to 27 June and October without the Easter break; daily periods as in script 41
  od <- od[daytype == "monthu" & month %in% c(3:6, 10) & o != "outside"]
  cat("Term-time Monday-Thursday boardings assigned to destinations:", round(sum(od$flow)), "\n")

  C   <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  xy  <- data.table(stop_name = C$stop_name, X = st_coordinates(C)[, 1], Y = st_coordinates(C)[, 2])
  cv  <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  td  <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")[, .(stop_name, b1 = ilr_res_vs_nonres, landuse_mix, wc_tree)]
  mun <- fread(file.path(s2, "municipal_catchments_network.csv"), encoding = "UTF-8")[, .(stop_name, muni_share)]
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  A <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(xy, cv[, .(stop_name, dist_km, dep, jobs_res, share_industry, osiedle_dom = osiedle)], td, mun))
  A[, osiedle := fifelse(!is.na(osiedle_dom), osiedle_dom, osi$name[st_nearest_feature(C, osi)][match(stop_name, C$stop_name)])]

  ## ── (1) destination profiles ────────────────────────────────────────────
  dest <- A[, .(d = stop_name, dX = X, dY = Y, d_centre = dist_km, d_b1 = b1, d_jobs = jobs_res, d_ind = share_industry, d_dep = dep, d_muni = muni_share, d_os = osiedle)]
  f <- merge(od, dest, by = "d", all.x = TRUE)
  f <- merge(f, A[, .(o = stop_name, oX = X, oY = Y, o_os = osiedle)], by = "o")
  f[, outside := d == "outside"]
  f[, trip_km := fifelse(outside, NA_real_, sqrt((dX - oX)^2 + (dY - oY)^2) / 1000)]
  prof <- function(ff) ff[, .(boardings = sum(flow_c),   # counted support, for the threshold
    trip_km = weighted.mean(trip_km, flow, na.rm = TRUE), dest_centre_km = weighted.mean(d_centre, flow, na.rm = TRUE),
    to_centre = sum(flow[!outside & d_centre <= 1.5]) / sum(flow), to_outside = sum(flow[outside]) / sum(flow),
    same_osiedle = sum(flow[!outside & d_os == o_os]) / sum(flow[!outside]),
    dest_b1 = weighted.mean(d_b1, flow, na.rm = TRUE), dest_jobs = weighted.mean(d_jobs, flow, na.rm = TRUE),
    dest_industry = weighted.mean(d_ind, flow, na.rm = TRUE), dest_dep = weighted.mean(d_dep, flow, na.rm = TRUE),
    dest_muni = weighted.mean(d_muni, flow, na.rm = TRUE)), by = o]
  pr <- rbind(prof(f)[, window := "all day"], prof(f[period == "am_peak"])[, window := "morning peak"])
  setnames(pr, "o", "stop_name")
  cv[, `:=`(lres = log1p(as.numeric(residents)), larea = log(catchment_area_m2))]
  pr <- merge(pr, cv[, .(stop_name, dep, dist_km, residents, flagged, share_0_18, share_61plus, lres, larea, osiedle)], by = "stop_name")
  pr <- merge(pr, td, by = "stop_name")
  pr <- pr[residents > 0 & flagged == FALSE & boardings >= fifelse(window == "all day", 2000, 500)]
  fwrite(pr, file.path(s2, "dest_profiles.csv"))
  cat("Origins analysed: all day", pr[window == "all day", .N], "| morning peak", pr[window == "morning peak", .N], "\n")

  Y <- c(trip_km = "Trip distance (km)", dest_centre_km = "Destination distance to centre (km)", to_centre = "Share to the centre (<1.5 km)",
         to_outside = "Share leaving the city", same_osiedle = "Share staying in the origin osiedle", dest_b1 = "Destination residential vs non-residential floor",
         dest_jobs = "Destination jobs per resident (log)", dest_industry = "Destination industry floor share", dest_dep = "Destination deprivation",
         dest_muni = "Destination municipal-housing share")
  pr[, tercile := cut(frank(dep, ties.method = "first"), quantile(seq_len(.N), c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                      labels = c("Least deprived third", "Middle third", "Most deprived third")), by = window]
  med <- pr[, lapply(.SD, median, na.rm = TRUE), by = .(window, tercile), .SDcols = names(Y)]
  cat("\nMedian destination profile by origin deprivation third:\n"); print(dcast(melt(med, id.vars = c("window", "tercile")), window + variable ~ tercile), digits = 3)

  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  BASE <- c("share_0_18", "share_61plus", "lres", "larea", "ns(dist_km, 3)")
  sets <- list(base = BASE, land_use = c(BASE, "b1", "landuse_mix", "wc_tree"), within_osiedle = c(BASE, "b1", "landuse_mix", "wc_tree", "factor(osiedle)"))
  fit <- function(D, y, X) {
    D <- copy(D[is.finite(get(y))]); for (v in c(y, "dep", "share_0_18", "share_61plus", "lres", "larea", "b1", "landuse_mix", "wc_tree")) D[, (v) := zs(get(v))]
    tryCatch(coef(lm(reformulate(c("dep", X), y), D))[["dep"]], error = function(e) NA_real_)
  }
  res <- rbindlist(lapply(c("all day", "morning peak"), function(w) {
    dw <- pr[window == w]; ix <- split(seq_len(nrow(dw)), dw$osiedle); on <- names(ix)
    rbindlist(lapply(names(Y), function(y) rbindlist(lapply(names(sets), function(s) {
      e <- fit(dw, y, sets[[s]])
      bt <- replicate(999, { k <- sample(on, replace = TRUE); D <- dw[unlist(ix[k])]
        if (s == "within_osiedle") D[, osiedle := paste(osiedle, rep(seq_along(k), lengths(ix[k])))]; fit(D, y, sets[[s]]) })
      data.table(window = w, measure = Y[[y]], controls = s, est = e, lo = quantile(bt, .025, na.rm = TRUE), hi = quantile(bt, .975, na.rm = TRUE),
                 p = min(1, 2 * min(mean(bt <= 0, na.rm = TRUE), mean(bt >= 0, na.rm = TRUE))))
    }))))
  }))
  res[, q := p.adjust(p, "BH"), by = .(window, controls)]
  fwrite(res, file.path(s2, "dest_tests.csv"))
  res[, s := sprintf("%5.2f [%5.2f,%5.2f]%s", est, lo, hi, fifelse(q < 0.05, " **", fifelse(p < 0.05, " *", "")))]
  cat("\nSlope of the destination profile on origin deprivation (SD units; * p < .05; ** BH q < .05)\n")
  print(dcast(res, window + measure ~ factor(controls, levels = names(sets)), value.var = "s"), width = 250)

  ## ── (2) density maps: morning-peak alightings by origin third ───────────
  am <- merge(f[period == "am_peak" & !outside], pr[window == "morning peak", .(o = stop_name, tercile)], by = "o")
  fl <- am[tercile != "Middle third", .(flow = sum(flow)), by = .(tercile = fifelse(tercile == "Most deprived third", "most", "least"), d, X = dX, Y = dY)]
  cat(sprintf("\nMorning-peak alightings: from the least deprived third %.0f, from the most deprived third %.0f (counted trips, term time)\n",
              fl[tercile == "least", sum(flow)], fl[tercile == "most", sum(flow)]))
  city <- st_union(osi)
  water <- st_read(file.path(out_dir, "output/stage1/tissue_inputs.gpkg"), "green_blue", quiet = TRUE)
  water <- st_transform(water[water$cover == "water", ], st_crs(C))
  r <- rast(ext(st_bbox(city)), resolution = 100, crs = st_crs(C)$wkt)
  g <- as.data.table(xyFromCell(r, seq_len(ncell(r))))
  g <- g[lengths(st_intersects(st_as_sf(g, coords = c("x", "y"), crs = st_crs(C)), city)) > 0]
  kde <- function(pts, bw) {
    w <- pts$flow / sum(pts$flow); v <- numeric(nrow(g))
    for (j in seq_len(nrow(pts))) v <- v + w[j] * exp(-((g$x - pts$X[j])^2 + (g$y - pts$Y[j])^2) / (2 * bw^2))
    100 * v / sum(v) / 0.01                                         # % of the group's alightings per km2 (cells of 0.01 km2)
  }
  S <- rbindlist(lapply(c(300, 500, 800), function(bw)
    copy(g)[, `:=`(least = kde(fl[tercile == "least"], bw), most = kde(fl[tercile == "most"], bw))][, `:=`(diff = most - least, bw = bw)]))

  thm <- theme_void(base_size = 10) + theme(plot.background = element_rect(fill = "#ffffff", colour = NA), plot.title = element_text(face = "bold", size = 10),
    legend.position = "bottom", legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.25, "cm"), legend.title.position = "top")
  ctx <- list(geom_sf(data = osi, fill = NA, colour = "#ffffff", linewidth = 0.25, inherit.aes = FALSE),
              geom_sf(data = water, fill = "#9fc3e6", colour = NA, inherit.aes = FALSE),
              geom_sf(data = city, fill = NA, colour = "#6b6a64", linewidth = 0.3, inherit.aes = FALSE))
  s5 <- S[bw == 500]; lim <- c(0, max(s5$least, s5$most)); m <- max(abs(s5$diff))
  panel <- function(v, title) ggplot(s5) + geom_raster(aes(x, y, fill = .data[[v]])) + ctx + coord_sf(datum = NA) +
    scale_fill_gradient(low = "#f6f4ee", high = "#4a2f96", limits = lim, name = "Share of alightings (% per km²)") + labs(title = title) + thm
  p_diff <- ggplot(s5) + geom_raster(aes(x, y, fill = diff)) + ctx + coord_sf(datum = NA) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f6f4ee", high = "#e34948", midpoint = 0, limits = c(-m, m), name = "Most minus least deprived (pp per km²)") +
    labs(title = "Difference") + thm
  fig <- (panel("least", "Riders from the least deprived third") | panel("most", "Riders from the most deprived third") | p_diff) + plot_annotation(
    title = "Where morning-peak riders alight: origins in the least vs the most deprived stop catchments",
    caption = "Term-time Monday-Thursday 06-09, March-June and October 2025. Destinations estimated by proportional downstream assignment on counted trips; Gaussian kernel, bandwidth 500 m, 100 m grid.",
    theme = theme(plot.title = element_text(face = "bold"), plot.caption = element_text(colour = "#898781", size = 7.5)))
  ggsave(file.path(out_dir, "output/figures/Fig_destination_density.png"), fig, width = 13, height = 5.6, dpi = 220, bg = "#ffffff")
  pb <- ggplot(S) + geom_raster(aes(x, y, fill = diff)) + facet_wrap(~ paste("Bandwidth", bw, "m")) + ctx + coord_sf(datum = NA) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f6f4ee", high = "#e34948", midpoint = 0, name = "Most minus least deprived (pp per km²)") + thm +
    labs(title = "Difference surface at three kernel bandwidths")
  ggsave(file.path(out_dir, "output/figures/Fig_destination_density_bandwidths.png"), pb, width = 13, height = 5.4, dpi = 200, bg = "#ffffff")
  cat("Correlation of the difference surface between bandwidths (300 vs 500, 500 vs 800):",
      round(cor(S[bw == 300, diff], S[bw == 500, diff]), 2), round(cor(S[bw == 500, diff], S[bw == 800, diff]), 2), "\n")
  cat("Saved dest_profiles.csv, dest_tests.csv, Fig_destination_density.png, Fig_destination_density_bandwidths.png\n")
})

## ── 46 · Destinations: network, timing and the lines riders use ────────────────────────────────────────────
section("46", "Destinations: network, timing and the lines riders use", {
  ## Why do riders from more and less deprived catchments end up in different
  ## places, and do they ride further? (Exploratory; flows by line and direction
  ## from script 40; Monday-Thursday, Saturdays, Sundays and holidays of March-June
  ## and October 2025.)
  ##
  ## With proportional downstream assignment, everyone boarding the same trip at
  ## the same stop shares one destination distribution. The only origin-specific
  ## information is how many board each line x direction in each period. Each
  ## origin's destination profile is therefore decomposed as
  ##   observed = network + timing + lines and directions used, where
  ##   network  = destinations if riders used each line x direction x period in
  ##              proportion to its departures at the stop (what the network
  ##              offers from where people live);
  ##   timing   = shift from the observed split of boardings across periods;
  ##   lines and directions used = the rest (which lines and directions riders
  ##              board; includes the stop's position along each line).
  ## The slope of each component on origin deprivation (common adjustments, CR2
  ## by osiedle, t with G - 2 df; everything in SD of the observed measure) adds
  ## up to the observed slope.
  ## Destination measures (flow-weighted over destination catchments):
  ##   ride_km (distance along the route), line_km (straight line), to_outside,
  ##   dest_dep, centre_km, to_centre (<= 1.5 km), same_osiedle; floor shares at
  ##   the destination: office, education and culture, university, health,
  ##   retail and services, industry and storage, housing; floor-area density;
  ##   zoned mixed residential-service, services and services-production land; jobs per resident
  ##   (workplace grid, log); paid-parking share; share alighting at
  ##   interchanges/terminals (where many riders transfer); share alighting in
  ##   catchments with a marketplace, a hospital or clinic, a welfare office
  ##   (social facility or employment agency) or a government office (OSM).
  ## Windows: weekday (Monday-Thursday) all day, early (04-06), morning peak (06-09),
  ## midday (09-15) and evening (19-); Saturday; Sunday and holidays.
  ## Checks: within osiedla (osiedle fixed effects); transfers at hubs (a share
  ## tau = 0.5 of alightings at interchanges and terminals continues as a
  ## typical boarding there); origin-side accessibility to office and university
  ## floor within 2.5 km (do less deprived residents live nearer white-collar
  ## workplaces?).
  ## Outputs: output/stage2/dest_mech_profiles.csv, dest_mech_tests.csv,
  ##          dest_mech_access.csv; output/figures/Fig_destination_mechanisms.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines); library(ggplot2)})
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R"))   # flows expanded from counted to all operated trips
  s2 <- file.path(out_dir, "output/stage2")
  MON <- c(3:6, 10)
  od <- rbindlist(lapply(sprintf("%02d", 1:10), function(m) expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))))))
  od <- od[daytype %in% c("monthu", "sat", "sunhol") & month %in% MON & o != "outside"]
  od <- od[, .(flow = sum(flow), fkm = sum(fkm), flow_c = sum(flow_c)), by = .(o, d, NR_LINI, dir, daytype, period)]
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))[daytype %in% c("monthu", "sat", "sunhol") & month %in% MON]
  dep <- cells[, .(dep = sum(departures)), by = .(o = stop_name, NR_LINI, dir, daytype, period = as.character(period))]
  rm(cells); gc()

  ## ── destination attributes ───────────────────────────────────────────────
  C   <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  xy  <- data.table(stop_name = C$stop_name, X = st_coordinates(C)[, 1], Y = st_coordinates(C)[, 2])
  cv  <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  td  <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")
  mun <- fread(file.path(s2, "municipal_catchments_network.csv"), encoding = "UTF-8")[, .(stop_name, muni_share)]
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  cat_poly <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  am <- readRDS(file.path(out_dir, "output/stage1/_osm_amenities_raw.rds"))
  am <- st_transform(am, st_crs(C))
  am_type <- with(am, fcase(amenity %in% "marketplace", "market",
                            amenity %in% c("hospital", "clinic"), "health_fac",
                            amenity %in% "social_facility" | office %in% "employment_agency", "welfare",
                            office %in% "government" | amenity %in% c("townhall", "courthouse"), "government", default = NA_character_))
  am <- am[!is.na(am_type), ]; am$type <- am_type[!is.na(am_type)]
  hit <- st_intersects(cat_poly, am)
  amc <- rbindlist(lapply(seq_along(hit), function(i) data.table(stop_name = cat_poly$stop_name[i], type = am$type[hit[[i]]])), fill = TRUE)[!is.na(type)]
  amc <- dcast(amc[, .N, by = .(stop_name, type)], stop_name ~ type, value.var = "N", fill = 0)
  cat("Catchments with a marketplace / hospital-clinic / welfare office / government office:",
      paste(sapply(c("market", "health_fac", "welfare", "government"), function(k) sum(amc[[k]] > 0)), collapse = " / "), "\n")

  A <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(xy,
    cv[, .(stop_name, dist_km, dep, jobs_res, spp_share, hub = terminal %in% TRUE | interchange %in% TRUE, osiedle_cv = osiedle)],
    td[, .(stop_name, office = share_office, educ = share_education_culture, univ = share_university, health = share_health,
           retail = share_retail_services, industry = share_industry_storage, housing = share_res_multi + share_res_single,
           floor_total, fsi, z_mixed, z_services, z_production, b1 = ilr_res_vs_nonres, landuse_mix, wc_tree)], mun, amc))
  for (k in c("market", "health_fac", "welfare", "government")) A[is.na(get(k)), (k) := 0]
  A[, osiedle := fifelse(!is.na(osiedle_cv) & osiedle_cv != "", osiedle_cv, osi$name[st_nearest_feature(C, osi)][match(stop_name, C$stop_name)])]
  DM <- c("ride_km", "line_km", "to_outside", "dest_dep", "centre_km", "to_centre", "same_osiedle", "office", "educ", "univ", "health", "retail",
          "industry", "housing", "fsi", "z_mixed", "z_services", "z_production", "jobs_res", "spp_share", "hub", "market", "health_fac", "welfare", "government")
  LAB <- c(ride_km = "Ride length along the route (km)", line_km = "Straight-line ride length (km)", to_outside = "Share leaving the city",
           dest_dep = "Destination deprivation", centre_km = "Destination distance to centre (km)", to_centre = "Share to the centre (<1.5 km)",
           same_osiedle = "Share staying in the origin osiedle", office = "Destination office floor share", educ = "Destination education & culture floor share",
           univ = "Destination university floor share", health = "Destination health floor share", retail = "Destination retail & services floor share",
           industry = "Destination industry & storage floor share", housing = "Destination housing floor share", fsi = "Destination floor-area density", z_mixed = "Destination zoned mixed residential-service land", z_services = "Destination zoned services land", z_production = "Destination zoned services-production land", jobs_res = "Destination jobs per resident (log)",
           spp_share = "Destination paid-parking share", hub = "Share alighting at interchanges/terminals", market = "Share alighting near a marketplace",
           health_fac = "Share alighting near a hospital/clinic", welfare = "Share alighting near a welfare office", government = "Share alighting near a government office")

  ## per-pair destination values (origin-specific: line_km, same_osiedle; others from the destination)
  D <- A[, .(d = stop_name, dX = X, dY = Y, centre_km = dist_km, dest_dep = dep, office, educ, univ, health, retail, industry, housing, fsi, z_mixed, z_services, z_production, jobs_res, spp_share,
             hub = as.numeric(hub), market = as.numeric(market > 0), health_fac = as.numeric(health_fac > 0), welfare = as.numeric(welfare > 0),
             government = as.numeric(government > 0), d_os = osiedle)]
  f <- merge(od, D, by = "d", all.x = TRUE)
  f <- merge(f, A[, .(o = stop_name, oX = X, oY = Y, o_os = osiedle)], by = "o")
  f[, outside := d == "outside"]
  f[, `:=`(ride_km = fkm / flow, line_km = fifelse(outside, NA_real_, sqrt((dX - oX)^2 + (dY - oY)^2) / 1000), to_outside = as.numeric(outside),
           to_centre = fifelse(outside, NA_real_, as.numeric(centre_km <= 1.5)), same_osiedle = fifelse(outside, NA_real_, as.numeric(d_os == o_os)))]

  ## ── transfer check: tau of alightings at hubs continue as typical boardings there ──
  tau <- 0.5
  onward <- f[, c(lapply(.SD, function(v) weighted.mean(v, flow, na.rm = TRUE)), list(w_on = sum(flow))), by = .(o, daytype, period), .SDcols = setdiff(DM, c("line_km", "same_osiedle"))]
  setnames(onward, c("o", setdiff(DM, c("line_km", "same_osiedle"))), c("d", paste0("on_", setdiff(DM, c("line_km", "same_osiedle")))))
  ft <- merge(f, onward[, !"w_on"], by = c("d", "daytype", "period"), all.x = TRUE)
  for (k in setdiff(DM, c("line_km", "same_osiedle", "hub"))) {
    if (k == "ride_km") ft[hub == 1 & !is.na(on_ride_km), ride_km := ride_km + tau * on_ride_km]
    else ft[hub == 1 & !is.na(get(paste0("on_", k))), (k) := (1 - tau) * get(k) + tau * get(paste0("on_", k))]
  }

  ## ── origin profiles: observed, network (departures) and timing-kept ─────
  ALLP <- c("night", "early", "am_peak", "midday", "pm_peak", "evening")
  WIN <- list("weekday all day" = list("monthu", ALLP), "weekday early (04-06)" = list("monthu", "early"),
              "weekday morning peak (06-09)" = list("monthu", "am_peak"), "weekday midday (09-15)" = list("monthu", "midday"),
              "weekday evening (19-)" = list("monthu", "evening"), "Saturday" = list("sat", ALLP), "Sunday and holidays" = list("sunhol", ALLP))
  profiles <- function(ff, w, label) {
    g <- ff[daytype %in% WIN[[w]][[1]] & period %in% WIN[[w]][[2]]]
    ## segment means of each measure
    seg <- g[, c(list(B = sum(flow), Bc = sum(flow_c)), lapply(.SD, function(v) { ok <- !is.na(v); if (any(ok)) sum(v[ok] * flow[ok]) / sum(flow[ok]) else NA_real_ })),
             by = .(o, NR_LINI, dir, daytype, period), .SDcols = DM]
    seg <- merge(seg, dep, by = c("o", "NR_LINI", "dir", "daytype", "period"), all.x = TRUE)[!is.na(dep) & dep > 0]
    seg[, Bp := sum(B), by = .(o, daytype, period)][, wB := dep / sum(dep) * Bp, by = .(o, daytype, period)]   # timing kept, line/direction by departures
    wm <- function(v, w) { ok <- !is.na(v); if (any(ok)) sum(v[ok] * w[ok]) / sum(w[ok]) else NA_real_ }
    ## boardings = counted support (for the threshold); profiles weight the expanded flows
    rbind(seg[, c(list(boardings = sum(Bc)), lapply(.SD, wm, w = B)), by = o, .SDcols = DM][, comp := "observed"],
          seg[, c(list(boardings = sum(Bc)), lapply(.SD, wm, w = dep)), by = o, .SDcols = DM][, comp := "network"],
          seg[, c(list(boardings = sum(Bc)), lapply(.SD, wm, w = wB)), by = o, .SDcols = DM][, comp := "timing_kept"])[, `:=`(window = w, transfers = label)]
  }
  PR <- rbindlist(c(lapply(names(WIN), function(w) profiles(f, w, "none")), lapply(names(WIN), function(w) profiles(ft, w, "tau = 0.5 at hubs"))))
  setnames(PR, "o", "stop_name")
  cv[, `:=`(lres = log1p(as.numeric(residents)), larea = log(catchment_area_m2))]
  PR <- merge(PR, cv[, .(stop_name, dep, dist_km, residents, flagged, share_0_18, share_61plus, lres, larea)], by = "stop_name")
  PR <- merge(PR, A[, .(stop_name, osiedle, b1, landuse_mix, wc_tree)], by = "stop_name")
  MINB <- c("weekday all day" = 2000, "weekday early (04-06)" = 100, "weekday morning peak (06-09)" = 500, "weekday midday (09-15)" = 500,
            "weekday evening (19-)" = 300, "Saturday" = 500, "Sunday and holidays" = 300)
  PR <- PR[residents > 0 & flagged == FALSE & boardings >= MINB[window]]
  fwrite(PR, file.path(s2, vtag("dest_mech_profiles.csv")))
  print(PR[comp == "observed" & transfers == "none", .(origins = .N, osiedla = uniqueN(osiedle), boardings = round(sum(boardings))), by = window])

  ## ── medians by thirds (observed vs network) ─────────────────────────────
  PR[, third := cut(frank(dep, ties.method = "first"), quantile(seq_len(.N), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("least", "middle", "most")),
     by = .(window, comp, transfers)]
  MED <- PR[transfers == "none" & comp %in% c("observed", "network"), lapply(.SD, median, na.rm = TRUE), by = .(window, comp, third), .SDcols = DM]
  fwrite(MED, file.path(s2, vtag("dest_mech_thirds.csv")))

  ## ── slopes on origin deprivation, decomposed ────────────────────────────
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  cr2 <- function(fit, cl) {
    X <- model.matrix(fit); e <- resid(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]
    B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) {
      i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); Am <- ev$vectors %*% diag(1 / sqrt(pmax(ev$values, 1e-10)), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (Am %*% e[i]); meat <- meat + u %*% t(u)
    }
    setNames(sqrt(diag(B %*% meat %*% B)), colnames(X))
  }
  BASE <- c("share_0_18", "share_61plus", "lres", "larea", "ns(dist_km, 3)")
  SETS <- list(base = BASE, land_use = c(BASE, "b1", "landuse_mix", "wc_tree"), within_osiedle = c(BASE, "b1", "landuse_mix", "wc_tree", "factor(osiedle)"))
  slope <- function(dd, y, X) {
    f <- lm(reformulate(c("dep", X), y), dd); used <- as.integer(rownames(model.matrix(f)))
    se <- cr2(f, dd$osiedle[used])[["dep"]]; G <- uniqueN(dd$osiedle[used])
    data.table(est = coef(f)[["dep"]], se = se, p = 2 * pt(-abs(coef(f)[["dep"]] / se), G - 2), n = length(used))
  }
  TS <- rbindlist(lapply(unique(PR$transfers), function(tr) rbindlist(lapply(names(WIN), function(w) rbindlist(lapply(DM, function(y) {
    W <- dcast(PR[window == w & transfers == tr], stop_name + dep + dist_km + share_0_18 + share_61plus + lres + larea + osiedle + b1 + landuse_mix + wc_tree ~ comp,
               value.var = y)
    W <- W[is.finite(observed) & is.finite(network) & is.finite(timing_kept)]
    if (nrow(W) < 30 || sd(W$observed) == 0) return(NULL)
    s <- sd(W$observed)
    W[, `:=`(obs = observed / s, net = network / s, timing = (timing_kept - network) / s, linedir = (observed - timing_kept) / s)]
    for (v in c("dep", "share_0_18", "share_61plus", "lres", "larea", "b1", "landuse_mix", "wc_tree")) W[, (v) := zs(get(v))]
    rbind(rbindlist(lapply(c("obs", "net", "timing", "linedir"), function(cmp) slope(W, cmp, SETS$base)[, `:=`(component = cmp, controls = "base")])),
          slope(W, "obs", SETS$land_use)[, `:=`(component = "obs", controls = "land_use")],
          slope(W, "obs", SETS$within_osiedle)[, `:=`(component = "obs", controls = "within_osiedle")],
          slope(W, "linedir", SETS$within_osiedle)[, `:=`(component = "linedir", controls = "within_osiedle")])[, `:=`(window = w, measure = y, transfers = tr)]
  }))))))
  TS[, q := p.adjust(p, "BH"), by = .(window, component, controls, transfers)]
  fwrite(TS, file.path(s2, vtag("dest_mech_tests.csv")))

  fmt <- function(e, p) sprintf("%5.2f%s", e, fifelse(p < 0.01, "**", fifelse(p < 0.05, "* ", "  ")))
  for (w in names(WIN)) {
    cat("\n==", w, "== slope on origin deprivation (SD of observed measure per SD), base adjustments, CR2; * p<.05 ** p<.01\n")
    X <- TS[window == w & transfers == "none" & controls == "base"]
    X <- dcast(X[, .(measure, component, s = fmt(est, p))], measure ~ factor(component, levels = c("obs", "net", "timing", "linedir")), value.var = "s")
    X <- merge(X, TS[window == w & transfers == "none" & component == "obs" & controls != "base", .(measure, controls, s = fmt(est, p))][
      , dcast(.SD, measure ~ controls, value.var = "s")], by = "measure")
    X <- merge(X, TS[window == w & transfers == "tau = 0.5 at hubs" & component == "obs" & controls == "base", .(measure, tau05 = fmt(est, p))], by = "measure", all.x = TRUE)
    X <- merge(X, MED[window == w & comp == "observed", .(measure = DM, least = unlist(.SD[third == "least", DM, with = FALSE]),
                                                           most = unlist(.SD[third == "most", DM, with = FALSE]))], by = "measure")
    X[, measure := factor(measure, levels = DM)]; setorder(X, measure)
    print(X[, .(measure, obs, net, timing, linedir, land_use, within_osiedle, tau05, least = signif(least, 3), most = signif(most, 3))], row.names = FALSE)
  }

  ## ── origin-side access: do less deprived residents live nearer white-collar floor? ──
  O <- merge(A[, .(stop_name, X, Y, office_m2 = floor_total * office, univ_m2 = floor_total * univ, educ_m2 = floor_total * educ, health_m2 = floor_total * health,
                   industry_m2 = floor_total * industry)], cv[, .(stop_name, dep, dist_km, residents, osiedle_cv = osiedle)], by = "stop_name")
  Dm <- as.matrix(dist(O[, .(X, Y)])) / 1000
  acc <- function(v, r = 2.5) as.vector((Dm <= r) %*% fcoalesce(v, 0))
  O[, `:=`(acc_office = acc(office_m2), acc_univ = acc(univ_m2), acc_educ = acc(educ_m2), acc_health = acc(health_m2), acc_industry = acc(industry_m2))]
  O <- O[residents > 0]
  AC <- rbindlist(lapply(c("acc_office", "acc_univ", "acc_educ", "acc_health", "acc_industry"), function(k) {
    O[, lk := log1p(get(k))]
    data.table(measure = k, rho = cor(O$dep, O$lk, method = "spearman"),
               net_of_distance = coef(lm(zs(lk) ~ zs(dep) + ns(dist_km, 3), O))[2],
               least_third_median = O[dep <= quantile(dep, 1/3), median(get(k))], most_third_median = O[dep >= quantile(dep, 2/3), median(get(k))])
  }))
  fwrite(AC, file.path(s2, vtag("dest_mech_access.csv")))
  cat("\nOrigin-side access: floor area (m2) within 2.5 km of the stop, by deprivation (stops with residents):\n"); print(AC, digits = 3)

  ## ── figure: observed slope and its parts, morning peak and midday ───────
  FW <- c("weekday morning peak (06-09)", "weekday midday (09-15)", "Saturday", "Sunday and holidays")
  FG <- TS[transfers == "none" & controls == "base" & window %in% FW &
             measure %in% c("ride_km", "dest_dep", "same_osiedle", "to_centre", "office", "univ", "health", "retail", "industry", "spp_share", "hub", "market")]
  FG <- FG[component != "timing"]
  FG[, component := factor(component, levels = c("obs", "net", "linedir"),
                           labels = c("Observed", "Network from where riders live", "Lines and directions used"))]
  FG[, measure := factor(LAB[measure], levels = rev(LAB[c("ride_km", "same_osiedle", "to_centre", "dest_dep", "office", "univ", "health", "retail", "industry",
                                                            "spp_share", "market", "hub")]))]
  FG[, window := factor(window, levels = FW)]
  pal <- c("Observed" = "#1f1f1f", "Network from where riders live" = "#2a78d6", "Lines and directions used" = "#e34948")
  g <- ggplot(FG, aes(est, measure, colour = component)) + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), position = position_dodge(width = 0.7), size = 0.25, linewidth = 0.45) +
    facet_wrap(~ window, ncol = 2) + scale_colour_manual(values = pal, name = NULL) +
    labs(x = "Slope on origin deprivation (SD of the observed destination measure per SD; 95% interval, CR2)", y = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#ffffff", colour = NA),
                                          plot.caption = element_text(colour = "#898781", size = 7.5), strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "output/figures", vtag("Fig_destination_mechanisms.png")), g, width = 10.5, height = 10, dpi = 200, bg = "#ffffff")
  cat("\nSaved dest_mech_profiles.csv, dest_mech_thirds.csv, dest_mech_tests.csv, dest_mech_access.csv, Fig_destination_mechanisms.png\n")
})

## ── 47 · Interchanges and transfers ────────────────────────────────────────────────────────────────────────
section("47", "Interchanges and transfers", {
  ## Which interchanges do riders from more and less deprived catchments use,
  ## and when do they connect to other lines there? (Exploratory; Monday-Thursday
  ## of March-June and October 2025; flows and interchange visits from script 40.)
  ##
  ## (1) Hub use. Riders from the least and the most deprived thirds of origins
  ## (complexes with residents, not flagged, at least 2,000 assigned weekday
  ## boardings) are followed to their alighting stop by proportional downstream
  ## assignment. For each of the 64 interchanges: each group's share of all its
  ## alightings, observed and under the network baseline of script 46 (boardings
  ## spread over the lines and directions at the origin in proportion to their
  ## departures), with 95% intervals from 499 bootstrap resamples of osiedla.
  ## (2) Transfer rates. The counts do not link a rider's alighting to the next
  ## boarding. At an interchange, riders who transfer raise the boardings of a
  ## departure in proportion to the riders who alighted from other lines in the
  ## minutes before it. For every counted departure (not the trip's last stop):
  ##   B_d = beta_p * X_d + f(headway) + line x direction x 10-minute-slot effects
  ##         + date x hour effects,
  ## where X_d = alightings from other lines arriving within w = 10 minutes before
  ## the departure (counted visits; uncounted visits get the mean of their line,
  ## direction and hour) and headway = minutes since the previous visit of the
  ## same line and direction. beta_p (by period of the day) is the share of those
  ## alightings boarding that departure; summed over all departures it gives the
  ## transfer rate T = share of alightings at the hub that board another line
  ## within w minutes. Hub-and-period rates are shrunk towards the city-wide rate
  ## of the period (empirical Bayes). Placebo: alightings arriving in the w
  ## minutes AFTER the departure (should give ~0). With line x direction x hour
  ## and date effects the placebo failed (0.08-0.21: demand peaks within the hour
  ## raise arrivals and boardings together); comparing the same 10-minute slot
  ## across days, net of date x hour shocks, brings it to ~0 except in the evening
  ## (about -0.08 on the busiest hubs). Checks: w = 5 and 15.
  ## Standard errors clustered by date. Transfers to stops of another name are
  ## not counted.
  ## (3) Group connections: each group's alightings at each interchange and hour
  ## times the hub's transfer rate in that period; connection opportunity = mean
  ## number of departures of other lines within 10 minutes after arriving.
  ## Outputs: output/stage2/hub_use.csv, hub_transfer_rates.csv, hub_connections.csv;
  ##          output/figures/Fig_hub_use.png, Fig_hub_connections.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines); library(sandwich); library(ggplot2); library(patchwork)})
  set.seed(20261008)
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R"))   # flows expanded from counted to all operated trips
  s2 <- file.path(out_dir, "output/stage2")
  MON <- c(3:6, 10)
  cv  <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  HUBS <- cv[interchange %in% TRUE, stop_name]
  PER <- function(h) as.character(cut(h, c(-1, 3, 5, 8, 14, 18, 99), labels = c("night", "early", "am_peak", "midday", "pm_peak", "evening")))

  ## ── origins and groups (as in script 46, weekday all day) ───────────────
  od <- rbindlist(lapply(sprintf("%02d", MON), function(m) expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))))))[daytype == "monthu" & o != "outside"]
  cv[, lres := log1p(as.numeric(residents))]
  B_o <- od[, .(boardings = sum(flow_c)), by = o]   # counted support, for the threshold
  orig <- merge(cv[residents > 0 & flagged == FALSE, .(o = stop_name, dep, osiedle)], B_o, by = "o")[boardings >= 2000]
  orig[, third := cut(frank(dep, ties.method = "first"), quantile(seq_len(.N), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("least", "middle", "most"))]
  orig <- orig[third != "middle"]
  cat("Origins: least deprived third", orig[third == "least", .N], "| most deprived third", orig[third == "most", .N], "\n")
  g <- merge(od, orig[, .(o, third, osiedle)], by = "o")

  ## ── (1) hub use: observed and network baseline ──────────────────────────
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))[daytype == "monthu" & month %in% MON]
  depw <- cells[, .(dep = sum(departures)), by = .(o = stop_name, NR_LINI, dir, period = as.character(period))]; rm(cells); gc()
  seg <- g[, .(flow = sum(flow)), by = .(o, third, osiedle, NR_LINI, dir, period, d)]
  seg[, Bs := sum(flow), by = .(o, NR_LINI, dir, period)]
  seg <- merge(seg, depw, by = c("o", "NR_LINI", "dir", "period"))
  seg[, Bo := sum(flow), by = o]
  seg[, ws := dep / sum(dep[!duplicated(paste(NR_LINI, dir, period))]), by = o]      # departure share of each segment at the origin
  seg[, null_flow := Bo * ws * flow / Bs]                                               # origin boardings spread over segments by departures
  per_o <- seg[, .(obs = sum(flow), net = sum(null_flow)), by = .(o, third, osiedle, d)]
  per_o[, `:=`(obs_tot = sum(obs), net_tot = sum(net)), by = o]
  use_fun <- function(P) {
    tot <- P[, .(obs_tot = sum(obs), net_tot = sum(net)), by = third]
    h <- P[d %in% HUBS, .(obs = sum(obs), net = sum(net)), by = .(third, d)]
    h <- merge(h, tot, by = "third"); h[, `:=`(obs_share = 100 * obs / obs_tot, net_share = 100 * net / net_tot)]
    w <- dcast(h, d ~ third, value.var = c("obs_share", "net_share"), fill = 0)
    w[, `:=`(diff_obs = obs_share_most - obs_share_least, diff_net = net_share_most - net_share_least)]
    w
  }
  U <- use_fun(per_o)
  os_list <- unique(orig$osiedle)
  BT <- rbindlist(lapply(1:499, function(b) {
    k <- sample(os_list, replace = TRUE)
    P <- rbindlist(lapply(seq_along(k), function(i) per_o[osiedle == k[i]]))
    use_fun(P)[, .(d, diff_obs, diff_net)]
  }))
  CI <- BT[, .(lo = quantile(diff_obs, .025), hi = quantile(diff_obs, .975), lo_net = quantile(diff_net, .025), hi_net = quantile(diff_net, .975)), by = d]
  U <- merge(U, CI, by = "d", all.x = TRUE)
  setorder(U, -diff_obs)
  fwrite(U, file.path(s2, "hub_use.csv"))
  cat("\nInterchanges most over-used by riders from the most deprived third (pp of all their alightings; observed vs network baseline):\n")
  print(U[1:12, .(hub = d, least = round(obs_share_least, 2), most = round(obs_share_most, 2), diff = round(diff_obs, 2), lo = round(lo, 2), hi = round(hi, 2),
                 net_diff = round(diff_net, 2))])
  cat("\n... and by riders from the least deprived third:\n")
  print(U[order(diff_obs)][1:12, .(hub = d, least = round(obs_share_least, 2), most = round(obs_share_most, 2), diff = round(diff_obs, 2), lo = round(lo, 2), hi = round(hi, 2),
                                   net_diff = round(diff_net, 2))])

  ## ── (2) transfer rates at each interchange ──────────────────────────────
  hol <- as.IDate(c("2025-01-01", "2025-01-06", "2025-04-20", "2025-04-21", "2025-05-01", "2025-05-03", "2025-06-08", "2025-06-19", "2025-08-15"))
  bridge <- as.IDate(c("2025-01-02", "2025-01-03", "2025-05-02", "2025-06-20")); dst <- as.IDate(c("2025-03-30", "2025-10-26"))
  HV <- rbindlist(lapply(sprintf("%02d", MON), function(m) readRDS(file.path(s2, "hubv", sprintf("hubvisits_2025_%s.rds", m)))))
  HV <- HV[!is.na(cx) & !DATA %in% c(hol, bridge, dst) & as.integer(format(as.Date(DATA), "%u")) <= 4 & !grepl("^2[0-9][0-9]$", NR_LINI)]
  HV[, period := PER(sched_hour)]
  ## imputed alightings for uncounted arrivals: mean of counted visits of the same hub, line, direction and hour (then hub, line, hour; then 0)
  HV[, A_imp := fcoalesce(A, mean(A, na.rm = TRUE)), by = .(cx, NR_LINI, dir, sched_hour)]
  HV[, A_imp := fcoalesce(A_imp, mean(A, na.rm = TRUE)), by = .(cx, NR_LINI, sched_hour)]
  HV[is.na(A_imp) | is.nan(A_imp) | first == TRUE, A_imp := 0]
  setorder(HV, cx, DATA, NR_LINI, dir, sec)
  HV[, headway := (sec - shift(sec)) / 60, by = .(cx, DATA, NR_LINI, dir)]
  cat("\nInterchange visits (Mon-Thu, term months):", nrow(HV), "| counted departures:", HV[counted == TRUE & last == FALSE, .N], "\n")

  window_sum <- function(H, w, after = FALSE) {
    arr <- H[A_imp > 0, .(DATA, NR_LINI_a = NR_LINI, sa = sec, A_imp)]
    dep <- H[last == FALSE, .(id, DATA, NR_LINI, lo = if (after) sec + 10L else sec - w * 60L, hi = if (after) sec + w * 60L else sec - 10L)]
    j <- arr[dep, on = .(DATA, sa >= lo, sa <= hi), .(id = i.id, NR_LINI = i.NR_LINI, NR_LINI_a = x.NR_LINI_a, A_imp = x.A_imp), nomatch = 0]
    j[NR_LINI_a != NR_LINI, .(X = sum(A_imp), n_arr = .N), by = id]
  }
  dm2 <- function(D, cols, f1, f2, it = 8) {                                            # two-way within transformation
    for (k in seq_len(it)) for (v in cols) { D[, (v) := get(v) - mean(get(v)), by = f1]; D[, (v) := get(v) - mean(get(v)), by = f2] }
    D
  }
  fit_hub <- function(h, w = 10) {
    H <- HV[cx == h]; H[, id := .I]
    X1 <- window_sum(H, w); X0 <- window_sum(H, w, after = TRUE); setnames(X0, c("X", "n_arr"), c("Xpost", "n_post"))
    H <- merge(merge(H, X1, by = "id", all.x = TRUE), X0, by = "id", all.x = TRUE)
    for (v in c("X", "Xpost", "n_arr", "n_post")) H[is.na(get(v)), (v) := 0]
    ## departures of other lines within w minutes after each arrival (connection opportunity)
    Dp <- H[last == FALSE, .(DATA, NR_LINI_d = NR_LINI, sd = sec)]
    Ar <- H[A_imp > 0, .(aid = id, DATA, NR_LINI, lo = sec + 10L, hi = sec + w * 60L, sched_hour, A_imp)]
    opp <- Dp[Ar, on = .(DATA, sd >= lo, sd <= hi), .(aid = i.aid, NR_LINI = i.NR_LINI, NR_LINI_d = x.NR_LINI_d), nomatch = NA][, .(n_conn = sum(!is.na(NR_LINI_d) & NR_LINI_d != NR_LINI)), by = aid]
    H[opp, n_conn := i.n_conn, on = c(id = "aid")]
    M <- H[counted == TRUE & last == FALSE & !is.na(headway) & headway > 0 & headway <= 60 & !is.na(B) & period %in% c("early", "am_peak", "midday", "pm_peak", "evening")]
    if (M[, uniqueN(DATA)] < 20 || nrow(M) < 500) return(NULL)
    M[, hw := pmin(headway, 30)]; hb <- tryCatch(ns(M$hw, 3), error = function(e) cbind(M$hw, log(M$hw), sqrt(M$hw))); M[, `:=`(h1 = hb[, 1], h2 = hb[, 2], h3 = hb[, 3])]
    M[, `:=`(stratum = paste(NR_LINI, dir, sec %/% 600), dh = paste(DATA, sched_hour))]; M <- M[, if (.N >= 3) .SD, by = stratum]
    for (p in unique(M$period)) M[, (paste0("X_", p)) := fifelse(period == p, X, 0)][, (paste0("P_", p)) := fifelse(period == p, Xpost, 0)]
    xc <- grep("^(X|P)_", names(M), value = TRUE)
    M <- dm2(M, c("B", xc, "h1", "h2", "h3"), "stratum", "dh", it = 10)
    f <- lm(reformulate(c(xc, "h1", "h2", "h3", "0"), "B"), M)
    V <- sandwich::vcovCL(f, cluster = M$DATA, type = "HC1")
    est <- coef(f); se <- sqrt(diag(V))
    ## transfer rate: beta x (sum of X over all departures) / (sum of alightings), by period
    tot <- H[period %in% c("early", "am_peak", "midday", "pm_peak", "evening"),
             .(sumX = sum(X[last == FALSE]), sumXpost = sum(Xpost[last == FALSE]), sumA = sum(A_imp), arrivals_A = sum(A_imp), opp = weighted.mean(n_conn, A_imp, na.rm = TRUE),
               departures = sum(last == FALSE), counted_dep = sum(counted & !last)), by = period]
    tot[, `:=`(beta = est[paste0("X_", period)], se = se[paste0("X_", period)], beta_post = est[paste0("P_", period)], se_post = se[paste0("P_", period)])]
    tot[, `:=`(T = beta * sumX / sumA, T_se = se * sumX / sumA, T_post = beta_post * sumXpost / sumA)]
    ## hourly connection opportunity for weighting group arrivals
    oh <- H[A_imp > 0, .(opp = weighted.mean(n_conn, A_imp, na.rm = TRUE), A = sum(A_imp)), by = sched_hour]
    list(rates = tot[, hub := h][, w := w], opp = oh[, hub := h])
  }
  hubs_here <- intersect(HUBS, unique(HV$cx))
  R10 <- lapply(hubs_here, fit_hub, w = 10); names(R10) <- hubs_here
  RT <- rbindlist(lapply(R10, `[[`, "rates")); OPP <- rbindlist(lapply(R10, `[[`, "opp"))
  RT5 <- rbindlist(lapply(hubs_here, function(h) fit_hub(h, 5)$rates)); RT15 <- rbindlist(lapply(hubs_here, function(h) fit_hub(h, 15)$rates))
  ## empirical-Bayes shrinkage of hub rates towards the city-wide rate of each period
  RT[, w_ := 1 / T_se^2]
  pool <- RT[is.finite(T) & is.finite(T_se) & T_se > 0, .(T_pool = sum(T * w_) / sum(w_), tau2 = max(0, var(T) - mean(T_se^2))), by = period]
  RT <- merge(RT, pool, by = "period")
  RT[, T_eb := fifelse(is.finite(T_se) & T_se > 0, (T / T_se^2 + T_pool / tau2) / (1 / T_se^2 + 1 / tau2), T_pool)]
  RT[tau2 == 0, T_eb := T_pool]
  RT[, T_eb := pmin(pmax(T_eb, 0), 1)]
  fwrite(rbind(RT[, !c("w_")], RT5, RT15, fill = TRUE), file.path(s2, "hub_transfer_rates.csv"))
  cat("\nCity-wide transfer rate by period (share of alightings at an interchange boarding another line within 10 min; arrivals-weighted mean of hub rates):\n")
  print(RT[, .(hubs = .N, T_raw = weighted.mean(T, sumA), T_eb = weighted.mean(T_eb, sumA), placebo = weighted.mean(T_post, sumA)), by = period][order(match(period, c("early", "am_peak", "midday", "pm_peak", "evening")))], digits = 3)
  cat("Window check (arrivals-weighted mean T): w = 5:", round(RT5[, weighted.mean(T, sumA, na.rm = TRUE)], 3), "| w = 10:", round(RT[, weighted.mean(T, sumA, na.rm = TRUE)], 3),
      "| w = 15:", round(RT15[, weighted.mean(T, sumA, na.rm = TRUE)], 3), "\n")
  cat("\nLargest interchanges (alightings, Mon-Thu): transfer rate by period (EB-shrunk), placebo in brackets\n")
  big <- RT[, .(A = sum(sumA)), by = hub][order(-A)][1:15, hub]
  print(dcast(RT[hub %in% big], hub ~ factor(period, levels = c("early", "am_peak", "midday", "pm_peak", "evening")),
              value.var = "T_eb", fun.aggregate = function(v) round(v, 2)))

  ## ── (3) group connections ───────────────────────────────────────────────
  gh <- g[d %in% hubs_here, .(flow = sum(flow)), by = .(third, d, hour_d)][, period := PER(hour_d)]
  gh <- merge(gh, RT[, .(d = hub, period, T_eb)], by = c("d", "period"), all.x = TRUE)
  gh <- merge(gh, OPP[, .(d = hub, hour_d = sched_hour, opp)], by = c("d", "hour_d"), all.x = TRUE)
  gh[is.na(T_eb), T_eb := 0]
  gh[, conn := flow * T_eb]
  tot_g <- g[, .(all_flow = sum(flow)), by = third]
  CG <- gh[, .(hub_alightings = sum(flow), connections = sum(conn), opp = weighted.mean(opp, flow, na.rm = TRUE)), by = third]
  CG <- merge(CG, tot_g, by = "third")[, `:=`(hub_share = 100 * hub_alightings / all_flow, conn_per_100_rides = 100 * connections / all_flow,
                                               conn_rate_at_hubs = 100 * connections / hub_alightings)]
  cat("\nEstimated connections at interchanges, by origin group (Mon-Thu):\n"); print(CG, digits = 3)
  CH <- gh[, .(alight = sum(flow), conn = sum(conn), opp = weighted.mean(opp, flow, na.rm = TRUE)), by = .(third, hour_d)]
  CH <- merge(CH, tot_g, by = "third")[, `:=`(conn_per_100 = 100 * conn / all_flow, conn_share_of_day = 100 * conn / sum(conn)), by = third]
  CHUB <- gh[, .(conn = sum(conn), alight = sum(flow)), by = .(third, d)][, conn_share := 100 * conn / sum(conn), by = third]
  fwrite(rbind(CG[, level := "group"], CH[, level := "hour"], CHUB[, level := "hub"], fill = TRUE), file.path(s2, "hub_connections.csv"))
  cat("\nWhen connections happen (share of each group's estimated connections by hour of arrival at the interchange):\n")
  print(dcast(CH, hour_d ~ third, value.var = c("conn_share_of_day", "opp"), fun.aggregate = function(v) round(v, 1)))
  cat("\nWhere they connect (top interchanges by share of each group's estimated connections):\n")
  print(dcast(CHUB, d ~ third, value.var = "conn_share", fill = 0)[order(-(least + most))][1:15], digits = 3)

  ## ── figures ──────────────────────────────────────────────────────────────
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  xy <- data.table(d = C$stop_name, X = st_coordinates(C)[, 1], Y = st_coordinates(C)[, 2])
  MU <- merge(U, xy, by = "d")[, size := obs_share_least + obs_share_most]
  lab <- MU[order(-abs(diff_obs))][1:12]
  m <- max(abs(MU$diff_obs))
  p1 <- ggplot() + geom_sf(data = osi, fill = "#f6f4ee", colour = "#ffffff", linewidth = 0.3) +
    geom_point(data = MU, aes(X, Y, size = size, fill = diff_obs), shape = 21, colour = "#3d3d3a", stroke = 0.25)
  p1 <- p1 + geom_text(data = lab, aes(X, Y, label = d), size = 2.6, vjust = -1.1, colour = "#2b2b28") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f6f4ee", high = "#e34948", midpoint = 0, limits = c(-m, m), name = "Most minus least\ndeprived (pp)") +
    scale_size_area(max_size = 9, name = "Share of the two\ngroups' alightings (%)") + coord_sf(datum = NA) + theme_void(base_size = 10) +
    labs(title = "Interchanges used by riders from the most and the least deprived catchments") +
    theme(plot.title = element_text(face = "bold", size = 10.5), plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(out_dir, "output/figures/Fig_hub_use.png"), p1, width = 9, height = 8, dpi = 220, bg = "#ffffff")
  CHp <- copy(CH)[, third := factor(third, levels = c("least", "most"), labels = c("Least deprived third", "Most deprived third"))]
  p2 <- ggplot(CHp, aes(hour_d + 0.5, conn_per_100, colour = third)) + geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
    scale_colour_manual(values = c("Least deprived third" = "#2a78d6", "Most deprived third" = "#e34948"), name = NULL) +
    scale_x_continuous(breaks = seq(4, 24, 2), limits = c(4, 24)) +
    labs(x = "Hour of arrival at the interchange", y = "Estimated connections per 100 rides", title = "When riders connect to other lines") +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5))
  p3 <- ggplot(CHp, aes(hour_d + 0.5, opp, colour = third)) + geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
    scale_colour_manual(values = c("Least deprived third" = "#2a78d6", "Most deprived third" = "#e34948"), name = NULL) +
    scale_x_continuous(breaks = seq(4, 24, 2), limits = c(4, 24)) +
    labs(x = "Hour of arrival at the interchange", y = "Departures of other lines within 10 min", title = "Connection opportunities where and when they arrive") +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5))
  ggsave(file.path(out_dir, "output/figures/Fig_hub_connections.png"), p2 / p3, width = 8, height = 8, dpi = 220, bg = "#ffffff")
  cat("\nSaved hub_use.csv, hub_transfer_rates.csv, hub_connections.csv, Fig_hub_use.png, Fig_hub_connections.png\n")
})

## ── 50 · Destination and interchange maps, and their reasons ───────────────────────────────────────────────
section("50", "Destination and interchange maps, and their reasons", {
  ## Where do the destination and interchange differences between riders from more
  ## and less deprived catchments come from, and why? (Exploratory; Monday-Thursday
  ## of March-June and October 2025; flows from script 40, profiles from script 46,
  ## transfer rates from script 47.)
  ##
  ## Maps
  ##   Fig_hubsheds.png   (a) deprivation of the stop catchments; (b) each origin
  ##                      catchment linked to the interchange where most of its
  ##                      riders alight (hub-sheds), line coloured by origin deprivation
  ##   Fig_flows_groups.png  morning-peak flows from the origins of the least and the
  ##                      most deprived thirds to the stops where they alight
  ##   Fig_network_vs_use_map.png  per origin catchment: destination deprivation and
  ##                      ride length offered by the network (riders spread over the
  ##                      lines and directions at the stop in proportion to departures),
  ##                      observed, and the difference (lines and directions used)
  ##   Fig_connections_map.png  estimated connections at each interchange by group
  ## Reasons: decomposition (Gelbach) of the deprivation gradient in each origin's
  ## ride length, destination deprivation, share alighting at interchanges, share
  ## carried to the centre and estimated connections per ride (weekday all day and
  ## morning peak) into the parts carried by
  ##   location    distance to the centre (spline)
  ##   service     log lines and log departures at the stop, tram share of
  ##               departures, share of departures on lines that reach the centre
  ##               (within 1.5 km) further along
  ##   population  shares aged 0-18, 19-25 and 61+, log residents, log catchment area
  ##   land use    residential vs non-residential balance, land-use mix, tree cover
  ## and a remainder; intervals from 499 bootstrap resamples of osiedla.
  ## Also: which origin osiedla account for each group's use of the main interchanges.
  ## Outputs: output/stage2/reasons_gelbach.csv, hub_origin_osiedla.csv; figures above.

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines); library(ggplot2); library(patchwork)})
  set.seed(20261008)
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R"))   # flows expanded from counted to all operated trips
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  MON <- c(3:6, 10)
  PER <- function(h) as.character(cut(h, c(-1, 3, 5, 8, 14, 18, 99), labels = c("night", "early", "am_peak", "midday", "pm_peak", "evening")))
  cv  <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  td  <- fread(file.path(s2, "catchment_tissue_descriptors_exclusive_400m_network.csv"), encoding = "UTF-8")[, .(stop_name, b1 = ilr_res_vs_nonres, landuse_mix, wc_tree)]
  HUBS <- cv[interchange %in% TRUE, stop_name]
  C   <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  CP  <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  city <- st_union(osi)
  xy  <- data.table(stop_name = C$stop_name, X = st_coordinates(C)[, 1], Y = st_coordinates(C)[, 2])

  ## ── flows, origins, transfer rates ──────────────────────────────────────
  od <- rbindlist(lapply(sprintf("%02d", MON), function(m) expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))))))[daytype == "monthu" & o != "outside"]
  RT <- fread(file.path(s2, "hub_transfer_rates.csv"))[w == 10 & !is.na(T_eb), .(d = hub, period, T_eb)]
  od[, p_d := PER(hour_d)]
  od <- merge(od, RT, by.x = c("d", "p_d"), by.y = c("d", "period"), all.x = TRUE)[is.na(T_eb), T_eb := 0]
  od <- merge(od, cv[, .(d = stop_name, centre_km = dist_km)], by = "d", all.x = TRUE)
  B_o <- od[, .(boardings = sum(flow_c)), by = o]   # counted support, for the threshold
  orig <- merge(cv[residents > 0 & flagged == FALSE, .(o = stop_name, dep, osiedle, dist_km, share_0_18, share_19_25, share_61plus,
                                                        lres = log1p(as.numeric(residents)), larea = log(catchment_area_m2))], B_o, by = "o")[boardings >= 2000]
  orig <- merge(orig, td[, .(o = stop_name, b1, landuse_mix, wc_tree)], by = "o")
  orig[, third := cut(frank(dep, ties.method = "first"), quantile(seq_len(.N), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("least", "middle", "most"))]

  ## service at the origin stop
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))[daytype == "monthu" & month %in% MON & night == FALSE]
  sv <- cells[, .(n_lines = uniqueN(NR_LINI), departures = sum(departures), tram_share = sum(departures[mode == "tram"]) / sum(departures)), by = .(o = stop_name)]
  segdep <- cells[, .(dep = sum(departures)), by = .(o = stop_name, NR_LINI, dir)]; rm(cells); gc()
  reach <- od[, .(reach_centre = any(centre_km <= 1.5, na.rm = TRUE)), by = .(o, NR_LINI, dir)]
  reach <- merge(reach, segdep, by = c("o", "NR_LINI", "dir"))[, .(direct_centre = sum(dep * reach_centre) / sum(dep)), by = o]
  orig <- merge(merge(orig, sv, by = "o"), reach, by = "o")
  orig[, `:=`(l_lines = log(n_lines), l_dep = log(departures))]

  ## origin outcomes, weekday all day and morning peak
  prof <- function(ff) ff[, .(hub_share = 100 * sum(flow[d %in% HUBS]) / sum(flow), conn_100 = 100 * sum(flow * T_eb) / sum(flow),
                              to_centre = 100 * sum(flow[!is.na(centre_km) & centre_km <= 1.5]) / sum(flow)), by = o]
  mp <- fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & comp == "observed" &
         window %in% c("weekday all day", "weekday morning peak (06-09)"), .(o = stop_name, window, ride_km, dest_dep)]
  OUT <- rbind(merge(prof(od), mp[window == "weekday all day"], by = "o")[, win := "weekday all day"],
               merge(prof(od[period == "am_peak"]), mp[window == "weekday morning peak (06-09)"], by = "o")[, win := "morning peak (06-09)"])
  D <- merge(orig, OUT[, !"window"], by = "o")
  cat("Origins:", D[win == "weekday all day", .N], "| osiedla:", D[win == "weekday all day", uniqueN(osiedle)], "\n")

  ## ── reasons: Gelbach decomposition of the deprivation gradient ──────────
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  GR <- list(location = c("dc1", "dc2", "dc3"), service = c("l_lines", "l_dep", "tram_share", "direct_centre"),
             population = c("share_0_18", "share_19_25", "share_61plus", "lres", "larea"), land_use = c("b1", "landuse_mix", "wc_tree"))
  gel <- function(dd, y) {
    dd <- copy(dd[is.finite(get(y))])
    bs <- ns(dd$dist_km, 3); dd[, `:=`(dc1 = bs[, 1], dc2 = bs[, 2], dc3 = bs[, 3])]
    for (v in c(y, "dep", unname(unlist(GR)))) dd[, (v) := zs(get(v))]
    base <- coef(lm(reformulate("dep", y), dd))[["dep"]]
    full <- coef(lm(reformulate(c("dep", unname(unlist(GR))), y), dd))
    gam <- vapply(unname(unlist(GR)), function(k) coef(lm(reformulate("dep", k), dd))[["dep"]], numeric(1))
    c(gradient = base, vapply(GR, function(g) sum(gam[g] * full[g]), numeric(1)), remainder = full[["dep"]])
  }
  YS <- c(ride_km = "Ride length", dest_dep = "Destination deprivation", hub_share = "Alighting at interchanges", to_centre = "Carried to the centre",
          conn_100 = "Connections per ride")
  GL <- rbindlist(lapply(unique(D$win), function(w) rbindlist(lapply(names(YS), function(y) {
    dw <- D[win == w]; ix <- split(seq_len(nrow(dw)), dw$osiedle); on <- names(ix)
    e <- gel(dw, y)
    bt <- t(replicate(499, { k <- sample(on, replace = TRUE); gel(dw[unlist(ix[k])], y) }))
    data.table(window = w, outcome = y, part = names(e), est = e, lo = apply(bt, 2, quantile, .025, na.rm = TRUE), hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
  }))))
  fwrite(GL, file.path(s2, "reasons_gelbach.csv"))
  GL[, s := sprintf("%5.2f%s", est, fifelse(lo > 0 | hi < 0, "*", " "))]
  cat("\nWhy riders from more deprived catchments differ: deprivation gradient (SD per SD, no controls) split into parts (* 95% interval excludes 0)\n")
  print(dcast(GL, window + outcome ~ factor(part, levels = c("gradient", names(GR), "remainder")), value.var = "s"), width = 200)
  cat("\nOrigin characteristics, least vs most deprived third (medians):\n")
  print(D[win == "weekday all day" & third != "middle", .(dist_km = median(dist_km), lines = median(n_lines), dep_per_day = median(departures) / 85,
                                                        tram_share = median(tram_share), direct_centre = median(direct_centre), share_19_25 = median(share_19_25),
                                                        ride_km = median(ride_km), hub_share = median(hub_share), conn_100 = median(conn_100)), by = third], digits = 3)

  ## ── which origin osiedla feed each group's use of the main interchanges ─
  G <- merge(od, orig[third != "middle", .(o, third, osiedle)], by = "o")
  gt <- G[, .(tot = sum(flow)), by = third]
  HO <- G[d %in% HUBS, .(flow = sum(flow)), by = .(third, d, osiedle)]
  HO <- merge(HO, gt, by = "third")[, share := 100 * flow / tot]
  U <- fread(file.path(s2, "hub_use.csv"), encoding = "UTF-8")
  top_h <- c(U[order(-diff_obs)][1:4, d], U[order(diff_obs)][1:4, d])
  HOt <- HO[d %in% top_h][order(d, third, -share)][, .SD[1:3], by = .(d, third)][!is.na(osiedle)]
  fwrite(HO, file.path(s2, "hub_origin_osiedla.csv"))
  cat("\nWhere each group's use of the main interchanges comes from (origin osiedla; % of the group's alightings):\n")
  print(HOt[, .(origins = paste(sprintf("%s %.2f", osiedle, share), collapse = "; ")), by = .(hub = d, third)], width = 200)

  ## ── maps ────────────────────────────────────────────────────────────────
  DIV <- c("#2a78d6", "#f0efec", "#e34948")
  base_theme <- theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10.5), plot.background = element_rect(fill = "#ffffff", colour = NA),
    legend.position = "bottom", legend.key.width = unit(1.3, "cm"), legend.key.height = unit(0.25, "cm"), legend.title.position = "top",
    plot.caption = element_text(colour = "#6b6a64", size = 7.5, hjust = 0))
  source(file.path(REPO, "R/map_style.R")); ML <- ms_layers(out_dir)   # shared map style: land, water, outline
  city_layer <- ms_base(ML, borders = FALSE)
  outline <- list(ms_water_top(ML), geom_sf(data = osi, fill = NA, colour = "#ffffff", linewidth = 0.3, inherit.aes = FALSE), ms_outline(ML))
  CPd <- merge(CP, cv[residents > 0, .(stop_name, dep)], by = "stop_name")
  mdep <- max(abs(CPd$dep), na.rm = TRUE)
  hub_xy <- merge(data.table(d = HUBS), xy, by.x = "d", by.y = "stop_name")
  lab_layer <- function(dd, size = 2.6) {
    if (requireNamespace("ggrepel", quietly = TRUE)) ggrepel::geom_text_repel(data = dd, aes(X, Y, label = d), size = size, colour = "#2b2b28", min.segment.length = 0.2,
                                                                               seed = 1, box.padding = 0.3, max.overlaps = 30)
    else geom_text(data = dd, aes(X, Y, label = d), size = size, vjust = -1, colour = "#2b2b28")
  }

  ## (1) hub-sheds
  main_hub <- od[d %in% HUBS, .(flow = sum(flow)), by = .(o, d)][order(o, -flow)][, .SD[1], by = o]
  main_hub <- merge(main_hub, orig[, .(o, dep)], by = "o")
  main_hub <- merge(merge(main_hub, xy, by.x = "o", by.y = "stop_name"), xy, by.x = "d", by.y = "stop_name", suffixes = c("", "_h"))
  hub_n <- main_hub[, .(origins = .N, mean_dep = mean(dep)), by = d]
  hub_n <- merge(hub_n, xy, by.x = "d", by.y = "stop_name")
  pa <- ggplot() + geom_sf(data = osi, fill = "#f7f6f2", colour = NA) + geom_sf(data = CPd, aes(fill = dep), colour = NA) + outline +
    scale_fill_gradient2(low = DIV[1], mid = DIV[2], high = DIV[3], midpoint = 0, limits = c(-mdep, mdep), name = "Catchment deprivation (z)") +
    coord_sf(datum = NA) + labs(title = "(a) Deprivation of stop catchments") + base_theme
  pb <- ggplot() + city_layer + outline +
    geom_segment(data = main_hub, aes(X, Y, xend = X_h, yend = Y_h, colour = dep), linewidth = 0.45, alpha = 0.85) +
    geom_point(data = main_hub, aes(X, Y, colour = dep), size = 0.9) +
    geom_point(data = hub_n, aes(X, Y, size = origins), shape = 21, fill = "#ffffff", colour = "#2b2b28", stroke = 0.6) +
    lab_layer(hub_n[origins >= 5]) +
    scale_colour_gradient2(low = DIV[1], mid = "#bdbbb3", high = DIV[3], midpoint = 0, limits = c(-mdep, mdep), name = "Origin deprivation (z)") +
    scale_size_area(max_size = 6, name = "Origins whose main interchange it is") + coord_sf(datum = NA) +
    labs(title = "(b) Hub-sheds: the interchange where most riders from each origin alight") + base_theme
  ggsave(file.path(fig, "Fig_hubsheds.png"), pa + pb, width = 14, height = 7.6, dpi = 220, bg = "#ffffff")   # caption in the manuscript

  ## (2) morning-peak flows by group
  FL <- G[period == "am_peak" & d != "outside", .(flow = sum(flow)), by = .(third, o, d)]
  FL <- merge(merge(FL, xy, by.x = "o", by.y = "stop_name"), xy, by.x = "d", by.y = "stop_name", suffixes = c("", "_d"))
  FL[, share := flow / sum(flow), by = third]
  FL <- FL[order(third, -share)][, cum := cumsum(share), by = third][cum <= 0.6]               # the flows carrying 60% of each group's riders
  RCV <- FL[, .(recv = sum(share)), by = .(third, d, X = X_d, Y = Y_d)]
  FLp <- copy(FL)[!(X == X_d & Y == Y_d)][, third := factor(third, levels = c("least", "most"), labels = c("Riders from the least deprived third", "Riders from the most deprived third"))]
  RCVp <- copy(RCV)[, third := factor(third, levels = c("least", "most"), labels = levels(FLp$third))]
  ORp <- merge(orig[third != "middle", .(o, third)], xy, by.x = "o", by.y = "stop_name")[, third := factor(third, levels = c("least", "most"), labels = levels(FLp$third))]
  pc <- ggplot() + city_layer + outline +
    geom_curve(data = FLp, aes(X, Y, xend = X_d, yend = Y_d, linewidth = share, colour = third), curvature = 0.18, alpha = 0.55, lineend = "round") +
    geom_point(data = ORp, aes(X, Y, colour = third), size = 0.9, shape = 16) +
    geom_point(data = RCVp, aes(X, Y, size = recv), shape = 21, fill = "#ffffff", colour = "#2b2b28", stroke = 0.5) +
    lab_layer(RCVp[order(-recv)][, .SD[1:6], by = third], size = 2.4) +
    facet_wrap(~ third) + scale_colour_manual(values = c("#2a78d6", "#e34948"), guide = "none") +
    scale_linewidth(range = c(0.15, 2.2), name = "Share of the group's morning-peak riders", labels = scales::percent) +
    scale_size_area(max_size = 6, name = "Share alighting at the stop", labels = scales::percent) + coord_sf(datum = NA) +
    base_theme + theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig, "Fig_flows_groups.png"), pc, width = 14, height = 7.6, dpi = 220, bg = "#ffffff")

  ## (3) network vs use, per origin catchment
  NV <- fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & window == "weekday all day" & comp %in% c("observed", "network"),
          .(stop_name, comp, dest_dep, ride_km)]
  NV <- dcast(NV, stop_name ~ comp, value.var = c("dest_dep", "ride_km"))
  NV[, `:=`(dest_dep_use = dest_dep_observed - dest_dep_network, ride_km_use = ride_km_observed - ride_km_network)]
  NVp <- merge(CP, NV, by = "stop_name")
  md <- max(abs(c(NVp$dest_dep_network, NVp$dest_dep_observed)), na.rm = TRUE); mu <- max(abs(NVp$dest_dep_use), na.rm = TRUE)
  mr <- range(c(NVp$ride_km_network, NVp$ride_km_observed), na.rm = TRUE); mru <- max(abs(NVp$ride_km_use), na.rm = TRUE)
  mk <- function(v, title, lim, div = TRUE, name) {
    g <- ggplot() + city_layer + geom_sf(data = NVp, aes(fill = .data[[v]]), colour = NA) + outline + coord_sf(datum = NA) + labs(title = title) + base_theme
    if (div) g + scale_fill_gradient2(low = DIV[1], mid = DIV[2], high = DIV[3], midpoint = 0, limits = c(-lim, lim), name = name)
    else g + scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = lim, name = name)
  }
  pn <- (mk("dest_dep_network", "Destination deprivation: offered by the network", md, TRUE, "Mean deprivation of destinations (z)") |
         mk("dest_dep_observed", "... reached by riders", md, TRUE, "Mean deprivation of destinations (z)") |
         mk("dest_dep_use", "... difference made by the lines and directions used", mu, TRUE, "Reached minus offered")) /
        (mk("ride_km_network", "Ride length: offered by the network", mr, FALSE, "km along the route") |
         mk("ride_km_observed", "... ridden", mr, FALSE, "km along the route") |
         mk("ride_km_use", "... difference made by the lines and directions used", mru, TRUE, "Ridden minus offered (km)"))
  ggsave(file.path(fig, "Fig_network_vs_use_map.png"), pn, width = 15, height = 11, dpi = 200, bg = "#ffffff")
  ## main-text version: destination deprivation only, three panels with tags; (a) and (b) share one legend
  sub3 <- function(v, s, lim, name) ggplot() + city_layer + geom_sf(data = NVp, aes(fill = .data[[v]]), colour = NA) + outline + coord_sf(datum = NA) +
    scale_fill_gradient2(low = DIV[1], mid = DIV[2], high = DIV[3], midpoint = 0, limits = c(-lim, lim), name = name) + labs(subtitle = s) + base_theme +
    theme(plot.subtitle = element_text(size = 9, hjust = 0.5), plot.tag = element_text(face = "bold"), legend.position = "bottom")
  p3 <- (sub3("dest_dep_network", "Offered by the network", md, "Deprivation of destinations (z)") |
         sub3("dest_dep_observed", "Reached by riders", md, "Deprivation of destinations (z)") |
         sub3("dest_dep_use", "Reached minus offered", mu, "Reached minus offered (z)")) +
    plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
  ggsave(file.path(fig, "Fig_network_vs_use_dest.png"), p3, width = 13, height = 5.6, dpi = 250, bg = "#ffffff")

  ## (4) estimated connections by interchange and group
  CN <- G[d %in% HUBS, .(conn = sum(flow * T_eb)), by = .(third, d)][, conn := 100 * conn / sum(conn), by = third]     # % of the group's estimated connections
  CN <- merge(CN, xy, by.x = "d", by.y = "stop_name")
  AMR <- RT[period == "am_peak", .(d, T_am = T_eb)]
  CN <- merge(CN, AMR, by = "d", all.x = TRUE)[, third := factor(third, levels = c("least", "most"), labels = c("Riders from the least deprived third", "Riders from the most deprived third"))]
  pd <- ggplot() + city_layer + outline +
    geom_point(data = CN, aes(X, Y, size = conn, fill = T_am), shape = 21, colour = "#2b2b28", stroke = 0.4) +
    lab_layer(CN[order(-conn)][, .SD[1:7], by = third], size = 2.4) + facet_wrap(~ third) +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", name = "Morning-peak transfer rate", labels = scales::percent) +
    scale_size_area(max_size = 9, name = "Share of the group's estimated connections (%)") + coord_sf(datum = NA) +
    base_theme + theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig, "Fig_connections_map.png"), pd, width = 14, height = 7.6, dpi = 220, bg = "#ffffff")
  ## (5) which interchanges each group uses: share of the group's alightings at each interchange, most minus least
  ##     deprived third, observed and offered by the network (hub_use.csv of script 47; bootstrap interval of osiedla)
  HU <- fread(file.path(s2, "hub_use.csv"), encoding = "UTF-8")
  HU <- merge(C[, "stop_name"], HU, by.x = "stop_name", by.y = "d")
  HU$size <- (HU$obs_share_least + HU$obs_share_most) / 2
  sigf <- function(lo, hi) factor(lo > 0 | hi < 0, levels = c(TRUE, FALSE), labels = c("Interval excludes zero", "Interval includes zero"))
  HU$sig_obs <- sigf(HU$lo, HU$hi); HU$sig_net <- sigf(HU$lo_net, HU$hi_net)   # each panel outlined by its own interval
  HU <- HU[order(HU$size), ]
  lim <- max(abs(c(HU$diff_obs, HU$diff_net)), na.rm = TRUE)
  lab <- HU[order(-abs(HU$diff_obs)), ][1:8, ]
  hb <- as.numeric(st_bbox(HU))   # map extent: the interchanges and 1.5 km around them
  hm <- function(v, s, sg, labels = FALSE) {
    g <- ggplot() + city_layer + outline +
      geom_sf(data = HU, aes(size = size, fill = .data[[v]], colour = .data[[sg]]), shape = 21, stroke = 0.5) +
      scale_fill_gradient2(low = DIV[1], mid = DIV[2], high = DIV[3], midpoint = 0, limits = c(-lim, lim),
                           name = "Most minus least deprived third\n(percentage points of alightings)") +
      scale_colour_manual(values = c("Interval excludes zero" = "#2b2a27", "Interval includes zero" = "#ffffff"), name = NULL, drop = FALSE,
                          guide = guide_legend(override.aes = list(size = 4, fill = "#d9d7d0", stroke = 0.8))) +
      scale_size_area(max_size = 9, breaks = c(1, 2, 3), name = "Share of alightings (%)") +
      coord_sf(xlim = hb[c(1, 3)] + c(-1500, 1500), ylim = hb[c(2, 4)] + c(-1500, 1500), datum = NA) + labs(subtitle = s) + base_theme +
      theme(plot.subtitle = element_text(size = 9.5, hjust = 0.5), plot.tag = element_text(face = "bold"), legend.position = "bottom")
    if (labels) g <- g + ggrepel::geom_text_repel(data = cbind(st_drop_geometry(lab), st_coordinates(st_centroid(st_geometry(lab)))),
                                                   aes(x = X, y = Y, label = stop_name), size = 2.6, colour = "#2b2a27",
                                                   min.segment.length = 0, box.padding = 0.4, seed = 1)
    g
  }
  ph <- (hm("diff_obs", "Observed", "sig_obs", labels = TRUE) | hm("diff_net", "Offered by the network", "sig_net")) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, "Fig_hub_use_map.png"), ph, width = 12, height = 6.8, dpi = 250, bg = "#ffffff")
  cat("\nSaved reasons_gelbach.csv, hub_origin_osiedla.csv, Fig_hubsheds.png, Fig_flows_groups.png, Fig_network_vs_use_map.png, Fig_network_vs_use_dest.png, Fig_connections_map.png, Fig_hub_use_map.png\n")
})

## ── 51 · Route heterogeneity ───────────────────────────────────────────────────────────────────────────────
section("51", "Route heterogeneity", {
  ## Route heterogeneity: are the trips of riders from more deprived catchments
  ## more or less concentrated on a few destinations, lines and directions, and
  ## does that change at weekends, on holidays and in the summer break?
  ## (Exploratory; flows by line and direction from script 40.)
  ##
  ## Per origin stop and window, from the destination distribution of its riders
  ## (proportional downstream assignment on counted trips; "outside" counts as
  ## one destination):
  ##   eff_dest   effective number of destinations, exp(Shannon entropy)
  ##   top3       share of riders going to the three largest destinations
  ##   eff_lines  effective number of line x direction combinations boarded
  ##   spread_km  standard distance of the destinations (km; flow-weighted)
  ##   dir_conc   directional concentration: mean resultant length of the
  ##              origin-to-destination bearings (0 = all directions, 1 = one)
  ## Each measure is also computed for the network baseline of script 46
  ## (boardings spread over the lines and directions at the stop in proportion
  ## to their departures); use = observed minus network.
  ## Windows: Monday-Thursday (term months March-June and October) whole day,
  ## morning peak and midday; Saturday; Sunday and public holidays; Monday-
  ## Thursday of the summer break (July-August).
  ## Tests: medians for the least and most deprived thirds (groups fixed on the
  ## weekday definition); slope on origin deprivation with the common adjustments
  ## (CR2 by osiedle, t with G - 2 df) for observed, network and use; change of
  ## the gradient from weekdays to each other window (per-origin difference).
  ## Group-level concentration: number of origin-destination pairs that carry
  ## half of each group's riders, share of the 20 largest pairs, share of the 5
  ## largest line-directions (descriptive: a bootstrap over origins duplicates
  ## origins and inflates concentration); largest
  ## routes of each group.
  ## Outputs: output/stage2/routes_origin.csv, routes_tests.csv, routes_group.csv,
  ##          routes_top.csv; output/figures/Fig_route_concentration.png,
  ##          Fig_route_heterogeneity_map.png, Fig_route_corridors.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(splines); library(ggplot2); library(patchwork)})
  set.seed(20261008)
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R"))   # flows expanded from counted to all operated trips
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  MON <- c(3:6, 10); SUM <- 7:8
  od <- rbindlist(lapply(sprintf("%02d", c(MON, SUM)), function(m)
    expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))))[daytype %in% c("monthu", "sat", "sunhol") & o != "outside",
                                                                .(flow = sum(flow), flow_c = sum(flow_c)), by = .(o, d, NR_LINI, dir, month, daytype, period)]))
  cells <- readRDS(file.path(s2, "stage2_line_stop_cells.rds"))[daytype %in% c("monthu", "sat", "sunhol") & month %in% c(MON, SUM)]
  depw <- cells[, .(dep = sum(departures)), by = .(o = stop_name, NR_LINI, dir, month, daytype, period = as.character(period))]; rm(cells); gc()
  WIN <- list("Mon-Thu"                  = list(day = "monthu", mon = MON, per = NULL, minb = 2000),
              "Mon-Thu morning peak"     = list(day = "monthu", mon = MON, per = "am_peak", minb = 500),
              "Mon-Thu midday"           = list(day = "monthu", mon = MON, per = "midday", minb = 500),
              "Saturday"                 = list(day = "sat", mon = MON, per = NULL, minb = 500),
              "Sunday and holidays"      = list(day = "sunhol", mon = MON, per = NULL, minb = 300),
              "Mon-Thu summer break"     = list(day = "monthu", mon = SUM, per = NULL, minb = 800))

  C  <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  xy <- data.table(stop_name = C$stop_name, X = st_coordinates(C)[, 1], Y = st_coordinates(C)[, 2])
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  cv[, `:=`(lres = log1p(as.numeric(residents)), larea = log(catchment_area_m2))]
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))

  ## ── distribution measures ───────────────────────────────────────────────
  measures <- function(D) {                                       # D: o, d, f (one origin-destination distribution per origin)
    D <- D[f > 0]                                                 # zero flows (riders cannot pass a stop where everyone alights)
    D <- merge(D, xy[, .(o = stop_name, oX = X, oY = Y)], by = "o")
    D <- merge(D, xy[, .(d = stop_name, dX = X, dY = Y)], by = "d", all.x = TRUE)
    D[, p := f / sum(f), by = o]
    D[, .(eff_dest = exp(-sum(p * log(p))), top3 = sum(head(sort(p, decreasing = TRUE), 3)),
          spread_km = { q <- p[!is.na(dX)] / sum(p[!is.na(dX)]); mx <- sum(q * dX[!is.na(dX)]); my <- sum(q * dY[!is.na(dX)])
                        sqrt(sum(q * ((dX[!is.na(dX)] - mx)^2 + (dY[!is.na(dX)] - my)^2))) / 1000 },
          dir_conc = { k <- !is.na(dX) & (dX != oX | dY != oY); q <- p[k] / sum(p[k]); th <- atan2(dY[k] - oY[k], dX[k] - oX[k])
                       sqrt(sum(q * cos(th))^2 + sum(q * sin(th))^2) }), by = o]
  }
  ent <- function(v) { p <- v[v > 0] / sum(v[v > 0]); exp(-sum(p * log(p))) }
  window_profiles <- function(w) {
    W <- WIN[[w]]
    g <- od[daytype == W$day & month %in% W$mon & (is.null(W$per) | period %in% W$per)]
    dp <- depw[daytype == W$day & month %in% W$mon & (is.null(W$per) | period %in% W$per), .(dep = sum(dep)), by = .(o, NR_LINI, dir, period)]
    seg <- g[, .(flow = sum(flow), flow_c = sum(flow_c)), by = .(o, NR_LINI, dir, period, d)]
    sB <- seg[, .(B = sum(flow), Bc = sum(flow_c)), by = .(o, NR_LINI, dir, period)]
    sB <- merge(sB, dp, by = c("o", "NR_LINI", "dir", "period"))[dep > 0 & B > 0]
    sB[, Bo := sum(B), by = o][, ws := dep / sum(dep), by = o]
    seg <- merge(seg, sB, by = c("o", "NR_LINI", "dir", "period"))
    obs <- seg[, .(f = sum(flow)), by = .(o, d)]
    net <- seg[, .(f = sum(Bo * ws * flow / B)), by = .(o, d)]
    lo <- sB[, .(eff_lines = ent(B), boardings = sum(Bc)), by = o]   # boardings = counted support, for the threshold
    ln <- sB[, .(eff_lines = ent(dep)), by = o]
    rbind(merge(measures(obs), lo, by = "o")[, comp := "observed"],
          merge(measures(net), ln[, boardings := NA_real_], by = "o")[, comp := "network"])[, window := w]
  }
  PR <- rbindlist(lapply(names(WIN), window_profiles))
  PR[, boardings := boardings[comp == "observed"][1], by = .(o, window)]
  PR <- merge(PR, cv[residents > 0 & flagged == FALSE, .(o = stop_name, dep, dist_km, osiedle, share_0_18, share_61plus, lres, larea)], by = "o")
  PR <- PR[boardings >= sapply(WIN, `[[`, "minb")[window]]
  ## groups fixed on the weekday definition
  G0 <- PR[window == "Mon-Thu" & comp == "observed", .(o, dep)]
  G0[, third := cut(frank(dep, ties.method = "first"), quantile(seq_len(.N), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("least", "middle", "most"))]
  PR <- merge(PR, G0[, .(o, third)], by = "o")
  fwrite(PR, file.path(s2, vtag("routes_origin.csv")))
  print(PR[comp == "observed", .(origins = .N), by = window])

  ## ── tests ───────────────────────────────────────────────────────────────
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  cr2 <- function(fit, cl) {
    X <- model.matrix(fit); e <- resid(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]
    B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) { i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); Am <- ev$vectors %*% diag(1 / sqrt(pmax(ev$values, 1e-10)), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (Am %*% e[i]); meat <- meat + u %*% t(u) }
    setNames(sqrt(diag(B %*% meat %*% B)), colnames(X))
  }
  COV <- c("share_0_18", "share_61plus", "lres", "larea", "ns(dist_km, 3)")
  slope <- function(dd, y, X = COV) {
    dd <- copy(dd[is.finite(get(y))]); for (v in c("dep", "share_0_18", "share_61plus", "lres", "larea")) dd[, (v) := zs(get(v))]
    f <- lm(reformulate(c("dep", X), y), dd); used <- as.integer(rownames(model.matrix(f))); se <- cr2(f, dd$osiedle[used])[["dep"]]; G <- uniqueN(dd$osiedle[used])
    data.table(est = coef(f)[["dep"]], se = se, p = 2 * pt(-abs(coef(f)[["dep"]] / se), G - 2), n = length(used))
  }
  M <- c("eff_dest", "top3", "eff_lines", "spread_km", "dir_conc")
  TS <- rbindlist(lapply(names(WIN), function(w) rbindlist(lapply(M, function(y) {
    W <- dcast(PR[window == w], o + dep + dist_km + osiedle + share_0_18 + share_61plus + lres + larea ~ comp, value.var = y)[is.finite(observed) & is.finite(network)]
    s <- sd(W$observed); W[, `:=`(obs = observed / s, net = network / s, use = (observed - network) / s)]
    rbind(slope(W, "obs")[, component := "observed"], slope(W, "net")[, component := "network"], slope(W, "use")[, component := "use"],
          slope(W, "obs", c(COV, "factor(osiedle)"))[, component := "observed, within osiedla"])[, `:=`(window = w, measure = y)]
  }))))
  TS[, q := p.adjust(p, "BH"), by = .(window, component)]
  ## change of the gradient from weekdays: per-origin difference, in SD of the weekday measure
  CHG <- rbindlist(lapply(setdiff(names(WIN), "Mon-Thu"), function(w) rbindlist(lapply(M, function(y) {
    W <- merge(PR[window == "Mon-Thu" & comp == "observed", c("o", "dep", "dist_km", "osiedle", "share_0_18", "share_61plus", "lres", "larea", y), with = FALSE],
               PR[window == w & comp == "observed", c("o", y), with = FALSE], by = "o", suffixes = c("_wd", "_w"))
    s <- sd(W[[paste0(y, "_wd")]]); W[, delta := (get(paste0(y, "_w")) - get(paste0(y, "_wd"))) / s]
    slope(W, "delta")[, `:=`(window = w, measure = y, mean_change = mean(W$delta), component = "change from Mon-Thu")]
  }))))
  fwrite(rbind(TS, CHG, fill = TRUE), file.path(s2, vtag("routes_tests.csv")))

  MED <- PR[comp == "observed" & third != "middle", lapply(.SD, median), by = .(window, third), .SDcols = M]
  fmt <- function(e, p) sprintf("%5.2f%s", e, fifelse(p < 0.01, "**", fifelse(p < 0.05, "* ", "  ")))
  for (y in M) {
    cat("\n==", y, "== medians (least / most deprived third) and slope on deprivation (SD per SD; * p<.05 ** p<.01)\n")
    X <- dcast(MED[, .(window, third, v = get(y))], window ~ third, value.var = "v")
    X <- merge(X, dcast(TS[measure == y, .(window, component, s = fmt(est, p))], window ~ component, value.var = "s"), by = "window")
    X <- merge(X, CHG[measure == y, .(window, change_slope = fmt(est, p), mean_change = round(mean_change, 2))], by = "window", all.x = TRUE)
    X[, window := factor(window, levels = names(WIN))]; setorder(X, window)
    print(X[, .(window, least = signif(least, 3), most = signif(most, 3), observed, network, use, `observed, within osiedla`, change_slope, mean_change)], row.names = FALSE)
  }

  ## ── group-level concentration and largest routes ────────────────────────
  conc <- function(g) {
    p <- sort(g$f / sum(g$f), decreasing = TRUE)
    data.table(pairs = length(p), n50 = which(cumsum(p) >= 0.5)[1], top20 = sum(head(p, 20)))
  }
  AGG <- lapply(names(WIN), function(w) {                         # per window: origin-destination and origin-line tables (small)
    W <- WIN[[w]]; g <- od[daytype == W$day & month %in% W$mon & (is.null(W$per) | period %in% W$per)]
    list(od = g[, .(f = sum(flow)), by = .(o, d)], ln = g[, .(f = sum(flow)), by = .(o, NR_LINI, dir)])
  }); names(AGG) <- names(WIN)
  grp <- function(ids, w) {
    a <- AGG[[w]]; idt <- data.table(o = ids)                     # ids may repeat (bootstrap): repeated origins count again
    g1 <- a$od[idt, on = "o", allow.cartesian = TRUE, nomatch = 0][, .(f = sum(f)), by = .(o, d)]
    g2 <- a$ln[idt, on = "o", allow.cartesian = TRUE, nomatch = 0][, .(f = sum(f)), by = .(NR_LINI, dir)]
    cbind(conc(g1), top5_lines = g2[order(-f)][, sum(head(f, 5)) / sum(f)])
  }
  GA <- rbindlist(lapply(names(WIN), function(w) rbindlist(lapply(c("least", "most"), function(t) {
    ids <- PR[window == w & comp == "observed" & third == t, o]
    grp(ids, w)[, `:=`(window = w, third = t, origins = length(ids))]
  }))))
  ## (no bootstrap here: resampling origins with replacement duplicates origins and inflates concentration; inference rests on the origin-level tests)
  fwrite(GA, file.path(s2, vtag("routes_group.csv")))
  cat("\nGroup-level concentration (origins fixed on the weekday thirds; descriptive):\n")
  print(GA[, .(window, third, origins, pairs, n50, top20 = round(top20, 3), top5_lines = round(top5_lines, 3))], row.names = FALSE)
  gs <- unique(fread(file.path(GTFS_DIR, "stops.txt"), encoding = "UTF-8")[, .(code = sub("[0-9]+$", "", stop_code), stop_name)])
  gs <- gs[!duplicated(code)]; TERM_code <- gs$code; TERM <- gs$stop_name                 # terminal cluster code -> stop name
  TOP <- rbindlist(lapply(c("Mon-Thu", "Saturday", "Sunday and holidays"), function(w) rbindlist(lapply(c("least", "most"), function(t) {
    W <- WIN[[w]]; ids <- PR[window == w & comp == "observed" & third == t, o]
    g <- od[o %in% ids & daytype == W$day & month %in% W$mon]
    r <- g[, .(f = sum(flow)), by = .(o, d)][, share := 100 * f / sum(f)][order(-share)][1:8]
    l <- g[, .(f = sum(flow)), by = .(NR_LINI, dir)][, share := 100 * f / sum(f)][order(-share)][1:5]
    l[, term := TERM[match(dir, TERM_code)]][is.na(term), term := dir]
    rbind(r[, .(type = "pair", route = paste(o, "->", d), share)], l[, .(type = "line", route = paste("line", NR_LINI, "to", term), share)])[, `:=`(window = w, third = t)]
  }))))
  fwrite(TOP, file.path(s2, vtag("routes_top.csv")))
  cat("\nLargest routes (share of the group's riders, %):\n")
  print(TOP[, .(routes = paste(sprintf("%s %.1f", route, share), collapse = "; ")), by = .(window, third, type)], width = 250)

  ## ── figures ──────────────────────────────────────────────────────────────
  ## (1) concentration curves
  CC <- rbindlist(lapply(names(WIN)[c(1, 4, 5, 6)], function(w) rbindlist(lapply(c("least", "most"), function(t) {
    W <- WIN[[w]]; ids <- PR[window == w & comp == "observed" & third == t, o]
    g <- od[o %in% ids & daytype == W$day & month %in% W$mon][, .(f = sum(flow)), by = .(o, d)][order(-f)]
    g[, .(rank = seq_len(.N), cum = cumsum(f) / sum(f))][rank <= 400][, `:=`(window = w, third = t)]
  }))))
  CC[, `:=`(window = factor(window, levels = names(WIN)), third = factor(third, levels = c("least", "most"), labels = c("Least deprived third", "Most deprived third")))]
  p1 <- ggplot(CC, aes(rank, cum, colour = third)) + geom_line(linewidth = 0.8) + facet_wrap(~ window, nrow = 1) +
    geom_hline(yintercept = 0.5, colour = "#9a9893", linewidth = 0.3, linetype = 2) +
    scale_colour_manual(values = c("Least deprived third" = "#2a78d6", "Most deprived third" = "#e34948"), name = NULL) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Origin-destination pairs, largest first", y = "Cumulative share of the group's riders") +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
                                          plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, vtag("Fig_route_concentration.png")), p1, width = 12, height = 4.2, dpi = 220, bg = "#ffffff")
  ## (2) map: effective number of destinations, weekday and Sunday, and change
  CP <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  MW <- dcast(PR[comp == "observed" & window %in% c("Mon-Thu", "Sunday and holidays"), .(stop_name = o, window, eff_dest)], stop_name ~ window, value.var = "eff_dest")
  setnames(MW, c("stop_name", "wd", "sun")); MW[, chg := sun - wd]
  MW <- merge(MW, PR[comp == "observed" & window == "Mon-Thu", .(stop_name = o, dirc = dir_conc)], by = "stop_name", all.x = TRUE)
  MP <- merge(CP, MW, by = "stop_name")
  lim <- range(MP$wd, na.rm = TRUE); m <- max(abs(MP$chg), na.rm = TRUE)
  base_theme <- theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10.5), plot.background = element_rect(fill = "#ffffff", colour = NA),
    legend.position = "bottom", legend.key.width = unit(1.3, "cm"), legend.key.height = unit(0.25, "cm"), legend.title.position = "top")
  lay <- list(geom_sf(data = osi, fill = "#f7f6f2", colour = "#d8d6cf", linewidth = 0.25, inherit.aes = FALSE))
  out <- geom_sf(data = osi, fill = NA, colour = "#ffffff", linewidth = 0.3, inherit.aes = FALSE)
  q1 <- ggplot() + lay + geom_sf(data = MP, aes(fill = wd), colour = NA) + out + coord_sf(datum = NA) +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = lim, name = "Effective number of destinations") + labs(title = "(a) Monday-Thursday") + base_theme
  q2 <- ggplot() + lay + geom_sf(data = MP[!is.na(MP$dirc), ], aes(fill = dirc), colour = NA) + out + coord_sf(datum = NA) +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = c(0, 1), name = "Directional concentration (0 = all directions, 1 = one)") +
    labs(title = "(b) Monday-Thursday: one direction or many") + base_theme
  q3 <- ggplot() + lay + geom_sf(data = MP[!is.na(MP$chg), ], aes(fill = chg), colour = NA) + out + coord_sf(datum = NA) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-m, m), name = "Sunday minus Monday-Thursday") +
    labs(title = "(c) Effective number of destinations: Sundays and holidays") + base_theme
  ggsave(file.path(fig, vtag("Fig_route_heterogeneity_map.png")), q1 | q2 | q3, width = 15, height = 6.2, dpi = 200, bg = "#ffffff")
  ## (3) the most condensed corridors: flows carrying the first 25% of each group's riders, Monday-Thursday and Sunday
  FC <- rbindlist(lapply(c("Mon-Thu", "Sunday and holidays"), function(w) rbindlist(lapply(c("least", "most"), function(t) {
    W <- WIN[[w]]; ids <- PR[window == w & comp == "observed" & third == t, o]
    g <- AGG[[w]]$od[o %in% ids & d != "outside"][, share := f / sum(f)][order(-share)][, cum := cumsum(share)]
    g[cum <= 0.25][, `:=`(window = w, third = t, n = .N)]
  }))))
  FC <- merge(merge(FC, xy, by.x = "o", by.y = "stop_name"), xy, by.x = "d", by.y = "stop_name", suffixes = c("", "_d"))
  FC[, panel := factor(paste0(fifelse(third == "least", "Least", "Most"), " deprived third, ", fifelse(window == "Mon-Thu", "Monday-Thursday", "Sundays and holidays"),
                              " (", n, " pairs)"))]
  FC[, panel := factor(panel, levels = unique(FC[order(third, window != "Mon-Thu"), panel]))]
  pc <- ggplot() + lay + out +
    geom_curve(data = FC, aes(X, Y, xend = X_d, yend = Y_d, linewidth = share, colour = third), curvature = 0.18, alpha = 0.6, lineend = "round") +
    geom_point(data = FC, aes(X, Y, colour = third), size = 0.8) + facet_wrap(~ panel, ncol = 2) +
    scale_colour_manual(values = c(least = "#2a78d6", most = "#e34948"), guide = "none") +
    scale_linewidth(range = c(0.2, 2), name = "Share of the group's riders", labels = scales::percent) + coord_sf(datum = NA) +
    base_theme + theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig, vtag("Fig_route_corridors.png")), pc, width = 11, height = 11, dpi = 200, bg = "#ffffff")
  cat("\nSaved routes_origin.csv, routes_tests.csv, routes_group.csv, routes_top.csv, Fig_route_concentration.png, Fig_route_heterogeneity_map.png\n")
})

## ── 53 · Route networks traced along the lines ─────────────────────────────────────────────────────────────
section("53", "Route networks traced along the lines", {
  ## The most common routes of riders from the most and the least deprived
  ## areas, on weekdays and at weekends, drawn along the network. (Exploratory;
  ## March-June and October 2025; flows by line and direction from script 40.)
  ##
  ## Each estimated ride (line, direction, boarding complex, alighting complex)
  ## is traced along the stop sequence the line actually runs (the three most
  ## common sequences of each line and direction in the vehicle records) and its
  ## riders are loaded onto every stop-to-stop segment between boarding and
  ## alighting; rides leaving the city are traced to the city limit. Segment
  ## loads are expressed as the share of the group's riders passing along the
  ## segment (both directions together).
  ## Groups: origins of the least and the most deprived thirds (fixed on the
  ## Monday-Thursday definition of script 51). Day windows: Monday-Thursday;
  ## weekend = Saturdays, Sundays and public holidays.
  ## Comparisons: overlap of the two groups' ride distributions over segments
  ## (sum of the smaller shares, 0-1), weekday against weekend; concentration
  ## (share of ride-segments on the busiest 5% of segments); the most common
  ## routes (line x boarding x alighting stop) of each group.
  ## Outputs: output/stage2/route_segments.csv, route_network_summary.csv,
  ##          route_top_rides.csv; output/figures/Fig_route_networks.png,
  ##          Fig_route_networks_difference.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  source(file.path(REPO, "R/od_expand.R"))   # flows expanded from counted to all operated trips
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  MON <- c(3:6, 10); t0 <- Sys.time()
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  gs <- fread(file.path(GTFS_DIR, "stops.txt"), encoding = "UTF-8")[, code := sub("[0-9]+$", "", stop_code)]
  gt <- unique(gs[, .(code, stop_name)])[stop_name %in% C$stop_name][!duplicated(code)]
  cc <- gs[, .(lon = mean(stop_lon), lat = mean(stop_lat)), by = code]
  xy <- st_coordinates(st_transform(st_as_sf(cc, coords = c("lon", "lat"), crs = 4326), st_crs(C))); cc[, `:=`(X = xy[, 1], Y = xy[, 2])]
  city_codes <- gt$code
  cx2code <- setNames(gt$code, gt$stop_name)

  ## ── stop sequences actually run (three most common per line x direction) ─
  SEQ <- rbindlist(lapply(sprintf("%02d", MON), function(m) {
    x <- readRDS(sprintf(file.path(PEKA_TABLES, "peka_2025_%s.rds"), m)); setDT(x)
    x <- unique(x[!is.na(SYM_PRZYST) & !is.na(NR_LINI) & NR_LINI != "", .(DT_KARTY, ID_KURSU, LP_PRZYST, NR_LINI, SYM_PRZYST)], by = c("DT_KARTY", "ID_KURSU", "LP_PRZYST"))
    setorder(x, DT_KARTY, ID_KURSU, LP_PRZYST)
    tr <- x[, .(NR_LINI = NR_LINI[1], dir = SYM_PRZYST[.N], seq = paste(SYM_PRZYST, collapse = ">")), by = .(DT_KARTY, ID_KURSU)]
    rm(x); gc(); tr[, .N, by = .(NR_LINI, dir, seq)]
  }))
  SEQ <- SEQ[, .(N = sum(N)), by = .(NR_LINI, dir, seq)][order(NR_LINI, dir, -N)][, rk := seq_len(.N), by = .(NR_LINI, dir)][rk <= 3]
  cat("Line-direction stop sequences:", SEQ[rk == 1, .N], "(", nrow(SEQ), "with variants )\n")

  ## ── groups and flows ─────────────────────────────────────────────────────
  GR <- fread(file.path(s2, "routes_origin.csv"), encoding = "UTF-8")[window == "Mon-Thu" & comp == "observed" & third != "middle", .(o, third)]
  od <- rbindlist(lapply(sprintf("%02d", MON), function(m) expand_od(readRDS(file.path(s2, "od", sprintf("od_2025_%s.rds", m))))[o %in% GR$o & daytype %in% c("monthu", "sat", "sunhol"),
          .(flow = sum(flow)), by = .(o, d, NR_LINI, dir, daytype)]))
  od[, win := fifelse(daytype == "monthu", "Monday-Thursday", "Weekend and holidays")]
  od <- merge(od[, .(flow = sum(flow)), by = .(o, d, NR_LINI, dir, win)], GR, by = "o")
  od[, `:=`(oc = cx2code[o], dc = fifelse(d == "outside", NA_character_, cx2code[d]))]

  ## ── trace each ride along its line's sequence ───────────────────────────
  trace_one <- function(ld, rides) {                                       # rides of one line x direction
    sq <- SEQ[NR_LINI == ld$NR_LINI & dir == ld$dir]
    rides[, `:=`(i = NA_integer_, j = NA_integer_, v = NA_integer_)]
    for (k in seq_len(nrow(sq))) {
      s <- strsplit(sq$seq[k], ">", fixed = TRUE)[[1]]; inside <- s %in% city_codes
      todo <- rides[is.na(i), which = TRUE]; if (!length(todo)) break
      pi <- match(rides$oc[todo], s)
      pj <- vapply(seq_along(todo), function(r) {
        a <- pi[r]; if (is.na(a)) return(NA_integer_)
        if (is.na(rides$dc[todo[r]])) { w <- which(!inside & seq_along(s) > a); return(if (length(w)) as.integer(w[1]) else NA_integer_) }   # leaving the city
        w <- which(s == rides$dc[todo[r]] & seq_along(s) > a); if (length(w)) as.integer(w[1]) else NA_integer_ }, integer(1))
      ok <- !is.na(pi) & !is.na(pj)
      rides[todo[ok], `:=`(i = pi[ok], j = pj[ok], v = k)]
    }
    TRACED <<- rbind(TRACED, rides[!is.na(i), .(traced = sum(flow)), by = .(third, win)])
    rides[!is.na(i)][, {
      s <- strsplit(sq$seq[v[1]], ">", fixed = TRUE)[[1]]
      idx <- unlist(Map(function(a, b) a:(b - 1L), i, j)); f <- rep(flow, j - i); g <- rep(third, j - i); w <- rep(win, j - i)
      data.table(a = s[idx], b = s[idx + 1L], flow = f, third = g, win = w)
    }, by = v][, v := NULL]
  }
  LD <- unique(od[, .(NR_LINI, dir)])
  TRACED <- data.table()
  SEG <- rbindlist(lapply(seq_len(nrow(LD)), function(r) {
    rides <- od[NR_LINI == LD$NR_LINI[r] & dir == LD$dir[r]]
    if (!nrow(SEQ[NR_LINI == LD$NR_LINI[r] & dir == LD$dir[r]])) return(NULL)
    trace_one(LD[r], copy(rides))
  }))
  tot <- od[, .(total = sum(flow)), by = .(third, win)]
  cat("Share of each group's rides traced along the network:\n")
  print(merge(TRACED[, .(traced = sum(traced)), by = .(third, win)], tot, by = c("third", "win"))[, share := round(traced / total, 3)][])
  ## segments: undirected key, city part only (both ends inside the city or the boundary segment)
  SEG[, `:=`(k1 = pmin(a, b), k2 = pmax(a, b))]
  SG <- SEG[, .(flow = sum(flow)), by = .(k1, k2, third, win)]
  SG <- merge(SG, tot, by = c("third", "win"))[, share := flow / total]                  # share of the group's riders passing the segment
  SG <- merge(SG, cc[, .(k1 = code, X1 = X, Y1 = Y)], by = "k1")
  SG <- merge(SG, cc[, .(k2 = code, X2 = X, Y2 = Y)], by = "k2")
  SG[, len_km := sqrt((X2 - X1)^2 + (Y2 - Y1)^2) / 1000]
  SG <- SG[len_km < 3]                                                                    # drop the rare jumps between non-adjacent stops
  fwrite(SG, file.path(s2, "route_segments.csv"))

  ## ── comparisons ──────────────────────────────────────────────────────────
  dist_of <- function(t, w) { s <- SG[third == t & win == w, .(key = paste(k1, k2), r = flow * len_km)]; s[, r := r / sum(r)] }   # ride-km distribution
  ov <- function(A, B) { m <- merge(A, B, by = "key", all = TRUE); m[is.na(m)] <- 0; sum(pmin(m$r.x, m$r.y)) }
  conc <- function(t, w) { s <- SG[third == t & win == w][order(-flow * len_km)]; k <- ceiling(0.05 * nrow(s)); sum(head(s$flow * s$len_km, k)) / sum(s$flow * s$len_km) }
  SUM <- rbind(
    data.table(comparison = "least vs most deprived", window = c("Monday-Thursday", "Weekend and holidays"),
               overlap = c(ov(dist_of("least", "Monday-Thursday"), dist_of("most", "Monday-Thursday")), ov(dist_of("least", "Weekend and holidays"), dist_of("most", "Weekend and holidays")))),
    data.table(comparison = "weekday vs weekend", window = c("least deprived third", "most deprived third"),
               overlap = c(ov(dist_of("least", "Monday-Thursday"), dist_of("least", "Weekend and holidays")), ov(dist_of("most", "Monday-Thursday"), dist_of("most", "Weekend and holidays")))))
  CONC <- CJ(third = c("least", "most"), win = c("Monday-Thursday", "Weekend and holidays"))[, top5pct := mapply(conc, third, win)]
  CONC <- merge(CONC, SG[, .(segments = .N, ride_km = sum(flow * len_km), mean_km = sum(flow * len_km) / total[1]), by = .(third, win)], by = c("third", "win"))
  fwrite(rbind(SUM, CONC, fill = TRUE), file.path(s2, "route_network_summary.csv"))
  cat("\nOverlap of ride-km distributions over segments (0 = disjoint, 1 = identical):\n"); print(SUM, digits = 3)
  cat("\nConcentration: share of ride-km on the busiest 5% of segments; mean km per ride traced:\n"); print(CONC, digits = 3)

  ## the most common routes (line x boarding x alighting stop)
  TOPR <- od[, .(flow = sum(flow)), by = .(third, win, NR_LINI, dir, o, d)][, share := 100 * flow / sum(flow), by = .(third, win)][order(third, win, -share)][, .SD[1:10], by = .(third, win)]
  TOPR[, term := gs[!duplicated(code)][match(dir, code), stop_name]]
  fwrite(TOPR, file.path(s2, "route_top_rides.csv"))
  cat("\nMost common routes (% of the group's rides):\n")
  print(TOPR[, .(routes = paste(sprintf("%s %s->%s %.2f", NR_LINI, o, d, share), collapse = "; ")), by = .(third, win)], width = 250)

  ## ── maps ────────────────────────────────────────────────────────────────
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  bb <- st_bbox(osi)
  base_theme <- theme_void(base_size = 10) + theme(plot.background = element_rect(fill = "#ffffff", colour = NA), strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom", legend.key.width = unit(1.3, "cm"), legend.key.height = unit(0.25, "cm"), legend.title.position = "top")
  lay <- list(geom_sf(data = osi, fill = "#f7f6f2", colour = "#e1dfd8", linewidth = 0.25, inherit.aes = FALSE))
  SGp <- copy(SG)[, panel := factor(paste0(fifelse(third == "least", "Least", "Most"), " deprived third: ", win),
                                     levels = c("Least deprived third: Monday-Thursday", "Least deprived third: Weekend and holidays",
                                                "Most deprived third: Monday-Thursday", "Most deprived third: Weekend and holidays"))]
  setorder(SGp, share)
  mx <- quantile(SGp$share, 0.995)
  p1 <- ggplot() + lay +
    geom_segment(data = SGp[share >= 0.002], aes(X1, Y1, xend = X2, yend = Y2, linewidth = pmin(share, mx), colour = third, alpha = pmin(share, mx)), lineend = "round") +
    facet_wrap(~ panel, ncol = 2) + scale_colour_manual(values = c(least = "#2a78d6", most = "#e34948"), guide = "none") +
    scale_linewidth(range = c(0.15, 3.2), name = "Share of the group's riders passing the segment", labels = scales::percent) + scale_alpha(range = c(0.25, 0.95), guide = "none") +
    coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], datum = NA) + base_theme
  ggsave(file.path(fig, "Fig_route_networks.png"), p1, width = 11, height = 11.5, dpi = 220, bg = "#ffffff")
  DF <- dcast(SG, k1 + k2 + X1 + Y1 + X2 + Y2 + win ~ third, value.var = "share", fill = 0)[, diff := most - least]
  m <- quantile(abs(DF$diff), 0.995)
  DF <- DF[abs(diff) >= 0.002][order(abs(diff))]
  p2 <- ggplot() + lay +
    geom_segment(data = DF, aes(X1, Y1, xend = X2, yend = Y2, colour = pmax(pmin(diff, m), -m), linewidth = pmin(abs(diff), m)), lineend = "round") +
    facet_wrap(~ win, ncol = 2) +
    scale_colour_gradient2(low = "#2a78d6", mid = "#d9d7d0", high = "#e34948", midpoint = 0, limits = c(-m, m), labels = scales::percent,
                           name = "Most minus least deprived third (share of riders)") +
    scale_linewidth(range = c(0.2, 3), guide = "none") + coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], datum = NA) + base_theme
  ggsave(file.path(fig, "Fig_route_networks_difference.png"), p2, width = 12, height = 6.6, dpi = 220, bg = "#ffffff")
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 54 · Access by the timetable ───────────────────────────────────────────────────────────────────────────
section("54", "Access by the timetable", {
  ## Does the place of deprived residents in the network make destinations
  ## harder to reach? Timetable-based accessibility from every stop complex
  ## (October 2025 GTFS; exploratory).
  ##
  ## Frequency-based transit graph per time window: stop nodes (stop clusters,
  ## i.e. complexes in the city) and line-pattern nodes (one per stop of each
  ## distinct stop sequence of a route and direction). Boarding costs half the
  ## headway of the pattern at that stop in the window (capped at 15 min); riding
  ## costs the median scheduled time between consecutive stops; alighting is free;
  ## walking links join stops within 500 m (straight line x 1.3, 80 m/min), used only
  ## as a transfer between rides or as the last leg (never chained; trips by transit).
  ## From every stop complex: generalised travel time to every other complex
  ## (wait + in-vehicle + walk), boardings on the fastest path (tie-break by a tiny
  ## cost per boarding) and the minimum number of boardings on any path.
  ## Windows: Tuesday-Thursday 05-06 (early), 07-09 (morning peak), 20-22
  ## (evening); Sunday 10-14.
  ## Per origin: city jobs reachable within 30 and 45 min (workplace grid jobs in
  ## the destination catchments), job-weighted mean travel time and transfers,
  ## share of jobs reachable without a transfer, and time and transfers to the
  ## centre (complexes within 1 km of Stary Rynek), the nearest hospital,
  ## university building and shopping centre (OpenStreetMap).
  ## Comparison: medians for the least and most deprived thirds of stops with
  ## residents; slope on catchment deprivation without and with distance to the
  ## centre (CR2 by osiedle).
  ## Outputs: output/stage2/access_origin.csv, access_tests.csv;
  ##          output/figures/Fig_access_map.png, Fig_access_windows.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(igraph); library(splines); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  G <- GTFS_DIR
  t0 <- Sys.time()
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  stops <- fread(file.path(G, "stops.txt"), encoding = "UTF-8")[, code := sub("[0-9]+$", "", stop_code)]
  gt <- unique(stops[, .(code, stop_name)])[stop_name %in% C$stop_name][!duplicated(code)]
  stops[gt, cx := i.stop_name, on = "code"]; stops[is.na(cx), cx := paste0("ext:", code)]       # node id: complex name in the city, cluster code outside
  nodexy <- stops[, .(lon = mean(stop_lon), lat = mean(stop_lat)), by = cx]
  xyp <- st_coordinates(st_transform(st_as_sf(nodexy, coords = c("lon", "lat"), crs = 4326), st_crs(C))); nodexy[, `:=`(X = xyp[, 1], Y = xyp[, 2])]
  trips <- fread(file.path(G, "trips.txt"), encoding = "UTF-8")
  st <- fread(file.path(G, "stop_times.txt"), encoding = "UTF-8", select = c("trip_id", "arrival_time", "departure_time", "stop_id", "stop_sequence"))
  tsec <- function(x) { p <- tstrsplit(x, ":", fixed = TRUE); as.integer(p[[1]]) * 3600L + as.integer(p[[2]]) * 60L + as.integer(p[[3]]) }
  st[, `:=`(arr = tsec(arrival_time), dep = tsec(departure_time))][, c("arrival_time", "departure_time") := NULL]
  st[stops, cx := i.cx, on = "stop_id"]
  setorder(st, trip_id, stop_sequence)
  ## pattern = route + direction + exact node sequence
  pat <- st[, .(seq = paste(cx, collapse = ">")), by = trip_id]
  pat <- merge(pat, trips[, .(trip_id, route_id, direction_id, service_id)], by = "trip_id")
  pat[, pattern := .GRP, by = .(route_id, direction_id, seq)]
  st <- merge(st, pat[, .(trip_id, pattern, service_id)], by = "trip_id")
  st[, pos := seq_len(.N), by = trip_id]
  cat("Trips:", nrow(pat), "| patterns:", uniqueN(pat$pattern), "| nodes:", uniqueN(st$cx), "\n")

  ## walking links between nodes within 500 m
  D <- as.matrix(dist(nodexy[, .(X, Y)]))
  wl <- which(D > 0 & D <= 500, arr.ind = TRUE)
  WALK <- data.table(from = nodexy$cx[wl[, 1]], to = nodexy$cx[wl[, 2]], w = D[wl] * 1.3 / 80, board = 0)

  WIN <- list("Weekday 05-06 (early)" = list(svc = 1, a = 5, b = 6), "Weekday 07-09 (morning peak)" = list(svc = 1, a = 7, b = 9),
              "Weekday 20-22 (evening)" = list(svc = 1, a = 20, b = 22), "Sunday 10-14" = list(svc = 4, a = 10, b = 14))
  build_graph <- function(W) {
    s <- st[service_id == W$svc]
    s[, `:=`(next_arr = shift(arr, -1), next_pos = shift(pos, -1)), by = trip_id]
    win <- s[dep >= W$a * 3600 & dep < W$b * 3600]
    hours <- W$b - W$a
    ## headway at each pattern stop in the window; ride times between consecutive stops
    hw <- win[!is.na(next_arr), .(n = .N), by = .(pattern, pos, cx)][, wait := pmin(60 * hours / n / 2, 15)]
    ride <- win[!is.na(next_arr), .(w = median((next_arr - dep) / 60)), by = .(pattern, pos)]
    ride[, w := pmax(w, 0.25)]
    pn <- function(p, i) paste0("p", p, "_", i)
    ## three layers per stop: D = ready to board, A = just alighted, F = arrived (terminal).
    ## Walking only as a transfer (A -> D of a nearby stop) or as the last leg (A -> F), so walks never chain.
    al <- s[pos > 1 & paste(pattern, pos - 1L) %in% ride[, paste(pattern, pos)], .(from = pn(pattern, pos), to = paste0("A|", cx))][!duplicated(paste(from, to))]
    nodes <- unique(st$cx)
    E <- rbind(hw[, .(from = paste0("D|", cx), to = pn(pattern, pos), w = wait, board = 1)],
               ride[, .(from = pn(pattern, pos), to = pn(pattern, pos + 1L), w, board = 0)],
               al[, .(from, to, w = 0, board = 0)],
               data.table(from = paste0("A|", nodes), to = paste0("D|", nodes), w = 0, board = 0),
               data.table(from = paste0("A|", nodes), to = paste0("F|", nodes), w = 0, board = 0),
               data.table(from = paste0("D|", nodes), to = paste0("F|", nodes), w = 0, board = 0),           # the origin itself
               WALK[, .(from = paste0("A|", from), to = paste0("D|", to), w, board = 0)],
               WALK[, .(from = paste0("A|", from), to = paste0("F|", to), w, board = 0)])
    E[w < 0, w := 0]
    graph_from_data_frame(E, directed = TRUE)
  }

  ## destinations and origins
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")
  city <- cv$stop_name
  am <- readRDS(file.path(out_dir, "output/stage1/_osm_amenities_raw.rds")); am <- st_transform(am, st_crs(C))
  CP <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  poi_cx <- function(sel) { p <- am[sel, ]; h <- st_intersects(p, CP); k <- lengths(h) > 0
    unique(c(CP$stop_name[unlist(h[k])], C$stop_name[st_nearest_feature(p[!k, ], C)])) }
  POI <- list(centre = cv[dist_km <= 1, stop_name],
              hospital = poi_cx(which(am$amenity %in% "hospital")),
              university = poi_cx(which(am$amenity %in% "university")),
              mall = poi_cx(which(am$shop %in% "mall")))
  cat("Destinations: centre", length(POI$centre), "| hospitals", length(POI$hospital), "| universities", length(POI$university), "| malls", length(POI$mall), "\n")
  jobs <- cv[, .(stop_name, jobs = fcoalesce(jobs, 0))]; JT <- sum(jobs$jobs)
  ORIG <- cv[residents > 0, stop_name]

  delta <- 1e-4; M <- 1000
  AC <- rbindlist(lapply(names(WIN), function(wn) {
    g <- build_graph(WIN[[wn]])
    V <- V(g)$name; o <- intersect(ORIG, sub("^D\\|", "", grep("^D\\|", V, value = TRUE))); dset <- intersect(city, sub("^F\\|", "", grep("^F\\|", V, value = TRUE)))
    w <- E(g)$w; b <- E(g)$board
    vo <- paste0("D|", o); vd <- paste0("F|", dset)
    T  <- distances(g, v = vo, to = vd, weights = w, mode = "out")
    Tb <- distances(g, v = vo, to = vd, weights = w + delta * b, mode = "out")
    Bm <- distances(g, v = vo, to = vd, weights = w + M * b, mode = "out")
    Bf <- round((Tb - T) / delta); Bmin <- floor(Bm / M)
    jw <- jobs$jobs[match(dset, jobs$stop_name)]; jw[is.na(jw)] <- 0
    res <- rbindlist(lapply(seq_along(o), function(i) {
      t <- T[i, ]; bf <- pmax(Bf[i, ] - 1, 0); bm <- pmax(Bmin[i, ] - 1, 0); ok <- is.finite(t)
      pt <- function(set) { k <- dset %in% set & ok; if (any(k)) { j <- which(k)[which.min(t[k])]; c(t[j], bf[j]) } else c(NA, NA) }
      ce <- pt(POI$centre); ho <- pt(POI$hospital); un <- pt(POI$university); ma <- pt(POI$mall)
      data.table(stop_name = o[i], jobs30 = sum(jw[ok & t <= 30]) / JT, jobs45 = sum(jw[ok & t <= 45]) / JT,
                 time_jobs = sum(jw[ok] * t[ok]) / sum(jw[ok]), transfers_jobs = sum(jw[ok] * bf[ok]) / sum(jw[ok]),
                 direct_jobs = sum(jw[ok & bm == 0]) / JT, reached = mean(ok),
                 t_centre = ce[1], x_centre = ce[2], t_hospital = ho[1], x_hospital = ho[2], t_univ = un[1], x_univ = un[2], t_mall = ma[1], x_mall = ma[2])
    }))
    cat(wn, ": graph", vcount(g), "nodes,", ecount(g), "edges; origins", length(o), "\n")
    res[, window := wn]
  }))

  ## ── comparison ──────────────────────────────────────────────────────────
  d <- merge(AC, cv[residents > 0 & flagged == FALSE, .(stop_name, dep, dist_km, osiedle, share_0_18, share_61plus, lres = log1p(as.numeric(residents)))], by = "stop_name")
  d[, third := cut(frank(dep, ties.method = "first"), quantile(seq_len(.N), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("least", "middle", "most")), by = window]
  fwrite(d, file.path(s2, "access_origin.csv"))
  zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  cr2 <- function(fit, cl) {
    X <- model.matrix(fit); e <- resid(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]
    B <- solve(crossprod(X)); meat <- matrix(0, ncol(X), ncol(X))
    for (g in unique(cl)) { i <- which(cl == g); Xg <- X[i, , drop = FALSE]; H <- Xg %*% B %*% t(Xg)
      ev <- eigen(diag(length(i)) - H, symmetric = TRUE); Am <- ev$vectors %*% diag(1 / sqrt(pmax(ev$values, 1e-10)), length(i)) %*% t(ev$vectors)
      u <- t(Xg) %*% (Am %*% e[i]); meat <- meat + u %*% t(u) }
    setNames(sqrt(diag(B %*% meat %*% B)), colnames(X))
  }
  MS <- c(jobs30 = "Jobs within 30 min (share of city)", jobs45 = "Jobs within 45 min (share of city)", time_jobs = "Mean travel time to jobs (min)",
          transfers_jobs = "Mean transfers to jobs (fastest path)", direct_jobs = "Jobs reachable without a transfer (share of city)",
          t_centre = "Time to the centre (min)", x_centre = "Transfers to the centre", t_hospital = "Time to the nearest hospital (min)",
          t_univ = "Time to the nearest university building (min)", t_mall = "Time to the nearest shopping centre (min)")
  TS <- rbindlist(lapply(unique(d$window), function(wn) rbindlist(lapply(names(MS), function(y) {
    dd <- d[window == wn & is.finite(get(y))]; dd[, yy := zs(get(y))]; dd[, dz := zs(dep)]
    f1 <- lm(yy ~ dz, dd); f2 <- lm(yy ~ dz + ns(dist_km, 3), dd)
    G1 <- uniqueN(dd$osiedle)
    data.table(window = wn, measure = y, least = dd[third == "least", median(get(y))], most = dd[third == "most", median(get(y))],
               raw = coef(f1)[["dz"]], p_raw = 2 * pt(-abs(coef(f1)[["dz"]] / cr2(f1, dd$osiedle)[["dz"]]), G1 - 2),
               net_dist = coef(f2)[["dz"]], p_net = 2 * pt(-abs(coef(f2)[["dz"]] / cr2(f2, dd$osiedle)[["dz"]]), G1 - 2))
  }))))
  fwrite(TS, file.path(s2, "access_tests.csv"))
  for (wn in unique(TS$window)) {
    cat("\n==", wn, "== medians (least / most deprived third); slope on deprivation raw and net of distance to the centre (SD per SD)\n")
    print(TS[window == wn, .(measure = MS[measure], least = signif(least, 3), most = signif(most, 3), raw = sprintf("%.2f (p=%.3f)", raw, p_raw),
                             net_of_distance = sprintf("%.2f (p=%.3f)", net_dist, p_net))], row.names = FALSE)
  }

  ## ── figures ──────────────────────────────────────────────────────────────
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), st_crs(C))
  base_theme <- theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10.5), plot.background = element_rect(fill = "#ffffff", colour = NA),
    legend.position = "bottom", legend.key.width = unit(1.3, "cm"), legend.key.height = unit(0.25, "cm"), legend.title.position = "top")
  lay <- list(geom_sf(data = osi, fill = "#f7f6f2", colour = "#e1dfd8", linewidth = 0.25, inherit.aes = FALSE))
  out <- geom_sf(data = osi, fill = NA, colour = "#ffffff", linewidth = 0.3, inherit.aes = FALSE)
  AMd <- merge(CP, d[window == "Weekday 07-09 (morning peak)"], by = "stop_name")
  CPd <- merge(CP, cv[residents > 0, .(stop_name, dep)], by = "stop_name"); mdep <- max(abs(CPd$dep))
  q1 <- ggplot() + lay + geom_sf(data = CPd, aes(fill = dep), colour = NA) + out + coord_sf(datum = NA) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-mdep, mdep), name = "Catchment deprivation (z)") +
    labs(title = "(a) Deprivation") + base_theme
  q2 <- ggplot() + lay + geom_sf(data = AMd, aes(fill = 100 * jobs30), colour = NA) + out + coord_sf(datum = NA) +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", name = "City jobs reachable within 30 min (%)") + labs(title = "(b) Jobs within 30 minutes, morning peak") + base_theme
  q3 <- ggplot() + lay + geom_sf(data = AMd, aes(fill = transfers_jobs), colour = NA) + out + coord_sf(datum = NA) +
    scale_fill_gradient(low = "#cde2fb", high = "#104281", name = "Mean transfers to the city's jobs") + labs(title = "(c) Transfers needed, morning peak") + base_theme
  ggsave(file.path(fig, "Fig_access_map.png"), q1 | q2 | q3, width = 15, height = 6.2, dpi = 200, bg = "#ffffff")
  FW <- d[third != "middle", .(jobs30 = 100 * median(jobs30), transfers = median(transfers_jobs), t_centre = median(t_centre, na.rm = TRUE)), by = .(window, third)]
  FW <- melt(FW, id.vars = c("window", "third"))
  FW[, `:=`(variable = factor(variable, levels = c("jobs30", "transfers", "t_centre"), labels = c("City jobs within 30 min (%)", "Mean transfers to jobs", "Time to the centre (min)")),
            window = factor(window, levels = names(WIN)), third = factor(third, levels = c("least", "most"), labels = c("Least deprived third", "Most deprived third")))]
  p4 <- ggplot(FW, aes(window, value, fill = third)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) + facet_wrap(~ variable, scales = "free_y") +
    scale_fill_manual(values = c("Least deprived third" = "#2a78d6", "Most deprived third" = "#e34948"), name = NULL) + labs(x = NULL, y = "Median over origin stops") +
    theme_minimal(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1),
                                          strip.text = element_text(face = "bold"), plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, "Fig_access_windows.png"), p4, width = 12, height = 4.8, dpi = 220, bg = "#ffffff")
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

