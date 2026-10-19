## DualPulse | stage 9: Main-text figures and graphical abstract
##
## Figures 2, 3 and 8 of the main text, the route networks and other maps for the supplement and the graphical abstract (Figure 1 is drawn in LaTeX; Figure 4 by section 75, Figure 5 by section 74 and Figures 6 and 7 by section 81 of stage 8).
## Supplementary figures are drawn by the sections that compute them.
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   57      Figure 2: how deprivation is assigned to stop catchments
##   60      Where the more and less deprived areas lie (supplement); network panel of Figure 3
##   64      Figure 8: the gaps and their parts
##   71      Deprivation of destinations (supplement)
##   70      Route networks of each third (supplementary S7)
##   69      Figure 3: what the network offers
##   72      Graphical abstract
##
## Run from the repository root: Rscript run.R 9   (one section: Rscript run.R 9:72)

if (!exists("section")) source("R/setup.R")

## ── 57 · Figure 2: how deprivation is assigned to stop catchments ──────────────────────────────────────────
section("57", "Figure 2: how deprivation is assigned to stop catchments", {
  ## How the deprivation measure and the walking catchments are built, and how
  ## the two are combined into the deprivation of each stop's catchment
  ## (manuscript §2.2c-d, §2.4-2.5). Reproduces the objects of scripts 20, 34 and
  ## 39; nothing is re-estimated.
  ##
  ## Main version (Figure 2), four panels: (a) the deprivation score by osiedle; in one
  ##   window, around Niezłomnych, where Jeżyce, Stare Miasto, Św. Łazarz and Wilda meet:
  ##   (b) the walking catchments on the pedestrian network, (c) catchment deprivation,
  ##   the residents' mean; (d) catchment deprivation of every stop, with the rings of
  ##   distance and the sectors that define location.
  ## Full version (supplementary Figure S1), twelve panels:
  ## Row 1, deprivation by osiedle: (a) registered unemployment and (b) social-
  ##   assistance poverty, each the city's standardised score averaged over the
  ##   2019 and 2021 monitoring rounds and signed so that higher values denote
  ##   greater deprivation; (c) their mean, the deprivation score.
  ## Row 2, people and catchments: (d) registered residents at geocoded addresses;
  ##   in a window over Winogrady, where an estate of high deprivation (Nowe
  ##   Winogrady Wschód) borders one of low deprivation (Nowe Winogrady Północ):
  ##   (e) the OpenStreetMap pedestrian network and the stop poles; (f) the
  ##   mutually exclusive walking catchments (each location assigned to the
  ##   complex whose nearest pole it reaches soonest on foot, up to 400 m).
  ## Row 3, combination: (g) residents in the window, coloured by the score of
  ##   their osiedle; (h) catchment deprivation, the population-weighted mean of
  ##   those scores; (i) catchment deprivation of every stop complex.
  ## Row 4, municipal housing in the square holding the most municipal residents.
  ## Output: output/figures/Fig_construction.png, Fig_construction_full.png (copied to manuscript/figures)

  suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  s1 <- file.path(out_dir, "output/stage1"); s2 <- file.path(out_dir, "output/stage2")
  fig <- file.path(out_dir, "output/figures"); ms_fig <- file.path(out_dir, "manuscript/figures")
  gtfs <- GTFS_DIR
  surf <- "#ffffff"; ink <- "#2b2a27"; ink2 <- "#5f5e5a"

  ## ── deprivation components (as in script 20) ────────────────────────────
  names_lp <- fread(file.path(out_dir, "data/diagnoza2021_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, osiedle)]
  g <- merge(fread(file.path(out_dir, "data/gpr_indicators_ABC_osiedla.csv")), names_lp, by = "lp")
  bz <- rbind(fread(file.path(out_dir, "data/diagnoza2019_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, round = 2019L, B_z = indicator_B_std_zscore)],
              fread(file.path(out_dir, "data/diagnoza2021_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, round = 2021L, B_z = indicator_B_std_zscore)])
  g <- merge(g, bz, by = c("lp", "round"))
  D <- g[round %in% c(2019L, 2021L), .(unemp = mean(-A_z), poverty = mean(-B_z)), by = osiedle][, dep := (unemp + poverty) / 2]
  P <- fread(file.path(s1, "s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, dep_check = dep)]
  chk <- merge(D, P, by = "osiedle"); stopifnot(nrow(chk) == 42, max(abs(chk$dep - chk$dep_check)) < 1e-9)
  cat("Deprivation score reproduced for all 42 osiedla\n")
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  osi <- merge(osi[, c("name")], D, by.x = "name", by.y = "osiedle")
  city <- st_union(osi)

  ## ── residents, catchments, poles ────────────────────────────────────────
  pp <- st_transform(st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE), 2180)
  pp$pop <- pp[[grep("ZAMELDOWANYCH", names(pp), value = TRUE)[1]]]
  pp <- pp[pp$pop > 0, "pop"]
  pp$osiedle <- osi$name[st_nearest_feature(pp, osi)]
  pp$dep <- osi$dep[match(pp$osiedle, osi$name)]
  cat(sprintf("Addresses with residents: %d; residents: %d\n", nrow(pp), sum(pp$pop)))
  vc <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[, .(stop_name, residents, cdep = dep)]
  vc <- merge(vc, cv, by = "stop_name", all.x = TRUE)
  vc$cdep[!(vc$residents > 0) | is.na(vc$residents)] <- NA
  cat(sprintf("Catchments: %d; with residents: %d; residents covered: %.1f%%\n", nrow(vc), sum(!is.na(vc$cdep)), 100 * sum(cv$residents, na.rm = TRUE) / sum(pp$pop)))
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  st <- fread(file.path(gtfs, "stops.txt"), encoding = "UTF-8")
  poles <- st_transform(st_as_sf(st, coords = c("stop_lon", "stop_lat"), crs = 4326), 2180)
  poles <- poles[poles$stop_name %in% C$stop_name & lengths(st_intersects(poles, city)) > 0, ]

  ## ── municipal housing (as in script 37) ─────────────────────────────────
  ap_crs <- st_crs(st_read(file.path(out_dir, "output/geopoz_address_points.gpkg"), quiet = TRUE))
  mb <- fread(file.path(s2, "municipal_buildings.csv"), encoding = "UTF-8")[pop > 0]
  mb <- st_transform(st_as_sf(mb, coords = c("x", "y"), crs = ap_crs), 2180)
  mc <- fread(file.path(s2, "municipal_catchments_network.csv"), encoding = "UTF-8")[residents > 0, .(stop_name, osiedle, muni_share)]
  mc <- merge(mc, D[, .(osiedle, os_dep = dep)], by = "osiedle")
  vc <- merge(vc, mc[, .(stop_name, muni_share)], by = "stop_name", all.x = TRUE)
  cat(sprintf("Municipal buildings with residents: %d; residents: %d; catchments with any: %d of %d; Spearman with catchment deprivation: %.2f\n",
              nrow(mb), sum(mb$pop), sum(mc$muni_share > 0), nrow(mc), cor(vc$muni_share, vc$cdep, method = "spearman", use = "complete.obs")))

  ## ── windows: 1.8 km squares centred on a stop complex, and what lies in them ──
  w0 <- readRDS(file.path(out_dir, "data/raw/osm/poznan_walk_ways.rds"))
  for (k in c("foot", "access")) if (!k %in% names(w0)) w0[[k]] <- NA_character_
  w0 <- w0[!(w0$foot %in% "no") & !(w0$access %in% c("no", "private") & !(w0$foot %in% c("yes", "designated", "permissive"))), ]
  w0 <- st_geometry(st_transform(w0, 2180))
  osi_lines <- st_line_merge(st_union(st_cast(st_geometry(osi), "MULTILINESTRING")))   # shared edges once, so dashes stay dashed
  window <- function(centre, lab_min = 1.5e5) {
    win <- st_as_sfc(st_bbox(st_buffer(st_centroid(st_geometry(C[C$stop_name == centre, ])), 900)))
    clip <- function(x) suppressWarnings(st_intersection(x, win))
    vw <- clip(vc); o <- clip(osi)
    ol <- suppressWarnings(st_point_on_surface(o)); ol <- ol[as.numeric(st_area(o))[match(ol$name, o$name)] > lab_min, ]
    xy <- st_coordinates(ol); bb <- st_bbox(win); m <- 330                # keep osiedle names inside the window
    xy[, 1] <- pmin(pmax(xy[, 1], bb[["xmin"]] + m), bb[["xmax"]] - m); xy[, 2] <- pmin(pmax(xy[, 2], bb[["ymin"]] + m), bb[["ymax"]] - m)
    st_geometry(ol) <- st_sfc(lapply(seq_len(nrow(xy)), function(i) st_point(xy[i, 1:2])), crs = st_crs(ol))
    W <- list(win = win, b = st_bbox(win), line = st_cast(st_geometry(win), "LINESTRING"), ways = clip(w0[lengths(st_intersects(w0, win)) > 0]),
              osi = clip(osi_lines), vc = vw, poles = poles[lengths(st_intersects(poles, win)) > 0, ], pp = pp[lengths(st_intersects(pp, win)) > 0, ],
              lab = suppressWarnings(st_point_on_surface(vw[!is.na(vw$cdep) & as.numeric(st_area(vw)) > 4e4, ])), os_lab = ol,
              mb = mb[lengths(st_intersects(mb, win)) > 0, ])
    cat(sprintf("Window around %s: %d catchments (%d with residents), %d poles, %d residents, %d municipal buildings with %d residents; osiedla: %s\n",
                centre, nrow(vw), sum(!is.na(vw$cdep)), nrow(W$poles), sum(W$pp$pop), nrow(W$mb), sum(W$mb$pop),
                paste(sprintf("%s %.2f", ol$name, ol$dep), collapse = "; ")))
    W
  }
  ## full version: Winogrady (W1) and the square holding the most municipal residents (W2); main version: one window (W3)
  cc <- data.table(stop_name = C$stop_name, st_coordinates(st_centroid(st_geometry(C))))
  mxy <- data.table(st_coordinates(mb), pop = mb$pop)
  cc[, mres := sapply(seq_len(.N), function(i) mxy[abs(X - cc$X[i]) <= 900 & abs(Y - cc$Y[i]) <= 900, sum(pop)])]
  W1 <- window("Os. Kosmonautów"); W2 <- window(cc[which.max(mres), stop_name], lab_min = 2e5); W3 <- window("Niezłomnych")

  ## ── panels ───────────────────────────────────────────────────────────────
  L <- ceiling(10 * max(abs(c(D$unemp, D$poverty)))) / 10
  dscale <- scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-L, L), na.value = "#dcdad2",
                                 name = "Deprivation (z; + = more deprived)", guide = guide_colourbar(barwidth = 12, barheight = 0.5, title.position = "top"), aesthetics = c("fill", "colour"))
  mscale <- scale_fill_gradient(low = "#f0efec", high = "#4a3aa7", limits = c(0, 0.4), oob = scales::squish, na.value = "#dcdad2", labels = function(x) paste0(100 * x, "%"),
                                name = "Residents in municipal buildings", guide = guide_colourbar(barwidth = 10, barheight = 0.5, title.position = "top"))
  mth <- theme_void(base_size = 9.5) + theme(plot.background = element_rect(fill = surf, colour = NA), legend.position = "bottom",
    legend.title = element_text(colour = ink, size = 8.5), legend.text = element_text(colour = ink2, size = 8), plot.tag = element_text(face = "bold"),
    plot.subtitle = element_text(colour = ink, size = 8.5, hjust = 0.5, margin = margin(b = 3)))
  sscale <- scale_size_area(max_size = 3.2, limits = c(0, max(pp$pop)), breaks = c(10, 100, 500), name = "Residents per address (g)", guide = guide_legend(title.position = "top", override.aes = list(fill = "#8f8d86")))
  fmt1 <- function(x) { v <- sprintf("%.1f", x); v[v == "-0.0"] <- "0.0"; sub("-", "−", v) }
  city_map <- function(v, sub) ggplot(osi) + geom_sf(aes(fill = .data[[v]]), colour = "#ffffff", linewidth = 0.25) + dscale + labs(subtitle = sub) + mth
  ## zoomed panels, for any window
  zoom <- function(W, p) p + coord_sf(xlim = W$b[c(1, 3)], ylim = W$b[c(2, 4)], expand = FALSE, datum = NA) +
    theme(panel.border = element_rect(fill = NA, colour = ink, linewidth = 0.5))
  os_lines <- function(W) geom_sf(data = W$osi, colour = ink, linewidth = 0.55, linetype = "22")
  os_names <- function(W, size = 2.1) geom_sf_label(data = W$os_lab, aes(label = gsub(" ", "\n", name)), size = size, colour = ink, lineheight = 0.85, fontface = "italic",
                                                     fill = "#ffffff", alpha = 0.75, linewidth = 0, label.padding = unit(0.06, "lines"))
  pole_pts <- function(W) list(geom_sf(data = W$poles, aes(shape = "Stops (poles in b and d, complexes in e)"), fill = ink, colour = "#ffffff", size = 1.5, stroke = 0.4),
                               scale_shape_manual(values = c("Stops (poles in b and d, complexes in e)" = 21), name = NULL))
  p_net <- function(W) zoom(W, ggplot() + geom_sf(data = W$ways, colour = "#8f8d86", linewidth = 0.22) + os_lines(W) + pole_pts(W) + os_names(W) +
    labs(subtitle = "Pedestrian network and stop poles") + mth)
  p_catch <- function(W) zoom(W, ggplot() + geom_sf(data = W$vc, fill = "#dcdad2", colour = "#ffffff", linewidth = 1) +
    geom_sf(data = W$ways, colour = "#b8b6ae", linewidth = 0.15) + os_lines(W) + pole_pts(W) +
    labs(subtitle = "Walking catchments (nearest pole, up to 400 m)") + mth)
  p_res <- function(W) zoom(W, ggplot() + geom_sf(data = W$vc, fill = "#f6f5f2", colour = "#c9c7c0", linewidth = 0.3) + os_lines(W) +
    geom_sf(data = W$pp[order(W$pp$pop), ], aes(size = pop, fill = dep), colour = "#9d9b94", stroke = 0.15, shape = 21) + sscale + dscale +
    labs(subtitle = "Residents, coloured by their neighbourhood's score") + mth)
  p_cdep <- function(W, lab = 2) zoom(W, ggplot() + geom_sf(data = W$vc, aes(fill = cdep), colour = "#ffffff", linewidth = 0.5) + dscale + os_lines(W) + pole_pts(W) +
    geom_sf_label(data = W$lab, aes(label = fmt1(cdep)), size = lab, colour = ink, fill = "#ffffff", alpha = 0.8, linewidth = 0, label.padding = unit(0.08, "lines")) +
    labs(subtitle = "Catchment deprivation: residents' mean") + mth)
  p_muni <- function(W, lab = 2.1) zoom(W, ggplot() + geom_sf(data = W$vc, aes(fill = muni_share), colour = "#ffffff", linewidth = 0.5) + mscale + os_lines(W) +
    geom_sf(data = W$mb, aes(size = pop), shape = 22, fill = ink, colour = "#ffffff", stroke = 0.2) + scale_size_area(max_size = 2.4, guide = "none") + os_names(W, lab) +
    labs(subtitle = "Municipal share of catchment residents") + mth)
  all_stops <- function() ggplot() + geom_sf(data = city, fill = "#f3f2ee", colour = "#c3c2b7", linewidth = 0.3) +
    geom_sf(data = vc[!is.na(vc$cdep), ], aes(fill = cdep), colour = NA) + geom_sf(data = vc[is.na(vc$cdep), ], fill = "#c9c7c0", colour = NA) +
    geom_sf(data = st_cast(st_geometry(osi), "MULTILINESTRING"), colour = "#ffffff", linewidth = 0.25) + dscale

  ## full version, twelve panels (supplement)
  pa <- city_map("unemp", "Registered unemployment")
  pb <- city_map("poverty", "Social-assistance poverty")
  pc <- city_map("dep", "Deprivation score: mean of a and b") + geom_sf(data = W1$line, colour = ink, linewidth = 0.5)
  pd <- ggplot() + geom_sf(data = city, fill = "#f3f2ee", colour = "#c3c2b7", linewidth = 0.3) +
    geom_sf(data = st_cast(st_geometry(osi), "MULTILINESTRING"), colour = "#ffffff", linewidth = 0.3) +
    geom_sf(data = pp, aes(size = pop), colour = ink, alpha = 0.3, stroke = 0, shape = 16) +
    scale_size_area(max_size = 1.3, guide = "none") +
    geom_sf(data = W1$line, colour = ink, linewidth = 0.5) + labs(subtitle = "Registered residents at addresses") + mth
  pe <- p_net(W1); pf <- p_catch(W1); pg <- p_res(W1); ph <- p_cdep(W1)
  pi <- all_stops() + geom_sf(data = W1$line, colour = ink, linewidth = 0.5) + labs(subtitle = "Catchment deprivation of every stop") + mth
  pj <- ggplot() + geom_sf(data = city, fill = "#f3f2ee", colour = "#c3c2b7", linewidth = 0.3) +
    geom_sf(data = st_cast(st_geometry(osi), "MULTILINESTRING"), colour = "#ffffff", linewidth = 0.3) +
    geom_sf(data = mb[order(mb$pop), ], aes(size = pop), shape = 21, fill = "#4a3aa7", colour = "#ffffff", stroke = 0.15, alpha = 0.85) +
    scale_size_area(max_size = 2.2, guide = "none") + geom_sf(data = W2$line, colour = ink, linewidth = 0.5) +
    labs(subtitle = "Municipal housing: buildings with residents") + mth
  pk <- p_muni(W2)
  pl <- ggplot(mc, aes(os_dep, 100 * muni_share)) + geom_hline(yintercept = 0, colour = "#c3c2b7") +
    geom_point(shape = 16, size = 1.1, alpha = 0.45, colour = "#4a3aa7", position = position_jitter(width = 0.03, height = 0, seed = 1)) +
    scale_y_sqrt(breaks = c(0, 1, 5, 10, 25, 50, 100)) + scale_x_continuous(breaks = -1:2) +
    labs(x = "Deprivation score of the catchment's main neighbourhood", y = "Residents in municipal\nbuildings (%, square-root scale)", subtitle = "Variation within neighbourhoods") +
    theme_minimal(base_size = 9.5) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(colour = ink, size = 8.5, hjust = 0.5, margin = margin(b = 3)),
                                           axis.title = element_text(colour = ink2, size = 8), axis.text = element_text(colour = ink2, size = 7.5), plot.tag = element_text(face = "bold"))
  f <- (pa + pb + pc) / (pd + pe + pf) / (pg + ph + pi) / (pj + pk + pl) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom", plot.background = element_rect(fill = surf, colour = NA)) &
    scale_shape_manual(values = c("Stops (poles in b and d, complexes in e)" = 21), labels = "Stop poles", name = NULL)
  ggsave(file.path(fig, "Fig_construction_full.png"), f, width = 11, height = 16.2, dpi = 300, bg = surf)

  ## main version, four panels: (a) the score; in one window (W3), (b) the catchments and (c) catchment deprivation;
  ## (d) catchment deprivation of every stop with the rings of distance and the sectors that define location.
  ## Residents and municipal housing are shown in the full version.
  source(file.path(REPO, "R/map_style.R"))
  mark <- list(geom_sf(data = W3$line, colour = ink, linewidth = 0.5),
               annotate("label", x = W3$b[["xmax"]] + 250, y = W3$b[["ymax"]], label = "b, c", hjust = 0, vjust = 1, size = 4.4, colour = ink, fill = "#ffffff",
                        linewidth = 0, label.padding = unit(0.1, "lines")))
  qa <- city_map("dep", "Deprivation score by neighbourhood") + mark
  qb <- p_catch(W3) + os_names(W3, 3.9) + labs(subtitle = "Walking catchments (up to 400 m)")
  qc <- p_cdep(W3, 3.3)
  qd <- all_stops() + ms_guides(1.75) + mark + labs(subtitle = "All stops, with rings and sectors") + mth
  f4 <- (qa + qb) / (qc + qd) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom", plot.background = element_rect(fill = surf, colour = NA), plot.subtitle = element_text(colour = ink, size = 12.5, hjust = 0.5, margin = margin(b = 3)),
          legend.title = element_text(colour = ink, size = 12), legend.text = element_text(colour = ink2, size = 11.5), plot.tag = element_text(face = "bold", size = 13.5)) &
    scale_shape_manual(values = c("Stops (poles in b and d, complexes in e)" = 21), labels = "Stop poles", name = NULL)
  ggsave(file.path(fig, "Fig_construction.png"), f4, width = 7.8, height = 8.4, dpi = 320, bg = surf)
  file.copy(file.path(fig, c("Fig_construction.png", "Fig_construction_full.png")), ms_fig, overwrite = TRUE)
  cat("Saved Fig_construction.png (four panels) and Fig_construction_full.png (twelve panels)\n")
})

## ── 60 · Where the more and less deprived areas lie (supplement); network panel of Figure 3 ────────────────
section("60", "Where the more and less deprived areas lie; the network panel of Figure 3", {
  ## Where the more and less deprived areas lie, and why (manuscript §3.1). The two-panel figure is kept for the
  ## supplement; panel (b), the network, is also saved (panel_network.rds) for Figure 3 (section 69), and the
  ## rings and sectors of panel (a) are drawn on Figure 2d (section 57).
  ## (a) Stop catchments of the least and the most deprived thirds (the comparison
  ##     groups of Table 2), with rings of distance (2, 4, 6 and 8 km) and the eight
  ##     sectors around the Old Market Square.
  ## (b) The transit network on an average Tuesday-Thursday (October 2025 GTFS):
  ##     segments between consecutive stop complexes, width proportional to the
  ##     scheduled departures along them, trams and buses apart; over the compact
  ##     mixed-use (pre-war) blocks of the land-use typology (script 16).
  ## Inputs: gap_thirds_stops.csv (script 56); stop_complexes.gpkg; GTFS;
  ##         output/stage1/cell_tissue_types.gpkg; osiedla boundaries.
  ## Output: output/figures/Fig_deprived_areas_map.png (copied to manuscript/figures)

  suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  s1 <- file.path(out_dir, "output/stage1"); s2 <- file.path(out_dir, "output/stage2")
  fig <- file.path(out_dir, "output/figures"); ms_fig <- file.path(out_dir, "manuscript/figures")
  gtfs <- GTFS_DIR
  surf <- "#ffffff"; ink <- "#2b2a27"; ink2 <- "#5f5e5a"

  osi  <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  city <- st_union(osi)
  C    <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  vc   <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  th   <- fread(file.path(s2, "gap_thirds_stops.csv"), encoding = "UTF-8")
  vc   <- merge(vc, th[, .(stop_name, third)], by = "stop_name", all.x = TRUE)
  vc$third <- factor(fifelse(is.na(vc$third), "other", vc$third), levels = c("least", "middle", "most", "other"))
  cat("Catchments by group:", paste(names(table(vc$third)), table(vc$third), collapse = "; "), "\n")

  ## transit network: consecutive stop complexes on Tuesday-Thursday trips
  st <- fread(file.path(gtfs, "stops.txt"), encoding = "UTF-8")[, .(stop_id, stop_name)]
  tr <- fread(file.path(gtfs, "trips.txt"), encoding = "UTF-8")[service_id == 1, .(trip_id, route_id)]
  rt <- fread(file.path(gtfs, "routes.txt"), encoding = "UTF-8")[, .(route_id, route_short_name, mode = fifelse(route_type == 0, "Tram", "Bus"))]
  tt <- fread(file.path(gtfs, "stop_times.txt"), encoding = "UTF-8")[trip_id %in% tr$trip_id, .(trip_id, stop_id, stop_sequence)]
  tt <- merge(merge(tt, st, by = "stop_id"), merge(tr, rt, by = "route_id"), by = "trip_id")
  tt <- tt[!grepl("^2[0-9][0-9]$", route_short_name)]                       # day lines
  setorder(tt, trip_id, stop_sequence)
  tt[, nxt := shift(stop_name, -1), by = trip_id]
  seg <- tt[!is.na(nxt) & nxt != stop_name & stop_name %in% C$stop_name & nxt %in% C$stop_name]
  seg[, `:=`(a = pmin(stop_name, nxt), b = pmax(stop_name, nxt))]
  SG <- seg[, .(departures = .N), by = .(a, b, mode)]
  xy <- data.table(stop_name = C$stop_name, X = st_coordinates(C)[, 1], Y = st_coordinates(C)[, 2])
  SG <- merge(merge(SG, xy[, .(a = stop_name, xa = X, ya = Y)], by = "a"), xy[, .(b = stop_name, xb = X, yb = Y)], by = "b")
  SG <- SG[sqrt((xa - xb)^2 + (ya - yb)^2) < 3000]                           # drop non-stop runs between distant stops
  segs <- st_sf(SG[, .(mode, departures)], geometry = st_sfc(lapply(seq_len(nrow(SG)), function(i) st_linestring(rbind(c(SG$xa[i], SG$ya[i]), c(SG$xb[i], SG$yb[i])))), crs = 2180))
  segs <- segs[order(segs$departures), ]
  cat(sprintf("Network segments: %d (tram %d, bus %d); departures per segment median %.0f, max %.0f\n", nrow(segs), sum(segs$mode == "Tram"), sum(segs$mode == "Bus"),
              median(segs$departures), max(segs$departures)))

  ## compact mixed-use (pre-war) blocks of the land-use typology
  cells <- st_read(file.path(s1, "cell_tissue_types.gpkg"), quiet = TRUE)
  compact <- st_union(cells[cells$type == 1, ])

  ## ── panels (shared map style, script map_style.R) ───────────────────────
  source(file.path(REPO, "R/map_style.R")); L <- ms_layers(out_dir)
  guides_ <- ms_guides()                                                     # rings of distance and the eight sectors
  pa <- ggplot() + ms_base(L) +
    geom_sf(data = vc[vc$third != "other", ], aes(fill = third), colour = NA) +
    scale_fill_manual(values = c(least = MSTY$least, middle = MSTY$middle, most = MSTY$most),
                      labels = c(least = "Least deprived third", middle = "Middle third", most = "Most deprived third"), name = NULL) +
    ms_water_top(L) + geom_sf(data = st_cast(st_geometry(osi), "MULTILINESTRING"), colour = "#ffffff", linewidth = 0.25) + ms_outline(L) + guides_ +
    ms_labels(L, c("Jeżyce", "Łazarz", "Wilda", "Górczyn", "Piątkowo", "Winogrady", "Rataje", "Chartowo", "Główna", "Krzesiny"), size = 2.3) +
    ms_scalebar(L, pos = c(0.0, 0.02)) + coord_sf(datum = NA) + labs(title = "Stop catchments by deprivation") + theme_ms()
  net_panel <- function(gs, base, title) ggplot() + ms_base(L) +
    geom_sf(data = compact, fill = "#e0c9a6", colour = NA) + ms_water_top(L) +
    geom_sf(data = st_cast(st_geometry(osi), "MULTILINESTRING"), colour = "#ffffff", linewidth = 0.25) + ms_outline(L) +
    geom_sf(data = segs, aes(linewidth = departures, colour = mode), lineend = "round") +
    scale_linewidth_continuous(range = c(0.12, 2.4), breaks = c(100, 300, 600), name = "Weekday departures") +
    scale_colour_manual(values = c(Tram = ink, Bus = "#8f8d86"), name = NULL) + gs +
    coord_sf(datum = NA) + labs(title = title) + theme_ms(base) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.2))) + theme(legend.box = "vertical", legend.spacing.y = unit(0.05, "cm"))
  pb <- net_panel(guides_, 9.5, "Transit network over compact pre-war blocks (beige)")
  f <- (pa | pb) + plot_annotation(tag_levels = "a") & theme(plot.background = element_rect(fill = surf, colour = NA))
  ggsave(file.path(fig, "Fig_deprived_areas_map.png"), f, width = 11, height = 6.4, dpi = 300, bg = surf)
  ## Figure 3a: the same panel with text sized for a figure saved 8.6 in wide
  saveRDS(net_panel(ms_guides(1.6), 13, "Network over pre-war blocks (beige)"), file.path(fig, "panel_network.rds"))
  cat("Saved Fig_deprived_areas_map.png (supplement) and panel_network.rds (Figure 3a)\n")
})

## ── 64 · Figure 8: the gaps and their parts ────────────────────────────────────────────────────────────────
section("64", "Figure 8: the gaps and their parts", {
  ## Figure replacing the 22-row table of the gaps (main text): for every measure, the raw slope on
  ## deprivation (open circle), the parts of it accounted for by location, service at the stop, land use and
  ## population (stacked bars; Gelbach decomposition of script 56, primary specification), and the
  ## remainder with its 95% CR2-Satterthwaite interval (black). The full table, with the medians of the least
  ## and most deprived thirds, moves to the supplement.
  ## Inputs: output/stage2/gap_decomposition.csv, gap_ladder.csv (script 56)
  ## Output: output/figures/Fig_decomposition.png

  suppressPackageStartupMessages({library(data.table); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  GB <- fread(file.path(s2, "gap_decomposition.csv")); RES <- fread(file.path(s2, "gap_ladder.csv"), encoding = "UTF-8")
  M <- unique(RES[step == "full", .(msr, dim, label, least, most, est, lo, hi, p)])
  lev <- c("volume", "tidal", "onset_b", "conc_b", "sat_b", "summer", "ride_km", "dest_dep", "to_centre", "same_osiedle", "eff_dest", "dir_conc", "hub_share", "conn_100",
           "jobs30", "transfers_jobs", "t_centre", "bridge", "irregularity", "dependence", "weather_rain", "night_vs_day")
  M <- M[match(lev, msr)]
  M[, label := sub("own osiedle", "own neighbourhood", label)]   # gap_ladder.csv written before the label change in section 56
  wrap <- c("How much and when" = "How much\nand when", "Routes and transfers" = "Routes and\ntransfers", "Access by the timetable" = "Access by\nthe timetable")
  M[dim %in% names(wrap), dim := wrap[dim]]
  M[, label := factor(label, levels = rev(label))][, dim := factor(dim, levels = unique(dim))]
  P <- merge(GB[part %in% c("location", "service", "landuse", "population")], M[, .(msr, label, dim)], by = "msr")
  P[, part := factor(part, levels = c("location", "service", "landuse", "population"), labels = c("Location", "Service at the stop", "Land use", "Population"))]
  R <- merge(GB[part == "raw", .(msr, raw = est)], M, by = "msr")
  pal <- c("Location" = "#2a78d6", "Service at the stop" = "#eb6834", "Land use" = "#1baf7a", "Population" = "#a9a7a0")
  xr <- range(c(P[, .(s = sum(pmax(est, 0))), by = msr]$s, P[, .(s = sum(pmin(est, 0))), by = msr]$s, R$raw, R$lo, R$hi)) * 1.05
  g <- ggplot() + geom_vline(xintercept = 0, colour = "#9a9893", linewidth = 0.35) +
    geom_col(data = P, aes(est, label, fill = part), width = 0.62, colour = "#ffffff", linewidth = 0.25) +
    geom_errorbarh(data = R, aes(xmin = lo, xmax = hi, y = label), height = 0, linewidth = 0.7, colour = "#1f1e1c") +
    geom_point(data = R, aes(est, label, shape = "Remainder (95% interval)"), size = 2.7, colour = "#1f1e1c", fill = "#1f1e1c") +
    geom_point(data = R, aes(raw, label, shape = "Raw gap"), size = 2.9, colour = "#1f1e1c", fill = "#ffffff", stroke = 0.6) +
    scale_shape_manual(values = c("Raw gap" = 21, "Remainder (95% interval)" = 18), name = NULL) +
    scale_fill_manual(values = pal, name = "Accounted for by") +
    scale_x_continuous(limits = xr, breaks = seq(-1, 1, 0.2)[seq(-1, 1, 0.2) >= xr[1] - 0.05 & seq(-1, 1, 0.2) <= xr[2] + 0.05]) +
    facet_grid(dim ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Slope on deprivation (SD per SD)", y = NULL) +
    guides(fill = guide_legend(order = 1, nrow = 2), shape = guide_legend(order = 2)) +
    theme_minimal(base_size = 12.5) +
    theme(legend.position = "top", legend.location = "plot", legend.box = "vertical", plot.margin = margin(6, 14, 6, 6), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1), axis.title.x = element_text(hjust = 0.5, size = 12), axis.text.y = element_text(size = 12, colour = "#2b2b28"),
          panel.spacing.y = unit(0.5, "lines"), plot.background = element_rect(fill = "#ffffff", colour = NA))
  ggsave(file.path(fig, "Fig_decomposition.png"), g, width = 8, height = 10.2, dpi = 280, bg = "#ffffff")
  cat("Saved Fig_decomposition.png\n")
})

## ── 71 · Deprivation of destinations (supplement) ──────────────────────────────────────────────────────────
section("71", "Deprivation of destinations (supplement)", {
  ## Where riders are carried (main text): deprivation of the destinations of each origin catchment's riders,
  ## weekdays, (a) offered by the network (riders assigned to every line and direction at their stop in
  ## proportion to its departures), (b) reached by riders, (c) reached minus offered; (d) destination against
  ## origin deprivation, with least-squares lines for the offered and the reached destinations.
  ## Input: dest_mech_profiles.csv (script 46). Output: output/figures/Fig_destinations_main.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  source(file.path(REPO, "R/map_style.R"))
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  L <- ms_layers(out_dir)
  CP <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  DP <- fread(file.path(s2, "dest_mech_profiles.csv"), encoding = "UTF-8")[transfers == "none" & window == "weekday all day" & comp %in% c("observed", "network"),
          .(stop_name, comp, dest_dep, dep)]
  W <- dcast(DP, stop_name + dep ~ comp, value.var = "dest_dep")[, diff := observed - network]
  lim <- max(abs(c(W$network, W$observed)), na.rm = TRUE)
  mp <- function(v, title, scale) { P <- merge(CP[, "stop_name"], W[, c("stop_name", v), with = FALSE], by = "stop_name"); names(P)[2] <- "val"
    ggplot() + ms_base(L) + geom_sf(data = P, aes(fill = val), colour = NA) + scale + ms_water_top(L) + ms_outline(L) +
      ms_labels(L, c("Piątkowo", "Winogrady", "Jeżyce", "Łazarz", "Wilda", "Rataje", "Chartowo", "Główna"), size = 2.2) + ms_centre(L) + ms_coord(L) +
      labs(title = title) + theme_ms() }
  sc_dep <- DEPSCALE(lim, name = "Deprivation of destinations (z)")
  pa <- mp("network", "Offered by the network", sc_dep) + ms_scalebar(L)
  pb <- mp("observed", "Reached by riders", sc_dep)
  dl <- max(abs(W$diff), na.rm = TRUE)
  pc <- mp("diff", "Reached minus offered", scale_fill_gradient2(low = "#2f8f83", mid = "#f0efec", high = "#a3527d", midpoint = 0, limits = c(-dl, dl),
                                                                 name = "Reached minus offered (z)"))
  S <- melt(W, id.vars = c("stop_name", "dep"), measure.vars = c("network", "observed"), variable.name = "comp", value.name = "dd")
  S[, comp := factor(comp, levels = c("network", "observed"), labels = c("Offered by the network", "Reached by riders"))]
  co <- S[, as.list(coef(lm(dd ~ dep))), by = comp]; setnames(co, c("comp", "a", "b"))
  pd <- ggplot(S, aes(dep, dd, colour = comp)) + geom_hline(yintercept = 0, colour = "#d9d7d0") + geom_vline(xintercept = 0, colour = "#d9d7d0") +
    geom_point(size = 0.9, alpha = 0.55, stroke = 0) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, formula = y ~ x) +
    scale_colour_manual(values = c("Offered by the network" = "#8f8d86", "Reached by riders" = "#4a3aa7"), name = NULL) +
    labs(title = "Destination against origin deprivation", x = "Deprivation of the origin catchment (z)", y = "Deprivation of destinations (z)") +
    theme_minimal(base_size = 9) + theme(plot.title = element_text(size = 9.5, hjust = 0.5, colour = MSTY$ink), legend.position = "bottom", legend.direction = "vertical",
                                         panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#ffffff", colour = NA), plot.tag = element_text(face = "bold", size = 11.5),
                                         axis.title = element_text(size = 8, colour = MSTY$ink))
  g <- ((pa | pb) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) / (pc | pd) + plot_annotation(tag_levels = "a")
  ggsave(file.path(fig, "Fig_destinations_main.png"), g, width = 10.5, height = 10.5, dpi = 250, bg = "#ffffff")
  cat(sprintf("Origins: %d; slopes of destination on origin deprivation (raw, unstandardised): offered %.2f, reached %.2f\n", nrow(W), co$b[1], co$b[2]))
  cat("Saved Fig_destinations_main.png\n")
})

## ── 70 · Route networks of each third (supplementary S7) ───────────────────────────────────────────────────
section("70", "Route networks of each third (supplementary S7)", {
  ## Route networks of riders from the least and the most deprived thirds (main text), Monday-Thursday:
  ## (a) least deprived third, (b) most deprived third: share of the group's riders passing each segment;
  ## (c) most minus least: segments used relatively more by riders from the most deprived third (red) or the
  ## least deprived third (blue), i.e. where a closure would fall more on one group than on the other.
  ## Shared map style (map_style.R). Weekend panels: supplementary figure of script 53.
  ## Input: output/stage2/route_segments.csv (script 53)
  ## Output: output/figures/Fig_route_networks_main.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  source(file.path(REPO, "R/map_style.R"))
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  L <- ms_layers(out_dir)
  SG <- fread(file.path(s2, "route_segments.csv"), encoding = "UTF-8")[win == "Monday-Thursday"]
  mx <- quantile(SG$share, 0.995)
  ## zoom to the extent of the segments carrying at least 0.2% of either group, with a 1.2 km margin
  ZB <- SG[share >= 0.002, .(x0 = min(pmin(X1, X2)) - 1200, x1 = max(pmax(X1, X2)) + 1200, y0 = min(pmin(Y1, Y2)) - 1200, y1 = max(pmax(Y1, Y2)) + 1200)]
  zoom <- coord_sf(datum = NA, expand = FALSE, xlim = c(ZB$x0, ZB$x1), ylim = c(ZB$y0, ZB$y1))
  Lz <- L; Lz$bbox <- st_bbox(c(xmin = ZB$x0, ymin = ZB$y0, xmax = ZB$x1, ymax = ZB$y1), crs = st_crs(2180))
  seg_layer <- function(d, col) list(geom_segment(data = d[share >= 0.002][order(share)], aes(X1, Y1, xend = X2, yend = Y2, linewidth = pmin(share, mx), alpha = pmin(share, mx)),
                                                  colour = col, lineend = "round"),
                                     scale_linewidth(range = c(0.15, 3.0), limits = c(0.002, mx), breaks = c(0.01, 0.03, 0.06), labels = scales::percent,
                                                     name = "Share of the group's riders passing the segment"),
                                     scale_alpha(range = c(0.3, 0.95), limits = c(0.002, mx), guide = "none"))
  LAB <- c("Piątkowo", "Winogrady", "Jeżyce", "Ogrody", "Łazarz", "Wilda", "Rataje", "Chartowo", "Śródka", "Górczyn")
  pa <- ggplot() + ms_base(L) + ms_water_top(L) + ms_outline(L) + seg_layer(SG[third == "least"], MSTY$least) + ms_labels(L, LAB, size = 2.3) + ms_centre(L) +
    ms_scalebar(Lz) + zoom + labs(title = "Least deprived third") + theme_ms()
  pb <- ggplot() + ms_base(L) + ms_water_top(L) + ms_outline(L) + seg_layer(SG[third == "most"], MSTY$most) + ms_labels(L, LAB, size = 2.3) + ms_centre(L) +
    zoom + labs(title = "Most deprived third") + theme_ms() + guides(linewidth = "none")
  DF <- dcast(SG, k1 + k2 + X1 + Y1 + X2 + Y2 ~ third, value.var = "share", fill = 0)[, diff := most - least]
  m <- quantile(abs(DF$diff), 0.995)
  pc <- ggplot() + ms_base(L) + ms_water_top(L) + ms_outline(L) +
    geom_segment(data = DF[abs(diff) >= 0.002][order(abs(diff))], aes(X1, Y1, xend = X2, yend = Y2, colour = pmax(pmin(diff, m), -m), linewidth = pmin(abs(diff), m)), lineend = "round") +
    scale_colour_gradient2(low = MSTY$least, mid = "#d9d7d0", high = MSTY$most, midpoint = 0, limits = c(-m, m), labels = function(v) paste0(ifelse(v > 0, "+", ""), round(100 * v), " pp"),
                           name = "Most minus least deprived third") +
    scale_linewidth(range = c(0.25, 3), guide = "none") + ms_labels(L, LAB, size = 2.3) + ms_centre(L) + zoom +
    labs(title = "Used relatively more by each group") + theme_ms()
  g <- (pa | pb | pc) + plot_annotation(tag_levels = "a") + plot_layout(guides = "keep")
  ggsave(file.path(fig, "Fig_route_networks_main.png"), g, width = 14, height = 6.4, dpi = 250, bg = "#ffffff")
  ## corridor summary for the text: segments where one group's share exceeds the other's by at least 1 point
  cat(sprintf("Segments with |most - least| >= 1 pp: %d favouring the least deprived third, %d favouring the most deprived third\n", DF[diff <= -0.01, .N], DF[diff >= 0.01, .N]))
  cat("Saved Fig_route_networks_main.png\n")
})

## ── 69 · Figure 3: what the network offers ─────────────────────────────────────────────────────────────────
section("69", "Figure 3: what the network offers", {
  ## What the network offers (main text): (a) the transit network over the compact pre-war blocks (drawn by
  ## section 60); (b) bivariate map of stop catchments, thirds of deprivation x thirds of morning-peak access
  ## (share of the city's jobs within 30 minutes); the legend gives the share of residents in each cell;
  ## (c) concentration curves of service, access and ride length over residents ranked from the most to the
  ## least deprived catchment (script 65 computes the indices).
  ## Bivariate colours: bilinear mix in CIELab of four corners (warm grey: less deprived, poor access; red:
  ## more deprived, poor access; blue: less deprived, good access; dark violet: more deprived, good access).
  ## Inputs: access_origin.csv (script 54), gap_thirds_stops.csv and gap_stops.rds (56), catchments (39)
  ## Output: output/figures/Fig_access_bivariate.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  source(file.path(REPO, "R/map_style.R"))
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  L <- ms_layers(out_dir)
  CP <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  TH <- fread(file.path(s2, "gap_thirds_stops.csv"), encoding = "UTF-8")
  AC <- fread(file.path(s2, "access_origin.csv"), encoding = "UTF-8")[window == "Weekday 07-09 (morning peak)", .(stop_name, jobs30 = 100 * jobs30, lres)]
  D <- merge(TH, AC, by = "stop_name")
  qa <- quantile(D$jobs30, c(1/3, 2/3))
  D[, `:=`(i = match(third, c("least", "middle", "most")), j = fcase(jobs30 <= qa[1], 1L, jobs30 <= qa[2], 2L, default = 3L), residents = expm1(lres))]

  ## bivariate palette
  lab <- function(h) grDevices::convertColor(t(grDevices::col2rgb(h)) / 255, from = "sRGB", to = "Lab")
  cn <- list(c11 = lab("#cfcbc2"), c31 = lab("#d8434b"), c13 = lab("#3a7bc8"), c33 = lab("#47265a"))
  BV <- CJ(i = 1:3, j = 1:3)[, col := mapply(function(i, j) { t <- (i - 1) / 2; s <- (j - 1) / 2
    v <- (1 - t) * (1 - s) * cn$c11 + t * (1 - s) * cn$c31 + (1 - t) * s * cn$c13 + t * s * cn$c33
    grDevices::rgb(pmin(pmax(grDevices::convertColor(v, from = "Lab", to = "sRGB"), 0), 1)) }, i, j)]
  D <- merge(D, BV, by = c("i", "j"))
  P <- merge(CP[, "stop_name"], D[, .(stop_name, col)], by = "stop_name")

  pm <- ggplot() + ms_base(L) + geom_sf(data = P, aes(fill = col), colour = NA) + scale_fill_identity() + ms_water_top(L) + ms_outline(L) +
    ms_labels(L, c("Stare Miasto", "Jeżyce", "Łazarz", "Wilda", "Piątkowo", "Winogrady", "Rataje", "Główna", "Krzesiny", "Głuszyna"), size = 3.9) +
    ms_centre(L) + ms_scalebar(L, km = 4, text = 1.6) + ms_coord(L) + theme_ms(13) + labs(title = "Deprivation and access to jobs")
  ## legend: 3 x 3 tiles with the share of residents in each cell
  LG <- merge(BV, D[, .(res = sum(residents)), by = .(i, j)], by = c("i", "j"), all.x = TRUE)[is.na(res), res := 0][, share := 100 * res / sum(res)]
  pl <- ggplot(LG, aes(j, i)) + geom_tile(aes(fill = col), colour = "#ffffff", linewidth = 0.8) + scale_fill_identity() +
    geom_text(aes(label = sprintf("%.0f%%", share), colour = ifelse(i + j >= 5 | (i == 3 & j == 1), "#ffffff", MSTY$ink)), size = 3.3) + scale_colour_identity() +
    scale_x_continuous(breaks = 1:3, labels = c("low", "", "high"), expand = c(0, 0)) + scale_y_continuous(breaks = 1:3, labels = c("least", "", "most"), expand = c(0, 0)) +
    labs(x = "Jobs within 30 min \u2192", y = "Deprivation \u2192") + coord_equal() +
    theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), axis.title = element_text(size = 11, colour = MSTY$ink), axis.text = element_text(size = 10, colour = MSTY$muted),
                                            plot.background = element_rect(fill = NA, colour = NA))
  pa <- pm                                                                  # the 3 x 3 legend (pl) sits beside the curves, hiding no catchments

  ## concentration curves (needs gap_stops.rds of script 56)
  X <- readRDS(file.path(s2, "gap_stops.rds")); X[, `:=`(residents = expm1(lres), departures = exp(l_dep))]
  X[, dep_per_1000 := 1000 * departures / residents]
  CV <- c(departures = "Departures at the stop", dep_per_1000 = "Departures per 1,000 residents", jobs30 = "Jobs within 30 min", ride_km = "Ride length")
  CC <- rbindlist(lapply(names(CV), function(v) { Z <- X[is.finite(get(v)) & residents > 0][order(-dep)]
    data.table(var = CV[[v]], x = c(0, cumsum(Z$residents) / sum(Z$residents)), y = c(0, cumsum(Z$residents * Z[[v]]) / sum(Z$residents * Z[[v]]))) }))
  CC[, var := factor(var, levels = CV)]
  pal <- c("Departures at the stop" = "#eb6834", "Departures per 1,000 residents" = "#8f8d86", "Jobs within 30 min" = MSTY$least, "Ride length" = "#1baf7a")
  ends <- CC[, .SD[which.min(abs(x - 0.33))], by = var]
  pb <- ggplot(CC, aes(x, y, colour = var)) + geom_abline(slope = 1, intercept = 0, colour = "#b9b7b0", linetype = "dashed", linewidth = 0.4) +
    geom_line(linewidth = 0.8) + scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous(labels = scales::percent, expand = c(0.01, 0)) + scale_y_continuous(labels = scales::percent, expand = c(0.01, 0)) + coord_equal() +
    labs(title = "Concentration over residents, most deprived first", x = "Residents, cumulative, from the most deprived catchment", y = "Cumulative share") +
    theme_minimal(base_size = 13) + theme(plot.title = element_text(size = 13, hjust = 0.5, colour = MSTY$ink), legend.position = c(0.64, 0.2), legend.text = element_text(size = 12),
                                         legend.key.height = unit(0.5, "cm"), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#ffffff", colour = NA),
                                         axis.title = element_text(size = 12.5, colour = MSTY$ink), axis.text = element_text(size = 11.5), plot.tag = element_text(face = "bold", size = 14))
  pn <- readRDS(file.path(fig, "panel_network.rds"))                        # the network over the pre-war blocks (section 60)
  g <- (pn | pa) / (plot_spacer() + pb + pl + plot_layout(widths = c(0.1, 1, 0.55))) + plot_layout(heights = c(1.05, 1)) +
    plot_annotation(tag_levels = list(c("a", "b", "c", ""))) & theme(plot.tag = element_text(face = "bold", size = 14))
  ggsave(file.path(fig, "Fig_access_bivariate.png"), g, width = 8.6, height = 9.6, dpi = 300, bg = "#ffffff")
  cat(sprintf("Access thirds: jobs within 30 min <= %.1f%% (low), <= %.1f%% (middle), above (high)\n", qa[1], qa[2]))
  print(LG[order(-i, j), .(deprivation = c("least", "middle", "most")[i], access = c("low", "middle", "high")[j], residents = round(res), share = round(share, 1))])
  cat("Saved Fig_access_bivariate.png\n")
})

## ── 72 · Graphical abstract ────────────────────────────────────────────────────────────────────────────────
section("72", "Graphical abstract", {
  ## Graphical abstract (Journal of Transport Geography; at least 1328 x 531 px, w x h): (1) where the most and
  ## least deprived stop catchments lie on the radial network; (2) for the headline measures, the slope on
  ## deprivation raw and net of position on the network (gap ladder of script 56, primary specification);
  ## (3) the main findings in short.
  ## Output: output/figures/Graphical_abstract.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
  out_dir <- ROOT
  source(file.path(REPO, "R/map_style.R"))
  s2 <- file.path(out_dir, "output/stage2"); fig <- file.path(out_dir, "output/figures")
  L <- ms_layers(out_dir)
  vc <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  th <- fread(file.path(s2, "gap_thirds_stops.csv"), encoding = "UTF-8")
  vc <- merge(vc, th[third != "middle", .(stop_name, third)], by = "stop_name")
  SG <- fread(file.path(s2, "route_segments.csv"), encoding = "UTF-8")[win == "Monday-Thursday", .(share = sum(share)), by = .(k1, k2, X1, Y1, X2, Y2)]
  p1 <- ggplot() + ms_base(L, borders = FALSE) + geom_sf(data = vc, aes(fill = third), colour = NA) +
    scale_fill_manual(values = c(least = MSTY$least, most = MSTY$most), labels = c(least = "Least deprived", most = "Most deprived"), name = NULL) +
    geom_segment(data = SG[share >= 0.01], aes(X1, Y1, xend = X2, yend = Y2, linewidth = share), colour = "#2b2b28", alpha = 0.55, lineend = "round") +
    scale_linewidth(range = c(0.2, 1.6), guide = "none") + ms_water_top(L) + ms_outline(L) + ms_centre(L) + ms_coord(L) +
    labs(title = "Deprived areas sit at the centre\nof a radial network") + theme_ms(10) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  R <- fread(file.path(s2, "gap_ladder.csv"), encoding = "UTF-8")[step %in% c("raw", "position")]
  K <- c(ride_km = "Ride length", conn_100 = "Connections per ride", jobs30 = "Jobs within 30 min", t_centre = "Time to the centre", dest_dep = "Deprivation of destinations")
  R <- R[msr %in% names(K)][, `:=`(lab = factor(K[msr], levels = rev(K)), step = factor(step, levels = c("raw", "position"), labels = c("Raw gap", "Net of network position")))]
  p2 <- ggplot(R, aes(est, lab, fill = step)) + geom_vline(xintercept = 0, colour = "#9a9893") +
    geom_col(aes(group = factor(step, levels = rev(levels(step)))), position = position_dodge(width = 0.72), width = 0.66) + scale_fill_manual(values = c("Raw gap" = "#8f8d86", "Net of network position" = "#eb6834"), name = NULL) +
    labs(title = "Position on the network accounts\nfor most of the gaps", x = "Most vs least deprived (SD per SD)", y = NULL) +
    theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5, colour = MSTY$ink), legend.position = "bottom", panel.grid.minor = element_blank(),
                                          panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = MSTY$ink, size = 9), plot.background = element_rect(fill = "#ffffff", colour = NA))
  txt <- data.table(y = c(5, 4.1, 3.5, 2.6, 2.0, 1.1, 0.5), t = c("Places, not people", "Shorter rides, fewer transfers and", "better job access follow position", "Destinations: lines from deprived areas run", "through deprived districts; riders narrow it", "Usage beyond position differs little,", "but a coarse measure limits what is ruled out"),
                    f = c("bold", "plain", "plain", "plain", "plain", "plain", "plain"))
  p3 <- ggplot(txt, aes(0, y, label = t, fontface = f)) + geom_text(hjust = 0, size = 3.5, colour = MSTY$ink) + xlim(0, 1) + ylim(0, 5.4) +
    theme_void(base_size = 10) + theme(plot.background = element_rect(fill = "#ffffff", colour = NA))
  g <- (p1 | p2 | p3) + plot_layout(widths = c(1, 1.15, 0.8))
  ggsave(file.path(fig, "Graphical_abstract.png"), g, width = 13.28, height = 5.31, dpi = 200, bg = "#ffffff")
  cat("Saved Graphical_abstract.png (2656 x 1062 px)\n")
})

