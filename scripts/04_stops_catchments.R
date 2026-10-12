## DualPulse | stage 4: Stop complexes, walking catchments and catchment covariates
##
## 607 stop complexes, mutually exclusive 400 m walking catchments, their land use, residents, deprivation
## and municipal housing (main text 2.2, supplementary S2).
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   21      Stop complexes, flags and catchment overlap
##   39      Mutually exclusive 400 m walking catchments
##   15b     Land-use descriptors of the walking catchments
##   34      Catchment covariates
##   37      Municipal housing
##   62      Catchment deprivation by round and catchment-osiedle weights
##
## Run from the repository root: Rscript run.R 4   (one section: Rscript run.R 4:62)

if (!exists("section")) source("R/setup.R")

## ── 21 · Stop complexes, flags and catchment overlap ───────────────────────────────────────────────────────
section("21", "Stop complexes, flags and catchment overlap", {
  ## Stop complexes and catchment overlap (manuscript §2.5).
  ##
  ## Merges the PEKA/ZTM stop poles within city limits into stop complexes by
  ## stop name, measures the spacing of complexes, and compares overlapping
  ## circular catchments (200/300/400/600 m) with mutually exclusive
  ## catchments (each location assigned to its nearest complex, up to 400 m):
  ## share of each catchment covered by others, residents covered, residents
  ## in two or more catchments, catchments per resident, and the rank
  ## correlation of catchment population between nearest-neighbour complexes
  ## (how similar the exposures of neighbouring stops are). Distances here are
  ## straight-line; the primary Stage 2 catchments use walking distance on the
  ## pedestrian network.
  ##
  ## Inputs: October 2025 GTFS (PEKA project), osiedle boundaries, geocoded
  ## address population. Outputs: output/stage2/stop_complexes.gpkg (complexes
  ## and exclusive straight-line catchments), output/stage2/catchment_overlap.csv

  suppressPackageStartupMessages({library(sf); library(data.table)})
  gtfs <- GTFS_DIR
  out_dir <- ROOT
  s2_dir <- file.path(out_dir, "output/stage2"); dir.create(s2_dir, showWarnings = FALSE)
  city <- st_union(st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180))

  ## stop poles in the city, mode served
  st <- fread(file.path(gtfs, "stops.txt"), encoding = "UTF-8")
  tr <- fread(file.path(gtfs, "trips.txt"), select = c("route_id", "trip_id"))
  ro <- fread(file.path(gtfs, "routes.txt"), select = c("route_id", "route_type"))
  stt <- unique(fread(file.path(gtfs, "stop_times.txt"), select = c("trip_id", "stop_id")))
  md <- unique(merge(merge(stt, tr, by = "trip_id"), ro, by = "route_id")[, .(stop_id, route_type)])
  st[, tram := stop_id %in% md[route_type == 0, stop_id]]
  s <- st_transform(st_as_sf(st, coords = c("stop_lon", "stop_lat"), crs = 4326), 2180)
  s <- s[lengths(st_intersects(s, city)) > 0, ]

  ## complexes: poles sharing a stop name
  cx <- as.data.table(cbind(st_drop_geometry(s), st_coordinates(s)))[, .(x = mean(X), y = mean(Y), poles = .N,
          spread = max(sqrt((X - mean(X))^2 + (Y - mean(Y))^2)), tram = any(tram)), by = stop_name]
  ## the AVL/APC records identify stops by cluster code (GTFS stop_code without the pole number): check it matches the name grouping
  pc <- as.data.table(st_drop_geometry(s))[, code := sub("[0-9]+$", "", stop_code)]
  cat("Cluster codes:", uniqueN(pc$code), "| names with >1 code:", pc[, uniqueN(code), by = stop_name][V1 > 1, .N],
      "| codes with >1 name:", pc[, uniqueN(stop_name), by = code][V1 > 1, .N], "\n")

  ## service structure per complex: lines, modes, share of calling trips that start or end there
  st_seq <- fread(file.path(gtfs, "stop_times.txt"), select = c("trip_id", "stop_id", "stop_sequence"))
  st_seq[, `:=`(first = stop_sequence == min(stop_sequence), last = stop_sequence == max(stop_sequence)), by = trip_id]
  calls <- merge(merge(st_seq, tr, by = "trip_id"), ro, by = "route_id")
  calls <- merge(calls, pc[, .(stop_id, stop_name)], by = "stop_id")
  svc <- calls[, .(lines = uniqueN(route_id), tram_bus = any(route_type == 0) & any(route_type == 3),
                   end_share = mean(first | last)), by = stop_name]
  cx <- merge(cx, svc, by = "stop_name", all.x = TRUE)
  cat("Lines per complex:", paste(names(quantile(cx$lines, c(.5, .9, .95), na.rm = TRUE)), round(quantile(cx$lines, c(.5, .9, .95), na.rm = TRUE)), collapse = " "), "\n")
  cat("Share of calling trips that start or end at the complex:", paste(c("50%", "90%", "95%"), round(quantile(cx$end_share, c(.5, .9, .95), na.rm = TRUE), 2), collapse = " "), "\n")
  cat("Stop poles in city:", nrow(s), "| complexes:", nrow(cx), "| tram-served complexes:", sum(cx$tram), "\n")
  cat("Distance of poles from their complex centre (m): median", round(median(cx$spread)), "| 95th pct", round(quantile(cx$spread, .95)),
      "| max", round(max(cx$spread)), "\n")
  cat("Complexes whose poles spread more than 200 m (check for same-name stops):\n"); print(cx[spread > 200, .(stop_name, poles, spread = round(spread))])
  C <- st_as_sf(cx, coords = c("x", "y"), crs = 2180)

  ## spacing
  D <- st_distance(C); diag(D) <- Inf
  nn <- apply(D, 1, which.min); nnd <- as.numeric(apply(D, 1, min))
  cat("\nNearest-neighbour distance between complexes (m):\n"); print(round(quantile(nnd, c(.1, .25, .5, .75, .9))))
  cat("Share with another complex within 200/300/400 m:", round(c(mean(nnd < 200), mean(nnd < 300), mean(nnd < 400)), 2), "\n")
  cat("Median nearest-neighbour distance, tram-served vs bus-only:", round(median(nnd[cx$tram])), round(median(nnd[!cx$tram])), "\n")

  ## residents and a 25 m lattice for area overlap
  p <- st_transform(st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE), 2180)
  pop <- p$LICZBA_ZAMELDOWANYCH_OGÓŁEM
  lat <- st_as_sf(st_make_grid(city, cellsize = 25, what = "centers")); lat <- lat[lengths(st_intersects(lat, city)) > 0, ]
  invert <- function(w, n) { inv <- vector("list", n); for (k in seq_along(w)) for (j in w[[k]]) inv[[j]] <- c(inv[[j]], k); inv }
  circles <- rbindlist(lapply(c(200, 300, 400, 600), function(r) {
    wl <- st_is_within_distance(lat, C, dist = r); nl <- lengths(wl)
    shared <- sapply(invert(wl, nrow(C)), function(ks) if (length(ks)) mean(nl[ks] >= 2) else NA)
    wp <- st_is_within_distance(p, C, dist = r); np <- lengths(wp)
    cpop <- sapply(invert(wp, nrow(C)), function(ks) sum(pop[ks]))
    data.table(catchment = paste0(r, " m circle"), median_share_shared = median(shared, na.rm = TRUE),
               residents_covered = sum(pop[np >= 1]) / sum(pop), residents_in_2plus = sum(pop[np >= 2]) / sum(pop),
               catchments_per_resident = sum(pop * np) / sum(pop), nn_rho_population = cor(cpop, cpop[nn], method = "spearman"))
  }))

  ## mutually exclusive: Voronoi cells of complexes, capped at 400 m
  v <- st_collection_extract(st_voronoi(st_union(C), envelope = st_buffer(city, 2000)))
  v <- v[unlist(st_intersects(C, st_sf(geometry = v)))]
  vc <- st_sf(stop_name = cx$stop_name, geometry = do.call(c, lapply(seq_len(nrow(C)), function(i) st_intersection(v[i], st_buffer(st_geometry(C)[i], 400)))))
  vpop <- sapply(st_intersects(vc, p), function(k) sum(pop[k]))
  np400 <- lengths(st_is_within_distance(p, C, dist = 400))
  excl <- data.table(catchment = "400 m exclusive (nearest complex)", median_share_shared = 0, residents_covered = sum(pop[np400 >= 1]) / sum(pop),
                     residents_in_2plus = 0, catchments_per_resident = sum(pop[np400 >= 1]) / sum(pop), nn_rho_population = cor(vpop, vpop[nn], method = "spearman"))
  res <- rbind(circles, excl)
  cat("\nCatchment overlap (straight-line distances):\n"); print(res, digits = 2)
  cat("Exclusive 400 m catchments: median area", round(median(as.numeric(st_area(vc))) / 1e6, 3), "km2 vs", round(pi * 0.4^2, 3), "km2 for a full circle\n")
  fwrite(res, file.path(s2_dir, "catchment_overlap.csv"))

  ## flags (criteria as in manuscript §2.5): terminal = at least a quarter of calling trips start or end there;
  ## interchange = served by at least 12 lines (top tenth); boundary = within 400 m of the city limit (catchment cut by the data edge)
  cx[, dist_boundary := as.numeric(st_distance(st_geometry(C), st_cast(city, "MULTILINESTRING")))]
  cx[, `:=`(terminal = end_share >= 0.25, interchange = lines >= 12, boundary = dist_boundary <= 400)]
  cat("\nFlags: terminal", sum(cx$terminal, na.rm = TRUE), "| interchange", sum(cx$interchange, na.rm = TRUE), "| boundary", sum(cx$boundary),
      "| any", sum(cx$terminal | cx$interchange | cx$boundary, na.rm = TRUE), "of", nrow(cx), "\n")

  ## residents and osiedla in the exclusive catchments
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  p$osiedle <- osi$name[st_nearest_feature(p, osi)]; p$pop <- pop
  jn <- as.data.table(st_drop_geometry(st_join(p[, c("pop", "osiedle")], vc[, "stop_name"], join = st_within)))[!is.na(stop_name)]
  by_os <- jn[, .(pop = sum(pop)), by = .(stop_name, osiedle)][, share := pop / sum(pop), by = stop_name]
  dom <- by_os[, .(residents = sum(pop), top = max(share), n_osiedla = .N), by = stop_name]
  cx <- merge(cx, dom, by = "stop_name", all.x = TRUE); cx[is.na(residents), residents := 0]
  cat("Complexes with no registered residents in their exclusive catchment:", sum(cx$residents == 0), "of", nrow(cx), "\n")
  cat("Catchments with residents from 2+ osiedla:", round(mean(cx[residents > 0]$n_osiedla > 1), 2),
      "| largest osiedle <90% of residents:", round(mean(cx[residents > 0]$top < 0.9), 2), "\n")
  fwrite(cx, file.path(s2_dir, "stop_complex_attributes.csv"))
  st_write(C, file.path(s2_dir, "stop_complexes.gpkg"), layer = "complexes", delete_dsn = TRUE, quiet = TRUE)
  st_write(vc, file.path(s2_dir, "stop_complexes.gpkg"), layer = "exclusive_400m_straightline", append = TRUE, quiet = TRUE)
  cat("Saved output/stage2/catchment_overlap.csv and stop_complexes.gpkg\n")
})

## ── 39 · Mutually exclusive 400 m walking catchments ───────────────────────────────────────────────────────
section("39", "Mutually exclusive 400 m walking catchments", {
  ## Mutually exclusive walking catchments (manuscript §2.5, primary catchment):
  ## each location is assigned to the stop complex whose nearest pole it can
  ## reach soonest on the pedestrian network, up to 400 m.
  ##
  ## Network: OpenStreetMap ways usable on foot in and around Poznań (Overpass,
  ## lz4.overpass-api.de first -- it answers when the main server is overloaded;
  ## cached in data/raw/osm/): footways, paths, pedestrian streets, steps,
  ## cycleways, tracks, living streets and all roads except motorways; ways
  ## tagged foot=no, or access=no/private without a foot permission, are
  ## dropped. Fallback when no server answers: BDOT10k road axes (OT_SKDR_L,
  ## motorways and expressways excluded) and pedestrian and cycle paths
  ## (OT_SKRP_L). Line ends that stop short of the network are joined to another
  ## line within 5 m. Lines are split into segments of at most 20 m, so that
  ## snapping a location to the nearest node stays close to the line.
  ## Allocation: one multi-source shortest-path search (igraph Dijkstra) from all
  ## stop poles of the 607 complexes at once (a super-source linked to the node
  ## nearest each pole, with the pole's snapping distance as link length). Every
  ## node receives its walking distance to the nearest pole and that pole's
  ## complex (followed back through the shortest-path tree).
  ## Catchments: a 20 m lattice over the city; a lattice point joins the complex
  ## of its nearest network node when the straight access to that node plus the
  ## network distance is at most 400 m; lattice cells are polygonised per complex.
  ## Outputs: layer exclusive_400m_network in output/stage2/stop_complexes.gpkg,
  ##          output/stage2/catchment_network_summary.csv,
  ##          output/figures/Fig_catchments_network_vs_straight.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(osmdata); library(igraph); library(nabor); library(terra); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")
  gtfs <- GTFS_DIR
  t0 <- Sys.time()
  city <- st_union(st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180))

  ## ── pedestrian network: OpenStreetMap, BDOT10k as fallback ───────────────
  walk_values <- c("footway", "path", "pedestrian", "steps", "cycleway", "track", "living_street", "residential", "service",
                   "unclassified", "road", "tertiary", "tertiary_link", "secondary", "secondary_link", "primary", "primary_link",
                   "trunk", "trunk_link", "corridor", "bridleway")
  osm_dir <- file.path(out_dir, "data/raw/osm"); dir.create(osm_dir, showWarnings = FALSE, recursive = TRUE)
  cache <- file.path(osm_dir, "poznan_walk_ways.rds")
  if (!file.exists(cache)) {
    bb <- st_bbox(st_transform(st_buffer(city, 1500), 4326))
    for (u in c("https://lz4.overpass-api.de/api/interpreter", "https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter")) {
      set_overpass_url(u)
      raw <- tryCatch(osmdata_sf(opq(bbox = unname(bb[c("xmin", "ymin", "xmax", "ymax")]), timeout = 900) %>%
                                   add_osm_feature(key = "highway", value = walk_values))$osm_lines,
                      error = function(e) { cat("Overpass failed at", u, ":", conditionMessage(e), "\n"); NULL })
      if (!is.null(raw) && nrow(raw) > 1000) {
        keep <- intersect(c("osm_id", "highway", "foot", "access", "geometry"), names(raw))
        saveRDS(raw[, keep], cache); cat("Overpass server:", u, "| ways fetched:", nrow(raw), "| downloaded", format(Sys.Date()), "\n"); break
      }
    }
  }
  if (file.exists(cache)) {
    w0 <- readRDS(cache)
    for (k in c("foot", "access")) if (!k %in% names(w0)) w0[[k]] <- NA_character_
    w0 <- w0[!(w0$foot %in% "no") & !(w0$access %in% c("no", "private") & !(w0$foot %in% c("yes", "designated", "permissive"))), ]
    w <- st_transform(st_geometry(w0), 2180); net_source <- "OpenStreetMap"
    cat("Network: OpenStreetMap,", length(w), "walkable ways\n")
  } else {
    bdot <- file.path(out_dir, "data/raw/bdot10k/PL.PZGiK.308.3064/BDOT10k")
    rd <- function(layer) st_read(file.path(bdot, paste0("PL.PZGiK.308.3064__", layer, ".xml")), quiet = TRUE)
    rdw <- rd("OT_SKDR_L"); rdw <- rdw[!rdw$klasaDrogi %in% c("autostrada", "droga ekspresowa"), ]
    pth <- rd("OT_SKRP_L"); pth <- pth[pth$kategoriaIstnienia == "eksploatowany", ]
    w <- st_transform(c(st_geometry(rdw), st_geometry(pth)), 2180); net_source <- "BDOT10k"
    cat(sprintf("Network: BDOT10k fallback, %d roads + %d pedestrian and cycle paths\n", nrow(rdw), nrow(pth)))
  }

  ## ── graph: nodes from segment end points (ways split at <= 20 m) ─────────
  w <- st_segmentize(st_cast(st_geometry(w), "LINESTRING"), dfMaxLength = 20)
  xy <- as.data.table(st_coordinates(w))[, .(X, Y, L1)]
  xy[, key := paste(round(X, 1), round(Y, 1))]
  nodes <- unique(xy[, .(key, X, Y)], by = "key")[, id := .I]
  xy[nodes, id := i.id, on = "key"]
  xy[, `:=`(nid = shift(id, -1), nX = shift(X, -1), nY = shift(Y, -1)), by = L1]
  E <- xy[!is.na(nid) & id != nid, .(from = id, to = nid, weight = sqrt((nX - X)^2 + (nY - Y)^2))]
  ## path and road ends that stop short of the network (digitising gaps) are joined to the nearest node of
  ## another line within 5 m
  deg <- tabulate(c(E$from, E$to), nbins = nrow(nodes))
  ends <- which(deg == 1)
  line_of <- xy[!duplicated(id), .(id, L1)][order(id)]$L1
  kn <- knn(as.matrix(nodes[, .(X, Y)]), as.matrix(nodes[ends, .(X, Y)]), k = 6)
  cand <- data.table(from = rep(ends, 6), to = as.vector(kn$nn.idx), weight = as.vector(kn$nn.dists))
  cand <- cand[from != to & weight <= 5 & line_of[from] != line_of[to]][order(from, weight)][!duplicated(from)]
  E <- rbind(E, cand)
  cat(sprintf("Dangling ends: %d; joined to another line within 5 m: %d\n", length(ends), nrow(cand)))
  g <- graph_from_data_frame(E, directed = FALSE, vertices = data.table(name = nodes$id))
  cmp <- components(g)
  giant <- which(cmp$membership == which.max(cmp$csize))
  cat(sprintf("Network: %d nodes, %d edges; largest connected part %.1f%% of nodes\n", vcount(g), ecount(g), 100 * length(giant) / vcount(g)))
  gx <- as.matrix(nodes[giant, .(X, Y)])

  ## ── stop poles and multi-source search ───────────────────────────────────
  st <- fread(file.path(gtfs, "stops.txt"), encoding = "UTF-8")
  sp <- st_transform(st_as_sf(st, coords = c("stop_lon", "stop_lat"), crs = 4326), 2180)
  sp <- sp[lengths(st_intersects(sp, city)) > 0, ]
  C <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  sp <- sp[sp$stop_name %in% C$stop_name, ]
  pk <- knn(gx, st_coordinates(sp), k = 1)
  poles <- data.table(stop_name = sp$stop_name, node = giant[pk$nn.idx[, 1]], snap = pk$nn.dists[, 1])
  cat(sprintf("Poles: %d of %d complexes; snapping distance median %.0f m, 95th pct %.0f m, max %.0f m\n",
              nrow(poles), uniqueN(poles$stop_name), median(poles$snap), quantile(poles$snap, .95), max(poles$snap)))
  pn <- poles[order(snap)][!duplicated(node)]                       # one complex per node: the closest pole
  S <- vcount(g) + 1L
  g2 <- add_vertices(g, 1, name = "S")
  g2 <- add_edges(g2, as.vector(rbind(S, pn$node)), weight = pn$snap)
  dist_n <- as.numeric(distances(g2, v = S, weights = E(g2)$weight))
  ## shortest-path tree from the distances alone (asking igraph for paths would build one path per node):
  ## a node's predecessor is a neighbour whose distance plus the edge length equals the node's distance
  E2 <- rbind(E[, .(a = from, b = to, w = weight)], E[, .(a = to, b = from, w = weight)])
  E2 <- E2[is.finite(dist_n[a]) & abs(dist_n[a] + w - dist_n[b]) < 1e-4][order(b, -dist_n[a])][!duplicated(b)]
  pred <- rep(NA_integer_, S); pred[E2$b] <- E2$a
  self <- pn$node[abs(dist_n[pn$node] - pn$snap) < 1e-4]           # poles whose node is reached directly from the super-source
  ## follow the tree back to the pole node that each node hangs from (pointer jumping)
  root <- pred; root[self] <- self; root[S] <- S
  repeat { nr <- root; ok <- !is.na(nr); nr[ok] <- root[nr[ok]]; if (identical(nr, root)) break; root <- nr }
  node_cx <- pn$stop_name[match(root, pn$node)]
  cat(sprintf("Nodes reached: %.1f%%; nodes within 400 m of a pole: %.1f%%\n", 100 * mean(is.finite(dist_n[-S])), 100 * mean(dist_n[-S] <= 400, na.rm = TRUE)))

  ## ── 20 m lattice, assignment and polygons ────────────────────────────────
  r <- rast(ext(st_bbox(city)), resolution = 20, crs = "EPSG:2180")
  lat <- as.data.table(xyFromCell(r, seq_len(ncell(r))))
  inside <- lengths(st_intersects(st_as_sf(lat, coords = c("x", "y"), crs = 2180), city)) > 0
  lk <- knn(gx, as.matrix(lat[, .(x, y)]), k = 1)
  nd <- giant[lk$nn.idx[, 1]]
  tot <- lk$nn.dists[, 1] + dist_n[nd]
  cx_names <- sort(unique(C$stop_name))
  val <- match(node_cx[nd], cx_names); val[!inside | !is.finite(tot) | tot > 400] <- NA
  values(r) <- val
  pol <- st_as_sf(as.polygons(r, dissolve = TRUE))
  names(pol)[1] <- "cx_id"
  pol$stop_name <- cx_names[pol$cx_id]
  pol <- st_make_valid(pol[, "stop_name"])
  cat(sprintf("Network catchments: %d of %d complexes; median area %.3f km2 (straight-line exclusive: see script 21)\n",
              nrow(pol), nrow(C), median(as.numeric(st_area(pol))) / 1e6))
  st_write(pol, file.path(s2, "stop_complexes.gpkg"), layer = "exclusive_400m_network", delete_layer = TRUE, quiet = TRUE)

  ## ── comparison with straight-line exclusive catchments ──────────────────
  vc <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_straightline", quiet = TRUE)
  p <- st_transform(st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE), 2180)
  pop <- p[[grep("ZAMELDOWANYCH", names(p), value = TRUE)[1]]]
  in_net <- lengths(st_intersects(p, pol)) > 0; in_str <- lengths(st_intersects(p, vc)) > 0
  cpop <- function(layer) { j <- st_intersects(layer, p); setNames(sapply(j, function(k) sum(pop[k])), layer$stop_name) }
  pn_pop <- cpop(pol); ps_pop <- cpop(vc)
  D <- st_distance(C); diag(D) <- Inf; nn <- C$stop_name[apply(D, 1, which.min)]
  nnr <- function(v) { a <- v[C$stop_name]; b <- v[nn]; a[is.na(a)] <- 0; b[is.na(b)] <- 0; cor(a, b, method = "spearman") }
  same <- as.data.table(st_drop_geometry(st_join(p[in_net & in_str, ], pol, join = st_within)))[, .(net = stop_name)]
  same$str <- as.data.table(st_drop_geometry(st_join(p[in_net & in_str, ], vc, join = st_within)))$stop_name
  summ <- data.table(catchment = c("400 m exclusive, straight line", "400 m exclusive, walking network"),
                     complexes = c(nrow(vc), nrow(pol)),
                     median_area_km2 = c(median(as.numeric(st_area(vc))), median(as.numeric(st_area(pol)))) / 1e6,
                     residents_covered = c(sum(pop[in_str]), sum(pop[in_net])) / sum(pop),
                     complexes_without_residents = c(sum(ps_pop[C$stop_name] %in% c(0, NA)), sum(pn_pop[C$stop_name] %in% c(0, NA))),
                     nn_rho_population = c(nnr(ps_pop), nnr(pn_pop)))
  print(summ, digits = 3)
  cat(sprintf("Residents covered by both: same complex in %.1f%% of cases (population-weighted)\n",
              100 * sum(pop[in_net & in_str][same$net == same$str], na.rm = TRUE) / sum(pop[in_net & in_str])))
  fwrite(summ, file.path(s2, "catchment_network_summary.csv"))

  ## ── figure: central window ───────────────────────────────────────────────
  ctr <- st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180)
  win <- st_as_sfc(st_bbox(st_buffer(ctr, 1800)))
  ways <- st_sf(geometry = w)[lengths(st_intersects(st_sf(geometry = w), win)) > 0, ]
  cmp_l <- rbind(st_sf(type = "Straight line (nearest complex centre)", geometry = st_geometry(st_intersection(vc, win))),
                 st_sf(type = "Walking network (nearest pole)", geometry = st_geometry(st_intersection(pol, win))))
  g <- ggplot() + geom_sf(data = ways, colour = "#d9d6cc", linewidth = 0.15) +
    geom_sf(data = cmp_l, fill = NA, colour = "#2a78d6", linewidth = 0.35) +
    geom_sf(data = st_intersection(st_geometry(C), win), colour = "#e34948", size = 0.9) + facet_wrap(~ type) +
    labs(title = "Exclusive 400 m stop catchments in central Poznań", caption = sprintf("Grey: %s pedestrian network; red: stop complexes.", net_source)) +
    theme_void(base_size = 10) + theme(plot.background = element_rect(fill = "#ffffff", colour = NA), plot.title = element_text(face = "bold"),
                                       strip.text = element_text(face = "bold"), plot.caption = element_text(colour = "#898781", size = 7.5))
  ggsave(file.path(out_dir, "output/figures/Fig_catchments_network_vs_straight.png"), g, width = 10, height = 5.4, dpi = 200, bg = "#ffffff")
  cat("Done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
})

## ── 15b · Land-use descriptors of the walking catchments ───────────────────────────────────────────────────
section("15b", "Land-use descriptors of the walking catchments", {
  source(file.path(REPO, "R/tissue_descriptors.R")); tissue_descriptors("exclusive_400m_network")
})

## ── 34 · Catchment covariates ──────────────────────────────────────────────────────────────────────────────
section("34", "Catchment covariates", {
  ## Catchment covariates for all 607 stop complexes (the pilot scripts 25-29
  ## built them for the 253 October pilot complexes only).
  ##
  ## Usage: Rscript 34_catchment_covariates_all.R [layer]
  ## Catchments: a layer of stop_complexes.gpkg -- exclusive_400m_straightline
  ## (default; the pilot catchments, output catchment_covariates_all.csv) or
  ## exclusive_400m_network (the manuscript's primary walking catchments, script
  ## 39; output catchment_covariates_network.csv).
  ## Per complex: registered residents, population-weighted catchment
  ## deprivation and dominant osiedle (script 25 method); distance to the centre;
  ## terminal / interchange / boundary flags; age shares 0-18, 19-25 and 61+;
  ## floor area by function, the
  ## residential vs non-residential balance, land-use mix, office, industry and
  ## retail floor shares, university and dormitory floor shares and WorldCover
  ## tree share (script 26 method); paid-parking-zone share (script 28 method); jobs per
  ## resident from the modelled workplace grid (script 29 method).
  ## Output: output/stage2/catchment_covariates_all.csv or catchment_covariates_network.csv

  suppressPackageStartupMessages({library(data.table); library(sf); library(terra)})
  out_dir <- ROOT
  s1 <- file.path(out_dir, "output/stage1"); s2 <- file.path(out_dir, "output/stage2")
  ## the walking catchments of script 39 (the straight-line pilot layer is retired)
  layer <- "exclusive_400m_network"
  out_csv <- file.path(s2, if (layer == "exclusive_400m_straightline") "catchment_covariates_all.csv" else "catchment_covariates_network.csv")
  vc  <- st_read(file.path(s2, "stop_complexes.gpkg"), layer, quiet = TRUE)
  cat("Catchment layer:", layer, "\n")
  C   <- st_read(file.path(s2, "stop_complexes.gpkg"), "complexes", quiet = TRUE)
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  P   <- fread(file.path(s1, "s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, dep)]
  att <- fread(file.path(s2, "stop_complex_attributes.csv"), encoding = "UTF-8")[, .(stop_name, terminal, interchange, boundary)]

  ## ── residents, deprivation, age ──────────────────────────────────────────
  pp <- st_transform(st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE), 2180)
  nm <- names(pp)
  pp$pop <- pp[[grep("ZAMELDOWANYCH", nm, value = TRUE)[1]]]
  pp$osiedle <- osi$name[st_nearest_feature(pp, osi)]
  jn <- as.data.table(st_drop_geometry(st_join(pp, vc[, "stop_name"], join = st_within)))[!is.na(stop_name) & pop > 0]
  jn <- merge(jn, P, by = "osiedle")
  cx <- jn[, .(residents = sum(pop), dep = sum(pop * dep) / sum(pop)), by = stop_name]
  dom <- jn[, .(pop = sum(pop)), by = .(stop_name, osiedle)][order(-pop)][!duplicated(stop_name), .(stop_name, osiedle)]
  cx <- merge(cx, dom, by = "stop_name")
  ## osiedle mix of catchment residents (manuscript §2.5)
  mix <- jn[, .(pop = sum(pop)), by = .(stop_name, osiedle)][, share := pop / sum(pop), by = stop_name][, .(n_os = .N, top = max(share)), by = stop_name]
  cat(sprintf("Catchments with residents: %d | residents from 2+ osiedla: %.0f%% | largest osiedle < 90%% of residents: %.0f%%\n", nrow(mix), 100 * mean(mix$n_os > 1), 100 * mean(mix$top < 0.9)))
  ## age shares from the register's age bands (0-2, 3-6, 7-14, 15-18, 19-25, 26-60, 61+)
  ag <- jn[, .(a0_18 = sum(WIEK_0_DO_2 + WIEK_3_DO_6 + WIEK_7_DO_14 + WIEK_15_DO_18, na.rm = TRUE),
               a19_25 = sum(WIEK_19_DO_25, na.rm = TRUE), a61 = sum(get(grep("^WIEK_61", names(jn), value = TRUE)), na.rm = TRUE),
               aall = sum(WIEK_0_DO_2 + WIEK_3_DO_6 + WIEK_7_DO_14 + WIEK_15_DO_18 + WIEK_19_DO_25 + WIEK_26_DO_60 +
                          get(grep("^WIEK_61", names(jn), value = TRUE)), na.rm = TRUE)), by = stop_name]
  ag <- ag[, .(stop_name, share_0_18 = a0_18 / aall, share_19_25 = a19_25 / aall, share_61plus = a61 / aall)]
  cx <- merge(cx, ag, by = "stop_name", all.x = TRUE)

  ## ── land use ─────────────────────────────────────────────────────────────
  b  <- st_read(file.path(s1, "tissue_inputs.gpkg"), "buildings", quiet = TRUE)
  b  <- b[b$fn != "auxiliary", ]
  bc <- as.data.table(st_drop_geometry(st_join(st_centroid(b[, c("fn", "floor", "university", "dormitory")]), vc[, "stop_name"], join = st_within)))[!is.na(stop_name)]
  fns <- c("res_multi", "res_single", "retail_services", "office", "industry_storage", "education_culture", "health", "other_nonres")
  fl <- dcast(bc[, .(floor = sum(floor)), by = .(stop_name, fn)], stop_name ~ fn, value.var = "floor", fill = 0)
  for (f in setdiff(fns, names(fl))) fl[, (f) := 0]
  fl[, tot := rowSums(.SD), .SDcols = fns]
  fl[, ilr_res_vs_nonres := sqrt(1/2) * log((res_multi + res_single + 100) / (tot - res_multi - res_single + 100))]
  fl[, landuse_mix := apply(as.matrix(.SD), 1, function(x) { if (sum(x) <= 0) return(NA_real_); p <- x[x > 0] / sum(x); -sum(p * log(p)) / log(length(x)) }), .SDcols = fns]
  fl[, `:=`(share_office = office / tot, share_industry = industry_storage / tot, share_retail = retail_services / tot)]
  ud <- bc[, .(share_university = sum(floor[university]) / sum(floor), share_dormitory = sum(floor[dormitory]) / sum(floor)), by = stop_name]
  wc <- rast(file.path(s1, "worldcover_poznan.tif"))
  tr <- terra::extract(wc == 10, vect(st_transform(vc, crs(wc))), fun = mean, na.rm = TRUE)
  lu <- merge(merge(fl[, .(stop_name, floor_total = tot, ilr_res_vs_nonres, landuse_mix, share_office, share_industry, share_retail)], ud, by = "stop_name", all = TRUE),
              data.table(stop_name = vc$stop_name, wc_tree = tr[[2]]), by = "stop_name", all = TRUE)

  ## ── parking zone and jobs ────────────────────────────────────────────────
  pm <- st_transform(st_read(file.path(out_dir, "data/raw/spp/parking_meters_2026-10-07.geojson"), quiet = TRUE), 2180)
  zone <- st_union(st_buffer(pm, 150))
  iz <- st_intersection(vc[, "stop_name"], zone)
  spp <- merge(data.table(stop_name = vc$stop_name, area = as.numeric(st_area(vc))), data.table(stop_name = iz$stop_name, inz = as.numeric(st_area(iz))), by = "stop_name", all.x = TRUE)
  spp[, spp_share := fifelse(is.na(inz), 0, inz / area)]
  wg <- st_transform(st_read(file.path(out_dir, "data/workplace_grid.gpkg"), quiet = TRUE), 2180); wg$cell_area <- as.numeric(st_area(wg))
  ix <- suppressWarnings(st_intersection(vc[, "stop_name"], wg[, c("workplaces", "cell_area")])); ix$f <- as.numeric(st_area(ix)) / ix$cell_area
  jobs <- as.data.table(st_drop_geometry(ix))[, .(jobs = sum(workplaces * f)), by = stop_name]

  ## ── assemble ─────────────────────────────────────────────────────────────
  rynek <- st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180)
  d <- data.table(stop_name = C$stop_name, dist_km = as.numeric(st_distance(C, rynek)) / 1000)
  d <- Reduce(function(a, b) merge(a, b, by = "stop_name", all.x = TRUE), list(d, att, cx, lu, spp[, .(stop_name, spp_share, catchment_area_m2 = area)], jobs))
  d[is.na(residents), residents := 0][is.na(jobs), jobs := 0]
  d[, jobs_res := log((jobs + 1) / (residents + 1))][, flagged := (terminal | interchange | boundary) %in% TRUE]
  fwrite(d, out_csv)
  cat(sprintf("Complexes: %d | with residents: %d | flagged: %d | with residents and not flagged: %d\n",
              nrow(d), d[residents > 0, .N], d[flagged == TRUE, .N], d[residents > 0 & flagged == FALSE, .N]))
})

## ── 37 · Municipal housing ─────────────────────────────────────────────────────────────────────────────────
section("37", "Municipal housing", {
  ## Municipal housing as a finer deprivation measure: the share of a stop
  ## catchment's registered residents who live in municipal buildings.
  ##
  ## Source: ZKZL (Zarząd Komunalnych Zasobów Lokalowych), "Wykaz budynków i
  ## garaży komunalnych", public HTML table of 710 entries (street and building
  ## numbers), saved 2026-10-07 to data/raw/zkzl/. The list covers buildings
  ## managed by ZKZL; flats the city owns in housing-community buildings
  ## (wspólnoty) are not itemised there.
  ## Parsing: corner buildings ("A 12 / B 12") are split; old addresses in
  ## brackets ("DAWNA ...") are dropped; number lists ("5, 5A, 6") are split;
  ## "BARAKI" entries are flagged. Garages and commercial units need no special
  ## handling: they have no registered residents.
  ## Geocoding: exact street + number match against the GEOPOZ address points;
  ## a "?" in the source stands for a lost "Ó"; barracks entries keep the
  ## building number; "24/25" falls back to "24"; last resort: closest street
  ## name (edit distance <= 3) among streets that have the number
  ## (as in script 04), with street names matched in three steps -- exact, after
  ## dropping a prefix (aleja, osiedle, plac, rondo), and by the GEOPOZ name that
  ## ends with the ZKZL name when exactly one such street has that number.
  ## Residents: registered population at the matched addresses (script 04).
  ## Usage: Rscript 37_municipal_housing.R [catchment layer]
  ## Outputs: output/stage2/municipal_buildings.csv, municipal_catchments.csv
  ##          (municipal_catchments_network.csv for the walking catchments),
  ##          output/figures/Fig_municipal_housing_map.png

  suppressPackageStartupMessages({library(data.table); library(sf); library(stringr); library(ggplot2)})
  out_dir <- ROOT
  s2 <- file.path(out_dir, "output/stage2")

  ## ── parse the ZKZL table ─────────────────────────────────────────────────
  h <- paste(readLines(file.path(out_dir, "data/raw/zkzl/wykaz_budynkow_komunalnych_2026-10-07.html"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
  rows <- regmatches(h, gregexpr("(?s)<tr[^>]*>.*?</tr>", h, perl = TRUE))[[1]]
  cells <- lapply(rows, function(r) str_squish(gsub("&nbsp;", " ", gsub("<[^>]+>", "", regmatches(r, gregexpr("(?s)<t[dh][^>]*>.*?</t[dh]>", r, perl = TRUE))[[1]]))))
  z <- data.table(lp = vapply(cells, `[`, "", 1), entry = vapply(cells, `[`, "", 2))[lp != "lp"]
  z[, `:=`(barracks = grepl("BARAK", entry), garage = grepl("GARA", entry))]
  cat("ZKZL entries:", nrow(z), "| barracks:", sum(z$barracks), "| mentioning garages:", sum(z$garage), "\n")
  clean <- function(e) {
    e <- gsub("?", "Ó", e, fixed = TRUE)                           # a lost 'Ó' in the source
    e <- gsub("([A-ZĄĆĘŁŃÓŚŹŻ])(\\d)", "\\1 \\2", e)               # "OPOLSKA81/83" -> "OPOLSKA 81/83"
    e <- sub(",? ?BARAK[^0-9]*[0-9 ,]*$", "", e)                  # "JARZĘBOWA 18, BARAK: 1, 2" -> "JARZĘBOWA 18"
    str_squish(gsub("\\+ ?GARA[^,/]*", "", gsub("\\([^)]*\\)", "", e)))
  }
  parts <- z[, .(part = str_squish(unlist(strsplit(clean(entry), "/(?= ?[A-ZĄĆĘŁŃÓŚŹŻ]{2})", perl = TRUE)))), by = .(lp, entry, barracks)]
  ## split street and numbers: numbers start at the first token that begins with a digit and is followed
  ## only by number-like tokens; "28 CZERWCA 1956 NR 140" keeps the number after NR
  split_addr <- function(p) {
    if (grepl(" NR ", p)) return(list(street = str_squish(sub(" NR .*$", "", p)), nums = str_squish(sub("^.* NR ", "", p))))
    m <- regexpr("\\s\\d+[A-Z]?(\\s?/\\s?\\d+[A-Z]?)?(\\s?(,|I|ORAZ)\\s?\\d+[A-Z]?(\\s?/\\s?\\d+[A-Z]?)?)*\\s*$", p, perl = TRUE)
    if (m < 0) return(list(street = p, nums = NA_character_))
    list(street = str_squish(substr(p, 1, m - 1)), nums = str_squish(substr(p, m, nchar(p))))
  }
  sp <- parts[, { s <- split_addr(part); .(street = s$street, nums = s$nums) }, by = .(lp, entry, barracks, part)]
  ad <- sp[!is.na(nums), .(nr = str_squish(unlist(strsplit(nums, "\\s*(,|\\sI\\s|ORAZ)\\s*")))), by = .(lp, entry, barracks, street)]
  ad[, nr := toupper(gsub("\\s", "", nr))]
  cat("Addresses after splitting:", nrow(ad), "from", uniqueN(ad$lp), "entries; entries without a number:", sp[is.na(nums), uniqueN(lp)], "\n")

  ## ── geocode against GEOPOZ address points ───────────────────────────────
  ap <- st_read(file.path(out_dir, "output/geopoz_address_points.gpkg"), quiet = TRUE)
  xy <- st_coordinates(ap)
  A <- data.table(street_g = toupper(str_squish(ap$NAZWY_UL_NAZWA)), typ = ap$TYPY_UL_TYP, nr = toupper(gsub("\\s", "", ap$NR)),
                  osiedle = ap$JEDN_POM, x = xy[, 1], y = xy[, 2])[!duplicated(paste(street_g, nr))]
  ## "Święty Marcin", "św. Trójcy", "ŚWIĘTEJ TRÓJCY", "ŚW. MARCIN" -> "ŚW MARCIN", "ŚW TRÓJCY"
  saint <- function(s) str_squish(gsub("(^|\\s)(ŚWIĘTEJ|ŚWIĘTEGO|ŚWIĘTY|ŚWIĘTA|ŚW\\.?)\\s+", "\\1ŚW ", toupper(s), perl = TRUE))
  A[, street_n := saint(street_g)]
  streets <- unique(A$street_n)
  strip <- function(s) str_squish(sub("^(ALEJA|ALEJE|AL\\.|OSIEDLE|OS\\.|PLAC|PL\\.|RONDO|RYNEK|ULICA|UL\\.)\\s+", "", s))
  ad[, s1 := saint(street)][, s2 := strip(s1)]
  key <- paste(A$street_n, A$nr)
  match_street <- function(s1, s2, nr) {
    if (s1 %in% streets) return(s1)
    if (s2 %in% streets) return(s2)
    ## GEOPOZ uses full names: "STANISŁAWA TACZAKA" for "TACZAKA", "GALLA ANONIMA" for "GALLA"
    cand <- streets[endsWith(streets, paste0(" ", s2)) | startsWith(streets, paste0(s2, " ")) | endsWith(s2, paste0(" ", streets))]
    cand <- cand[paste(cand, nr) %in% key]
    if (length(cand) == 1) return(cand)
    ## fuzzy: closest street name (edit distance <= 3) among streets that have this number
    have <- A$street_n[A$nr == nr]
    if (length(have)) { d <- adist(s2, have); if (min(d) <= 3 && sum(d == min(d)) == 1) return(have[which.min(d)]) }
    NA_character_
  }
  ad[, street_g := mapply(match_street, s1, s2, nr)]
  ## numbers such as "24/25" that GEOPOZ records under the first number only
  ad[is.na(street_g) & grepl("/", nr), nr := sub("/.*$", "", nr)]
  ad[is.na(street_g), street_g := mapply(match_street, s1, s2, nr)]
  setnames(ad, "street_g", "street_n")
  ad <- merge(ad, A, by = c("street_n", "nr"), all.x = TRUE)
  cat(sprintf("Geocoded addresses: %d of %d (%.1f%%); entries with at least one geocoded address: %d of %d\n",
              ad[!is.na(x), .N], nrow(ad), 100 * ad[, mean(!is.na(x))], ad[!is.na(x), uniqueN(lp)], nrow(z)))
  cat("Unmatched examples:\n"); print(head(unique(ad[is.na(x), .(entry, street, nr)]), 20))

  ## ── residents in municipal buildings ─────────────────────────────────────
  pp <- st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE)
  pp$pop <- pp[[grep("ZAMELDOWANYCH", names(pp), value = TRUE)[1]]]
  P <- data.table(street_g = toupper(str_squish(pp$ULICA)), nr = toupper(gsub("\\s", "", pp$NR_DOMU)), pop = pp$pop)[, .(pop = sum(pop)), by = .(street_g, nr)]
  mb <- unique(ad[!is.na(x), .(street_g, nr, barracks, osiedle, x, y)])[, barracks := any(barracks), by = .(street_g, nr)]
  mb <- unique(mb)
  mb <- merge(mb, P, by = c("street_g", "nr"), all.x = TRUE)[is.na(pop), pop := 0]
  fwrite(mb, file.path(s2, "municipal_buildings.csv"))
  cat(sprintf("Municipal addresses geocoded: %d; with registered residents: %d; residents in them: %d (city register total %d)\n",
              nrow(mb), mb[pop > 0, .N], sum(mb$pop), sum(pp$pop)))
  cat("Residents in barracks entries:", mb[barracks == TRUE, sum(pop)], "\n")

  ## ── catchments and osiedla ───────────────────────────────────────────────
  ## catchment layer as argument: exclusive_400m_straightline (default, pilot) or exclusive_400m_network (primary)
  ## the walking catchments of script 39 (the straight-line pilot layer is retired)
  layer <- "exclusive_400m_network"
  net <- layer != "exclusive_400m_straightline"
  cat("Catchment layer:", layer, "\n")
  vc <- st_read(file.path(s2, "stop_complexes.gpkg"), layer, quiet = TRUE)
  mbs <- st_transform(st_as_sf(mb, coords = c("x", "y"), crs = st_crs(ap)), st_crs(vc))
  mj <- as.data.table(st_drop_geometry(st_join(mbs, vc[, "stop_name"], join = st_within)))[!is.na(stop_name)]
  cv <- fread(file.path(s2, if (net) "catchment_covariates_network.csv" else "catchment_covariates_all.csv"), encoding = "UTF-8")
  mc <- merge(cv[, .(stop_name, osiedle, residents, dep)], mj[, .(muni_residents = sum(pop), muni_buildings = sum(pop > 0)), by = stop_name], by = "stop_name", all.x = TRUE)
  mc[is.na(muni_residents), `:=`(muni_residents = 0, muni_buildings = 0)]
  mc[, muni_share := fifelse(residents > 0, muni_residents / residents, NA_real_)]
  fwrite(mc, file.path(s2, if (net) "municipal_catchments_network.csv" else "municipal_catchments.csv"))
  cat(sprintf("\nCatchments with residents: %d; with any municipal residents: %d; share distribution (those with any): median %.3f, p90 %.3f, max %.3f\n",
              mc[residents > 0, .N], mc[muni_residents > 0, .N], mc[muni_residents > 0, median(muni_share)], mc[muni_residents > 0, quantile(muni_share, .9)], mc[residents > 0, max(muni_share)]))

  ## validation against the osiedle indicators
  G <- fread(file.path(out_dir, "output/stage1/s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, dep, B2021, B_pooled)]
  pq <- as.data.table(st_drop_geometry(st_join(st_transform(pp[, "pop"], st_crs(vc)), st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE)[, "name"], st_crs(vc)), join = st_nearest_feature)))
  os <- merge(pq[, .(residents = sum(pop)), by = .(osiedle = name)], mb[, .(muni_residents = sum(pop)), by = osiedle], by = "osiedle", all.x = TRUE)[is.na(muni_residents), muni_residents := 0]
  os[, muni_share := muni_residents / residents]
  os <- merge(os, G, by = "osiedle")
  cat(sprintf("Osiedle level (42): rho(municipal share, deprivation) %.2f | rho(municipal share, poverty 2021) %.2f | rho(municipal share, poverty pooled) %.2f\n",
              cor(os$muni_share, os$dep, method = "spearman"), cor(os$muni_share, os$B2021, method = "spearman"), cor(os$muni_share, os$B_pooled, method = "spearman")))
  print(os[order(-muni_share)][1:8, .(osiedle, muni_share = round(muni_share, 3), dep = round(dep, 2), B2021)])
  w <- mc[residents > 0]
  cat(sprintf("Catchment level: rho(municipal share, catchment deprivation) %.2f; share of variance within osiedla %.0f%%\n",
              cor(w$muni_share, w$dep, method = "spearman"), 100 * var(w$muni_share - ave(w$muni_share, w$osiedle)) / var(w$muni_share)))

  ## ── map ──────────────────────────────────────────────────────────────────
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  osi <- merge(osi, os[, .(name = osiedle, dep)], by = "name")
  mp <- st_transform(mbs[mbs$pop > 0, ], 2180)
  g <- ggplot() + geom_sf(data = osi, aes(fill = dep), colour = "#ffffff", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f2efe6", high = "#e34948", midpoint = 0, name = "Osiedle deprivation") +
    geom_sf(data = mp, aes(size = pop, shape = barracks), colour = "#0b0b0a", alpha = 0.7) +
    scale_size_area(max_size = 3, name = "Registered residents") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("Building", "Barracks"), name = NULL) +
    labs(title = "Municipal housing (ZKZL buildings) and osiedle deprivation, Poznań",
         caption = "ZKZL list of municipal buildings (2026), geocoded to GEOPOZ address points; residents from the municipal address register.") +
    theme_void(base_size = 10) + theme(plot.background = element_rect(fill = "#ffffff", colour = NA), plot.title = element_text(face = "bold"),
                                       plot.caption = element_text(colour = "#898781", size = 7.5))
  ggsave(file.path(out_dir, "output/figures/Fig_municipal_housing_map.png"), g, width = 9, height = 8, dpi = 200, bg = "#ffffff")
  cat("Saved municipal_buildings.csv, municipal_catchments.csv, Fig_municipal_housing_map.png\n")
})

## ── 62 · Catchment deprivation by round and catchment-osiedle weights ──────────────────────────────────────
section("62", "Catchment deprivation by round and catchment-osiedle weights", {
  ## Catchment deprivation by round, for the measurement-error checks of script 56.
  ## The score of the paper is the mean of the 2019 and 2021 rounds of the city's
  ## standardised unemployment and poverty indicators (script 20). Each catchment
  ## takes the mean osiedle score of its registered residents (as in script 34).
  ## Here the same resident weights give the score of each round (2015, 2019,
  ## 2021) and of the three-round mean, and each catchment's share of residents
  ## from its main osiedle (catchments drawing at least 90% from one osiedle are
  ## "pure": their deprivation is not a mix across osiedle boundaries).
  ## Output: output/stage2/catchment_dep_rounds.csv

  suppressPackageStartupMessages({library(sf); library(data.table)})
  out_dir <- ROOT
  s1 <- file.path(out_dir, "output/stage1"); s2 <- file.path(out_dir, "output/stage2")
  vc  <- st_read(file.path(s2, "stop_complexes.gpkg"), "exclusive_400m_network", quiet = TRUE)
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  P   <- fread(file.path(s1, "s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, dep, dep3, dep2015, dep2019, dep2021)]
  pp <- st_transform(st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE), 2180)
  pp$pop <- pp[[grep("ZAMELDOWANYCH", names(pp), value = TRUE)[1]]]
  pp$osiedle <- osi$name[st_nearest_feature(pp, osi)]
  jn <- as.data.table(st_drop_geometry(st_join(pp[, c("pop", "osiedle")], vc[, "stop_name"], join = st_within)))[!is.na(stop_name) & pop > 0]
  jn <- merge(jn, P, by = "osiedle")
  cx <- jn[, .(residents = sum(pop), dep = sum(pop * dep) / sum(pop), dep3 = sum(pop * dep3) / sum(pop), dep2015 = sum(pop * dep2015) / sum(pop),
               dep2019 = sum(pop * dep2019) / sum(pop), dep2021 = sum(pop * dep2021) / sum(pop)), by = stop_name]
  top <- jn[, .(pop = sum(pop)), by = .(stop_name, osiedle)][, .(top_share = max(pop) / sum(pop)), by = stop_name]
  cx <- merge(cx, top, by = "stop_name")
  cv <- fread(file.path(s2, "catchment_covariates_network.csv"), encoding = "UTF-8")[, .(stop_name, dep_cv = dep)]
  chk <- merge(cx, cv, by = "stop_name")
  cat(sprintf("Catchments: %d | check against catchment_covariates_network.csv: max |difference| = %.2e\n", nrow(cx), max(abs(chk$dep - chk$dep_cv))))
  cat(sprintf("Pure catchments (>= 90%% of residents from one osiedle): %.0f%%\n", 100 * mean(cx$top_share >= 0.9)))
  cat("Correlations between round scores at catchment level (Pearson):\n"); print(round(cor(cx[, .(dep, dep3, dep2015, dep2019, dep2021)]), 2))
  ## osiedle-level reliability of the two-round mean (Spearman-Brown on the between-round correlation)
  r <- P[, cor(dep2019, dep2021)]; rs <- P[, cor(dep2019, dep2021, method = "spearman")]
  cat(sprintf("Osiedla: 2019 vs 2021 r = %.2f (Spearman %.2f); reliability of the two-round mean 2r/(1+r) = %.2f (Spearman-based %.2f)\n", r, rs, 2 * r / (1 + r), 2 * rs / (1 + rs)))
  fwrite(cx, file.path(s2, "catchment_dep_rounds.csv"))
  ## residents of each catchment by osiedle (weights that map any osiedle-level field to catchments; script 66)
  fwrite(jn[, .(pop = sum(pop)), by = .(stop_name, osiedle)][, w := pop / sum(pop), by = stop_name], file.path(s2, "catchment_osiedle_weights.csv"))
  cat("Saved catchment_dep_rounds.csv, catchment_osiedle_weights.csv\n")
})

