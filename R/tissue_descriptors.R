## Stage 1: urban-tissue descriptors for every populated 200 m cell, measured
## in the 400 m radius around the cell centroid (the walking-catchment scale
## used for stops later). Buffers are clipped to the city boundary and every
## density/share uses the clipped area, because BDOT10k covers the city only.
##
## Families:
##   built form   : GSI (coverage), FSI (floor-area density), L = FSI/GSI
##                  (mean storeys), OSR = (1-GSI)/FSI, median footprint,
##                  building count density, residential population density
##   spacing/grain: median gap to the nearest neighbouring main building,
##                  share of main buildings touching another (shared walls),
##                  median urban-block area, official BDOT built-up character
##                  (compact / dense / loose)
##   function     : floor-area shares by function, ilr balances, overall
##                  land-use mix (entropy), university and student-dormitory
##                  floor shares, multi-family buildings with a registered
##                  commercial function, zoning shares (Studium 2023),
##                  amenity density and diversity (OSM)
##   street net   : street length density, intersection density, dead-end
##                  share, street length per intersection, main-road and
##                  internal-road length shares, orientation entropy
##   green/open   : woodland, grassland, water (BDOT); tree, grass, built-up
##                  cover (ESA WorldCover 10 m, satellite-derived)
##
## Run 14 first. Output: output/stage1/cell_tissue_descriptors.gpkg
## Called as tissue_descriptors("cells") in stage 2 and tissue_descriptors("exclusive_400m_network") in stage 4.

tissue_descriptors <- function(units_mode = "cells") {
  suppressPackageStartupMessages({library(sf); library(dplyr); library(terra); library(exactextractr); library(nabor)})

  out_dir <- ROOT
  s1_dir  <- file.path(out_dir, "output/stage1")
  gpkg    <- file.path(s1_dir, "tissue_inputs.gpkg")
  R_BUF   <- 400

  ## ── Units: populated cells, 400 m buffers clipped to the city ────────────
  ## Usage: Rscript 15_urban_tissue_descriptors.R [units]
  ##   units = "cells" (default, Stage 1) or a layer of output/stage2/stop_complexes.gpkg
  ##   (e.g. exclusive_400m_network): the same descriptors measured within each stop
  ##   catchment, written to output/stage2/catchment_tissue_descriptors_<layer>.csv


  city  <- st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE) %>%
    st_transform(2180) %>% st_union() %>% st_make_valid()
  if (units_mode == "cells") {
    cells <- st_read(file.path(out_dir, "output/population_weighted_deprivation_200m.gpkg"), quiet = TRUE) %>%
      filter(population > 0) %>% select(cell_id, osiedle, population, rate = indicator_B_poverty_per1000)
    cent  <- st_centroid(st_geometry(cells))
    buf   <- st_intersection(st_buffer(cent, R_BUF, nQuadSegs = 16), city)
    if (length(buf) != nrow(cells)) stop("buffer clipping changed the number of units")
  } else {
    cells <- st_read(file.path(out_dir, "output/stage2/stop_complexes.gpkg"), units_mode, quiet = TRUE) %>%
      st_transform(2180) %>% st_make_valid() %>% select(cell_id = stop_name)
    buf   <- st_geometry(cells)
  }
  bufsf <- st_sf(cell_id = cells$cell_id, geometry = buf)
  A_m2  <- as.numeric(st_area(bufsf)); A_km2 <- A_m2 / 1e6
  cat("Units:", nrow(cells), " median clipped buffer area (km2):", round(median(A_km2), 3), "\n")
  in_buf <- function(pts) st_intersects(bufsf, pts)   # list: for each buffer, indices of points
  sum_by <- function(idx, v) vapply(idx, function(i) sum(v[i]), numeric(1))

  ## ── Buildings: per-building spacing first ────────────────────────────────
  b <- st_read(gpkg, layer = "buildings", quiet = TRUE)
  main <- b$fn != "auxiliary"
  bm <- b[main, ]
  cm <- st_coordinates(st_centroid(st_geometry(bm)))
  nn <- nabor::knn(cm, cm, k = 9)$nn.idx[, -1]          # 8 nearest centroids, excluding self
  gap <- rep(NA_real_, nrow(bm))
  geo <- st_geometry(bm)
  for (j in seq_len(ncol(nn))) {
    d <- as.numeric(st_distance(geo, geo[nn[, j]], by_element = TRUE))
    gap <- pmin(gap, d, na.rm = TRUE)
  }
  bm$gap <- gap; bm$touch <- gap < 0.5
  cat("Main buildings:", nrow(bm), " median gap (m):", round(median(gap), 1), " touching share:", round(mean(bm$touch), 3), "\n")

  ## street interface: setback of each main building from the nearest
  ## carriageway edge = distance to the carriageway centreline (BDOT10k
  ## OT_SKJZ_L) minus half its recorded width. Perimeter blocks sit on the
  ## pavement line (a few metres); estate slabs stand back tens of metres.
  skjz <- st_read(file.path(out_dir, "data/raw/bdot10k/PL.PZGiK.308.3064/BDOT10k/PL.PZGiK.308.3064__OT_SKJZ_L.xml"), quiet = TRUE)
  skjz <- skjz[skjz$kategoriaIstnienia %in% "eksploatowany" & !skjz$klasaDrogi %in% c("autostrada", "droga ekspresowa") &
               skjz$polozenie %in% "na powierzchni gruntu", ]
  halfw <- suppressWarnings(as.numeric(as.character(skjz$szerokoscNawierzchni))) / 2
  halfw[is.na(halfw)] <- median(halfw, na.rm = TRUE)
  skjz <- st_transform(st_sf(halfw = halfw, geometry = st_geometry(skjz)), 2180)
  nf <- st_nearest_feature(geo, skjz)
  bm$setback <- pmax(0, as.numeric(st_distance(geo, st_geometry(skjz)[nf], by_element = TRUE)) - skjz$halfw[nf])
  cat("Median building setback from carriageway edge (m):", round(median(bm$setback), 1),
      " share within 5 m:", round(mean(bm$setback <= 5), 3), "\n")

  cb_all  <- st_centroid(st_geometry(b))
  cb_main <- st_centroid(geo)
  ia <- in_buf(cb_all); im <- in_buf(cb_main)

  fp_all  <- sum_by(ia, b$footprint); fl_all <- sum_by(ia, b$floor)
  fp_main <- sum_by(im, bm$footprint); fl_main <- sum_by(im, bm$floor)
  d <- data.frame(cell_id = cells$cell_id,
    gsi        = fp_all / A_m2,
    fsi        = fl_all / A_m2,
    storeys    = ifelse(fp_main > 0, fl_main / fp_main, NA),
    bldg_dens  = lengths(im) / A_km2,
    med_footprint = vapply(im, function(i) if (length(i)) median(bm$footprint[i]) else NA_real_, numeric(1)),
    med_gap    = vapply(im, function(i) if (length(i) >= 3) median(bm$gap[i]) else NA_real_, numeric(1)),
    touch_share= vapply(im, function(i) if (length(i) >= 3) mean(bm$touch[i]) else NA_real_, numeric(1)),
    med_setback = vapply(im, function(i) if (length(i) >= 3) median(bm$setback[i]) else NA_real_, numeric(1)),
    frontage_share = vapply(im, function(i) if (length(i) >= 3) mean(bm$setback[i] <= 5) else NA_real_, numeric(1)))
  d$osr <- ifelse(d$fsi > 0, (1 - d$gsi) / d$fsi, NA)

  ## floor-area shares by function (main buildings)
  fns <- c("res_multi","res_single","retail_services","office","industry_storage","education_culture","health","other_nonres")
  for (f in fns) d[[paste0("fl_", f)]] <- sum_by(im, bm$floor * (bm$fn == f))
  tot <- rowSums(d[paste0("fl_", fns)])
  for (f in fns) d[[paste0("share_", f)]] <- ifelse(tot > 0, d[[paste0("fl_", f)]] / tot, NA)
  res    <- d$fl_res_multi + d$fl_res_single
  serv   <- d$fl_retail_services + d$fl_office + d$fl_education_culture + d$fl_health
  prod   <- d$fl_industry_storage
  nonres <- tot - res
  pc <- 100  # pseudo-count (m2 of floor area) for zero parts in log-ratios
  d$ilr_res_vs_nonres   <- sqrt(1/2) * log((res + pc) / (nonres + pc))
  d$ilr_multi_vs_single <- sqrt(1/2) * log((d$fl_res_multi + pc) / (d$fl_res_single + pc))
  d$ilr_serv_vs_prod    <- sqrt(1/2) * log((serv + pc) / (prod + pc))
  ## remaining balances of the sequential binary partition (manuscript §2.7), amalgamated in the same way:
  ## services: retail and services vs workplaces; workplaces: office vs education, culture and health;
  ## production and other: industry and storage vs other non-residential
  work <- d$fl_office + d$fl_education_culture + d$fl_health
  d$ilr_retail_vs_work  <- sqrt(1/2) * log((d$fl_retail_services + pc) / (work + pc))
  d$ilr_office_vs_eduhealth <- sqrt(1/2) * log((d$fl_office + pc) / (d$fl_education_culture + d$fl_health + pc))
  d$ilr_industry_vs_other   <- sqrt(1/2) * log((d$fl_industry_storage + pc) / (d$fl_other_nonres + pc))
  d$floor_total <- tot
  p_nr <- as.matrix(d[paste0("fl_", setdiff(fns, c("res_multi","res_single")))])
  d$nonres_mix <- apply(p_nr, 1, function(x) { s <- sum(x); if (s <= 0) return(NA); p <- x[x > 0] / s; -sum(p * log(p)) / log(length(x)) })
  d$mixed_ground_share <- vapply(im, function(i) { k <- i[bm$fn[i] == "res_multi"]; if (length(k)) mean(bm$mixed_ground[k]) else NA_real_ }, numeric(1))
  ## overall land-use mix: normalised Shannon entropy of floor area across all eight function classes
  p_all <- as.matrix(d[paste0("fl_", fns)])
  d$landuse_mix <- apply(p_all, 1, function(x) { s <- sum(x); if (s <= 0) return(NA); p <- x[x > 0] / s; -sum(p * log(p)) / log(length(x)) })
  ## university buildings and student dormitories as shares of total floor area (seasonal student presence)
  d$share_university <- ifelse(tot > 0, sum_by(im, bm$floor * bm$university) / tot, NA)
  d$share_dormitory  <- ifelse(tot > 0, sum_by(im, bm$floor * bm$dormitory) / tot, NA)
  d[paste0("fl_", fns)] <- NULL

  ## ── Residential population density ───────────────────────────────────────
  pop <- st_read(file.path(out_dir, "output/pop_poznan_geocoded.gpkg"), quiet = TRUE)
  pcol <- grep("ZAMELDOWANYCH", names(pop), value = TRUE)[1]
  pop  <- st_transform(pop, 2180)
  d$res_density <- sum_by(in_buf(pop), pop[[pcol]]) / A_km2

  ## ── Amenities ────────────────────────────────────────────────────────────
  am <- st_read(gpkg, layer = "amenities", quiet = TRUE)
  iam <- in_buf(am)
  d$amenity_dens <- lengths(iam) / A_km2
  cats <- sort(unique(am$cat))
  d$amenity_div <- vapply(iam, function(i) { if (length(i) < 2) return(NA_real_); p <- table(factor(am$cat[i], levels = cats)) / length(i); p <- p[p > 0]; -sum(p * log(p)) / log(length(cats)) }, numeric(1))
  d$retail_food_dens <- vapply(iam, function(i) sum(am$cat[i] %in% c("retail", "food_drink")), numeric(1)) / A_km2

  ## ── Street network ───────────────────────────────────────────────────────
  rd <- st_read(gpkg, layer = "roads", quiet = TRUE) %>% filter(!klasa %in% c("autostrada", "droga ekspresowa"))
  rc <- st_coordinates(rd)                                  # X, Y, L1 (feature index)
  first <- !duplicated(rc[, "L1"]); last <- !duplicated(rc[, "L1"], fromLast = TRUE)
  ends <- rbind(rc[first, 1:2], rc[last, 1:2])
  key  <- paste(round(ends[, 1]), round(ends[, 2]))
  deg  <- table(key)
  nodes <- data.frame(key = names(deg), degree = as.integer(deg))
  xy <- do.call(rbind, strsplit(nodes$key, " ")); nodes$x <- as.numeric(xy[, 1]); nodes$y <- as.numeric(xy[, 2])
  nodes <- st_as_sf(nodes, coords = c("x", "y"), crs = 2180)
  inod <- in_buf(nodes)
  n_int  <- vapply(inod, function(i) sum(nodes$degree[i] >= 3), numeric(1))
  n_dead <- vapply(inod, function(i) sum(nodes$degree[i] == 1), numeric(1))
  pieces <- suppressWarnings(st_intersection(rd %>% select(klasa), bufsf))
  pieces <- suppressWarnings(st_collection_extract(pieces, "LINESTRING"))
  pieces <- st_cast(st_cast(pieces, "MULTILINESTRING", warn = FALSE), "LINESTRING", warn = FALSE)
  pieces$len <- as.numeric(st_length(pieces))
  pcrd <- st_coordinates(pieces)
  f1 <- !duplicated(pcrd[, "L1"]); f2 <- !duplicated(pcrd[, "L1"], fromLast = TRUE)
  pieces$bearing <- (atan2(pcrd[f2, "X"] - pcrd[f1, "X"], pcrd[f2, "Y"] - pcrd[f1, "Y"]) * 180 / pi) %% 180
  pp <- st_drop_geometry(pieces)
  agg <- pp %>% group_by(cell_id) %>% summarise(
    street_len = sum(len),
    main_len   = sum(len[klasa %in% c("droga główna", "droga główna ruchu przyśpieszonego", "droga zbiorcza")]),
    internal_len = sum(len[klasa == "droga wewnętrzna"]),
    orient_entropy = { h <- tapply(len, cut(bearing, seq(0, 180, by = 5), include.lowest = TRUE), sum); h <- h[!is.na(h) & h > 0]; p <- h / sum(h); -sum(p * log(p)) / log(36) },
    .groups = "drop")
  d <- d %>% left_join(agg, by = "cell_id")
  d$street_len_dens <- d$street_len / 1000 / A_km2          # km per km2
  d$intersection_dens <- n_int / A_km2
  d$deadend_share <- ifelse(n_int + n_dead > 0, n_dead / (n_int + n_dead), NA)
  d$len_per_intersection <- ifelse(n_int > 0, d$street_len / n_int, NA)
  d$main_road_share <- d$main_len / d$street_len
  d$internal_road_share <- d$internal_len / d$street_len
  d[c("street_len", "main_len", "internal_len")] <- NULL

  ## urban blocks = faces of the road network
  blk <- st_polygonize(st_union(st_geometry(rd)))
  blk <- st_collection_extract(blk, "POLYGON")
  blk <- st_sf(geometry = blk); blk$area <- as.numeric(st_area(blk))
  blk <- blk[blk$area >= 500, ]
  iblk <- in_buf(st_centroid(st_geometry(blk)))
  d$med_block_area_ha <- vapply(iblk, function(i) if (length(i)) median(blk$area[i]) / 1e4 else NA_real_, numeric(1))
  cat("Blocks:", nrow(blk), " median (ha):", round(median(blk$area) / 1e4, 2), "\n")

  ## ── Area shares via 5 m rasters (zoning, built-up character, BDOT cover) ─
  tmpl <- rast(ext(vect(city)), resolution = 5, crs = "EPSG:2180")
  as_rows <- function(f, n_out) {   # exact_extract returns a list or a matrix depending on version
    if (is.list(f) && !is.data.frame(f)) f <- do.call(rbind, f)
    f <- as.matrix(f)
    if (nrow(f) != nrow(bufsf) && ncol(f) == nrow(bufsf)) f <- t(f)
    stopifnot(nrow(f) == nrow(bufsf), ncol(f) == n_out); f
  }
  frac_of <- function(sfpoly, field, levels) {
    sfpoly$code <- match(sfpoly[[field]], levels)
    r <- rasterize(vect(sfpoly[!is.na(sfpoly$code), ]), tmpl, field = "code")
    f <- exact_extract(r, bufsf, function(v, cov) { s <- sum(cov); vapply(seq_along(levels), function(k) sum(cov[!is.na(v) & v == k]) / s, numeric(1)) }, progress = FALSE)
    f <- as_rows(f, length(levels)); colnames(f) <- levels; f
  }
  zon <- st_read(file.path(out_dir, "output/studium2023_landuse_zones.gpkg"), quiet = TRUE) %>% st_transform(2180)
  ## Studium symbols carry starred variants (ZO*, MN*, US1*, ...) of the same category: strip the star first
  zsym <- sub("\\*$", "", zon$symb_t_o)
  zon$zclass <- case_when(
    zsym == "MN" ~ "z_single", zsym == "MW" ~ "z_multi",
    zsym %in% c("MN/MW", "MW/MN", "MN/U", "MW/U") ~ "z_mixed",
    zsym %in% c("U", "US", "US1", "US2") ~ "z_services", zsym == "U/P" ~ "z_production",
    zsym %in% c("ZN", "ZO", "ZP", "ZD", "ZC", "ZS", "ZR", "UF") ~ "z_green",
    grepl("^(kd|kk|kt|IT|TKL)", zsym) ~ "z_transport_infra")
  if (any(is.na(zon$zclass))) stop("Unclassified Studium symbols: ", paste(unique(zon$symb_t_o[is.na(zon$zclass)]), collapse = ", "))
  zf <- frac_of(zon, "zclass", c("z_single", "z_multi", "z_mixed", "z_services", "z_production", "z_green", "z_transport_infra"))
  d <- cbind(d, as.data.frame(zf))

  zb <- st_read(gpkg, layer = "builtup", quiet = TRUE)
  cf <- frac_of(zb, "charakter", c("zwarta", "gęsta", "luźna"))
  bu <- rowSums(cf)
  d$builtup_compact <- ifelse(bu > 0, cf[, "zwarta"] / bu, NA)
  d$builtup_dense   <- ifelse(bu > 0, cf[, "gęsta"] / bu, NA)
  d$builtup_loose   <- ifelse(bu > 0, cf[, "luźna"] / bu, NA)

  gb <- st_read(gpkg, layer = "green_blue", quiet = TRUE)
  gf <- frac_of(gb, "cover", c("woodland", "grassland", "water", "arable"))
  d$woodland <- gf[, "woodland"]; d$grassland <- gf[, "grassland"]; d$water <- gf[, "water"]

  ## ── ESA WorldCover (satellite) ───────────────────────────────────────────
  wc <- rast(file.path(s1_dir, "worldcover_poznan.tif"))
  wf <- exact_extract(wc, bufsf, function(v, cov) { s <- sum(cov); c(tree = sum(cov[v %in% 10]), grass = sum(cov[v %in% 30]), built = sum(cov[v %in% 50])) / s }, progress = FALSE)
  wf <- as_rows(wf, 3)
  d$wc_tree <- wf[, 1]; d$wc_grass <- wf[, 2]; d$wc_built <- wf[, 3]

  ## ── Save ─────────────────────────────────────────────────────────────────
  d$area_km2 <- A_km2
  if (units_mode != "cells") {
    names(d)[names(d) == "cell_id"] <- "stop_name"
    write.csv(d, file.path(out_dir, "output/stage2", paste0("catchment_tissue_descriptors_", units_mode, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
    cat("Saved output/stage2/catchment_tissue_descriptors_", units_mode, ".csv (", nrow(d), " catchments)\n", sep = "")
    return(invisible(NULL))
  }
  out <- cells %>% left_join(d, by = "cell_id")
  st_write(out, file.path(s1_dir, "cell_tissue_descriptors.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  cat("\nDescriptor summary:\n"); print(summary(st_drop_geometry(out)[, setdiff(names(d), "cell_id")]))
  cat("\nSaved output/stage1/cell_tissue_descriptors.gpkg (", nrow(out), "cells )\n")
}
