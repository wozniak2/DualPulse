## Sharing the network (sections 80 and 82). The traced ride-km of section 79 as sparse matrices of cells
## (line x direction x segment x period) and of corridors (undirected segments) against origin catchments,
## and the indices of section 79 for any vector of catchment deprivation, so that a benchmark can recompute
## them for other allocations of the same scores over the same flows.

## thirds of catchments by deprivation (the rule of section 56)
seg_thirds <- function(d) { q <- quantile(d, c(1/3, 2/3)); list(m = as.numeric(d >= q[2]), l = as.numeric(d <= q[1])) }

## matrices of one window: CE = segregation_network_cells.rds (79), G = catchments with column o (stop name)
seg_setup <- function(CE, G, wn) {
  E <- CE[win == wn & o %in% G$o]
  E[, `:=`(vi = match(cell, unique(cell)), ci = match(corr, unique(corr)), oi = match(o, G$o))]
  Mv <- Matrix::sparseMatrix(i = E$vi, j = E$oi, x = E$w, dims = c(max(E$vi), nrow(G)))
  Mc <- Matrix::sparseMatrix(i = E$ci, j = E$oi, x = E$w, dims = c(max(E$ci), nrow(G)))
  K <- unique(E[, .(ci, th, dctr)], by = "ci")[order(ci)]
  list(Mv = Mv, Mc = Mc, vw = Matrix::rowSums(Mv), cw = Matrix::rowSums(Mc), ew = abs(cos(K$th)) >= abs(sin(K$th)), inner = K$dctr <= 2)
}

## indices for catchment deprivation d and thirds tt (m = most, l = least deprived, 0/1 per catchment)
seg_ind <- function(S, d, tt = seg_thirds(d)) {
  Mv <- S$Mv; Mc <- S$Mc; vw <- S$vw; cw <- S$cw
  vd <- as.vector(Mv %*% d); vd2 <- as.vector(Mv %*% d^2); vm <- as.vector(Mv %*% tt$m); vl <- as.vector(Mv %*% tt$l)
  cd <- as.vector(Mc %*% d); cm <- as.vector(Mc %*% tt$m); cl <- as.vector(Mc %*% tt$l)
  W <- sum(vw); mu <- sum(vd) / W; tv <- sum(vd2) / W - mu^2
  eta <- function(w, wd) { k <- w > 0; sum(w[k] * (wd[k] / w[k] - mu)^2) / W / tv }
  iso <- function(w, g) { k <- w > 0; sum(g[k]^2 / w[k]) / sum(g) }
  exl <- function(w, gm, gl) { k <- w > 0; sum(gl[k] * gm[k] / w[k]) / sum(gl) }
  sm <- cm / sum(cm); sl <- cl / sum(cl)
  Pm <- sum(vm) / W; Pl <- sum(vl) / W
  nis <- function(i, p) (i - p) / (1 - p)                                # isolation net of group size (0 = even mixing, 1 = complete)
  c(niso_most_vehicle = nis(iso(vw, vm), Pm), niso_least_vehicle = nis(iso(vw, vl), Pl),
    niso_most_corridor = nis(iso(cw, cm), Pm), niso_least_corridor = nis(iso(cw, cl), Pl),
    rel_exp_least_to_most_vehicle = exl(vw, vm, vl) / Pm, rel_exp_least_to_most_corridor = exl(cw, cm, cl) / Pm,
    eta2_vehicle = eta(vw, vd), eta2_corridor = eta(cw, cd), iso_most_vehicle = iso(vw, vm), iso_least_vehicle = iso(vw, vl),
    exp_least_to_most_vehicle = exl(vw, vm, vl), iso_most_corridor = iso(cw, cm), iso_least_corridor = iso(cw, cl),
    exp_least_to_most_corridor = exl(cw, cm, cl), overlap_corridor = sum(pmin(sm, sl)), overlap_beyond2km = sum(pmin(sm, sl)[!S$inner]),
    ew_share_most = sum(sm[S$ew]), ew_share_least = sum(sl[S$ew]), ew_most_minus_least = sum(sm[S$ew]) - sum(sl[S$ew]),
    share_most = Pm, share_least = Pl)
}

## share of benchmark draws (rows of B) at least as separated as the observed indices obs; the axis difference
## is two-sided, overlap and exposure are separated when low, shares have no direction (kept for reference)
seg_p <- function(B, obs) {
  hi <- !grepl("^overlap|^exp_|^rel_exp|^share", colnames(B))
  vapply(seq_len(ncol(B)), function(k) {
    m <- colnames(B)[k]
    if (m == "ew_most_minus_least") mean(abs(B[, k]) >= abs(obs[m])) else if (hi[k]) mean(B[, k] >= obs[m]) else mean(B[, k] <= obs[m])
  }, numeric(1))
}

## turned and mirrored maps of the osiedle scores (section 82): the map is turned around the Old Market Square in
## 15-degree steps (23 rotations) and mirrored across 24 axes through it (every 7.5 degrees). Each catchment of G
## (columns o = stop name, dep) takes the score that the transformed map places at its stop (that of the nearest
## osiedle where the transformed location falls outside the city), and its deprivation shifts by the difference
## between this score and the score of the osiedle at its own stop. Returns the maps (TR) and the deprivation
## vectors (D), with the share of stops placed outside the city and each map's correlation with the observed one.
seg_turned <- function(G, out_dir) {
  P <- data.table::fread(file.path(out_dir, "output/stage1/s1ms_osiedle_profiles.csv"), encoding = "UTF-8")[, .(osiedle, s = dep)]
  OS <- sf::st_transform(sf::st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)[, "name"]
  OS$s <- P$s[match(OS$name, P$osiedle)]
  if (anyNA(OS$s)) stop("osiedla without a score")
  city <- sf::st_union(OS)
  C <- sf::st_read(file.path(out_dir, "output/stage2/stop_complexes.gpkg"), "complexes", quiet = TRUE)
  XY <- sf::st_coordinates(C)[match(G$o, C$stop_name), , drop = FALSE]
  if (anyNA(XY)) stop("catchments without a stop location")
  ctr <- sf::st_coordinates(sf::st_transform(sf::st_sfc(sf::st_point(c(16.9335, 52.4083)), crs = 4326), 2180))[1, ]
  score_at <- function(xy) {
    pt <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = 2180)
    k <- sf::st_nearest_feature(pt, OS)
    list(s = OS$s[k], name = OS$name[k], outside = lengths(sf::st_intersects(pt, city)) == 0)
  }
  own <- score_at(XY)
  dx <- XY[, 1] - ctr[1]; dy <- XY[, 2] - ctr[2]
  inv <- function(type, a) {                                                # the score at a stop is the original score at its inverse image
    r <- a * pi / 180
    if (type == "rotation") cbind(ctr[1] + cos(r) * dx + sin(r) * dy, ctr[2] - sin(r) * dx + cos(r) * dy)
    else cbind(ctr[1] + cos(2 * r) * dx + sin(2 * r) * dy, ctr[2] + sin(2 * r) * dx - cos(2 * r) * dy)
  }
  TR <- rbind(data.table::data.table(type = "rotation", angle = seq(15, 345, 15)), data.table::data.table(type = "mirror", angle = seq(0, 172.5, 7.5)))
  Z <- lapply(seq_len(nrow(TR)), function(i) { z <- score_at(inv(TR$type[i], TR$angle[i])); list(d = G$dep + z$s - own$s, outside = mean(z$outside)) })
  TR[, `:=`(outside = vapply(Z, `[[`, numeric(1), "outside"), r_observed = vapply(Z, function(z) cor(z$d, G$dep), numeric(1)))]
  list(TR = TR, D = lapply(Z, `[[`, "d"), own_match = mean(own$name == G$osiedle))
}
