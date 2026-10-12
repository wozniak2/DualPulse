## DualPulse | stage 3: The deprivation score and where deprived areas sit in land use
##
## Land use around residents' homes aggregated to the 42 osiedla, the deprivation score (mean of the
## standardised unemployment and poverty scores, 2019 and 2021) and its land-use associations
## (main text 3.1, supplementary S2 and S4).
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   16      Core descriptors, osiedle profiles and land-use types
##   18      Land use and deprivation in each round
##   19      Area status indicators
##   20      Deprivation score and its land-use associations
##
## Run from the repository root: Rscript run.R 3   (one section: Rscript run.R 3:20)

if (!exists("section")) source("R/setup.R")

## ── 16 · Core descriptors, osiedle profiles and land-use types ─────────────────────────────────────────────
section("16", "Core descriptors, osiedle profiles and land-use types", {
  ## Stage 1 (RQ1): do the urban tissues of more and less deprived osiedla
  ## differ structurally, and which characteristics distinguish them?
  ##
  ## The 37 descriptors from 15 are narrowed to a CORE set before any test:
  ## one or two per dimension of urban tissue (intensity, height, spacing,
  ## street interface, grain, function, destinations, housing type, street
  ## connectivity, permeability, greenery), chosen on theory in a fixed
  ## priority order and pruned for redundancy (a candidate is dropped if its
  ## cell-level Spearman |rho| with an already-kept descriptor exceeds 0.8).
  ## The selection uses descriptor-descriptor correlations only, never the
  ## poverty rate. The full set is kept as supplementary material.
  ##
  ## 1. Osiedle profiles: each descriptor aggregated over the populated cells
  ##    of an osiedle with population weights (weighted medians for
  ##    per-building medians), i.e. the tissue residents actually live in.
  ## 2. Association with the 2021 poverty rate (n = 42): Spearman rho,
  ##    permutation p (9,999), bootstrap 95% CI, Benjamini-Hochberg FDR within
  ##    the core set; sensitivity without the two single-institution osiedla;
  ##    quadratic rank test for non-monotonic relationships.
  ## 3. Tercile profiles of the core set.
  ## 4. Typology: k-means on the standardised core descriptors (k chosen by
  ##    silhouette among 4-8); population share of each type per osiedle.
  ## 5. Separability: within-osiedle variance share; leave-one-out R^2 of the
  ##    osiedle rate on tissue principal components.
  ## 6. Maps: tissue types on one map; 2021 poverty rate (MOPR beneficiaries
  ##    per 1,000) by osiedle; both side by side.
  ##
  ## Run 14 and 15 first. Outputs in output/stage1/ and output/figures/.

  suppressPackageStartupMessages({library(sf); library(dplyr); library(data.table); library(ggplot2); library(cluster); library(patchwork); library(tidyr)})
  set.seed(20261005)

  out_dir <- ROOT
  s1_dir  <- file.path(out_dir, "output/stage1")
  fig_dir <- file.path(out_dir, "output/figures")

  cells <- st_read(file.path(s1_dir, "cell_tissue_descriptors.gpkg"), quiet = TRUE)
  gt    <- fread(file.path(out_dir, "data/diagnoza2021_osiedla_ground_truth.csv"), encoding = "UTF-8")
  osi   <- st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE) %>% st_transform(2180)
  INSTITUTION_OSIEDLA <- c("Główna", "Krzesiny-Pokrzywno-Garaszewo")

  ## ── Descriptor catalogue (full set) ──────────────────────────────────────
  cat_tbl <- tribble(
    ~var, ~label, ~family, ~agg,
    "gsi", "Building coverage (GSI)", "Built form", "mean",
    "fsi", "Floor-area density (FSI)", "Built form", "mean",
    "storeys", "Mean storeys", "Built form", "mean",
    "osr", "Open-space ratio (OSR)", "Built form", "median",
    "med_footprint", "Median building footprint", "Built form", "median",
    "bldg_dens", "Buildings per km2", "Built form", "mean",
    "res_density", "Residents per km2", "Built form", "mean",
    "med_gap", "Gap to nearest building", "Spacing and street interface", "median",
    "touch_share", "Buildings with shared walls", "Spacing and street interface", "mean",
    "med_setback", "Building setback from street", "Spacing and street interface", "median",
    "frontage_share", "Buildings within 5 m of street", "Spacing and street interface", "mean",
    "med_block_area_ha", "Median block area", "Spacing and street interface", "median",
    "builtup_compact", "Compact built-up land (BDOT)", "Spacing and street interface", "mean",
    "builtup_loose", "Loose built-up land (BDOT)", "Spacing and street interface", "mean",
    "ilr_res_vs_nonres", "Residential vs non-residential floor", "Function and amenities", "mean",
    "ilr_multi_vs_single", "Multi- vs single-family floor", "Function and amenities", "mean",
    "ilr_serv_vs_prod", "Services vs production floor", "Function and amenities", "mean",
    "nonres_mix", "Non-residential function mix", "Function and amenities", "mean",
    "mixed_ground_share", "Housing with commercial ground floor", "Function and amenities", "mean",
    "amenity_dens", "Amenities per km2", "Function and amenities", "mean",
    "amenity_div", "Amenity diversity", "Function and amenities", "mean",
    "z_single", "Zoned single-family", "Function and amenities", "mean",
    "z_multi", "Zoned multi-family", "Function and amenities", "mean",
    "z_mixed", "Zoned mixed residential-service", "Function and amenities", "mean",
    "z_production", "Zoned services-industry", "Function and amenities", "mean",
    "street_len_dens", "Street length density", "Street network", "mean",
    "intersection_dens", "Intersections per km2", "Street network", "mean",
    "deadend_share", "Dead-end share", "Street network", "mean",
    "len_per_intersection", "Street length per intersection", "Street network", "median",
    "main_road_share", "Main and collector road share", "Street network", "mean",
    "internal_road_share", "Internal road share", "Street network", "mean",
    "orient_entropy", "Street orientation entropy", "Street network", "mean",
    "woodland", "Woodland (BDOT)", "Green, open and sealed", "mean",
    "grassland", "Grassland (BDOT)", "Green, open and sealed", "mean",
    "z_green", "Zoned green space", "Green, open and sealed", "mean",
    "wc_tree", "Tree cover (WorldCover)", "Green, open and sealed", "mean",
    "wc_built", "Built-up cover (WorldCover)", "Green, open and sealed", "mean")
  fam_levels <- c("Built form", "Spacing and street interface", "Function and amenities", "Street network", "Green, open and sealed")

  ## core candidates in theory-driven priority order, with the tissue dimension each represents
  core_candidates <- c(fsi = "intensity", storeys = "height", med_gap = "spacing", med_setback = "street interface",
                       med_block_area_ha = "grain", ilr_res_vs_nonres = "function balance", amenity_dens = "destinations",
                       ilr_multi_vs_single = "housing type", intersection_dens = "connectivity", deadend_share = "permeability",
                       wc_tree = "greenery", res_density = "residential intensity", touch_share = "contiguity")
  cd <- st_drop_geometry(cells)
  R <- cor(cd[names(core_candidates)], method = "spearman", use = "pairwise.complete.obs")
  core <- character(0)
  for (v in names(core_candidates)) if (all(abs(R[v, core]) <= 0.8)) core <- c(core, v)
  cat("Cell-level Spearman correlations among core candidates:\n"); print(round(R, 2))
  cat("\nCore set (", length(core), "):", paste(core, collapse = ", "), "\n")
  cat("Dropped as redundant:", paste(setdiff(names(core_candidates), core), collapse = ", "), "\n")
  cat_tbl$core <- cat_tbl$var %in% core
  fwrite(as.data.table(R, keep.rownames = "var"), file.path(s1_dir, "core_candidate_correlations.csv"))

  wmedian <- function(x, w) { ok <- !is.na(x) & !is.na(w) & w > 0; x <- x[ok]; w <- w[ok]; if (!length(x)) return(NA_real_); o <- order(x); x <- x[o]; w <- cumsum(w[o]) / sum(w); x[which(w >= 0.5)[1]] }
  wmean   <- function(x, w) { ok <- !is.na(x) & w > 0; if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }

  prof <- cd %>% group_by(osiedle) %>% summarise(across(all_of(cat_tbl$var), ~ {
    a <- cat_tbl$agg[cat_tbl$var == cur_column()]
    if (a == "median") wmedian(.x, population) else wmean(.x, population) }), pop = sum(population), n_cells = n(), .groups = "drop") %>%
    left_join(gt %>% select(osiedle, rate = indicator_B_poverty_per1000), by = "osiedle")
  stopifnot(nrow(prof) == 42, !any(is.na(prof$rate)))
  fwrite(prof, file.path(s1_dir, "osiedle_tissue_profiles.csv"))

  ## ── 2. Associations with the poverty rate ───────────────────────────────
  spear <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
  assoc <- bind_rows(lapply(cat_tbl$var, function(v) {
    x <- prof[[v]]; y <- prof$rate; ok <- !is.na(x)
    rho <- spear(x[ok], y[ok])
    perm <- replicate(9999, spear(x[ok], sample(y[ok])))
    boot <- replicate(2000, { i <- sample(which(ok), replace = TRUE); spear(x[i], y[i]) })
    keep <- ok & !prof$osiedle %in% INSTITUTION_OSIEDLA
    rx <- rank(x[ok]); ry <- rank(y[ok]); m1 <- lm(ry ~ rx); m2 <- lm(ry ~ rx + I(rx^2))
    data.frame(var = v, n = sum(ok), rho = rho,
               p_perm = (1 + sum(abs(perm) >= abs(rho))) / (1 + length(perm)),
               ci_lo = quantile(boot, .025, na.rm = TRUE), ci_hi = quantile(boot, .975, na.rm = TRUE),
               rho_no_institution = spear(x[keep], y[keep]),
               p_quadratic = anova(m1, m2)$`Pr(>F)`[2], quad_sign = sign(coef(m2)[3]))
  })) %>% left_join(cat_tbl, by = "var")
  assoc$q_fdr_full <- p.adjust(assoc$p_perm, "BH")
  assoc$q_fdr <- NA_real_; assoc$q_fdr[assoc$core] <- p.adjust(assoc$p_perm[assoc$core], "BH")
  fwrite(assoc, file.path(s1_dir, "stage1_associations.csv"))
  cat("\n=== CORE: Spearman rho with the 2021 poverty rate (42 osiedla); FDR within core ===\n")
  print(assoc %>% filter(core) %>% arrange(desc(abs(rho))) %>% transmute(label, rho = round(rho, 2), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi),
        q = signif(q_fdr, 2), rho_noinst = round(rho_no_institution, 2), p_quad = signif(p_quadratic, 2)), row.names = FALSE)

  ## ── 3. Tercile profiles (core) ───────────────────────────────────────────
  prof$tercile <- cut(rank(prof$rate, ties.method = "first"), c(0, 14, 28, 42), labels = c("Least deprived", "Middle", "Most deprived"))
  zs <- prof %>% mutate(across(all_of(core), ~ as.numeric(scale(.x))))
  terc <- zs %>% group_by(tercile) %>% summarise(across(all_of(core), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  fwrite(terc, file.path(s1_dir, "tercile_profiles.csv"))

  ## ── 4. Typology on the core set ─────────────────────────────────────────
  tr <- function(v, x) switch(v, fsi = log(x + 0.01), med_gap = log1p(x), med_setback = log1p(x), med_block_area_ha = log(x),
                              amenity_dens = log1p(x), res_density = log1p(x), x)
  feat <- as.data.frame(lapply(setNames(core, core), function(v) tr(v, cd[[v]])))
  for (j in names(feat)) feat[[j]][is.na(feat[[j]])] <- median(feat[[j]], na.rm = TRUE)
  X <- scale(feat); D <- dist(X)
  sil <- sapply(4:8, function(k) { km <- kmeans(X, k, nstart = 50, iter.max = 100); mean(silhouette(km$cluster, D)[, 3]) })
  names(sil) <- 4:8; cat("\nSilhouette by k:\n"); print(round(sil, 3))
  ## parsimony: the smallest k whose silhouette is within 0.01 of the best (tissue is a continuum, so the curve is flat)
  k_best <- min(as.integer(names(sil)[sil >= max(sil) - 0.01]))
  cat("Chosen k:", k_best, "\n")
  km <- kmeans(X, k_best, nstart = 100, iter.max = 200)
  cl <- km$cluster
  raw_means <- cd %>% mutate(cl = cl) %>% group_by(cl) %>%
    summarise(cells = n(), residents = sum(population), fsi = mean(fsi), storeys = mean(storeys, na.rm = TRUE),
              med_gap = median(med_gap, na.rm = TRUE), setback = median(med_setback, na.rm = TRUE), touch = mean(touch_share, na.rm = TRUE),
              block_ha = median(med_block_area_ha, na.rm = TRUE), single_family = mean(share_res_single, na.rm = TRUE),
              amenity_dens = median(amenity_dens), int_dens = mean(intersection_dens), deadend = mean(deadend_share, na.rm = TRUE),
              tree = mean(wc_tree), mean_rate = weighted.mean(rate, population), .groups = "drop")
  ## name types from their profiles (k = 4): compact blocks = highest FSI; estates = lowest single-family share of
  ## the rest; loose peripheral = lowest FSI of the rest; the remaining one = single-family and terraced housing
  if (k_best == 4) {
    rest <- raw_means$cl; at <- function(v) raw_means[[v]][match(rest, raw_means$cl)]
    compact <- rest[which.max(at("fsi"))]; rest <- setdiff(rest, compact)
    estates <- rest[which.min(at("single_family"))]; rest <- setdiff(rest, estates)
    loose   <- rest[which.min(at("fsi"))]; rest <- setdiff(rest, loose)
    ord <- c(compact, estates, rest, loose)
    type_names <- c("Compact mixed-use blocks", "Modernist multi-family estates", "Single-family and terraced housing", "Loose peripheral tissue")
  } else {
    ord <- order(tapply(cd$fsi, cl, mean), decreasing = TRUE); type_names <- paste("Type", seq_len(k_best))
  }
  cells$type <- match(cl, ord)
  raw_means <- raw_means %>% mutate(type = match(cl, ord), name = type_names[type]) %>% arrange(type) %>% select(type, name, everything(), -cl)
  cat("\nType profiles (raw units):\n"); print(as.data.frame(raw_means), digits = 3)
  fwrite(raw_means, file.path(s1_dir, "typology_profiles.csv"))

  comp <- st_drop_geometry(cells) %>% group_by(osiedle, type) %>% summarise(pop = sum(population), .groups = "drop") %>%
    group_by(osiedle) %>% mutate(share = pop / sum(pop)) %>% ungroup() %>%
    complete(osiedle, type, fill = list(pop = 0, share = 0)) %>% left_join(gt %>% select(osiedle, rate = indicator_B_poverty_per1000), by = "osiedle")
  type_assoc <- comp %>% group_by(type) %>% summarise(rho = spear(share, rate),
    p_perm = { r0 <- spear(share, rate); (1 + sum(abs(replicate(9999, spear(share, sample(rate)))) >= abs(r0))) / 10000 }, .groups = "drop") %>%
    mutate(q_fdr = p.adjust(p_perm, "BH"), name = type_names[type])
  cat("\nType share vs poverty rate (Spearman across osiedla):\n"); print(as.data.frame(type_assoc), digits = 3)
  fwrite(comp, file.path(s1_dir, "osiedle_type_shares.csv")); fwrite(type_assoc, file.path(s1_dir, "type_share_associations.csv"))

  ## ── 5. Separability diagnostics (core) ───────────────────────────────────
  within_share <- sapply(core, function(v) {
    y <- cd[[v]]; ok <- !is.na(y)
    1 - summary(lm(y[ok] ~ factor(cd$osiedle[ok]), weights = cd$population[ok]))$r.squared })
  cat("\nShare of cell-level variance within osiedla:\n"); print(round(within_share, 2))
  P <- prof %>% select(all_of(core)) %>% mutate(across(everything(), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))
  yy <- log1p(prof$rate)
  loo <- function(ncomp) { pred <- sapply(1:42, function(i) {
    pci <- prcomp(P[-i, ], scale. = TRUE); s <- predict(pci, P[i, , drop = FALSE])[, 1:ncomp, drop = FALSE]
    m <- lm(yy[-i] ~ ., data = data.frame(pci$x[, 1:ncomp, drop = FALSE])); predict(m, newdata = data.frame(s)) })
    1 - sum((yy - pred)^2) / sum((yy - mean(yy))^2) }
  loo_r2 <- sapply(1:5, loo)
  cat("\nLOO R2 of log poverty rate on the first k tissue PCs:", paste0("k=", 1:5, ": ", round(loo_r2, 2), collapse = ", "), "\n")
  fwrite(data.frame(var = core, within_osiedle_share = within_share), file.path(s1_dir, "separability_within_share.csv"))
  fwrite(data.frame(k = 1:5, loo_r2 = loo_r2), file.path(s1_dir, "separability_loo_r2.csv"))
  st_write(cells %>% select(cell_id, osiedle, population, rate, type), file.path(s1_dir, "cell_tissue_types.gpkg"), delete_dsn = TRUE, quiet = TRUE)

  ## ── Figures (light surface; palette validated all-pairs for 4 slots) ─────
  ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; surf <- "#fcfcfb"; blue <- "#2a78d6"
  type_pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7", "#e34948", "#008300", "#eda100", "#e87ba4")[seq_len(k_best)]
  names(type_pal) <- type_names
  th <- theme_minimal(base_size = 10) + theme(
    plot.background = element_rect(fill = surf, colour = NA), panel.background = element_rect(fill = surf, colour = NA),
    panel.grid.major = element_line(colour = grid, linewidth = 0.25), panel.grid.minor = element_blank(),
    text = element_text(colour = ink), axis.text = element_text(colour = ink2), strip.text = element_text(colour = ink, face = "bold", hjust = 0),
    plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = ink2), plot.caption = element_text(colour = muted, size = 7.5),
    legend.position = "top", legend.title = element_blank())
  map_th <- theme_void(base_size = 10) + theme(plot.background = element_rect(fill = surf, colour = NA), plot.title = element_text(face = "bold", colour = ink),
    plot.subtitle = element_text(colour = ink2), plot.caption = element_text(colour = muted, size = 7.5), legend.title = element_text(colour = ink, size = 9),
    legend.text = element_text(colour = ink2), plot.margin = margin(8, 8, 8, 8))

  ## (a) core associations; (a2) full set, supplementary
  lollipop <- function(df, qcol, title, h) {
    df <- df %>% mutate(family = factor(family, fam_levels), sig = ifelse(.data[[qcol]] < 0.10, "FDR q < 0.10", "Not significant"), label = reorder(label, rho))
    ggplot(df, aes(x = rho, y = label)) + geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.4) +
      geom_segment(aes(x = ci_lo, xend = ci_hi, yend = label), colour = "#c3c2b7", linewidth = 0.6) +
      geom_point(aes(colour = sig), size = 2.4) + scale_colour_manual(values = c("FDR q < 0.10" = blue, "Not significant" = muted), drop = FALSE) +
      facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
      labs(x = "Spearman correlation with the 2021 poverty rate (42 osiedla), bootstrap 95% CI", y = NULL, title = title,
           caption = "Population-weighted osiedle profiles of 400 m neighbourhoods around populated 200 m cells. Permutation p-values, Benjamini-Hochberg FDR.") +
      th + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 8))
  }
  ggsave(file.path(fig_dir, "Fig_stage1_tissue_correlations.png"), lollipop(assoc %>% filter(core), "q_fdr", "Core urban-tissue indicators and deprivation", 6),
         width = 9, height = 5.5, dpi = 200, bg = surf)
  ggsave(file.path(fig_dir, "Fig_stage1_tissue_correlations_full.png"), lollipop(assoc, "q_fdr_full", "All urban-tissue descriptors and deprivation (supplementary)", 10),
         width = 9, height = 10.5, dpi = 200, bg = surf)

  ## (b) tercile heatmap (core)
  tl <- terc %>% pivot_longer(-tercile, names_to = "var", values_to = "z") %>% left_join(cat_tbl, by = "var") %>%
    mutate(label = factor(label, rev(cat_tbl$label[cat_tbl$var %in% core])))
  p3 <- ggplot(tl, aes(x = tercile, y = label, fill = z)) + geom_tile(colour = surf, linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f", z)), size = 2.8, colour = ink) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-1, 1), oob = scales::squish, name = "Mean z-score") +
    labs(x = NULL, y = NULL, title = "Core tissue profile by deprivation tercile (14 osiedla each)") +
    th + theme(legend.position = "right", legend.title = element_text(), panel.grid = element_blank())
  ggsave(file.path(fig_dir, "Fig_stage1_tercile_profiles.png"), p3, width = 7.5, height = 5, dpi = 200, bg = surf)

  ## (c) type composition of each osiedle, ordered by poverty rate
  cp <- comp %>% mutate(lbl = sprintf("%s (%.0f)", osiedle, rate), type = factor(type_names[type], type_names))
  cp$lbl <- reorder(cp$lbl, cp$rate)
  p5 <- ggplot(cp, aes(x = share, y = lbl, fill = type)) + geom_col(width = 0.8, colour = surf, linewidth = 0.3) +
    scale_fill_manual(values = type_pal) + scale_x_continuous(labels = scales::percent) + guides(fill = guide_legend(nrow = 2)) +
    labs(x = "Share of registered residents", y = "Osiedle (2021 poverty rate per 1,000), most deprived at top",
         title = "Urban-tissue composition of osiedla along the deprivation gradient") + th
  ggsave(file.path(fig_dir, "Fig_stage1_osiedle_composition.png"), p5, width = 9, height = 10, dpi = 200, bg = surf)

  ## (d) maps: tissue types on one map; 2021 poverty rate by osiedle; combined
  grid_all <- st_read(file.path(out_dir, "output/income_map_grid_200m.gpkg"), quiet = TRUE)
  tmap <- grid_all %>% inner_join(st_drop_geometry(cells) %>% select(cell_id, type), by = "cell_id") %>% mutate(type = factor(type_names[type], type_names))
  city_outline <- st_union(osi)
  m_types <- ggplot() +
    geom_sf(data = city_outline, fill = "#f3f2ee", colour = "#c3c2b7", linewidth = 0.3) +
    geom_sf(data = tmap, aes(fill = type), colour = NA) +
    geom_sf(data = osi, fill = NA, colour = "#ffffff", linewidth = 0.35) +
    scale_fill_manual(values = type_pal, name = "Urban-tissue type") + guides(fill = guide_legend(ncol = 1)) +
    labs(title = "Urban-tissue types", subtitle = "Populated 200 m cells; white lines = osiedle boundaries") + map_th +
    theme(legend.position = "bottom")
  osi_rate <- osi %>% left_join(gt %>% select(name = osiedle, rate = indicator_B_poverty_per1000), by = "name")
  lab_pts <- suppressWarnings(st_point_on_surface(osi_rate))
  m_rate <- ggplot(osi_rate) +
    geom_sf(aes(fill = rate), colour = "#ffffff", linewidth = 0.35) +
    geom_sf_text(data = lab_pts, aes(label = round(rate), colour = rate > 35), size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "#ffffff"), guide = "none") +
    scale_fill_gradientn(colours = c("#fbe7e6", "#f3a9a7", "#e34948", "#a32a29", "#5e1515"), name = "MOPR beneficiaries\nper 1,000 residents, 2021",
                         guide = guide_colourbar(barwidth = 12, barheight = 0.6, title.position = "top")) +
    labs(title = "Deprivation by osiedle (MOPR, 2021)", subtitle = "Social-assistance beneficiaries per 1,000 residents; GPR monitoring report") +
    map_th + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "Fig_stage1_typology_single_map.png"), m_types, width = 7.5, height = 8.5, dpi = 200, bg = surf)
  ggsave(file.path(fig_dir, "Fig_deprivation_osiedla_2021.png"), m_rate, width = 7.5, height = 8.5, dpi = 200, bg = surf)
  ggsave(file.path(fig_dir, "Fig_stage1_types_vs_deprivation_maps.png"), m_types + m_rate + plot_layout(ncol = 2) &
           theme(plot.background = element_rect(fill = surf, colour = NA)),
         width = 14, height = 8.5, dpi = 200, bg = surf)

  cat("\nSaved Stage 1 tables to output/stage1/ and figures to output/figures/\n")
})

## ── 18 · Land use and deprivation in each round ────────────────────────────────────────────────────────────
section("18", "Land use and deprivation in each round", {
  ## Land use and deprivation: which kinds of land use are more deprived
  ## osiedla made of, and what do their residents live next to?
  ##
  ## Two views per osiedle, each related to the poverty rate (MOPR
  ## beneficiaries per 1,000) of the 2015, 2019 and 2021 GPR rounds:
  ## A. Whole-osiedle composition: area shares of Studium 2023 zoning classes
  ##    and of BDOT10k built-up land types, and shares of building floor area
  ##    by registered function, within the osiedle boundary.
  ## B. Residents' surroundings (population-weighted over populated 200 m
  ##    cells): share of the 400 m neighbourhood that is industrial-storage
  ##    built-up land or rail land, share of residents living within 200 m of
  ##    each, and the floor-area function shares from 15.
  ## Spearman rho with permutation p (4,999 permutations), BH-FDR across the
  ## variables of this script within each round; "mean rank" = the rate
  ## averaged over the within-round ranks of the three rounds, which damps
  ## single-round anomalies (e.g. Jezyce 2021).
  ##
  ## Run 14 and 15 first. Outputs: output/stage1/landuse_osiedle.csv,
  ## output/stage1/landuse_associations.csv, output/figures/Fig_landuse_*.png

  suppressPackageStartupMessages({library(sf); library(dplyr); library(data.table); library(ggplot2); library(tidyr); library(patchwork)})
  set.seed(20261005)
  out_dir <- ROOT
  s1_dir  <- file.path(out_dir, "output/stage1"); fig_dir <- file.path(out_dir, "output/figures")
  gpkg <- file.path(s1_dir, "tissue_inputs.gpkg")

  osi <- st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE) %>% st_transform(2180) %>% select(osiedle = name)
  osi$area <- as.numeric(st_area(osi))
  rd <- function(f, y) fread(file.path(out_dir, "data", f), encoding = "UTF-8")[, .(osiedle, rate = indicator_B_poverty_per1000, round = y)]
  rates <- rbind(rd("diagnoza2016_osiedla_ground_truth.csv", "2015"), rd("diagnoza2019_osiedla_ground_truth.csv", "2019"),
                 rd("diagnoza2021_osiedla_ground_truth.csv", "2021"))
  rates[, rk := frank(rate), by = round]
  mr <- rates[, .(rate = mean(rk)), by = osiedle][, round := "mean rank"]
  rates_all <- rbind(rates[, .(osiedle, rate, round)], mr)
  stopifnot(setequal(osi$osiedle, mr$osiedle))

  area_shares <- function(poly, cls, prefix) {
    x <- suppressWarnings(st_intersection(st_make_valid(poly[, cls]), osi[, "osiedle"]))
    x$a <- as.numeric(st_area(x))
    st_drop_geometry(x) %>% group_by(osiedle, k = paste0(prefix, .data[[cls]])) %>% summarise(a = sum(a), .groups = "drop") %>%
      left_join(st_drop_geometry(osi), by = "osiedle") %>% mutate(v = a / area) %>% select(osiedle, k, v) %>%
      pivot_wider(names_from = k, values_from = v, values_fill = 0)
  }

  ## ── A. Whole-osiedle composition ────────────────────────────────────────
  zon <- st_read(file.path(out_dir, "output/studium2023_landuse_zones.gpkg"), quiet = TRUE) %>% st_transform(2180)
  zsym <- sub("\\*$", "", zon$symb_t_o)
  zon$cls <- case_when(
    zsym == "MN" ~ "single_family", zsym == "MW" ~ "multi_family", zsym %in% c("MN/MW", "MW/MN") ~ "mixed_housing",
    zsym %in% c("MN/U", "MW/U") ~ "housing_services", zsym %in% c("U", "US", "US1", "US2") ~ "services_sport",
    zsym == "U/P" ~ "services_production", zsym %in% c("ZN", "ZO", "ZP", "ZD", "ZC", "ZS", "ZR", "UF") ~ "green",
    zsym == "kk" ~ "rail", grepl("^kd", zsym) ~ "roads", grepl("^(kt|IT|TKL)", zsym) ~ "other_infrastructure")
  stopifnot(!any(is.na(zon$cls)))
  A_zon <- area_shares(zon, "cls", "zon_")

  bu <- st_read(gpkg, "builtup", quiet = TRUE)
  bu$cls <- recode(bu$rodzaj, "jednorodzinna" = "single_family", "wielorodzinna" = "multi_family", "przemysłowo-składowa" = "industrial",
                   "handlowo-usługowa" = "commercial", "pozostała zabudowa" = "other")
  A_bu <- area_shares(bu, "cls", "bdot_")

  bld <- st_read(gpkg, "buildings", quiet = TRUE) %>% filter(fn != "auxiliary")
  bc <- st_join(st_centroid(bld[, c("fn", "floor")]), osi[, "osiedle"], join = st_within) %>% st_drop_geometry() %>% filter(!is.na(osiedle))
  A_fl <- bc %>% group_by(osiedle, k = paste0("floor_", fn)) %>% summarise(f = sum(floor), .groups = "drop") %>%
    group_by(osiedle) %>% mutate(v = f / sum(f)) %>% select(-f) %>% pivot_wider(names_from = k, values_from = v, values_fill = 0) %>% ungroup()

  ## ── B. Residents' surroundings ───────────────────────────────────────────
  cells <- st_read(file.path(s1_dir, "cell_tissue_descriptors.gpkg"), quiet = TRUE)
  city <- st_union(osi)
  cen <- st_centroid(st_geometry(cells))
  buf <- st_intersection(st_buffer(cen, 400), city)
  ind <- st_union(st_make_valid(bu[bu$cls == "industrial", ]))
  rail <- st_union(st_make_valid(zon[zon$cls == "rail", ]))
  share_in <- function(target) { i <- st_intersection(st_sf(id = seq_along(buf), geometry = buf), target)
    out <- numeric(length(buf)); out[i$id] <- as.numeric(st_area(i)); out / as.numeric(st_area(buf)) }
  cells$nb_industrial <- share_in(ind)
  cells$nb_rail <- share_in(rail)
  cells$near_industrial <- as.numeric(st_distance(cen, ind)) <= 200
  cells$near_rail <- as.numeric(st_distance(cen, rail)) <= 200
  cat("Residents within 200 m of industrial land:", round(weighted.mean(cells$near_industrial, cells$population), 3),
      "; of rail land:", round(weighted.mean(cells$near_rail, cells$population), 3), "\n")

  wm <- function(x, w) { ok <- !is.na(x); sum(x[ok] * w[ok]) / sum(w[ok]) }
  B <- st_drop_geometry(cells) %>% group_by(osiedle) %>% summarise(
    res_nb_industrial = wm(nb_industrial, population), res_nb_rail = wm(nb_rail, population),
    res_near_industrial = wm(near_industrial, population), res_near_rail = wm(near_rail, population),
    res_floor_industry = wm(share_industry_storage, population), res_floor_retail = wm(share_retail_services, population),
    res_floor_office = wm(share_office, population), res_floor_educ = wm(share_education_culture, population),
    res_floor_multi = wm(share_res_multi, population), res_floor_single = wm(share_res_single, population),
    res_tree = wm(wc_tree, population), .groups = "drop")

  lu <- st_drop_geometry(osi) %>% select(osiedle) %>% left_join(A_zon, by = "osiedle") %>% left_join(A_bu, by = "osiedle") %>%
    left_join(A_fl, by = "osiedle") %>% left_join(B, by = "osiedle") %>% mutate(across(-osiedle, ~ coalesce(.x, 0)))
  fwrite(lu, file.path(s1_dir, "landuse_osiedle.csv"))

  ## ── Associations ──────────────────────────────────────────────────────────
  vars <- setdiff(names(lu), c("osiedle", "zon_NA"))
  lab <- c(zon_single_family = "Zoned single-family", zon_multi_family = "Zoned multi-family", zon_mixed_housing = "Zoned mixed housing",
           zon_housing_services = "Zoned housing with services", zon_services_sport = "Zoned services and sport",
           zon_services_production = "Zoned services-production", zon_green = "Zoned green", zon_rail = "Zoned rail",
           zon_roads = "Zoned roads", zon_other_infrastructure = "Zoned airport, tram, utilities",
           bdot_single_family = "Single-family built-up land", bdot_multi_family = "Multi-family built-up land",
           bdot_industrial = "Industrial-storage built-up land", bdot_commercial = "Commercial built-up land", bdot_other = "Other built-up land",
           floor_res_multi = "Floor: multi-family housing", floor_res_single = "Floor: single-family housing", floor_industry_storage = "Floor: industry and storage",
           floor_retail_services = "Floor: retail and services", floor_office = "Floor: office", floor_education_culture = "Floor: education and culture",
           floor_health = "Floor: health", floor_other_nonres = "Floor: other non-residential",
           res_nb_industrial = "Industrial land within 400 m", res_nb_rail = "Rail land within 400 m",
           res_near_industrial = "Residents within 200 m of industrial land", res_near_rail = "Residents within 200 m of rail land",
           res_floor_industry = "Nearby floor: industry", res_floor_retail = "Nearby floor: retail and services", res_floor_office = "Nearby floor: office",
           res_floor_educ = "Nearby floor: education and culture", res_floor_multi = "Nearby floor: multi-family", res_floor_single = "Nearby floor: single-family",
           res_tree = "Tree cover within 400 m")
  view_of <- function(v) ifelse(startsWith(v, "res_"), "B. Residents' surroundings", "A. Whole osiedle")
  spear <- function(x, y) cor(x, y, method = "spearman")
  res <- rbindlist(lapply(c("2015", "2019", "2021", "mean rank"), function(r) {
    p <- merge(as.data.table(lu), rates_all[round == r], by = "osiedle")
    rbindlist(lapply(vars, function(v) { x <- p[[v]]; rho <- spear(x, p$rate)
      perm <- replicate(4999, spear(x, sample(p$rate)))
      keep <- !p$osiedle %in% c("Główna", "Krzesiny-Pokrzywno-Garaszewo")
      data.table(round = r, var = v, rho = rho, p_perm = (1 + sum(abs(perm) >= abs(rho))) / 5000, rho_no_institution = spear(x[keep], p$rate[keep])) }))
  }))
  res[, q_fdr := p.adjust(p_perm, "BH"), by = round]
  res[, `:=`(label = lab[var], view = view_of(var))]
  fwrite(res, file.path(s1_dir, "landuse_associations.csv"))
  w <- dcast(res, view + label ~ round, value.var = "rho")
  q <- dcast(res, label ~ round, value.var = "q_fdr"); setnames(q, -1, paste0("q_", names(q)[-1]))
  w <- merge(w, q, by = "label")[order(-abs(`mean rank`))]
  cat("\nSpearman rho with the poverty rate (q = BH-FDR within round):\n")
  print(w[, .(view = substr(view, 1, 1), label, r2015 = round(`2015`, 2), r2019 = round(`2019`, 2), r2021 = round(`2021`, 2),
              mean_rank = round(`mean rank`, 2), q_mean = signif(`q_mean rank`, 2), q_2021 = signif(q_2021, 2))], nrows = 100)

  cat("\nLand use of the most and least deprived osiedla (mean rank over rounds):\n")
  show <- c("zon_single_family", "zon_multi_family", "zon_housing_services", "zon_services_production", "zon_rail", "zon_green",
            "bdot_industrial", "res_near_industrial", "res_near_rail", "res_tree")
  tab <- merge(as.data.table(lu)[, c("osiedle", show), with = FALSE], dcast(rates, osiedle ~ round, value.var = "rate"), by = "osiedle")
  tab <- merge(tab, mr[, .(osiedle, mean_rank = rate)], by = "osiedle")[order(-mean_rank)]
  print(tab[c(1:8, 35:42), c("osiedle", "2015", "2019", "2021", show), with = FALSE], digits = 2)
  fwrite(tab, file.path(s1_dir, "landuse_osiedle_ranked.csv"))

  ## centrality check: do the land-use signals reduce to distance from the Old Market Square?
  rynek <- st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180)
  cells$dist_km <- as.numeric(st_distance(cen, rynek)) / 1000
  dc <- st_drop_geometry(cells) %>% group_by(osiedle) %>% summarise(dist_km = weighted.mean(dist_km, population))
  pc <- function(x, y, z) cor(resid(lm(rank(x) ~ rank(z))), resid(lm(rank(y) ~ rank(z))))
  cp <- merge(merge(as.data.table(lu), as.data.table(dc), by = "osiedle"), rates_all, by = "osiedle", allow.cartesian = TRUE)
  cent <- cp[, .(rho_distance = spear(dist_km, rate),
                 office = spear(res_floor_office, rate), office_net_of_distance = pc(res_floor_office, rate, dist_km),
                 housing_services = spear(zon_housing_services, rate), housing_services_net = pc(zon_housing_services, rate, dist_km),
                 tree = spear(res_tree, rate), tree_net_of_distance = pc(res_tree, rate, dist_km)), by = round]
  cat("\nCentrality check (partial rank correlations controlling for residents' mean distance to the Old Market Square):\n"); print(cent, digits = 2)
  fwrite(cent, file.path(s1_dir, "landuse_centrality_check.csv"))

  ## figure: four land-use measures vs the deprivation rank (mean of the three rounds' ranks)
  suppressPackageStartupMessages(library(ggrepel))
  ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; surf <- "#fcfcfb"; blue <- "#2a78d6"
  panels <- list(c("res_floor_office", "Office floor near residents", "% of floor area within 400 m"),
                 c("zon_housing_services", "Zoned housing with services", "% of osiedle area (Studium MW/U, MN/U)"),
                 c("res_tree", "Tree cover near residents", "% of land within 400 m (WorldCover)"),
                 c("res_near_industrial", "Living next to industrial land", "% of residents within 200 m"))
  pd <- merge(as.data.table(lu), mr[, .(osiedle, mean_rank = rate)], by = "osiedle")
  plots <- lapply(panels, function(pp) {
    v <- pp[1]; rr <- res[var == v]
    sub <- sprintf("rho = %.2f / %.2f / %.2f (2015 / 2019 / 2021)", rr[round == "2015", rho], rr[round == "2019", rho], rr[round == "2021", rho])
    dd <- pd[, .(osiedle, x = 100 * get(v), y = mean_rank)]
    dd[, lab := ifelse(y >= 33 | y <= 6 | x >= quantile(x, 0.93), osiedle, "")]
    ggplot(dd, aes(x, y)) + geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#c3c2b7", linewidth = 0.6, span = 1) +
      geom_point(colour = blue, size = 2.2) +
      geom_text_repel(aes(label = lab), size = 2.5, colour = ink2, min.segment.length = 0.2, segment.colour = "#c3c2b7", max.overlaps = 30, seed = 1) +
      labs(title = pp[2], subtitle = sub, x = pp[3], y = NULL) +
      theme_minimal(base_size = 9.5) + theme(plot.background = element_rect(fill = surf, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = grid, linewidth = 0.25), axis.text = element_text(colour = ink2),
        plot.title = element_text(face = "bold", colour = ink), plot.subtitle = element_text(colour = ink2, size = 8))
  })
  fig <- wrap_plots(plots, ncol = 2) + plot_annotation(
    title = "Land use and deprivation across Poznań's 42 osiedla",
    subtitle = "Vertical axis: deprivation rank, mean of the 2015, 2019 and 2021 ranks of the MOPR poverty rate (42 = most deprived)",
    caption = "Residents' surroundings are population-weighted over populated 200 m cells. Grey line: loess smoother.",
    theme = theme(plot.background = element_rect(fill = surf, colour = NA), plot.title = element_text(face = "bold", colour = ink),
                  plot.subtitle = element_text(colour = ink2), plot.caption = element_text(colour = muted, size = 7.5))) &
    theme(plot.background = element_rect(fill = surf, colour = NA))
  ggsave(file.path(fig_dir, "Fig_landuse_vs_deprivation.png"), fig, width = 11, height = 9, dpi = 200, bg = surf)
  cat("\nSaved landuse_*.csv and Fig_landuse_vs_deprivation.png\n")
})

## ── 19 · Area status indicators ────────────────────────────────────────────────────────────────────────────
section("19", "Area status indicators", {
  ## A fuller picture of osiedle status: a second deprivation indicator and a
  ## housing-market indicator alongside the MOPR poverty rate.
  ##
  ## 1. GPR indicators (city's own tables, all 42 osiedla, three rounds,
  ##    transcribed from GPR2023_ocena_aktualnosci.pdf pp. 29, 31, 33 into
  ##    data/gpr_indicators_ABC_osiedla.csv and checked against the printed
  ##    z-scores): A registered unemployment per 1,000, B MOPR beneficiaries per
  ##    1,000, C election turnout. Economic deprivation = mean of the city's own
  ##    standardised A and B, signed so that higher = more deprived; pooled =
  ##    mean over the two most recent rounds, 2019 and 2021 (damps single-round
  ##    errors; the three-round mean is kept as dep_AB_3r). C is kept as
  ##    a separate civic-participation dimension (sensitivity: A+B+C).
  ## 2. Apartment price index per osiedle from RCN (live WFS layer, Poznan
  ##    city): residential units, free-market arm's-length sales, whole
  ##    ownership share; hedonic log price per m2 on osiedle + year + unit size
  ##    (spline) + rooms + floor class; index = osiedle effect relative to the
  ##    mean osiedle, in %. Windows: 2014-16 (2015 round), 2018-20 (2019),
  ##    2020-22 (2021), and pooled 2014-22. Market type (new build vs resale) is
  ##    not recorded for most sales, so new-build premiums stay in the index.
  ## 3. Typology: pooled deprivation x pooled price index (split at medians).
  ## 4. Do the Stage 1 tissue and land-use associations change when the
  ##    outcome is the composite or the price index instead of B alone?
  ##
  ## Prices partly capitalise urban form (centrality, greenery, building
  ## stock), so the price index is reported as a separate dimension of area
  ## status, never merged into the deprivation index used as an exposure.
  ##
  ## Run 15, 16 and 18 first. Outputs: output/stage1/area_status_*.csv,
  ## output/figures/Fig_area_status_*.png

  suppressPackageStartupMessages({library(sf); library(dplyr); library(data.table); library(ggplot2); library(splines); library(patchwork); library(ggrepel)})
  set.seed(20261005)
  out_dir <- ROOT
  s1_dir <- file.path(out_dir, "output/stage1"); fig_dir <- file.path(out_dir, "output/figures")
  rcn_gpkg <- RCN_GPKG

  osi <- st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE) %>% st_transform(2180) %>% select(osiedle = name)
  names_lp <- fread(file.path(out_dir, "data/diagnoza2021_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, osiedle)]
  pop <- fread(file.path(out_dir, "data/diagnoza2016_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, population)]

  ## ── 1. GPR indicators and the deprivation composite ─────────────────────
  g <- merge(fread(file.path(out_dir, "data/gpr_indicators_ABC_osiedla.csv")), names_lp, by = "lp")
  bz <- rbind(fread(file.path(out_dir, "data/diagnoza2016_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, round = 2015L, B_z = NA_real_)],
              fread(file.path(out_dir, "data/diagnoza2019_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, round = 2019L, B_z = indicator_B_std_zscore)],
              fread(file.path(out_dir, "data/diagnoza2021_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, round = 2021L, B_z = indicator_B_std_zscore)])
  g <- merge(g, bz, by = c("lp", "round"))
  ## the 2015 B z-scores are not in our CSV; standardise like the city does (mean of osiedla, sample SD, sign: high = better)
  g[round == 2015, B_z := -(B_poverty_per1000 - mean(B_poverty_per1000)) / sd(B_poverty_per1000)]
  g[, `:=`(dep_AB = -(A_z + B_z) / 2, dep_ABC = -(A_z + B_z + C_z) / 3)]
  cat("Spearman correlations between indicators, by round:\n")
  print(g[, .(A_B = cor(A_unemployment_per1000, B_poverty_per1000, method = "spearman"),
              A_C = cor(A_unemployment_per1000, C_turnout_pct, method = "spearman"),
              B_C = cor(B_poverty_per1000, C_turnout_pct, method = "spearman"),
              AB_vs_ABC = cor(dep_AB, dep_ABC, method = "spearman")), by = round], digits = 2)
  stab <- dcast(g, osiedle ~ round, value.var = c("A_unemployment_per1000", "B_poverty_per1000", "dep_AB"))
  cat("\nRank stability across rounds (2015-19, 2019-21, 2015-21):\n")
  for (v in c("A_unemployment_per1000", "B_poverty_per1000", "dep_AB")) { s <- stab[, paste0(v, "_", c(2015, 2019, 2021)), with = FALSE]
    cat(sprintf("  %-24s %.2f  %.2f  %.2f\n", v, cor(s[[1]], s[[2]], method = "spearman"), cor(s[[2]], s[[3]], method = "spearman"), cor(s[[1]], s[[3]], method = "spearman"))) }
  g[, B_rk := frank(B_poverty_per1000), by = round]
  r2 <- c(2019L, 2021L)  # primary measure: the two most recent rounds
  pooled <- g[, .(dep_AB = mean(dep_AB[round %in% r2]), dep_AB_3r = mean(dep_AB), dep_ABC = mean(dep_ABC[round %in% r2]), A_z = mean(A_z[round %in% r2]), B_z = mean(B_z[round %in% r2]), C_z = mean(C_z[round %in% r2]),
                  A_mean = mean(A_unemployment_per1000[round %in% r2]), B_rank = mean(B_rk[round %in% r2])), by = osiedle]

  ## ── 2. Apartment price index ─────────────────────────────────────────────
  lok <- st_read(rcn_gpkg, quiet = TRUE, query = paste(
    "SELECT tran_lokalny_id_iip, tran_rodzaj_trans, tran_cena_brutto, dok_data, nier_udzial, lok_funkcja, lok_liczba_izb, lok_nr_kond, lok_pow_uzyt, lok_cena_brutto, geom",
    "FROM lokale WHERE teryt = '3064'")) %>% st_transform(2180)
  lok <- st_join(lok, osi, join = st_within)
  L <- as.data.table(st_drop_geometry(lok))
  L[, `:=`(year = as.integer(substr(dok_data, 1, 4)), area = as.numeric(lok_pow_uzyt),
           price = as.numeric(fifelse(is.na(lok_cena_brutto), tran_cena_brutto, lok_cena_brutto)),
           rooms = suppressWarnings(as.integer(lok_liczba_izb)), floor = suppressWarnings(as.integer(lok_nr_kond)))]
  n0 <- nrow(L)
  L <- L[lok_funkcja == "mieszkalna" & tran_rodzaj_trans == "wolnyRynek" & (is.na(nier_udzial) | nier_udzial == "1/1") &
         !is.na(osiedle) & between(year, 2014, 2022) & between(area, 18, 250) & price > 0]
  L[, ppsm := price / area]
  L <- L[, .SD[between(ppsm, quantile(ppsm, 0.01), quantile(ppsm, 0.99))], by = year]
  L[, `:=`(rooms_c = factor(fifelse(is.na(rooms), "NA", fifelse(rooms >= 5, "5+", as.character(pmax(rooms, 1))))),
           floor_c = factor(fifelse(is.na(floor), "NA", fifelse(floor <= 0, "ground", fifelse(floor <= 3, "1-3", fifelse(floor <= 9, "4-9", "10+"))))))]
  cat("\nApartment sales used:", nrow(L), "of", n0, "city unit records (2014-2022, residential, free market, whole share, trimmed)\n")

  price_index <- function(d, label) {
    m <- lm(log(ppsm) ~ 0 + osiedle + factor(year) + ns(log(area), 3) + rooms_c + floor_c, data = d)
    b <- coef(m)[grepl("^osiedle", names(coef(m)))]; se <- sqrt(diag(vcov(m)))[names(b)]
    data.table(osiedle = sub("^osiedle", "", names(b)), b = b, se = se)[, `:=`(index = 100 * (exp(b - mean(b)) - 1), window = label)] %>%
      merge(d[, .(n = .N, median_ppsm = median(ppsm)), by = osiedle], by = "osiedle")
  }
  pi_all <- rbind(price_index(L[between(year, 2014, 2016)], "2015"), price_index(L[between(year, 2018, 2020)], "2019"),
                  price_index(L[between(year, 2020, 2022)], "2021"), price_index(L, "pooled"))
  cat("Osiedla with fewer than 30 sales in a window:\n"); print(pi_all[n < 30, .(window, osiedle, n)])
  pw <- dcast(pi_all, osiedle ~ window, value.var = "index")
  cat("Price index rank stability (2015-19, 2019-21, 2015-21):",
      round(c(cor(pw$`2015`, pw$`2019`, method = "spearman"), cor(pw$`2019`, pw$`2021`, method = "spearman"), cor(pw$`2015`, pw$`2021`, method = "spearman")), 2), "\n")
  fwrite(pi_all, file.path(s1_dir, "area_status_price_index.csv"))

  ## ── 3. Indicators together; typology ─────────────────────────────────────
  pr <- merge(g, pi_all[window != "pooled", .(osiedle, round = as.integer(window), price_index = index)], by = c("osiedle", "round"))
  cat("\nSpearman: price index vs deprivation indicators, by round:\n")
  print(pr[, .(price_vs_A = cor(price_index, A_unemployment_per1000, method = "spearman"), price_vs_B = cor(price_index, B_poverty_per1000, method = "spearman"),
               price_vs_depAB = cor(price_index, dep_AB, method = "spearman"), price_vs_turnout = cor(price_index, C_turnout_pct, method = "spearman")), by = round], digits = 2)
  st <- merge(pooled, pi_all[window == "pooled", .(osiedle, price_index = index, n_sales = n, median_ppsm)], by = "osiedle")
  st <- merge(st, merge(names_lp, pop, by = "lp")[, .(osiedle, population)], by = "osiedle")
  cat("Pooled: price vs economic deprivation rho =", round(cor(st$price_index, st$dep_AB, method = "spearman"), 2),
      "; vs turnout =", round(cor(st$price_index, st$C_z, method = "spearman"), 2), "\n")
  type_lv <- c("Deprived, lower prices", "Deprived, higher prices", "Less deprived, higher prices", "Less deprived, lower prices")
  st[, type := factor(fcase(dep_AB >= median(dep_AB) & price_index < median(price_index), type_lv[1],
                            dep_AB >= median(dep_AB), type_lv[2], price_index >= median(price_index), type_lv[3], default = type_lv[4]), type_lv)]
  cat("\nTypology (pooled deprivation x pooled price index, median splits):\n")
  for (t in type_lv) { s <- st[type == t][order(-dep_AB)]; cat(sprintf("  %s (%d osiedla, %.0fk residents): %s\n", t, nrow(s), sum(s$population) / 1000, paste(s$osiedle, collapse = ", "))) }
  fwrite(st, file.path(s1_dir, "area_status_osiedla.csv"))

  ## ── 4. Tissue and land use against the three outcomes ────────────────────
  prof <- fread(file.path(s1_dir, "osiedle_tissue_profiles.csv"), encoding = "UTF-8")
  lu <- fread(file.path(s1_dir, "landuse_osiedle.csv"), encoding = "UTF-8")
  assoc16 <- fread(file.path(s1_dir, "stage1_associations.csv"), encoding = "UTF-8")
  core <- assoc16[core == TRUE, .(var, label)]
  lu_vars <- data.table(var = c("res_floor_office", "zon_housing_services", "zon_single_family", "res_near_industrial", "res_near_rail", "bdot_industrial"),
                        label = c("Office floor near residents", "Zoned housing with services", "Zoned single-family", "Residents within 200 m of industrial land",
                                  "Residents within 200 m of rail land", "Industrial-storage built-up land"))
  X <- merge(merge(prof[, c("osiedle", core$var), with = FALSE], lu[, c("osiedle", lu_vars$var), with = FALSE], by = "osiedle"), st, by = "osiedle")
  outcomes <- c(B_rank = "Poverty rate only (B)", dep_AB = "Unemployment + poverty (A+B)", price_index = "Apartment price index")
  vv <- rbind(core[, block := "Core urban tissue"], lu_vars[, block := "Land use"])
  cmp <- rbindlist(lapply(names(outcomes), function(o) rbindlist(lapply(vv$var, function(v) {
    x <- X[[v]]; y <- X[[o]]; ok <- !is.na(x); rho <- cor(x[ok], y[ok], method = "spearman")
    perm <- replicate(4999, cor(x[ok], sample(y[ok]), method = "spearman"))
    data.table(outcome = outcomes[[o]], var = v, rho = rho, p = (1 + sum(abs(perm) >= abs(rho))) / 5000) }))))
  cmp <- merge(cmp, vv, by = "var"); cmp[, q := p.adjust(p, "BH"), by = .(outcome, block)]
  cat("\nSpearman rho (pooled over rounds) of tissue and land use with each outcome (q = BH within outcome x block):\n")
  print(dcast(cmp, block + label ~ outcome, value.var = "rho")[order(block)], digits = 2)
  cat("\nq < 0.10:\n"); print(cmp[q < 0.10, .(outcome, label, rho = round(rho, 2), q = signif(q, 2))][order(outcome)])
  fwrite(cmp, file.path(s1_dir, "area_status_tissue_associations.csv"))

  ## ── Figures ──────────────────────────────────────────────────────────────
  ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; surf <- "#fcfcfb"
  map_th <- theme_void(base_size = 9.5) + theme(plot.background = element_rect(fill = surf, colour = NA), plot.title = element_text(face = "bold", colour = ink),
    plot.subtitle = element_text(colour = ink2, size = 8), legend.position = "bottom", legend.title = element_text(colour = ink, size = 8),
    legend.text = element_text(colour = ink2, size = 7.5), plot.margin = margin(6, 6, 6, 6))
  om <- merge(osi, st, by = "osiedle")
  labpt <- suppressWarnings(st_point_on_surface(om))
  div <- function(lim, name) scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-lim, lim), oob = scales::squish, name = name,
                                                  guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.position = "top"))
  lab_layer <- function(v, fmt, dark) geom_sf_text(data = labpt, aes(label = sprintf(fmt, .data[[v]]), colour = abs(.data[[v]]) > dark), size = 2.1, show.legend = FALSE)
  lab_cols <- scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "#ffffff"), guide = "none")
  mA <- ggplot(om) + geom_sf(aes(fill = -A_z), colour = "#ffffff", linewidth = 0.3) + div(2.5, "z-score, + = more unemployment") +
    labs(title = "A. Registered unemployment", subtitle = "Mean of the 2019 and 2021 standardised values") + map_th
  mB <- ggplot(om) + geom_sf(aes(fill = -B_z), colour = "#ffffff", linewidth = 0.3) + div(2.5, "z-score, + = more MOPR beneficiaries") +
    labs(title = "B. Social-assistance poverty (MOPR)", subtitle = "Mean of the 2019 and 2021 standardised values") + map_th
  mD <- ggplot(om) + geom_sf(aes(fill = dep_AB), colour = "#ffffff", linewidth = 0.3) + lab_layer("dep_AB", "%.1f", 1.6) + lab_cols +
    div(2.5, "mean z of A and B, + = more deprived") + labs(title = "Economic deprivation (A + B)", subtitle = "Mean of the 2019 and 2021 GPR rounds") + map_th
  mP <- ggplot(om) + geom_sf(aes(fill = price_index), colour = "#ffffff", linewidth = 0.3) + lab_layer("price_index", "%+.0f", 25) + lab_cols +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-35, 35), oob = scales::squish,
                         name = "% vs mean osiedle", guide = guide_colourbar(barwidth = 9, barheight = 0.5, title.position = "top")) +
    labs(title = "Apartment price index", subtitle = "Hedonic price per m2, 2014-2022 sales (year, size, rooms, floor held constant)") + map_th
  f1 <- (mA + mB) / (mD + mP) + plot_annotation(theme = theme(plot.background = element_rect(fill = surf, colour = NA))) & theme(plot.background = element_rect(fill = surf, colour = NA))
  ggsave(file.path(fig_dir, "Fig_area_status_indicators.png"), f1, width = 11, height = 12, dpi = 200, bg = surf)

  type_pal <- setNames(c("#eb6834", "#4a3aa7", "#2a78d6", "#1baf7a"), type_lv)
  mT <- ggplot(om) + geom_sf(aes(fill = type), colour = "#ffffff", linewidth = 0.35) + scale_fill_manual(values = type_pal, name = NULL, drop = FALSE) +
    guides(fill = guide_legend(ncol = 2)) + labs(title = "Deprivation and housing-market status", subtitle = "Pooled economic deprivation x apartment price index, split at the medians") + map_th
  sc <- ggplot(st, aes(price_index, dep_AB)) +
    geom_hline(yintercept = median(st$dep_AB), colour = "#c3c2b7", linewidth = 0.4) + geom_vline(xintercept = median(st$price_index), colour = "#c3c2b7", linewidth = 0.4) +
    geom_point(aes(colour = type, size = population / 1000)) + scale_colour_manual(values = type_pal, guide = "none") + scale_size_area(max_size = 5, name = "Residents (k)") +
    geom_text_repel(aes(label = osiedle), size = 2.3, colour = ink2, max.overlaps = 15, segment.colour = "#c3c2b7", seed = 2) +
    labs(x = "Apartment price index (% vs mean osiedle)", y = "Economic deprivation (mean z of A and B, + = more deprived)",
         title = "Each osiedle", subtitle = sprintf("Spearman rho = %.2f", cor(st$price_index, st$dep_AB, method = "spearman"))) +
    theme_minimal(base_size = 9.5) + theme(plot.background = element_rect(fill = surf, colour = NA), panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid, linewidth = 0.25), axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = ink2), legend.position = "bottom")
  f2 <- (mT + sc + plot_layout(widths = c(1, 1.1))) & theme(plot.background = element_rect(fill = surf, colour = NA))
  ggsave(file.path(fig_dir, "Fig_area_status_typology.png"), f2, width = 14, height = 7.5, dpi = 200, bg = surf)

  oc_pal <- setNames(c("#898781", "#e34948", "#2a78d6"), outcomes)
  cp <- cmp[, .(block, label = reorder(label, rho * (outcome == outcomes[["dep_AB"]]), sum), outcome = factor(outcome, outcomes), rho, sig = q < 0.10)]
  f3 <- ggplot(cp, aes(rho, label)) + geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.4) +
    geom_line(aes(group = label), colour = grid, linewidth = 0.8) +
    geom_point(aes(colour = outcome, shape = sig), size = 2.4, stroke = 0.9) +
    scale_colour_manual(values = oc_pal, name = NULL) + scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("q >= 0.10", "FDR q < 0.10"), name = NULL) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_x_continuous(limits = c(-0.8, 0.8)) +
    labs(x = "Spearman correlation across the 42 osiedla (deprivation outcomes averaged over the 2019 and 2021 rounds)", y = NULL,
         title = "Urban tissue and land use against three measures of osiedle status") +
    theme_minimal(base_size = 9.5) + theme(plot.background = element_rect(fill = surf, colour = NA), panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid, linewidth = 0.25), axis.text = element_text(colour = ink2), strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 8), legend.position = "top", legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0), plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "Fig_area_status_tissue_associations.png"), f3, width = 9.5, height = 6.5, dpi = 200, bg = surf)
  cat("\nSaved area_status_*.csv and Fig_area_status_*.png\n")
})

## ── 20 · Deprivation score and its land-use associations ───────────────────────────────────────────────────
section("20", "Deprivation score and its land-use associations", {
  ## Stage 1 results for the manuscript (RQ1: how is deprivation embedded in
  ## land use?), against the manuscript's deprivation measure: the mean of the
  ## city's standardised unemployment (A) and poverty (B) scores, averaged over
  ## the two most recent GPR rounds, 2019 and 2021 (as in 19). Sensitivity:
  ## each round separately (2015, 2019, 2021), the three-round mean, and the
  ## poverty indicator alone (2019 and 2021).
  ##
  ## A. Core land-use descriptors (9): Spearman rho, bootstrap 95% CI,
  ##    permutation p, BH-FDR within the core; per round; poverty alone;
  ##    without the two single-institution osiedla; net of residents' distance
  ##    to the Old Market Square; quadratic rank test.
  ## B. Composition measures (residents' 400 m neighbourhoods,
  ##    population-weighted): same statistics, FDR within the set.
  ## C. Land-use types: residents' mean deprivation by type; Spearman of type
  ##    shares; compositional regression of deprivation on ilr balances of the
  ##    type shares (zero replacement 0.005), permutation p for the F statistic.
  ## D. Separability and overlap: within-osiedle variance share; leave-one-out
  ##    R2 of deprivation on land-use principal components; common support of
  ##    each core descriptor between the most and least deprived thirds of
  ##    residents.
  ##
  ## Run 15, 16, 18 and 19 first. Outputs: output/stage1/s1ms_*.csv,
  ## output/figures/Fig_s1ms_*.png (copied to manuscript/figures/).

  suppressPackageStartupMessages({library(sf); library(dplyr); library(data.table); library(ggplot2); library(patchwork)})
  set.seed(20261005)
  out_dir <- ROOT
  s1_dir <- file.path(out_dir, "output/stage1"); fig_dir <- file.path(out_dir, "output/figures")
  ms_fig <- file.path(out_dir, "manuscript/figures"); dir.create(ms_fig, showWarnings = FALSE)
  INSTITUTION_OSIEDLA <- c("Główna", "Krzesiny-Pokrzywno-Garaszewo")

  ## ── Deprivation measure (identical to 19) ───────────────────────────────
  names_lp <- fread(file.path(out_dir, "data/diagnoza2021_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, osiedle)]
  g <- merge(fread(file.path(out_dir, "data/gpr_indicators_ABC_osiedla.csv")), names_lp, by = "lp")
  bz <- rbind(data.table(lp = 1:42, round = 2015L, B_z = NA_real_),
              fread(file.path(out_dir, "data/diagnoza2019_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, round = 2019L, B_z = indicator_B_std_zscore)],
              fread(file.path(out_dir, "data/diagnoza2021_osiedla_ground_truth.csv"), encoding = "UTF-8")[, .(lp, round = 2021L, B_z = indicator_B_std_zscore)])
  g <- merge(g, bz, by = c("lp", "round"))
  g[round == 2015, B_z := -(B_poverty_per1000 - mean(B_poverty_per1000)) / sd(B_poverty_per1000)]
  g[, dep := -(A_z + B_z) / 2][, B_rk := frank(B_poverty_per1000), by = round]
  r2 <- c(2019L, 2021L)
  D <- g[, .(dep = mean(dep[round %in% r2]), dep3 = mean(dep), B_pooled = mean(B_rk[round %in% r2])), by = osiedle]
  D <- merge(D, dcast(g, osiedle ~ round, value.var = "dep"), by = "osiedle"); setnames(D, c("2015", "2019", "2021"), c("dep2015", "dep2019", "dep2021"))
  D <- merge(D, g[round == 2021, .(osiedle, B2021 = B_poverty_per1000)], by = "osiedle")

  ## ── Osiedle profiles of residents' surroundings ─────────────────────────
  cells <- st_read(file.path(s1_dir, "cell_tissue_descriptors.gpkg"), quiet = TRUE)
  rynek <- st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180)
  cells$dist_km <- as.numeric(st_distance(st_centroid(st_geometry(cells)), rynek)) / 1000
  cd <- as.data.table(st_drop_geometry(cells))
  assoc16 <- fread(file.path(s1_dir, "stage1_associations.csv"), encoding = "UTF-8")
  core <- assoc16[core == TRUE, .(var, label, family, agg)]
  comp <- data.table(var = c("share_res_multi", "share_res_single", "share_retail_services", "share_office", "share_education_culture", "share_health",
                             "share_industry_storage", "share_university", "share_dormitory", "landuse_mix", "mixed_ground_share", "amenity_dens", "amenity_div",
                             "z_single", "z_multi", "z_mixed", "z_services", "z_production", "z_green", "z_transport_infra"),
                     label = c("Floor: multi-family housing", "Floor: single-family housing", "Floor: retail and services", "Floor: office", "Floor: education and culture",
                               "Floor: health", "Floor: industry and storage", "Floor: university", "Floor: student dormitories", "Land-use mix",
                               "Multi-family with commercial use", "Amenity density", "Amenity diversity",
                               "Zoned single-family", "Zoned multi-family", "Zoned mixed residential-service", "Zoned services and sport", "Zoned services-production",
                               "Zoned green space", "Zoned transport and infrastructure"))
  wmed <- function(x, w) { ok <- !is.na(x) & w > 0; x <- x[ok]; w <- w[ok]; o <- order(x); x[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]] }
  wmn <- function(x, w) { ok <- !is.na(x) & w > 0; sum(x[ok] * w[ok]) / sum(w[ok]) }
  P <- cd[, c(lapply(setNames(core$var, core$var), function(v) if (core[var == v, agg] == "median") wmed(get(v), population) else wmn(get(v), population)),
              lapply(setNames(comp$var, comp$var), function(v) wmn(get(v), population)), .(dist_km = wmn(dist_km, population))), by = osiedle]
  lu <- fread(file.path(s1_dir, "landuse_osiedle.csv"), encoding = "UTF-8")[, .(osiedle, res_nb_industrial, res_nb_rail)]
  P <- merge(merge(P, lu, by = "osiedle"), D, by = "osiedle"); stopifnot(nrow(P) == 42)
  comp <- rbind(comp, data.table(var = c("res_nb_industrial", "res_nb_rail"), label = c("Industrial-storage land within 400 m", "Rail land within 400 m")))

  ## ── A/B. Associations ─────────────────────────────────────────────────────
  sp <- function(x, y) cor(x, y, method = "spearman", use = "complete.obs")
  pc <- function(x, y, z) { ok <- !is.na(x); cor(resid(lm(rank(y[ok]) ~ rank(z[ok]))), resid(lm(rank(x[ok]) ~ rank(z[ok])))) }
  assoc <- function(set, setname) rbindlist(lapply(seq_len(nrow(set)), function(i) {
    v <- set$var[i]; x <- P[[v]]; y <- P$dep; ok <- !is.na(x)
    rho <- sp(x, y); perm <- replicate(9999, sp(x[ok], sample(y[ok])))
    boot <- replicate(2000, { j <- sample(which(ok), replace = TRUE); sp(x[j], y[j]) })
    keep <- ok & !P$osiedle %in% INSTITUTION_OSIEDLA
    rx <- rank(x[ok]); ry <- rank(y[ok]); m1 <- lm(ry ~ rx); m2 <- lm(ry ~ rx + I(rx^2))
    data.table(set = setname, var = v, label = set$label[i], rho = rho, lo = quantile(boot, .025), hi = quantile(boot, .975),
               p = (1 + sum(abs(perm) >= abs(rho))) / 10000, r2015 = sp(x, P$dep2015), r2019 = sp(x, P$dep2019), r2021 = sp(x, P$dep2021),
               r_three = sp(x, P$dep3), r_poverty = sp(x, P$B_pooled), r_noinst = sp(x[keep], y[keep]), r_net_dist = pc(x, y, P$dist_km), p_quad = anova(m1, m2)$`Pr(>F)`[2])
  }))
  A <- rbind(assoc(core, "core"), assoc(comp, "composition"))
  A[, q := p.adjust(p, "BH"), by = set]
  fwrite(A, file.path(s1_dir, "s1ms_associations.csv"))
  fmt <- function(d) d[, .(label, rho = round(rho, 2), ci = sprintf("[%.2f, %.2f]", lo, hi), q = signif(q, 2), r2015 = round(r2015, 2), r2019 = round(r2019, 2),
                           r2021 = round(r2021, 2), three = round(r_three, 2), pov = round(r_poverty, 2), noinst = round(r_noinst, 2), net_dist = round(r_net_dist, 2), p_quad = signif(p_quad, 2))]
  cat("Distance to centre vs deprivation: pooled", round(sp(P$dist_km, P$dep), 2), "; by round",
      round(sp(P$dist_km, P$dep2015), 2), round(sp(P$dist_km, P$dep2019), 2), round(sp(P$dist_km, P$dep2021), 2), "\n")
  cat("\n=== A. Core descriptors vs pooled deprivation (FDR within core) ===\n"); print(fmt(A[set == "core"][order(-abs(rho))]), row.names = FALSE)
  cat("\n=== B. Composition (residents' surroundings) vs pooled deprivation (FDR within set) ===\n"); print(fmt(A[set == "composition"][order(-abs(rho))]), row.names = FALSE)
  ## shape of the two curvilinear associations: mean deprivation by tercile of the descriptor (14 osiedla each)
  for (v in c("fsi", "wc_tree")) { tt <- P[, .(dep = mean(dep), lo = min(get(v)), hi = max(get(v)), osiedla = paste(head(osiedle[order(-dep)], 3), collapse = ", ")),
                                     by = .(tercile = cut(rank(get(v), ties.method = "first"), c(0, 14, 28, 42), labels = c("low", "middle", "high")))][order(tercile)]
    cat("\nMean deprivation by tercile of", v, "\n"); print(tt, digits = 2) }

  ## ── C. Land-use types ─────────────────────────────────────────────────────
  ty <- as.data.table(st_drop_geometry(st_read(file.path(s1_dir, "cell_tissue_types.gpkg"), quiet = TRUE)))[, .(cell_id, osiedle, population, type)]
  tn <- fread(file.path(s1_dir, "typology_profiles.csv"), encoding = "UTF-8")[, .(type, name, residents)]
  ty <- merge(ty, D[, .(osiedle, dep)], by = "osiedle")
  bytype <- ty[, .(residents = sum(population), mean_dep = weighted.mean(dep, population), sd_dep = sqrt(sum(population * (dep - weighted.mean(dep, population))^2) / sum(population)),
                   osiedla_present = uniqueN(osiedle[population > 0])), by = type][order(type)]
  bytype <- merge(bytype, tn[, .(type, name)], by = "type")
  cat("\n=== C. Residents' mean deprivation by land-use type ===\n"); print(bytype, digits = 3)
  sh <- dcast(ty[, .(pop = sum(population)), by = .(osiedle, type)], osiedle ~ type, value.var = "pop", fill = 0)
  S <- as.matrix(sh[, -1]); S <- S / rowSums(S); S[S == 0] <- 0.005; S <- S / rowSums(S)
  ilr <- data.table(osiedle = sh$osiedle,
                    b1 = sqrt(4 / 4) * log(sqrt(S[, 1] * S[, 2]) / sqrt(S[, 3] * S[, 4])),   # multi-family types (compact, estates) vs low-rise (single-family, loose)
                    b2 = sqrt(1 / 2) * log(S[, 1] / S[, 2]),                               # compact blocks vs estates
                    b3 = sqrt(1 / 2) * log(S[, 3] / S[, 4]))                               # single-family vs loose peripheral
  R <- merge(ilr, D[, .(osiedle, dep)], by = "osiedle")
  fit <- lm(dep ~ b1 + b2 + b3, data = R); F0 <- summary(fit)$fstatistic[1]
  Fperm <- replicate(9999, summary(lm(sample(R$dep) ~ b1 + b2 + b3, data = R))$fstatistic[1])
  type_share_rho <- sapply(1:4, function(k) sp(S[match(R$osiedle, sh$osiedle), k], R$dep))
  cat("Compositional regression of deprivation on type balances: R2 =", round(summary(fit)$r.squared, 3), "adj R2 =", round(summary(fit)$adj.r.squared, 3),
      "permutation p =", signif((1 + sum(Fperm >= F0)) / 10000, 2), "\n"); print(round(coef(summary(fit)), 3))
  cat("Spearman of type share with deprivation (types 1-4):", round(type_share_rho, 2), "\n")
  fwrite(cbind(bytype, share_rho = type_share_rho[bytype$type]), file.path(s1_dir, "s1ms_types.csv"))

  ## ── D. Separability and overlap ──────────────────────────────────────────
  cd <- merge(cd, D[, .(osiedle, dep)], by = "osiedle")
  within <- sapply(core$var, function(v) { y <- cd[[v]]; ok <- !is.na(y); 1 - summary(lm(y[ok] ~ factor(cd$osiedle[ok]), weights = cd$population[ok]))$r.squared })
  Pm <- as.data.frame(P[, core$var, with = FALSE]); for (j in names(Pm)) Pm[[j]][is.na(Pm[[j]])] <- median(Pm[[j]], na.rm = TRUE)
  yy <- P$dep
  loo <- sapply(1:5, function(k) { pred <- sapply(1:42, function(i) { pca <- prcomp(Pm[-i, ], scale. = TRUE)
    m <- lm(yy[-i] ~ ., data = data.frame(pca$x[, 1:k, drop = FALSE])); predict(m, newdata = data.frame(predict(pca, Pm[i, , drop = FALSE])[, 1:k, drop = FALSE])) })
    1 - sum((yy - pred)^2) / sum((yy - mean(yy))^2) })
  wq <- function(x, w, p) { o <- order(x); x <- x[o]; w <- cumsum(w[o]) / sum(w); x[which(w >= p)[1]] }
  o <- order(cd$dep); cw <- cumsum(cd$population[o]) / sum(cd$population)
  cut_lo <- cd$dep[o][which(cw >= 1/3)[1]]; cut_hi <- cd$dep[o][which(cw >= 2/3)[1]]
  lo_t <- cd[dep <= cut_lo]; hi_t <- cd[dep > cut_hi]
  support <- sapply(core$var, function(v) { a <- lo_t[!is.na(get(v))]; b <- hi_t[!is.na(get(v))]
    rng <- c(wq(a[[v]], a$population, 0.05), wq(a[[v]], a$population, 0.95)); sum(b$population[between(b[[v]], rng[1], rng[2])]) / sum(b$population) })
  sep <- data.table(var = core$var, label = core$label, within_share = within, common_support = support)
  cat("\n=== D. Separability ===\n"); print(sep[, .(label, within = round(within_share, 2), support = round(common_support, 2))], row.names = FALSE)
  cat("LOO R2 of deprivation on first k land-use PCs:", paste0("k=", 1:5, ": ", round(loo, 2), collapse = ", "), "\n")
  fwrite(sep, file.path(s1_dir, "s1ms_separability.csv")); fwrite(data.table(k = 1:5, loo_r2 = loo), file.path(s1_dir, "s1ms_loo_r2.csv"))
  fwrite(P, file.path(s1_dir, "s1ms_osiedle_profiles.csv"))

  ## ── Numbers quoted in the manuscript text ────────────────────────────────
  sp2 <- function(a, b) round(cor(a, b, method = "spearman"), 2)
  cat("\n=== Text numbers ===\n")
  cat("Agreement of measures: 2019+2021 vs three-round", sp2(P$dep, P$dep3), "| vs 2021 alone", sp2(P$dep, P$dep2021), "| 2019 vs 2021", sp2(P$dep2019, P$dep2021), "\n")
  cat("Distance to centre: 2019+2021", sp2(P$dist_km, P$dep), "| three-round", sp2(P$dist_km, P$dep3), "\n")
  cat("Range of deprivation:", round(range(P$dep), 2), "\n")
  print(P[order(-dep)][c(1:6, 37:42), .(osiedle, dep = round(dep, 2))])
  print(P[osiedle %like% "Nowe Winogrady", .(osiedle, dep = round(dep, 2), fsi = round(fsi, 2))])
  tq <- function(v) P[, .(dep = round(mean(dep), 2), lo = round(min(get(v)), 2), hi = round(max(get(v)), 2)),
                       by = .(t = cut(rank(get(v), ties.method = "first"), c(0, 14, 28, 42), labels = c("low", "mid", "high")))][order(t)]
  cat("FSI thirds:\n"); print(tq("fsi")); cat("Tree-cover thirds:\n"); print(tq("wc_tree"))
  P[, fsi_t := cut(rank(fsi, ties.method = "first"), c(0, 14, 28, 42), labels = c("low", "mid", "high"))]
  cat("Least deprived 8 and their FSI third:\n"); print(P[order(dep)][1:8, .(osiedle, dep = round(dep, 2), fsi_t)])
  cat("Most vs least deprived third of residents: FSI", round(c(weighted.mean(hi_t$fsi, hi_t$population), weighted.mean(lo_t$fsi, lo_t$population)), 2),
      "| tree", round(c(weighted.mean(hi_t$wc_tree, hi_t$population), weighted.mean(lo_t$wc_tree, lo_t$population)), 2), "\n")
  lab_tex <- c(z_mixed = "Zoned mixed residential--service", share_office = "Office floor", amenity_div = "Amenity diversity", landuse_mix = "Land-use mix",
               ilr_res_vs_nonres = "Residential vs non-residential floor$^{a}$", z_single = "Zoned single-family", res_nb_industrial = "Industrial--storage land",
               res_nb_rail = "Rail land", fsi = "Floor-area density$^{a}$", med_gap = "Gap to nearest building$^{a}$", med_setback = "Setback from street$^{a}$",
               med_block_area_ha = "Block area$^{a}$", intersection_dens = "Intersection density$^{a}$", deadend_share = "Dead-end share$^{a}$", wc_tree = "Tree cover$^{a}$")
  ft <- function(x) { s <- sprintf("%.2f", x); ifelse(x < 0, paste0("$-$", substring(s, 2)), s) }
  fq <- function(q) ifelse(q < 0.01, "<0.01", sprintf("%.2f", q))
  tex <- vapply(names(lab_tex), function(v) { r <- A[var == v]
    sprintf("%s & %s (%s, %s) & %s & %s & %s & %s & %s & %s \\\\", lab_tex[[v]], ft(r$rho), ft(r$lo), ft(r$hi), fq(r$q), ft(r$r2015), ft(r$r2019), ft(r$r2021), ft(r$r_poverty), ft(r$r_net_dist)) }, "")
  writeLines(tex, file.path(s1_dir, "s1ms_table2_rows.tex"), useBytes = TRUE)
  cat("Table 2 rows written to output/stage1/s1ms_table2_rows.tex\n")
  cat("Poverty-only (2019+2021) correlations:\n"); print(A[var %in% names(lab_tex), .(label, pov = round(r_poverty, 2), noinst = round(r_noinst, 2))])

  ## ── Figures ──────────────────────────────────────────────────────────────
  ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid <- "#e1e0d9"; surf <- "#ffffff"; blue <- "#2a78d6"
  fam_map <- c("Built form" = "Intensity", "Spacing and street interface" = "Configuration", "Street network" = "Configuration",
               "Function and amenities" = "Composition", "Green, open and sealed" = "Green and open space")
  pd <- copy(A)[, family := fifelse(set == "composition", "Composition", fam_map[core[match(A$var, core$var), family]])]
  pd[set == "core" & var %in% c("ilr_res_vs_nonres", "ilr_multi_vs_single"), family := "Composition"]
  pd[, `:=`(family = factor(family, c("Composition", "Intensity", "Configuration", "Green and open space")),
            sig = factor(fifelse(q < 0.10, "FDR q < 0.10", "Not significant"), c("FDR q < 0.10", "Not significant")),
            label = fifelse(set == "core", paste0(label, " *"), label))]
  pd[, label := reorder(label, rho)]
  f1 <- ggplot(pd, aes(rho, label)) + geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, yend = label), colour = "#c3c2b7", linewidth = 0.6) +
    geom_point(aes(colour = sig), size = 2.2) + scale_colour_manual(values = c("FDR q < 0.10" = blue, "Not significant" = muted), name = NULL) +
    facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_x_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.75, 0.75, 0.25)) +
    labs(x = "Spearman correlation with deprivation (42 osiedla), bootstrap 95% interval", y = NULL) +
    theme_minimal(base_size = 9.5) + theme(plot.background = element_rect(fill = surf, colour = NA), panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid, linewidth = 0.25), axis.text = element_text(colour = ink2), strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5, colour = ink), legend.position = "top")
  ggsave(file.path(fig_dir, "Fig_s1ms_associations.png"), f1, width = 8, height = 8.2, dpi = 300, bg = surf)

  osi <- st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE) %>% st_transform(2180) %>% select(osiedle = name) %>% left_join(D, by = "osiedle")
  grid_all <- st_read(file.path(out_dir, "output/income_map_grid_200m.gpkg"), quiet = TRUE)
  type_names <- bytype[order(type), name]; type_pal <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"), type_names)
  tmap <- grid_all %>% inner_join(ty[, .(cell_id, type)], by = "cell_id") %>% mutate(type = factor(type_names[type], type_names))
  lab <- suppressWarnings(st_point_on_surface(osi))
  mth <- theme_void(base_size = 9.5) + theme(plot.background = element_rect(fill = surf, colour = NA), legend.position = "bottom",
    legend.title = element_text(colour = ink, size = 8.5), legend.text = element_text(colour = ink2, size = 8), plot.tag = element_text(face = "bold"))
  mA <- ggplot(osi) + geom_sf(aes(fill = dep), colour = "#ffffff", linewidth = 0.3) +
    geom_sf_text(data = lab, aes(label = sprintf("%.1f", dep), colour = abs(dep) > 1.6), size = 2.1, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "#ffffff"), guide = "none") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-2.6, 2.6), oob = scales::squish,
                         name = "Deprivation (mean z, + = more deprived)", guide = guide_colourbar(barwidth = 10, barheight = 0.5, title.position = "top")) + mth
  mB <- ggplot() + geom_sf(data = st_union(osi), fill = "#f3f2ee", colour = "#c3c2b7", linewidth = 0.3) + geom_sf(data = tmap, aes(fill = type), colour = NA) +
    geom_sf(data = osi, fill = NA, colour = "#ffffff", linewidth = 0.3) + scale_fill_manual(values = type_pal, name = NULL) + guides(fill = guide_legend(ncol = 2)) + mth
  f2 <- (mA + mB) + plot_annotation(tag_levels = "a") & theme(plot.background = element_rect(fill = surf, colour = NA))
  ggsave(file.path(fig_dir, "Fig_s1ms_maps.png"), f2, width = 11, height = 6.6, dpi = 300, bg = surf)
  file.copy(file.path(fig_dir, c("Fig_s1ms_associations.png", "Fig_s1ms_maps.png")), ms_fig, overwrite = TRUE)
  cat("\nSaved s1ms_*.csv and Fig_s1ms_*.png (also in manuscript/figures/)\n")
})

