## Shared cartographic style for the maps of the paper (sourced by the figure scripts).
## One basemap (city land, the Warta and lakes, white osiedle boundaries), one palette, reference labels for
## a few districts (white halo), a 2-km scale bar, and one theme with lowercase panel tags. North is up on every map.
## Colours: least deprived #2a78d6 (blue), most deprived #e34948 (red), neutral #d9d7d0; the diverging
## deprivation scale runs between them through a light neutral.

suppressPackageStartupMessages({library(sf); library(data.table); library(ggplot2)})
MSTY <- list(least = "#2a78d6", most = "#e34948", middle = "#d9d7d0", land = "#f3f2ee", edge = "#ffffff", water = "#c4daee",
             water_line = "#a9c7e3", ink = "#2b2b28", muted = "#6b6a64", line = "#8f8d86")
ms_layers <- function(out_dir = ROOT) {
  osi <- st_transform(st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE), 2180)
  city <- st_union(st_make_valid(osi))
  water <- st_read(file.path(out_dir, "output/stage1/tissue_inputs.gpkg"), "green_blue", quiet = TRUE)
  water <- st_transform(water[water$cover == "water", ], 2180)
  water <- suppressWarnings(st_intersection(st_make_valid(water), city))
  water <- water[as.numeric(st_area(water)) > 2000, ]          # the river, its branches and lakes; no ponds
  lab <- data.table(osiedle = c("Stare Miasto", "Jeżyce", "Św. Łazarz", "Wilda", "Górczyn", "Piątkowo", "Nowe Winogrady Północ", "Rataje",
                                "Chartowo", "Ostrów Tumski-Śródka-Zawady-Komandoria", "Główna", "Krzesiny-Pokrzywno-Garaszewo", "Ogrody", "Naramowice",
                                "Głuszyna", "Fabianowo-Kotowo"),
                    label = c("Stare Miasto", "Jeżyce", "Łazarz", "Wilda", "Górczyn", "Piątkowo", "Winogrady", "Rataje",
                              "Chartowo", "Śródka", "Główna", "Krzesiny", "Ogrody", "Naramowice", "Głuszyna", "Fabianowo-Kotowo"))
  pts <- st_point_on_surface(osi[osi$name %in% lab$osiedle, "name"])
  lab <- merge(lab, data.table(osiedle = pts$name, X = st_coordinates(pts)[, 1], Y = st_coordinates(pts)[, 2]), by = "osiedle")
  ctr <- st_coordinates(st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180))
  list(osi = osi, city = city, water = water, labels = lab, centre = data.table(X = ctr[1], Y = ctr[2]), bbox = st_bbox(osi))
}
## basemap: land, water, osiedle boundaries (drawn under the data); city outline on top via ms_outline()
ms_base <- function(L, borders = TRUE) list(
  geom_sf(data = L$city, fill = MSTY$land, colour = NA, inherit.aes = FALSE),
  geom_sf(data = L$water, fill = MSTY$water, colour = NA, inherit.aes = FALSE),
  if (borders) geom_sf(data = L$osi, fill = NA, colour = MSTY$edge, linewidth = 0.25, inherit.aes = FALSE))
ms_water_top <- function(L) geom_sf(data = L$water, fill = MSTY$water, colour = NA, alpha = 0.9, inherit.aes = FALSE)
ms_outline <- function(L) geom_sf(data = L$city, fill = NA, colour = "#b9b7b0", linewidth = 0.35, inherit.aes = FALSE)
## district names with a white halo; `which` selects labels by display name
ms_labels <- function(L, which = NULL, size = 2.5) {
  d <- if (is.null(which)) L$labels else L$labels[label %in% which]
  ggrepel::geom_text_repel(data = d, aes(X, Y, label = label), size = size, colour = MSTY$ink, fontface = "italic", bg.colour = "#ffffff", bg.r = 0.15,
                           min.segment.length = Inf, force = 0.3, seed = 7, box.padding = 0.05, inherit.aes = FALSE)
}
ms_centre <- function(L) geom_point(data = L$centre, aes(X, Y), shape = 3, size = 2.2, stroke = 0.7, colour = MSTY$ink, inherit.aes = FALSE)
## 2-km scale bar (0-1-2 km, alternating) in the bottom-left corner; `text` scales the labels
ms_scalebar <- function(L, km = 2, pos = c(0.05, 0.06), h = 0.010, text = 1) {
  b <- L$bbox; x0 <- b[["xmin"]] + pos[1] * (b[["xmax"]] - b[["xmin"]]); y0 <- b[["ymin"]] + pos[2] * (b[["ymax"]] - b[["ymin"]]); hh <- h * (b[["ymax"]] - b[["ymin"]])
  seg <- data.table(x1 = x0 + c(0, 500 * km), x2 = x0 + c(500 * km, 1000 * km), f = c(MSTY$ink, "#ffffff"))
  list(annotate("rect", xmin = seg$x1, xmax = seg$x2, ymin = y0, ymax = y0 + hh, fill = seg$f, colour = MSTY$ink, linewidth = 0.25),
       annotate("text", x = x0 + c(0, 500 * km, 1000 * km), y = y0 + 2.6 * hh, label = c("0", km / 2, paste(km, "km")), size = 2.3 * text, colour = MSTY$ink, vjust = 0))
}
ms_coord <- function(L, pad = 300) coord_sf(datum = NA, expand = FALSE, xlim = c(L$bbox[["xmin"]] - pad, L$bbox[["xmax"]] + pad), ylim = c(L$bbox[["ymin"]] - pad, L$bbox[["ymax"]] + pad))
## theme: void map, plain panel titles, legends at the bottom, bold lowercase tags (patchwork tag_levels = "a")
theme_ms <- function(base = 9.5) theme_void(base_size = base) + theme(
  plot.title = element_text(size = base, colour = MSTY$ink, hjust = 0.5, margin = margin(b = 3)),
  plot.tag = element_text(face = "bold", size = base + 2), legend.position = "bottom", legend.title = element_text(size = base - 1, colour = MSTY$ink),
  legend.text = element_text(size = base - 1.5, colour = MSTY$muted), legend.key.height = unit(0.22, "cm"), legend.key.width = unit(1.1, "cm"),
  legend.title.position = "top", plot.background = element_rect(fill = "#ffffff", colour = NA), plot.margin = margin(4, 4, 4, 4))
DEPSCALE <- function(lim, name = "Catchment deprivation (z)", ...) scale_fill_gradient2(low = MSTY$least, mid = "#f0efec", high = MSTY$most, midpoint = 0,
                                                                          limits = c(-lim, lim), oob = scales::squish, name = name, ...)
## rings of distance (2, 4, 6 and 8 km) and the eight sectors around the Old Market Square; `text` scales the labels
ms_guides <- function(text = 1) {
  ctr <- st_transform(st_sfc(st_point(c(16.9335, 52.4083)), crs = 4326), 2180); c0 <- st_coordinates(ctr)[1, 1:2]
  rings <- do.call(c, lapply(c(2, 4, 6, 8), function(r) st_cast(st_buffer(ctr, r * 1000, nQuadSegs = 60), "LINESTRING")))
  rl <- st_sf(r = paste0(c(2, 4, 6, 8), " km"), geometry = st_sfc(lapply(c(2, 4, 6, 8), function(r) st_point(c0 + c(r * 1000 * cos(pi / 4) + 150, r * 1000 * sin(pi / 4)))), crs = 2180))
  rays <- st_sfc(lapply(seq(22.5, 337.5, 45) * pi / 180, function(a) st_linestring(rbind(c0 + 1000 * c(cos(a), sin(a)), c0 + 9500 * c(cos(a), sin(a))))), crs = 2180)
  sect <- st_sf(lab = c("E", "NE", "N", "NW", "W", "SW", "S", "SE"),
                geometry = st_sfc(lapply(seq(0, 315, 45) * pi / 180, function(a) st_point(c0 + 10200 * c(cos(a), sin(a)))), crs = 2180))
  list(geom_sf(data = rings, colour = MSTY$ink, linewidth = 0.25, linetype = "22", alpha = 0.55, inherit.aes = FALSE),
       geom_sf(data = rays, colour = MSTY$ink, linewidth = 0.2, alpha = 0.3, inherit.aes = FALSE),
       geom_sf_text(data = rl, aes(label = r), size = 2.2 * text, colour = MSTY$muted, hjust = 0, inherit.aes = FALSE),
       geom_sf_text(data = sect, aes(label = lab), size = 2.6 * text, colour = MSTY$ink, fontface = "bold", inherit.aes = FALSE),
       geom_sf(data = ctr, shape = 3, size = 2, stroke = 0.6, colour = MSTY$ink, inherit.aes = FALSE))
}
