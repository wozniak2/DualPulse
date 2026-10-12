## DualPulse | stage 2: Land-use inputs and descriptors
##
## BDOT10k buildings, roads and land cover, OpenStreetMap amenities and ESA WorldCover, and the land-use
## descriptors of the 400 m neighbourhood of every populated cell (supplementary S1).
##
## Sections (named by the numbered scripts they replace, as in the project notes):
##   14      Land-use inputs: BDOT10k, OpenStreetMap, WorldCover
##   15a     Land-use descriptors of populated cells
##
## Run from the repository root: Rscript run.R 2   (one section: Rscript run.R 2:15a)

if (!exists("section")) source("R/setup.R")

## ── 14 · Land-use inputs: BDOT10k, OpenStreetMap, WorldCover ───────────────────────────────────────────────
section("14", "Land-use inputs: BDOT10k, OpenStreetMap, WorldCover", {
  ## Stage 1 inputs: assemble the layers that describe urban tissue.
  ##
  ## Sources (each used for what it observes directly):
  ##   - BDOT10k for Poznań city (powiat 3064), GUGiK open data
  ##     https://opendata.geoportal.gov.pl/bdot10k/schemat2021/30/3064_GML.zip
  ##     (downloaded to data/raw/bdot10k/ on 2026-10-05, 36.5 MB):
  ##       OT_BUBD_A  buildings: footprint, general + detailed function, storeys
  ##       OT_SKDR_L  road axes with road class (internal/local/collector/main/...)
  ##       OT_PTZB_A  built-up land with type and character (compact/dense/loose)
  ##       OT_PTLZ_A, OT_PTTR_A, OT_PTWP_A  woodland, grassland, water
  ##   - OpenStreetMap amenities and shops (Overpass, cached)
  ##   - ESA WorldCover 10 m 2021 (satellite-derived land cover), read as a
  ##     windowed COG stream from the public AWS bucket -- no full-tile download.
  ##
  ## Outputs: output/stage1/tissue_inputs.gpkg (layers) + worldcover_poznan.tif

  suppressPackageStartupMessages({library(sf); library(dplyr); library(osmdata); library(terra)})

  out_dir <- ROOT
  bdot    <- file.path(out_dir, "data/raw/bdot10k/PL.PZGiK.308.3064/BDOT10k")
  s1_dir  <- file.path(out_dir, "output/stage1"); dir.create(s1_dir, showWarnings = FALSE, recursive = TRUE)
  gpkg    <- file.path(s1_dir, "tissue_inputs.gpkg")
  rd <- function(layer) st_read(file.path(bdot, paste0("PL.PZGiK.308.3064__", layer, ".xml")), quiet = TRUE)
  flat <- function(x) vapply(x, function(v) paste(v, collapse = "|"), character(1))

  ## ── Buildings ────────────────────────────────────────────────────────────
  b <- rd("OT_BUBD_A") %>% filter(kategoriaIstnienia %in% c("eksploatowany", "nieczynny"))
  b$det <- flat(b$funkcjaSzczegolowaBudynku)
  ogol  <- b$funkcjaOgolnaBudynku
  b$fn <- case_when(
    grepl("wielorodzinny", b$det)                                    ~ "res_multi",
    grepl("jednorodzinny", b$det)                                    ~ "res_single",
    ogol == "budynki mieszkalne" & grepl("letniskowy", b$det)        ~ "auxiliary",
    ogol == "budynki mieszkalne"                                     ~ "res_multi",  # residential, detail missing
    ogol == "budynki biurowe"                                        ~ "office",
    ogol == "budynki handlowo-usługowe"                              ~ "retail_services",
    ogol %in% c("budynki przemysłowe", "zbiorniki, silosy i budynki magazynowe") ~ "industry_storage",
    ogol == "budynki oświaty, nauki i kultury oraz budynki sportowe" ~ "education_culture",
    ogol == "budynki szpitali i inne budynki opieki zdrowotnej"      ~ "health",
    grepl("garaż|gospodarczy|stacja transformatorowa|domek kempingowy|letniskowy", b$det) ~ "auxiliary",
    TRUE                                                             ~ "other_nonres"
  )
  b$storeys <- as.numeric(b$liczbaKondygnacji)
  med_st <- tapply(b$storeys, b$fn, median, na.rm = TRUE)
  b$storeys[is.na(b$storeys)] <- med_st[b$fn[is.na(b$storeys)]]
  b$footprint <- as.numeric(st_area(b))
  b$floor     <- b$footprint * b$storeys
  b$mixed_ground <- grepl("wielorodzinny", b$det) & grepl("handlowo|apteka|restauracja|usług", b$det)
  ## universities and student dormitories, identified from the detailed function; they stay in their
  ## general class (education and culture; housing) and are measured as separate shares in 15
  b$university <- grepl("szkoła wyższa", b$det)
  b$dormitory  <- grepl("dom studencki", b$det)
  cat("University buildings:", sum(b$university), "| student dormitories:", sum(b$dormitory), "\n")
  print(table(function_class = b$fn, university = b$university)[, "TRUE", drop = FALSE]); print(table(function_class = b$fn, dormitory = b$dormitory)[, "TRUE", drop = FALSE])
  b <- b %>% select(lokalnyId, fn, storeys, footprint, floor, mixed_ground, university, dormitory)
  cat("Buildings kept:", nrow(b), "\n"); print(table(b$fn))
  st_write(b, gpkg, layer = "buildings", delete_layer = TRUE, quiet = TRUE)

  ## ── Roads (axes) ─────────────────────────────────────────────────────────
  r <- rd("OT_SKDR_L") %>% filter(kategoriaIstnienia == "eksploatowany")
  r <- r %>% transmute(klasa = klasaDrogi, zarzadzanie = kategoriaZarzadzania, polozenie = polozenie,
                       jezdnie = as.numeric(liczbaJezdniDrogi))
  r <- st_cast(r, "LINESTRING", warn = FALSE)
  cat("\nRoad axes:", nrow(r), "\n"); print(table(r$klasa))
  st_write(r, gpkg, layer = "roads", delete_layer = TRUE, quiet = TRUE)

  ## ── Land cover (BDOT) ────────────────────────────────────────────────────
  zb <- rd("OT_PTZB_A") %>% transmute(rodzaj, charakter)
  st_write(zb, gpkg, layer = "builtup", delete_layer = TRUE, quiet = TRUE)
  green <- bind_rows(
    rd("OT_PTLZ_A") %>% transmute(cover = "woodland"),
    rd("OT_PTTR_A") %>% transmute(cover = ifelse(rodzaj == "uprawa na gruntach ornych", "arable", "grassland")),
    rd("OT_PTWP_A") %>% transmute(cover = "water")
  )
  st_write(green, gpkg, layer = "green_blue", delete_layer = TRUE, quiet = TRUE)
  cat("\nBuilt-up polygons:", nrow(zb), " green/blue polygons:", nrow(green), "\n")

  ## ── OSM amenities, from a local extract ──────────────────────────────────
  ## The public Overpass API returned repeated 504 timeouts on 2026-10-05, so
  ## amenities are read with GDAL's OSM driver from the August 2026 extract
  ## already held by the sibling PoznanABM project (clipped from Geofabrik's
  ## wielkopolskie-latest). Points carry their tags in other_tags; areas carry
  ## amenity/shop/leisure/tourism/office as columns.
  poz <- st_read(file.path(out_dir, "output/osiedla_42_boundaries.gpkg"), quiet = TRUE) %>% st_union()
  bb  <- st_bbox(st_transform(poz, 4326))
  osm_file <- OSM_EXTRACT
  keys <- c("shop", "amenity", "leisure", "tourism", "office")   # priority when a feature has several
  cache <- file.path(s1_dir, "_osm_amenities_raw.rds")
  if (!file.exists(cache)) {
    where_pts <- paste(sprintf("other_tags LIKE '%%\"%s\"=>%%'", keys), collapse = " OR ")
    p <- st_read(osm_file, layer = "points", query = paste("SELECT osm_id, other_tags FROM points WHERE", where_pts), quiet = TRUE)
    tagval <- function(tags, k) ifelse(grepl(sprintf('"%s"=>"', k), tags),
                                       sub(sprintf('.*"%s"=>"([^"]*)".*', k), "\\1", tags), NA_character_)
    for (k in keys) p[[k]] <- tagval(p$other_tags, k)
    where_pol <- paste(sprintf("%s IS NOT NULL", keys), collapse = " OR ")
    a <- st_read(osm_file, layer = "multipolygons",
                 query = paste("SELECT osm_id, osm_way_id, shop, amenity, leisure, tourism, office FROM multipolygons WHERE", where_pol), quiet = TRUE)
    a$osm_id <- ifelse(is.na(a$osm_id), paste0("w", a$osm_way_id), a$osm_id)
    a <- suppressWarnings(st_point_on_surface(st_make_valid(a)))
    raw <- bind_rows(p %>% select(osm_id, all_of(keys)), a %>% select(osm_id, all_of(keys)))
    saveRDS(raw, cache)
  }
  raw <- readRDS(cache)
  first_key <- apply(st_drop_geometry(raw)[keys], 1, function(r) keys[which(!is.na(r))[1]])
  raw$key <- first_key
  raw$value <- vapply(seq_len(nrow(raw)), function(i) raw[[raw$key[i]]][i], character(1))
  am <- raw %>% select(osm_id, key, value) %>% distinct(osm_id, .keep_all = TRUE) %>% st_transform(2180)
  cat("\nOSM features with amenity/shop/leisure/tourism/office tags in extract:", nrow(am), "\n")
  ## categories; street furniture and parking are not destinations
  cat_of <- function(key, v) case_when(
    key == "shop"                                                                        ~ "retail",
    key == "amenity" & v %in% c("restaurant","cafe","bar","pub","fast_food","ice_cream","food_court","biergarten") ~ "food_drink",
    key == "amenity" & v %in% c("school","kindergarten","university","college","library","language_school","music_school","childcare") ~ "education",
    key == "amenity" & v %in% c("hospital","clinic","doctors","dentist","pharmacy","veterinary","nursing_home","social_facility") ~ "health_care",
    key == "amenity" & v %in% c("bank","post_office","townhall","courthouse","police","marketplace","community_centre","place_of_worship","social_centre") ~ "civic_services",
    key == "amenity" & v %in% c("theatre","cinema","arts_centre","nightclub") ~ "culture",
    key %in% c("leisure") & v %in% c("sports_centre","fitness_centre","swimming_pool","stadium","sports_hall","playground","pitch") ~ "sport_play",
    key == "tourism" & v %in% c("hotel","hostel","museum","gallery","guest_house","attraction") ~ "culture",
    key == "office"                                                                      ~ "office_poi",
    TRUE ~ NA_character_)
  am$cat <- cat_of(am$key, am$value)
  am <- am %>% filter(!is.na(cat)) %>% select(osm_id, key, value, cat)
  am <- am[lengths(st_intersects(am, st_transform(poz, 2180))) > 0, ]
  cat("\nAmenities inside city:", nrow(am), "\n"); print(table(am$cat))
  st_write(am, gpkg, layer = "amenities", delete_layer = TRUE, quiet = TRUE)

  ## ── ESA WorldCover 10 m (satellite-derived), windowed read ───────────────
  wc_url <- "/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N51E015_Map.tif"
  wc <- rast(wc_url)
  ext_ll <- ext(bb["xmin"] - 0.01, bb["xmax"] + 0.01, bb["ymin"] - 0.01, bb["ymax"] + 0.01)
  wc_c <- crop(wc, ext_ll)
  wc_2180 <- project(wc_c, "EPSG:2180", method = "near", res = 10)
  writeRaster(wc_2180, file.path(s1_dir, "worldcover_poznan.tif"), overwrite = TRUE)
  cat("\nWorldCover classes (cells):\n"); print(freq(wc_2180))
  cat("\nSaved", gpkg, "\n")
})

## ── 15a · Land-use descriptors of populated cells ──────────────────────────────────────────────────────────
section("15a", "Land-use descriptors of populated cells", {
  source(file.path(REPO, "R/tissue_descriptors.R")); tissue_descriptors("cells")
})

