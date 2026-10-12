# Data

The pipeline reads its inputs from `ROOT/data/` (see `config_template.R`). Only small public inputs are tracked here.

## Tracked (public)

| File | What it is | Used by |
|---|---|---|
| `gpr_indicators_ABC_osiedla.csv` | Registered unemployment (A) and social-assistance beneficiaries (B) per 1,000 residents, and turnout (C), with the city's z-scores, for the 42 osiedla in 2015, 2019 and 2021. Transcribed from the city's revitalisation monitoring report (*Gminny Program Rewitalizacji*, 2023 assessment); every A and C value reproduces its printed z-score within 0.02. | sections 19, 20, 57 |
| `diagnoza2016_osiedla_ground_truth.csv`, `diagnoza2019_…`, `diagnoza2021_…` | Indicator B and its z-score per osiedle in each round, from the same reports (2016 also gives osiedle areas, used as a geometry check). | sections 01, 07, 19, 20, 57 |
| `raw/spp/parking_meters_2026-10-07.geojson` | Parking meters of the paid-parking zone (city open data). | section 34 |
| `raw/weather/s_d_330_2025.csv` | Daily synoptic records, Poznań-Ławica, 2025 (IMGW-PIB public data), with the format notes. | section 52 |

## Not tracked

| Input | Where it comes from |
|---|---|
| `pop_geocoded.csv` | Registered residents by address, age and sex. Not public. |
| PEKA vehicle location and door-count records, January–October 2025 | Not public. The ingest (stage 5) expects the monthly archives in `PEKA_RAW` and writes `peka_2025_MM.rds` to `PEKA_DIR/output/tables`. |
| GTFS timetable, October 2025 | ZTM Poznań open data (`PEKA_DIR/gtfs/extracted_202510`). |
| `raw/zkzl/wykaz_budynkow_komunalnych_2026-10-07.html` | The public list of municipal buildings of the city's housing company (ZKZL), saved from https://zkzl.poznan.pl/administracja-i-eksploatacja-budynkow/wykaz-budynkow-komunalnych/ on 7 October 2026 (section 37). Not tracked because the saved page embeds the site's own Google Maps key. |
| `raw/bdot10k/` | BDOT10k for powiat 3064 (GUGiK open data): `https://opendata.geoportal.gov.pl/bdot10k/schemat2021/30/3064_GML.zip`, unzipped. |
| `raw/osm/poznan_walk_ways.rds` | OpenStreetMap pedestrian network, downloaded by section 39 through Overpass and cached. |
| `workplace_grid.gpkg` | Modelled 200 m employment grid (jobs per cell). |
| `poz.gpkg` (`CITY_DIR`) | City boundary of Poznań. |
| `rcn_poznan_metro.gpkg` (`RCN_GPKG`) | Property transactions from the national price register (RCN), used for an apartment price index in section 19. |
| OpenStreetMap extract (`OSM_EXTRACT`) | Amenities for section 14; cached in `output/stage1/_osm_amenities_raw.rds` after the first run. |

Osiedle boundaries (section 01), address points (03) and the 2023 structure-plan zoning (06) are downloaded from OpenStreetMap and the city's GEOPOZ services on the first run and cached in `output/`.
