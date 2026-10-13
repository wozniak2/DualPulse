# The dual pulse of a city

**How deprivation shapes usage of public transport?**

Marcin Wozniak (Adam Mickiewicz University, Poznań) and Gabriele Filomena (University of Liverpool)

Code and figures for the paper *The dual pulse of a city. How deprivation shapes usage of public transport?* (manuscript in preparation).

![Graphical abstract](figures/Graphical_abstract.png)

## The question

Public transit usage carries a dual pulse: that of the places a network serves and of the people who live there. We ask how usage differs between more and less deprived areas of Poznań, Poland, and whether the differences lie in what the network offers them, in land use and population, or in the use riders make of the network.

```mermaid
flowchart LR
  D["Deprivation<br/>of the area"] -- "sorting by housing" --> P["The network's offer<br/>position, lines, where they lead"]
  D --> L["Land use<br/>what the catchment hosts"]
  D -. "who lives there" .-> R["Residents' circumstances<br/>work hours, car access, income"]
  P --> U["Usage<br/>how much and when · where to ·<br/>what the network offers · responses"]
  L --> U
  R -.-> U
```

Usage combines what the network offers at a stop with the use riders make of it. Dashed: the part that the offer, land use and population do not account for.

## Data

- **Door counts.** Ten months (January–October 2025) of PEKA vehicle location and passenger-count records: 67.9 million operated stop visits, 45% of them counted. Expanded to all operated visits, they give 148 million boardings, about 640,000 on an average term-time weekday.
- **Stops.** 607 stop complexes, each with a mutually exclusive 400 m walking catchment on the OpenStreetMap pedestrian network; 428 complexes with residents are compared.
- **Deprivation.** The city's revitalisation indicators of unemployment and social assistance for its 42 osiedla (2019 and 2021), assigned to catchments through address-level registered population. Municipal housing gives a finer indicator within osiedla.
- **Land use.** Building function and floor area (BDOT10k), the 2023 structure plan, OpenStreetMap amenities, ESA WorldCover tree cover and a modelled employment grid.

## Main findings

- **Door counts separate the offer from its use.** A baseline that spreads each stop's riders over its lines in proportion to departures gives what the network offers there: where its lines lead and how far. The observed measure minus the offer is the use riders make of it, the lines and directions they board.
- **The gaps lie in the offer.** Deprived areas sit near the centre of a radial network (median stop 3.3 against 5.4 km from the Old Market Square), which offers them shorter rides, fewer connections and 32% of the city's jobs within 30 minutes, against 9%. Measured directly, the offer captures 64–84% of the variation between origins in seven of nine destination and route measures, against 11–55% for distance and sector.
- **Destinations record the offer first.** The lines that serve deprived areas run through other deprived districts, so the network would carry their riders to places 0.82 SD more deprived per SD of origin deprivation; most of this (0.66) is geography, because deprived places lie near each other. The use riders make of their lines narrows the gap to 0.60. It is the one robust difference in use, and it comes from riders of the least deprived areas travelling into the inner city more than their lines offer.
- **Two streams that meet in the centre.** Traced along the network, riders from the least deprived third ride north–south, along the rapid tram from Piątkowo and Winogrady, and riders from the most deprived third west–east, through Ogrody, Jeżyce, Łazarz and Śródka. The two share 41% of their ride-kilometres, and the lines offer the same streams. Where riders get off they meet: arrivals overlap 0.73 within 1.5 km of the centre, against 0.41 beyond 5 km. Against random and rotated maps of deprivation, the direction of the streams follows from where deprivation clusters; only the isolation of the least deprived third on its own lines exceeds every map.
- **The advantage is contingent.** If the same 42 deprivation scores rose with distance from the centre, the same network would give the most deprived third 3.8% of the city's jobs within 30 minutes, against 42.5% for the least deprived.
- **Beyond the offer, usage differs little,** and it responds alike to leave days, service cuts and weather. With 42 osiedla and a noisy deprivation score, moderate differences cannot be ruled out.

## Figures

### Main text

| | |
|:---:|:---:|
| ![Construction](figures/Fig_construction.png) | ![The offer](figures/Fig_access_bivariate.png) |
| **Figure 2.** How deprivation is assigned to stop catchments: osiedle score; in one window where four osiedla meet, walking catchments and catchment deprivation; all stops with the rings and sectors that define location. | **Figure 3.** What the network offers: the network over the compact pre-war blocks, deprivation against jobs reachable within 30 minutes, and concentration curves. |
| ![Periphery](figures/Fig_periphery.png) | ![Offer and use](figures/Fig_offer_use.png) |
| **Figure 4.** What the network offers at the periphery: osiedla more than 5 km from the centre, the deprived ones in red. | **Figure 5.** Offer and use: each gap split into what the network offers and the use riders make of it, use with the offer held fixed, and where use changes the deprivation of destinations. |
| ![Two streams](figures/Fig_two_streams.png) | ![Sharing benchmarks](figures/Fig_sharing_benchmark.png) |
| **Figure 6.** Two streams on one network: riding (corridors), getting off (stops) and the direction of each third's ride-kilometres. | **Figure 7.** Sharing the network, on the lines and where riders get off, against deprivation scattered at random and the map turned or mirrored around the centre. |
| ![Decomposition](figures/Fig_decomposition.png) | |
| **Figure 8.** The 22 gaps between more and less deprived catchments and the parts accounted for by location, service, land use and population. | |

Figure 1 is the framework drawn above.

<details>
<summary><b>Supplementary figures (S1–S18)</b></summary>

| | |
|:---:|:---:|
| ![S1](figures/Fig_construction_full.png) | ![S2](figures/Fig_s1ms_associations.png) |
| **S1.** How deprivation is assigned to stop catchments, in full. | **S2.** Land use around residents' homes and deprivation across the 42 osiedla. |
| ![S3](figures/Fig_s1ms_maps.png) | ![S4](figures/Fig_sharpen_onset.png) |
| **S3.** Deprivation by osiedle and land-use types of populated cells. | **S4.** The early-boarding gradient across refined measures. |
| ![S5](figures/Fig_hubsheds.png) | ![S6](figures/Fig_flows_groups.png) |
| **S5.** Catchment deprivation and hub-sheds. | **S6.** Largest morning-peak flows from the least and the most deprived thirds. |
| ![S7](figures/Fig_destination_density.png) | ![S8](figures/Fig_destination_mechanisms.png) |
| **S7.** Where morning-peak riders alight (kernel density). | **S8.** Destination slopes split into the network part and the lines used. |
| ![S9](figures/Fig_network_vs_use_map.png) | ![S10](figures/Fig_connections_map.png) |
| **S9.** Destinations offered by the network and reached, by origin. | **S10.** Where riders connect to other lines. |
| ![S11](figures/Fig_route_networks_difference.png) | ![S12](figures/Fig_route_heterogeneity_map.png) |
| **S11.** Route networks: most minus least deprived third. | **S12.** Route heterogeneity by origin. |
| ![S13](figures/Fig_hub_use_map.png) | ![S14](figures/Fig_route_networks.png) |
| **S13.** Use of the interchanges, observed and offered. | **S14.** Route networks on weekdays and at weekends. |
| ![S15](figures/Fig_behaviour.png) | ![S16](figures/Fig_extremes.png) |
| **S15.** Responses to leave days, service cuts and weather. | **S16.** The extremes of deprivation. |
| ![S17](figures/Fig_gap_ladder.png) | ![S18](figures/Fig_gap_loo.png) |
| **S17.** Slopes raw, net of position and net of land use and population. | **S18.** Remainders with one osiedle left out. |

Drawn for the supplement, not yet numbered in it:

| | |
|:---:|:---:|
| ![Route networks](figures/Fig_route_networks_main.png) | ![Offer and weekends](figures/Fig_two_streams_offer_weekend.png) |
| Route networks of the least and the most deprived thirds, and their difference. | The two streams as the lines offer them on weekdays, and as ridden at weekends. |
| ![Rotation benchmark](figures/Fig_rotation_benchmark.png) | ![Meeting places](figures/Fig_meeting_places.png) |
| Each measure of sharing over the 47 turned and mirrored maps of deprivation. | Where the two thirds get off, and the share of each third's arrivals at the 24 busiest stops. |
| ![Deprived areas](figures/Fig_deprived_areas_map.png) | ![Destinations](figures/Fig_destinations_main.png) |
| Where the more and less deprived areas lie, and the network over the compact pre-war blocks (formerly Figure 3). | Deprivation of destinations offered, reached, and their difference (formerly Figure 5). |

</details>

## Pipeline

Nine stage scripts, run in order. Each is a sequence of sections; a section is one step of the analysis and keeps the number of the script it came from in the project notes (for example, section 40 is the trip-level cleaning).

```mermaid
flowchart LR
  S1["1 · Areas and<br/>population"] --> S2["2 · Land use"] --> S3["3 · Deprivation score<br/>and land use"]
  S1 --> S4["4 · Stops and<br/>catchments"]
  S3 --> S4
  K[("PEKA door counts<br/>Jan–Oct 2025")] --> S5["5 · Passenger counts"]
  S4 --> S5
  S5 --> S6["6 · Destinations,<br/>routes, access"]
  S5 --> S7["7 · Usage and<br/>responses"]
  S6 --> S8["8 · Gaps and<br/>inference"]
  S7 --> S8
  S8 --> S9["9 · Figures"]
```

| Stage | Sections | What it does |
|---|---|---|
| [`01_areas_population.R`](scripts/01_areas_population.R) | 01–07 | Osiedle boundaries, the 200 m grid, address points, registered residents geocoded to them, the 2023 zoning, populated cells. |
| [`02_land_use.R`](scripts/02_land_use.R) | 14, 15a | BDOT10k, OpenStreetMap and WorldCover inputs; land-use descriptors of populated cells. |
| [`03_deprivation_land_use.R`](scripts/03_deprivation_land_use.R) | 16, 18–20 | The deprivation score and where deprived areas sit in land use (supplementary S4). |
| [`04_stops_catchments.R`](scripts/04_stops_catchments.R) | 21, 39, 15b, 34, 37, 62 | Stop complexes, walking catchments, their land use, residents, deprivation and municipal housing. |
| [`05_passenger_counts.R`](scripts/05_passenger_counts.R) | ingest, 40–42, 22, 23, 58 | Ingest, trip-level cleaning, stop-to-stop flows by three estimators, expansion, same-trip indices, descriptive checks. |
| [`06_destinations_routes_access.R`](scripts/06_destinations_routes_access.R) | 45–47, 50, 51, 53, 54 | Destinations against a network baseline, interchanges and transfers, routes, timetable-based access. |
| [`07_usage_responses.R`](scripts/07_usage_responses.R) | 43, 48, 52, 61 | Tests stated before modelling, early boardings, responses to leave days, service cuts and weather, tram and bus. |
| [`08_gaps_inference.R`](scripts/08_gaps_inference.R) | 56, 55, 59, 63, 65–68, 73–83 | The 22 gaps raw, net of position and net of land use and population; decomposition; CR2-Satterthwaite inference, wild cluster bootstrap, measurement error, equivalence; spatial placebo; equity and counterfactual geographies; estimator checks; offer and use; peripheral deprivation with osiedle-level randomisation inference; the early offer and its use; what location stands for (jobs nearby, structure of service); segregation by destination (geography, network, use); sharing the network on corridors, lines and stops, against deprivation scattered at random and the map turned or mirrored; Figures 4–7. |
| [`09_figures.R`](scripts/09_figures.R) | 57, 60, 64, 71, 70, 69, 72 | Main-text figures 2, 3 and 8, maps for the supplement and the graphical abstract (Figures 4–7 are drawn by sections 75, 74 and 81). |

Helpers in [`R/`](R): `setup.R` (paths and the section runner), `od_expand.R` (expansion of counted flows to all operated trips), `cr_inference.R` (CR2 standard errors with Satterthwaite degrees of freedom, effective clusters, wild cluster bootstrap), `seg_index.R` (indices of sharing the network and the turned maps of deprivation), `map_style.R` (shared cartographic style) and `tissue_descriptors.R` (land-use descriptors for cells and catchments).

### Methods in brief

- **Same-trip indices.** Each usage measure compares a stop with the other stops served by the same trips, so the timetable is held fixed.
- **Destinations from door counts.** Counts do not link boardings to alightings. On every counted trip, the Tsygalnitsky recursion (the biproportional fit of a uniform seed) assigns alightings to everyone on board, reproducing boardings, alightings and loads. A distance-seeded fit, the earlier proportional rule and a synthetic test serve as checks. A network baseline assigns riders to every line and direction at their stop in proportion to departures.
- **Transfers.** Boardings on each departure at the 64 interchanges are regressed on alightings from other lines in the preceding 10 minutes, with later alightings as a placebo.
- **Offer and use.** For the nine destination and route measures, each observed slope on deprivation splits into the part offered by the network baseline and the use made of it; use is also estimated with the offer held fixed.
- **Sharing the network.** Each estimated ride is traced along its line's stop sequence and takes the deprivation of the catchment where it boards. Overlap, isolation (net of group size), exposure and the variance ratio of riders' deprivation are computed on corridors, on lines (line, direction, segment and period) and at the stops where riders get off, and compared with the same flows under 999 random allocations of the 42 osiedle scores and 47 turned or mirrored maps of them.
- **Comparisons.** Slopes of each standardised measure on standardised catchment deprivation, raw, net of position (distance spline, sector and, for measures other than the same-trip indices, service at the stop) and net of land use and population. Gelbach decomposition; CR2 standard errors by osiedle with Satterthwaite degrees of freedom; Benjamini–Hochberg and Romano–Wolf control; correction for the score's reliability; equivalence tests at ±0.2 SD.

## Running it

Requirements: R 4.4 with `data.table`, `sf`, `terra`, `exactextractr`, `nabor`, `igraph`, `osmdata`, `curl`, `jsonlite`, `dplyr`, `tidyr`, `stringr`, `readr`, `sandwich`, `cluster`, `Matrix`, `ggplot2`, `patchwork`, `scales` and `ggrepel`; 7-Zip for the ingest.

1. Copy `config_template.R` to `config_local.R` and set the paths. `ROOT` holds `data/`; results go to `ROOT/output/` and the paper's figures and tables to `ROOT/manuscript/`.
2. From the repository root:

```r
Rscript run.R            # all stages
Rscript run.R 6          # one stage
Rscript run.R 8:56,63    # sections 56 and 63 of stage 8
```

Each stage runs in its own R process, and each section in a fresh frame with only base packages attached, as if it were a script of its own.

## Data availability

The passenger-count records and the address-level register of residents are not public, so the pipeline cannot be rerun from this repository alone. The small public inputs are in [`data/`](data), with a list of everything else the pipeline reads and where it comes from.

## Repository structure

```
DualPulse/
├── run.R                 runs the stages
├── config_template.R     paths (copy to config_local.R)
├── R/                    setup and helpers
├── scripts/              the nine stages
├── figures/              figures of the paper
└── data/                 public inputs and a list of the rest
```
