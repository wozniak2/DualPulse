# The dual pulse of a city

**How does deprivation shape the use of public transport?**

Marcin Wozniak (Adam Mickiewicz University, Poznań) and Gabriele Filomena (University of Liverpool)

Code and figures for the paper *The dual pulse of a city. How does deprivation shape the use of public transport?* (manuscript in preparation).

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
- **Deprivation.** The city's revitalisation indicators of unemployment and social assistance for its 42 neighbourhoods (osiedla) (2019 and 2021), assigned to catchments through address-level registered population. Municipal housing gives a finer indicator within neighbourhoods.
- **Land use.** Building function and floor area (BDOT10k), the 2023 structure plan, OpenStreetMap amenities, ESA WorldCover tree cover and a modelled employment grid.

## Main findings

- **Door counts separate the offer from its use.** A baseline that spreads each stop's riders over its lines in proportion to departures gives what the network offers there: where its lines lead and how far. The observed measure minus the offer is the use riders make of it, the lines and directions they board.
- **The gaps lie in the offer.** Deprived areas sit near the centre of a radial network (median stop 3.3 against 5.4 km from the Old Market Square), which offers them shorter rides, fewer connections and 32% of the city's jobs within 30 minutes, against 9%. Measured directly, the offer captures 64–84% of the variation between origins in seven of nine destination and route measures, against 11–55% for distance and sector.
- **Destinations record the offer first.** The lines that serve deprived areas run through other deprived districts, so the network would carry their riders to places 0.82 SD more deprived per SD of origin deprivation; most of this (0.66) is geography, because deprived places lie near each other. The use riders make of their lines narrows the gap to 0.60. It is the one robust difference in use, and it comes from riders of the least deprived areas travelling into the inner city more than their lines offer.
- **Two streams that meet in the centre.** Traced along the network, riders from the least deprived third ride north–south, along the rapid tram from Piątkowo and Winogrady, and riders from the most deprived third west–east, through Ogrody, Jeżyce, Łazarz and Śródka. The two share 41% of their ride-kilometres, and the lines offer the same streams. Where riders get off they meet: arrivals overlap 0.73 within 1.5 km of the centre, against 0.41 beyond 5 km. Against random and rotated maps of deprivation, the direction of the streams follows from where deprivation clusters; only the isolation of the least deprived third on its own lines exceeds every map.
- **The advantage is contingent.** If the same 42 deprivation scores rose with distance from the centre, the same network would give the most deprived third 3.8% of the city's jobs within 30 minutes, against 42.5% for the least deprived.
- **Beyond the offer, usage differs little,** and it responds alike to leave days, service cuts and weather. With 42 neighbourhoods and a noisy deprivation score, moderate differences cannot be ruled out.

## Figures

### Major results

| |
|:---:|
| ![The offer](figures/Fig_access_bivariate.png) |
| **Figure 3. Deprived core, prime access.** The network over the compact pre-war blocks, catchments by deprivation and by jobs reachable within 30 minutes, and concentration curves: access is progressive, but by position. |
| ![Two streams](figures/Fig_two_streams.png) |
| **Figure 6. Two streams that cross at the centre.** Riding (corridors), getting off (stops) and the direction of each third's ride-kilometres. |
| ![Decomposition](figures/Fig_decomposition.png) |
| **Figure 8. Places, not people.** The 22 gaps between more and less deprived catchments and the parts accounted for by location, service, land use and population. |

<details>
<summary><b>The other figures of the paper</b></summary>

| | |
|:---:|:---:|
| ![Construction](figures/Fig_construction.png) | ![Periphery](figures/Fig_periphery.png) |
| **Figure 2.** How deprivation is assigned to stop catchments. | **Figure 4.** What the network offers at the periphery: the deprived neighbourhoods in red. |
| ![Offer and use](figures/Fig_offer_use.png) | ![Sharing benchmarks](figures/Fig_sharing_benchmark.png) |
| **Figure 5.** Offer and use in where riders are carried. | **Figure 7.** Sharing the network against deprivation scattered at random and the map turned or mirrored around the centre. |
| ![S1](figures/Fig_construction_full.png) | |
| **Supplementary Figure S1.** How deprivation is assigned to stop catchments, in full. | |

Figure 1 is the framework drawn above. Exploratory figures drawn by the pipeline but not used in the paper are in [`figures/`](figures/).

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
| [`01_areas_population.R`](scripts/01_areas_population.R) | 01–07 | Neighbourhood (neighbourhood) boundaries, the 200 m grid, address points, registered residents geocoded to them, the 2023 zoning, populated cells. |
| [`02_land_use.R`](scripts/02_land_use.R) | 14, 15a | BDOT10k, OpenStreetMap and WorldCover inputs; land-use descriptors of populated cells. |
| [`03_deprivation_land_use.R`](scripts/03_deprivation_land_use.R) | 16, 18–20 | The deprivation score and where deprived areas sit in land use (supplementary S4). |
| [`04_stops_catchments.R`](scripts/04_stops_catchments.R) | 21, 39, 15b, 34, 37, 62 | Stop complexes, walking catchments, their land use, residents, deprivation and municipal housing. |
| [`05_passenger_counts.R`](scripts/05_passenger_counts.R) | ingest, 40–42, 22, 23, 58 | Ingest, trip-level cleaning, stop-to-stop flows by three estimators, expansion, same-trip indices, descriptive checks. |
| [`06_destinations_routes_access.R`](scripts/06_destinations_routes_access.R) | 45–47, 50, 51, 53, 54 | Destinations against a network baseline, interchanges and transfers, routes, timetable-based access. |
| [`07_usage_responses.R`](scripts/07_usage_responses.R) | 43, 48, 52, 61 | Tests stated before modelling, early boardings, responses to leave days, service cuts and weather, tram and bus. |
| [`08_gaps_inference.R`](scripts/08_gaps_inference.R) | 56, 55, 59, 63, 65–68, 73–86 | The 22 gaps raw, net of position and net of land use and population; decomposition; CR2-Satterthwaite inference, wild cluster bootstrap, measurement error, equivalence; spatial placebo; equity and counterfactual geographies; estimator checks; offer and use; peripheral deprivation with neighbourhood-level randomisation inference; the early offer and its use; what location stands for (jobs nearby, structure of service); segregation by destination (geography, network, use); sharing the network on corridors, lines and stops, against deprivation scattered at random and the map turned or mirrored; usage and the link to the centre at the periphery; the offer weighted by departures or by boardings; the supplement tables of sharing, re-allocations and the periphery; Figures 4–7. |
| [`09_figures.R`](scripts/09_figures.R) | 57, 60, 64, 71, 70, 69, 72 | Main-text figures 2, 3 and 8, maps for the supplement and the graphical abstract (Figures 4–7 are drawn by sections 75, 74 and 81). |

Helpers in [`R/`](R): `setup.R` (paths and the section runner), `od_expand.R` (expansion of counted flows to all operated trips), `cr_inference.R` (CR2 standard errors with Satterthwaite degrees of freedom, effective clusters, wild cluster bootstrap), `seg_index.R` (indices of sharing the network and the turned maps of deprivation), `map_style.R` (shared cartographic style) and `tissue_descriptors.R` (land-use descriptors for cells and catchments).

### Methods in brief

- **Same-trip indices.** Each usage measure compares a stop with the other stops served by the same trips, so the timetable is held fixed.
- **Destinations from door counts.** Counts do not link boardings to alightings. On every counted trip, the Tsygalnitsky recursion (the biproportional fit of a uniform seed) assigns alightings to everyone on board, reproducing boardings, alightings and loads. A distance-seeded fit, the earlier proportional rule and a synthetic test serve as checks. A network baseline assigns riders to every line and direction at their stop in proportion to departures.
- **Transfers.** Boardings on each departure at the 64 interchanges are regressed on alightings from other lines in the preceding 10 minutes, with later alightings as a placebo.
- **Offer and use.** For the nine destination and route measures, each observed slope on deprivation splits into the part offered by the network baseline and the use made of it; use is also estimated with the offer held fixed.
- **Sharing the network.** Each estimated ride is traced along its line's stop sequence and takes the deprivation of the catchment where it boards. Overlap, isolation (net of group size), exposure and the variance ratio of riders' deprivation are computed on corridors, on lines (line, direction, segment and period) and at the stops where riders get off, and compared with the same flows under 999 random allocations of the 42 neighbourhood scores and 47 turned or mirrored maps of them.
- **Comparisons.** Slopes of each standardised measure on standardised catchment deprivation, raw, net of position (distance spline, sector and, for measures other than the same-trip indices, service at the stop) and net of land use and population. Gelbach decomposition; CR2 standard errors by neighbourhood with Satterthwaite degrees of freedom; Benjamini–Hochberg and Romano–Wolf control; correction for the score's reliability; equivalence tests at ±0.2 SD.

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
