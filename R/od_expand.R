## Helper, sourced by the scripts that use the stop-to-stop flows of script 40
## (45, 46, 47, 50, 51, 53, 56). The flows come from valid counted trips only,
## and counter coverage differs between lines, so every flow is expanded to all
## operated trips: it is multiplied by the expansion factor of its origin
## stratum, operated over validly counted visits of the boarding stop x line x
## direction x month x day type x period (table od_expansion_factors.rds,
## written by script 41). Fallbacks for strata without counted visits in that
## table (mainly stop-days excluded as diversions): the same without the month,
## then line x direction x month x day type x period over all stops, then 1.
## Both flow and fkm (flow x ride km) are expanded. The counted flow is kept as
## flow_c, so that minimum-boarding thresholds for origins can still refer to
## counted support (the samples of origins stay those of the counted analysis).
## OD estimator: the environment variable OD_VARIANT selects the flows of script 40:
## "tsy" (default; Tsygalnitsky recursion = null-seed biproportional fit), "dist"
## (biproportional fit with a distance seed) or "prop" (proportional to downstream
## alightings, the earlier rule). The chosen variant becomes flow/fkm.

## vtag(): output (and variant-specific input) file names carry the variant, e.g. dest_mech_profiles_prop.csv,
## so runs with OD_VARIANT = dist or prop never overwrite the main (tsy) results
vtag <- function(f, variant = Sys.getenv("OD_VARIANT", "tsy")) if (variant == "tsy") f else sub("(\\.[A-Za-z]+)$", paste0("_", variant, "\\1"), f)

expand_od <- function(x, s2 = file.path(ROOT, "output/stage2"),
                      variant = Sys.getenv("OD_VARIANT", "tsy")) {
  if (!exists(".od_expansion", envir = .GlobalEnv)) assign(".od_expansion", readRDS(file.path(s2, "od_expansion_factors.rds")), envir = .GlobalEnv)
  F <- get(".od_expansion", envir = .GlobalEnv)
  if (variant %in% c("dist", "prop")) x[, `:=`(flow = get(paste0("flow_", variant)), fkm = get(paste0("fkm_", variant)))]
  else if (variant != "tsy") stop("OD_VARIANT must be tsy, dist or prop")
  vc <- intersect(c("flow_prop", "fkm_prop", "flow_dist", "fkm_dist"), names(x)); if (length(vc)) x[, (vc) := NULL]
  x <- merge(x, F$e1, by = c("o", "NR_LINI", "dir", "month", "daytype", "period"), all.x = TRUE)
  x <- merge(x, F$e2, by = c("o", "NR_LINI", "dir", "daytype", "period"), all.x = TRUE)
  x <- merge(x, F$e3, by = c("NR_LINI", "dir", "month", "daytype", "period"), all.x = TRUE)
  x[, e := fcoalesce(e1, e2, e3, 1)]
  x[, flow_c := flow]
  x[, `:=`(flow = flow * e, fkm = fkm * e)]
  x[, c("e1", "e2", "e3", "e") := NULL]
  x[]
}
