





#### ============================================================
#### SCRIPT 01b (v2, LIGHTWEIGHT): BUILD REQUIRED DRE SPATIAL OBJECTS
#### ============================================================
#### Run this in the SAME R session where you already ran the AMBER
#### v2 prep script — it reuses rivers_sf, amber_basins_sf,
#### inverts_basins_sf, and basins_sf, which are already built there.
#### No files are reloaded. This only adds the one thing AMBER's
#### snap_to_river() never computed: position ALONG the matched river
#### segment (line_distance_m), needed for the DRE "same_segment"
#### network-distance calculation.
####
#### REQUIRES: rivers_sf, amber_basins_sf, inverts_basins_sf, basins_sf
#### already in your session (all built by 01_prepare_amber_data_v2.R).
#### ============================================================

library(purrr)   #### for map2_dbl() in the linear-referencing step below

required_from_amber_script <- c("rivers_sf", "amber_basins_sf", "inverts_basins_sf", "basins_sf")
missing_now <- required_from_amber_script[!sapply(required_from_amber_script, exists)]
if (length(missing_now) > 0) {
  stop("Missing from session: ", paste(missing_now, collapse = ", "),
       ". Run the AMBER v2 prep script (at least through Section 4) in this ",
       "same session first.")
}

SNAP_THRESHOLD_M <- 2000   #### matches the AMBER script's SNAP_THRESHOLD_M

#### No external package needed for linear referencing — lwgeom's
#### exported functions vary by version and st_linelocate_point isn't
#### available in yours. This does the same thing manually: project
#### the point onto each segment of the matched line, keep the
#### closest projection, and return the cumulative distance along the
#### line to that point. Standard point-to-polyline linear referencing.
get_line_distance_m <- function(line_geom, point_geom) {
  coords <- sf::st_coordinates(line_geom)[, c("X", "Y"), drop = FALSE]
  pt <- as.numeric(sf::st_coordinates(point_geom)[1, c("X", "Y")])
  
  seg_vec_x <- diff(coords[, "X"]); seg_vec_y <- diff(coords[, "Y"])
  seg_lengths <- sqrt(seg_vec_x^2 + seg_vec_y^2)
  cum_dist <- c(0, cumsum(seg_lengths))
  
  n_seg <- length(seg_lengths)
  best_dist <- Inf; best_pos <- 0
  
  for (i in seq_len(n_seg)) {
    Ax <- coords[i, "X"]; Ay <- coords[i, "Y"]
    Bx <- coords[i + 1, "X"]; By <- coords[i + 1, "Y"]
    ABx <- Bx - Ax; ABy <- By - Ay
    APx <- pt[1] - Ax; APy <- pt[2] - Ay
    denom <- ABx^2 + ABy^2
    t <- if (denom > 0) (APx * ABx + APy * ABy) / denom else 0
    t <- max(0, min(1, t))
    projx <- Ax + t * ABx; projy <- Ay + t * ABy
    d <- sqrt((pt[1] - projx)^2 + (pt[2] - projy)^2)
    if (d < best_dist) { best_dist <- d; best_pos <- cum_dist[i] + t * seg_lengths[i] }
  }
  best_pos
}

#### Snap-with-linear-referencing — same core logic as AMBER's
#### snap_to_river(), plus line_distance_m via the manual function
#### above. Assumes points_sf ALREADY has HYBAS_ID (true for
#### amber_basins_sf and inverts_basins_sf, both already joined to
#### basins_sf in the AMBER script's Section 4). NOTE: this loops row
#### by row for the linear-referencing step, so it's slower than a
#### vectorised operation — expect it to take a little while for
#### community_position_df given the larger row count, but each
#### individual line has few vertices so it should still be minutes,
#### not hours.
snap_to_river_with_position <- function(points_sf, rivers_sf, threshold_m = 2000) {
  pts_3035 <- st_transform(points_sf, 3035)
  riv_3035 <- st_transform(rivers_sf, 3035)
  
  nearest_idx <- st_nearest_feature(pts_3035, riv_3035)
  dist_m <- as.numeric(st_distance(pts_3035, riv_3035[nearest_idx, ], by_element = TRUE))
  keep <- dist_m <= threshold_m
  
  pts_keep <- pts_3035[keep, ]
  riv_matched <- riv_3035[nearest_idx[keep], ]
  
  cat("Computing linear referencing for", nrow(pts_keep), "points...\n")
  line_distance_m <- purrr::map2_dbl(
    st_geometry(riv_matched), st_geometry(pts_keep), get_line_distance_m
  )
  
  out <- bind_cols(
    st_drop_geometry(pts_keep),
    st_drop_geometry(riv_matched) %>% dplyr::select(HYRIV_ID, NEXT_DOWN, MAIN_RIV, LENGTH_KM, DIST_DN_KM)
  )
  out$dist_to_river_m <- dist_m[keep]
  out$line_distance_m <- line_distance_m
  out
}

#### 1. hydrorivers_sf — just an alias, already loaded #############
hydrorivers_sf <- rivers_sf

#### 2. amber_position_df — reuses amber_basins_sf (already has
#### AMBER_ID, barrier_type_reduced, HYBAS_ID from the AMBER script)
amber_position_df <- snap_to_river_with_position(amber_basins_sf, rivers_sf, SNAP_THRESHOLD_M) %>%
  mutate(amber_type = type, amber_type_reduced = barrier_type_reduced) %>%
  filter(!is.na(HYBAS_ID))

cat("\namber_position_df built:", nrow(amber_position_df), "barriers\n")

#### 3. community_position_df — reuses inverts_basins_sf, but that's
#### one row per TAXON RECORD (many per site-year); dedupe to one row
#### per site_id x sample_year BEFORE snapping, since they all share
#### the same coordinates and would otherwise be redundantly snapped
#### hundreds of thousands of times.
community_points_dedup <- inverts_basins_sf %>%
  distinct(site_id, sample_year, .keep_all = TRUE)

cat("\nDeduplicated to", nrow(community_points_dedup), "unique site-years (from",
    nrow(inverts_basins_sf), "taxon records) before snapping.\n")

community_position_df <- snap_to_river_with_position(community_points_dedup, rivers_sf, SNAP_THRESHOLD_M) %>%
  filter(!is.na(HYBAS_ID))

cat("community_position_df built:", nrow(community_position_df),
    "site-years,", n_distinct(community_position_df$site_id), "sites\n")

cat("\n============================================================\n")
cat("3 of 4 required objects built: hydrorivers_sf, amber_position_df,\n")
cat("community_position_df. dre_position_df still needs your raw DRE\n")
cat("removals file — see 01c_build_dre_position_df.R next (its Step C\n")
cat("now reuses THIS script's snap_to_river_with_position(), which\n")
cat("expects HYBAS_ID to already be joined — add an explicit\n")
cat("st_join(dre_raw_sf, basins_sf %>% select(HYBAS_ID), join = st_intersects)\n")
cat("step before calling it, since dre_raw_sf won't have HYBAS_ID yet\n")
cat("the way amber_basins_sf/inverts_basins_sf already did.\n")
cat("============================================================\n")





#### 1c ============================================================
#### STEP A — INSPECT THE DRE DATABASE FIRST
#### ============================================================
#### Run this on its own first. It just loads the file and shows you
#### the real column names/structure, since I don't know them yet.
#### ============================================================

library(readxl)
library(dplyr)

dre_raw_path <- "C:/Users/User/OneDrive - Queen's University Belfast/PhD/Chapters/Chapter 3 European Connectivity (Inverts)/UpToDateWork/EUConnectivity_Inverts/Data/Removed Barriers/Copy of DRE database_Ellen Dolan - this is what we used.xlsx"

#### If the file isn't actually .xlsx, this will error clearly and
#### tell you so — swap read_xlsx() for read_csv() (readr package) if
#### it turns out to be a .csv instead.
dre_raw_inspect <- read_xlsx(dre_raw_path)

cat("\nColumn names in your DRE database:\n")
print(names(dre_raw_inspect))

cat("\nFirst few rows:\n")
print(head(dre_raw_inspect))

cat("\nStructure (types):\n")
str(dre_raw_inspect)

#### ============================================================
#### STEP B — BUILD dre_position_df
#### ============================================================
#### ADJUST THE COLUMN NAMES BELOW to match what Step A printed.
#### My best guesses, based on the column names already referenced
#### throughout the DRE analysis code you'd previously shared
#### (removal_id, YearRemoved, removal_type), are used as defaults —
#### but the coordinate column names are a genuine guess and most
#### likely to need adjusting.
#### ============================================================

library(sf)
library(readxl)
library(dplyr)
library(stringr)
library(lwgeom)

select <- dplyr::select
filter <- dplyr::filter

sf_use_s2(FALSE)

clean_removal_type <- function(x) {
  x <- str_trim(tolower(as.character(x)))
  case_when(str_detect(x, "dam") ~ "dam", str_detect(x, "weir") ~ "weir",
            str_detect(x, "culvert") ~ "culvert", TRUE ~ "other")
}

parse_num <- function(x) {
  x <- str_trim(as.character(x)); x <- na_if(x, ""); x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

dre_raw <- read_xlsx(dre_raw_path)

#### ---- ADJUST THESE COLUMN NAMES if Step A showed different ones ----
colnames(dre_raw)
colnames(dre_raw)[4] <- "YearRemoved"
colnames(dre_raw)[6] <- "Latitude"
colnames(dre_raw)[7] <- "Longitude"

dre_raw_sf <- dre_raw %>%
  transmute(
    removal_id = row_number(),
    YearRemoved = suppressWarnings(as.integer(YearRemoved)),      #### <- adjust if named e.g. "Year", "Year_removed", "Year of removal"
    removal_type = as.character(Type),                             #### <- adjust if named e.g. "Barrier_type", "removal_type"
    removal_type_reduced = clean_removal_type(removal_type),
    Longitude_X = parse_num(Longitude),                             #### <- adjust if named e.g. "Long", "X", "Lon"
    Latitude_Y = parse_num(Latitude)                                #### <- adjust if named e.g. "Lat", "Y"
  ) %>%
  filter(!is.na(YearRemoved), !is.na(Longitude_X), !is.na(Latitude_Y)) %>%
  st_as_sf(coords = c("Longitude_X", "Latitude_Y"), crs = 4326, remove = FALSE)

cat("\nDRE removals with valid year + coordinates:", nrow(dre_raw_sf), "\n")
cat("Year range:", paste(range(dre_raw_sf$YearRemoved, na.rm = TRUE), collapse = "-"), "\n")



#### ============================================================
#### STEP C — SNAP TO HYDRORIVERS (reuses snap_to_river_with_position()
#### from 01b_build_dre_spatial_objects.R — run that script FIRST in
#### this same session, so hydrorivers_sf, basins_sf, and the helper
#### function already exist)
#### ============================================================

if (!exists("snap_to_river_with_position") || !exists("hydrorivers_sf") || !exists("basins_sf")) {
  stop("Run 01b_build_dre_spatial_objects.R FIRST in this same R session — ",
       "it defines snap_to_river_with_position(), hydrorivers_sf, and basins_sf.")
}

SNAP_THRESHOLD_M <- 2000

#### snap_to_river_with_position() expects HYBAS_ID already joined —
#### unlike amber_basins_sf/inverts_basins_sf (which got this from the
#### AMBER script's Section 4), dre_raw_sf is brand new and needs it
#### done explicitly here first.
dre_basins_sf <- st_join(dre_raw_sf, basins_sf %>% dplyr::select(HYBAS_ID), join = st_intersects) %>%
  filter(!is.na(HYBAS_ID))

cat("\nDRE removals retained after point-in-polygon join to HydroBASINS:",
    nrow(dre_basins_sf), "of", nrow(dre_raw_sf), "\n")

dre_position_df <- snap_to_river_with_position(dre_basins_sf, hydrorivers_sf, SNAP_THRESHOLD_M) %>%
  filter(!is.na(HYBAS_ID))

cat("\ndre_position_df built:", nrow(dre_position_df), "removals retained (within 2km snap threshold)\n")
cat("Removals dropped by the 2km snap:", nrow(dre_raw_sf) - nrow(dre_position_df), "\n")

cat("\nRemoval type breakdown:\n")
print(table(dre_position_df$removal_type_reduced))

cat("\n============================================================\n")
cat("All 4 required objects should now exist: hydrorivers_sf,\n")
cat("amber_position_df, community_position_df, dre_position_df.\n")
cat("Confirm with: sapply(c('hydrorivers_sf','amber_position_df',\n")
cat("  'community_position_df','dre_position_df'), exists)\n")
cat("Then run 02_prepare_dre_data_v2.R.\n")
cat("============================================================\n")



#### Quick structural sanity check before the big DRE prep run ####
cat("Required spatial objects present:\n")
print(sapply(c("hydrorivers_sf", "amber_position_df", "community_position_df", "dre_position_df"), exists))

cat("\nRequired AMBER intermediate objects present (needed for covariate joins):\n")
print(sapply(c("climate_wide", "river_len_df", "hyde_full", "hyde_1961", "site_river_position"), exists))

cat("\ndre_position_df — column check:\n")
print(names(dre_position_df))
cat("\ndre_position_df — required columns present:\n")
print(all(c("removal_id", "YearRemoved", "removal_type", "removal_type_reduced",
            "HYBAS_ID", "HYRIV_ID", "MAIN_RIV", "line_distance_m") %in% names(dre_position_df)))

cat("\ndre_position_df — year range and row count:\n")
cat("Rows:", nrow(dre_position_df), "\n")
cat("Year range:", paste(range(dre_position_df$YearRemoved, na.rm = TRUE), collapse = "-"), "\n")
cat("Removal type breakdown:\n")
print(table(dre_position_df$removal_type_reduced))

cat("\namber_position_df — required columns present:\n")
print(all(c("AMBER_ID", "HYRIV_ID", "amber_type", "amber_type_reduced", "line_distance_m") %in% names(amber_position_df)))

cat("\ncommunity_position_df — required columns present:\n")
print(all(c("site_id", "sample_year", "country", "method", "HYBAS_ID", "HYRIV_ID", "MAIN_RIV", "line_distance_m") %in% names(community_position_df)))
cat("Unique sites:", n_distinct(community_position_df$site_id), "\n")
cat("Unique site-years:", nrow(community_position_df), "\n")


### save important things #####
saveRDS(climate_wide, "Data/Processed/climate_wide.rds")
saveRDS(river_len_df, "Data/Processed/river_len_df.rds")
saveRDS(hyde_full, "Data/Processed/hyde_full.rds")
saveRDS(hyde_1961, "Data/Processed/hyde_1961.rds")
saveRDS(site_river_position, "Data/Processed/site_river_position.rds")

#### Also worth saving the 4 spatial objects you just spent time building —
#### rebuilding dre_position_df in particular took real effort.
saveRDS(hydrorivers_sf, "Data/Processed/hydrorivers_sf.rds")
saveRDS(amber_position_df, "Data/Processed/amber_position_df.rds")
saveRDS(community_position_df, "Data/Processed/community_position_df.rds")
saveRDS(dre_position_df, "Data/Processed/dre_position_df.rds")



#### ============================================================
#### SCRIPT 02: PREPARE DRE (BARRIER REMOVAL) DATA FOR FINAL MODELS ##########
#### v2 — REBUILT TO MATCH THE AMBER 4-MODEL STRUCTURE EXACTLY
####
#### Produces TWO model structures x FOUR responses = 8 model-ready
#### datasets, all sharing the SAME random-effect structure as AMBER
#### and the SAME response definitions as AMBER:
####
#### RESPONSES (identical definitions to AMBER A1i-A4i):
####   1. total_richness    — all unique taxa
####   2. native_richness   — native-confirmed taxa only
####   3. prop_non_native_total  — non-native / total_richness (bounded, BINOMIAL — see note)
####   4. prop_non_native_known — non-native / (native+non-native) (bounded, BINOMIAL — see note)
####
#### FAMILY NOTE: proportion responses use BINOMIAL, not beta-binomial,
#### per the explicit diagnosed decision from the confirmed glmmTMB
#### exploration — beta-binomial's dispersion parameter repeatedly hit
#### a numerical boundary (~1e7-1e8) and a plain binomial with the same
#### random-effect structure fit just as well without the boundary
#### problem. Re-checked here (Section 12e) against the EXPANDED
#### predictor/design, not assumed to transfer automatically.
####
#### STRUCTURES:
####   BASELINE  — original nearest-connected-removal design: one row
####               per post-removal community sample, pre-removal value
####               as a covariate, pre-value x time interaction
####               (confirmed important via AIC in the earlier glmmTMB
####               work — retained here).
####   YEAR0     — continuous per-site time series including EVERY
####               community sample (not just pre/post pairs) and
####               NON-REMOVAL CONTROL BASINS. removal_status x time
####               is the key test: does removal actually change the
####               trajectory, or would recovery/decline have happened
####               anyway?
####
#### RANDOM EFFECTS (SAME for all 8 models, matching AMBER exactly):
####   (1 | country) + (1 | HYBAS_ID) + (1 | site_id) +
####   (1 | method) + (1 | sample_year)
####
#### COVARIATES (SAME base set as AMBER, PLUS DRE-specific removal terms):
####   z_dist_up, z_dist_mouth, z_log_population_density,
####   z_mean_temp_annual, z_precip_annual, z_log_river_km
####   (+ z_network_distance_km, z_n_amber_barriers_on_path,
####     network_direction, z_additional_removals_since_event —
####     BASELINE structure only, since these describe a specific
####     removal event and don't apply to unremoved control basins)
####
#### REMOVAL-HISTORY COVARIATES (NEW — capture multiple-removal exposure):
####   z_n_removals_before_first_measurement — basin-level, both
####     structures. Historical removals predating monitoring.
####   z_n_prior_removals_at_focal_event — BASELINE only. Cumulative
####     removals before the specific removal being analysed in that row.
####   z_n_removals_to_date — YEAR0 only. Cumulative removals up to
####     and including each sample (per-sample, not per-event).
####   z_years_since_last_removal — YEAR0 only. Recency of the MOST
####     RECENT removal (not the first) — neutral-zero for control/
####     pre-removal rows so those aren't dropped via NA (see Section 9c-ii).
####
#### KEY DESIGN DECISION: community-level diversity metrics (total
#### richness, native richness, both proportions) are NOT recomputed
#### here from scratch. They are loaded directly from the AMBER v2
#### prep script's saved outputs (site_diversity_base.rds and
#### site_prop_base_non_native_taxa.rds), so DRE responses are
#### GUARANTEED identical in definition to AMBER responses — this is
#### what "naturally follows the 4 AMBER models" means concretely.
#### ============================================================

#### 0. Setup ######################################################

library(sf)
library(car)        # VIF
library(MASS)        # glm.nb for dispersion/zero-inflation checks — loaded BEFORE dplyr, same masking fix as AMBER
library(DHARMa)      # zero-inflation / dispersion diagnostics
library(performance) # check_collinearity, check_singularity cross-check
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(ggplot2)

#### Guard against namespace masking (same issue hit in the AMBER script)
select <- dplyr::select
filter <- dplyr::filter

sf_use_s2(FALSE)
options(scipen = 999)

z_safe <- function(x) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) < 2 || length(unique(na.omit(x))) <= 1) return(rep(0, length(x)))
  as.numeric(scale(x))
}

amber_out_dir <- "Data/Processed/Final_Invert_Models_v2"
dre_out_dir <- "Data/Processed/DRE_Final_Models_v2"
diag_dir <- file.path(dre_out_dir, "Diagnostics")
dir.create(dre_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

record_log <- tibble(step = character(), n_rows = integer(), n_site_years = integer(), note = character())
log_step <- function(df, step_name, note = "", id_cols = c("site_id", "sample_year")) {
  n_sy <- if (length(id_cols) > 0 && all(id_cols %in% names(df))) n_distinct(df[id_cols]) else NA_integer_
  record_log <<- bind_rows(record_log, tibble(step = step_name, n_rows = nrow(df), n_site_years = n_sy, note = note))
  invisible(df)
}

#### 0a. Required upstream objects ################################
#### This script ASSUMES the following already exist in the session
#### from your earlier DRE pipeline (Sections 1-25, not shown here):
####   community_position_df  — sf, one row per site x year, with
####                             HYBAS_ID, HYRIV_ID, MAIN_RIV, country,
####                             method, line_distance_m
####   dre_position_df        — sf, one row per DRE removal, with
####                             removal_id, YearRemoved, removal_type,
####                             removal_type_reduced, HYBAS_ID,
####                             HYRIV_ID, MAIN_RIV, line_distance_m
####   amber_position_df      — sf, one row per AMBER barrier, with
####                             AMBER_ID, HYRIV_ID, amber_type_reduced,
####                             line_distance_m
####   hydrorivers_sf         — the same HydroRIVERS layer used in AMBER
#### If any of these don't exist, STOP and rebuild them first (they
#### are the spatial join/snap step, analogous to AMBER Sections 4-8).
required_objects <- c("community_position_df", "dre_position_df", "amber_position_df", "hydrorivers_sf")
missing_objects <- required_objects[!sapply(required_objects, exists)]
if (length(missing_objects) > 0) {
  stop("Missing required upstream objects: ", paste(missing_objects, collapse = ", "),
       ". These must be built by the earlier DRE spatial-join pipeline before running this script.")
}

micro_out_dir <- "Data/Processed/DRE_ClosestBarrierEvents/MicrobasinNetwork"
dir.create(micro_out_dir, recursive = TRUE, showWarnings = FALSE)


#### 1. Load AMBER's canonical community diversity metrics ########
#### This is the key fix: DRE responses now come from the EXACT SAME
#### source as AMBER responses, guaranteeing identical definitions.
amber_diversity_path <- file.path(amber_out_dir, "site_diversity_base.rds")
amber_prop_path <- file.path(amber_out_dir, "site_prop_base_non_native_taxa.rds")

if (!file.exists(amber_diversity_path) || !file.exists(amber_prop_path)) {
  stop("AMBER outputs not found at ", amber_out_dir, ". Run 01_prepare_amber_data_v2.R first — ",
       "DRE responses are built directly from its outputs to guarantee identical definitions.")
}

amber_site_diversity <- readRDS(amber_diversity_path) %>%
  mutate(site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year)))

amber_site_prop <- readRDS(amber_prop_path) %>%
  mutate(site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year)))

#### Canonical community-level response table: one row per site x year.
#### prop_non_native_total  = non_native_taxa / total_richness   (matches A3i)
#### prop_non_native_known  = non_native_taxa / total_known_taxa (matches A4i)
canonical_community_responses <- amber_site_diversity %>%
  select(site_id, sample_year, total_richness, native_richness,
         n_native_taxa, n_non_native_taxa, total_known_taxa) %>%
  left_join(
    amber_site_prop %>% select(site_id, sample_year, prop_non_native_total = prop_non_native_taxa),
    by = c("site_id", "sample_year")
  ) %>%
  mutate(
    prop_non_native_known = if_else(total_known_taxa > 0, n_non_native_taxa / total_known_taxa, NA_real_)
  )

log_step(canonical_community_responses, "01_canonical_responses_from_amber",
         "loaded directly from AMBER v2 outputs — guarantees identical response definitions")

cat("\nCanonical community responses loaded from AMBER outputs:\n")
print(head(canonical_community_responses))


#### 2. Rebuild community_samples_network with canonical responses #
community_samples_network <- community_position_df %>%
  st_drop_geometry() %>%
  mutate(
    site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year)),
    HYBAS_ID = as.character(HYBAS_ID), HYRIV_ID = as.character(HYRIV_ID),
    community_id = paste(site_id, sample_year, sep = "_")
  ) %>%
  select(community_id, site_id, sample_year, country, method, HYBAS_ID, HYRIV_ID, MAIN_RIV,
         community_line_distance_m = line_distance_m) %>%
  inner_join(canonical_community_responses, by = c("site_id", "sample_year"))

log_step(community_samples_network, "02_community_samples_network",
         "spatial data joined to canonical AMBER-consistent responses")

cat("\nCommunities lost joining spatial data to canonical responses (should be 0 or small):\n")
cat("Spatial rows:", n_distinct(community_position_df$site_id, community_position_df$sample_year), "\n")
cat("Matched rows:", nrow(community_samples_network), "\n")

dre_microbasin <- dre_position_df %>%
  st_drop_geometry() %>%
  mutate(removal_id = as.integer(removal_id), YearRemoved = as.integer(YearRemoved),
         HYBAS_ID = as.character(HYBAS_ID), HYRIV_ID = as.character(HYRIV_ID)) %>%
  filter(!is.na(YearRemoved)) %>%
  select(removal_id, YearRemoved, removal_type, removal_type_reduced,
         removal_HYBAS_ID = HYBAS_ID, removal_HYRIV_ID = HYRIV_ID,
         removal_MAIN_RIV = MAIN_RIV, removal_line_distance_m = line_distance_m)

amber_network <- amber_position_df %>%
  st_drop_geometry() %>%
  mutate(AMBER_ID = as.integer(AMBER_ID), HYRIV_ID = as.character(HYRIV_ID)) %>%
  select(AMBER_ID, HYRIV_ID, amber_type, amber_type_reduced, amber_line_distance_m = line_distance_m)


#### 3. River topology + network distance functions (reused, unchanged) ####
river_topo <- hydrorivers_sf %>%
  st_drop_geometry() %>%
  transmute(HYRIV_ID = as.character(HYRIV_ID), NEXT_DOWN = as.character(NEXT_DOWN),
            MAIN_RIV = as.character(MAIN_RIV), LENGTH_M = as.numeric(LENGTH_KM) * 1000,
            DIST_DN_KM = as.numeric(DIST_DN_KM)) %>%
  mutate(NEXT_DOWN = if_else(NEXT_DOWN %in% c("0", "", "NA", NA_character_), NA_character_, NEXT_DOWN))

next_down_vec <- setNames(river_topo$NEXT_DOWN, river_topo$HYRIV_ID)
length_m_vec <- setNames(river_topo$LENGTH_M, river_topo$HYRIV_ID)

#### Memoized downstream-path traversal. The naive version (single
#### starting node, path built as a growing vector) has two problems
#### that compound badly on real river networks:
####   1. No memoization — most tributaries converge onto the same
####      main-stem trunk, so unrelated starting segments recompute
####      the SAME downstream trunk from scratch, over and over.
####   2. `nxt %in% path` cycle-checks against the ENTIRE path built
####      so far on every step — O(path length) per step, so a path
####      of length L costs O(L^2) just for cycle detection.
#### Fix: cache the downstream path for EVERY node visited along the
#### way (not just the requested starting node), using an environment
#### for O(1) lookups both for the cache and for cycle detection. Once
#### a walk hits any previously-cached node, it reuses that node's
#### already-known path and stops immediately — so each edge in the
#### network is only ever traversed once across the WHOLE computation,
#### not once per starting segment.
path_cache_env <- new.env(parent = emptyenv())

get_downstream_path <- function(start_hyriv, max_steps = 10000) {
  start_hyriv <- as.character(start_hyriv)
  if (is.na(start_hyriv) || !start_hyriv %in% names(next_down_vec)) return(character(0))
  if (exists(start_hyriv, envir = path_cache_env, inherits = FALSE)) {
    return(get(start_hyriv, envir = path_cache_env, inherits = FALSE))
  }
  
  chain <- character(0)
  visited_env <- new.env(parent = emptyenv())
  current <- start_hyriv
  
  for (i in seq_len(max_steps)) {
    if (exists(current, envir = path_cache_env, inherits = FALSE)) {
      #### Hit a node whose downstream path is already known — reuse
      #### it and stop. This is what avoids re-walking shared trunks.
      chain <- c(chain, get(current, envir = path_cache_env, inherits = FALSE))
      break
    }
    if (exists(current, envir = visited_env, inherits = FALSE)) break  #### O(1) cycle guard
    assign(current, TRUE, envir = visited_env)
    chain <- c(chain, current)
    
    nxt <- next_down_vec[[current]]
    if (is.null(nxt) || is.na(nxt) || !nxt %in% names(next_down_vec)) break
    current <- nxt
  }
  
  #### Cache the result for EVERY node on this walk, not just the
  #### starting node — a future call starting anywhere along this
  #### same chain will hit the cache immediately.
  for (j in seq_along(chain)) {
    node <- chain[j]
    if (!exists(node, envir = path_cache_env, inherits = FALSE)) {
      assign(node, chain[j:length(chain)], envir = path_cache_env)
    }
  }
  
  chain
}

count_amber_barriers_on_path <- function(path_string, amber_df) {
  if (is.na(path_string)) return(NA_integer_)
  path_ids <- unlist(strsplit(path_string, "\\|"))
  amber_df %>% filter(HYRIV_ID %in% path_ids) %>% summarise(n = n_distinct(AMBER_ID)) %>% pull(n)
}

get_network_connection_segments_cached <- function(community_hyriv, removal_hyriv, path_cache) {
  community_hyriv <- as.character(community_hyriv); removal_hyriv <- as.character(removal_hyriv)
  if (is.na(community_hyriv) || is.na(removal_hyriv) ||
      !community_hyriv %in% names(length_m_vec) || !removal_hyriv %in% names(length_m_vec)) {
    return(tibble(network_direction = NA_character_, segment_path_distance_m = NA_real_,
                  segment_path_distance_km = NA_real_, n_segment_steps = NA_integer_, path_hyriv_ids = NA_character_))
  }
  if (community_hyriv == removal_hyriv) {
    return(tibble(network_direction = "same_segment", segment_path_distance_m = 0,
                  segment_path_distance_km = 0, n_segment_steps = 0L, path_hyriv_ids = community_hyriv))
  }
  community_down_path <- if (community_hyriv %in% names(path_cache)) path_cache[[community_hyriv]] else character(0)
  removal_down_path <- if (removal_hyriv %in% names(path_cache)) path_cache[[removal_hyriv]] else character(0)
  
  if (removal_hyriv %in% community_down_path) {
    idx <- match(removal_hyriv, community_down_path); path_ids <- community_down_path[seq_len(idx)]
    d_m <- sum(length_m_vec[path_ids], na.rm = TRUE)
    return(tibble(network_direction = "downstream", segment_path_distance_m = d_m,
                  segment_path_distance_km = d_m / 1000, n_segment_steps = length(path_ids) - 1L,
                  path_hyriv_ids = paste(path_ids, collapse = "|")))
  }
  if (community_hyriv %in% removal_down_path) {
    idx <- match(community_hyriv, removal_down_path); path_ids <- removal_down_path[seq_len(idx)]
    d_m <- sum(length_m_vec[path_ids], na.rm = TRUE)
    return(tibble(network_direction = "upstream", segment_path_distance_m = d_m,
                  segment_path_distance_km = d_m / 1000, n_segment_steps = length(path_ids) - 1L,
                  path_hyriv_ids = paste(rev(path_ids), collapse = "|")))
  }
  common_ids <- intersect(community_down_path, removal_down_path)
  if (length(common_ids) > 0) {
    first_common <- common_ids[1]
    community_idx <- match(first_common, community_down_path); removal_idx <- match(first_common, removal_down_path)
    community_to_common <- community_down_path[seq_len(community_idx)]; removal_to_common <- removal_down_path[seq_len(removal_idx)]
    path_ids <- unique(c(community_to_common, removal_to_common))
    d_m <- sum(length_m_vec[path_ids], na.rm = TRUE)
    return(tibble(network_direction = "other_branch_same_microbasin", segment_path_distance_m = d_m,
                  segment_path_distance_km = d_m / 1000,
                  n_segment_steps = length(community_to_common) + length(removal_to_common) - 2L,
                  path_hyriv_ids = paste(path_ids, collapse = "|")))
  }
  tibble(network_direction = "not_connected_in_hydrorivers", segment_path_distance_m = NA_real_,
         segment_path_distance_km = NA_real_, n_segment_steps = NA_integer_, path_hyriv_ids = NA_character_)
}


#### 4. Bring in AMBER's site-level network-position covariates ####
#### (z_dist_up, z_dist_mouth, climate, population density, river km)
#### DRE sites are a subset of the same invertebrate dataset AMBER
#### uses, so these join directly by site_id + sample_year / HYBAS_ID.
site_river_position_path <- file.path(amber_out_dir, "..", "site_river_position.rds")  # adjust if saved elsewhere
amber_a1i <- readRDS(file.path(amber_out_dir, "A1i_TotalRichnesswithtype_data.rds")) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  distinct(site_id, sample_year, .keep_all = TRUE) %>%
  transmute(
    site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year)),
    dist_up_source = z_dist_up, dist_mouth_source = z_dist_mouth,
    log_population_density_source = z_log_population_density,
    mean_temp_annual_source = z_mean_temp_annual, precip_annual_source = z_precip_annual,
    log_river_km_source = z_log_river_km
  )

cat("\nNOTE: the AMBER covariate join above pulls ALREADY Z-SCORED AMBER\n")
cat("columns as a shortcut. This is fine for correlation/QA checks, but for\n")
cat("the FINAL DRE model data below, raw (unscaled) versions are re-derived\n")
cat("and re-scored within the DRE sample so z-scores reflect the DRE\n")
cat("sample's own mean/SD, not AMBER's. If you'd rather share AMBER's exact\n")
cat("scaling, swap the z_safe() calls in Sections 7-8 for the *_source columns above.\n")


#### ============================================================
#### 5. MICROBASIN REMOVAL EXPOSURE — EFFECT-SIZE EVENTS
#### (adapted from the original Section 26-27; unchanged logic,
#### now operating on community_samples_network with canonical
#### AMBER-consistent responses)
#### ============================================================

community_intervals_micro <- community_samples_network %>%
  arrange(site_id, sample_year) %>%
  group_by(site_id) %>%
  mutate(
    pre_community_id = lag(community_id), pre_year = lag(sample_year),
    pre_HYBAS_ID = lag(HYBAS_ID), pre_HYRIV_ID = lag(HYRIV_ID), pre_MAIN_RIV = lag(MAIN_RIV),
    pre_line_distance_m = lag(community_line_distance_m),
    pre_total_richness = lag(total_richness), pre_native_richness = lag(native_richness),
    pre_prop_non_native_total = lag(prop_non_native_total), pre_prop_non_native_known = lag(prop_non_native_known),
    pre_total_known_taxa = lag(total_known_taxa), pre_n_non_native_taxa = lag(n_non_native_taxa),
    pre_n_native_taxa = lag(n_native_taxa),
    post_community_id = community_id, post_year = sample_year,
    post_HYBAS_ID = HYBAS_ID, post_HYRIV_ID = HYRIV_ID, post_MAIN_RIV = MAIN_RIV,
    post_line_distance_m = community_line_distance_m,
    post_total_richness = total_richness, post_native_richness = native_richness,
    post_prop_non_native_total = prop_non_native_total, post_prop_non_native_known = prop_non_native_known,
    post_total_known_taxa = total_known_taxa, post_n_non_native_taxa = n_non_native_taxa,
    post_n_native_taxa = n_native_taxa
  ) %>%
  ungroup() %>%
  filter(!is.na(pre_year), post_year > pre_year, pre_HYBAS_ID == post_HYBAS_ID) %>%
  select(site_id, country, method, HYBAS_ID = post_HYBAS_ID,
         pre_community_id, post_community_id, pre_year, post_year,
         pre_HYRIV_ID, post_HYRIV_ID, pre_MAIN_RIV, post_MAIN_RIV,
         pre_line_distance_m, post_line_distance_m,
         pre_total_richness, post_total_richness, pre_native_richness, post_native_richness,
         pre_prop_non_native_total, post_prop_non_native_total,
         pre_prop_non_native_known, post_prop_non_native_known,
         pre_total_known_taxa, post_total_known_taxa,
         pre_n_non_native_taxa, post_n_non_native_taxa, pre_n_native_taxa, post_n_native_taxa)

log_step(community_intervals_micro, "05_community_intervals", "consecutive-sample intervals within same basin")


#### ============================================================
#### 6. BASELINE STRUCTURE: nearest-connected-removal recovery dataset
#### (adapted from Sections 28-31 of the original — same core logic:
#### find last pre-removal sample, all post-removal samples, network
#### distance, de-duplicate to nearest connected removal per
#### post-observation)
#### ============================================================

site_removal_candidates <- community_samples_network %>%
  distinct(site_id, country, method, HYBAS_ID) %>%
  inner_join(dre_microbasin, by = c("HYBAS_ID" = "removal_HYBAS_ID"), relationship = "many-to-many")

pre_baselines <- site_removal_candidates %>%
  left_join(community_samples_network, by = c("site_id", "HYBAS_ID"),
            relationship = "many-to-many", suffix = c("_removal", "_sample")) %>%
  filter(sample_year < YearRemoved) %>%
  group_by(site_id, removal_id, YearRemoved) %>%
  slice_max(order_by = sample_year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    site_id, country = country_sample, method = method_sample, HYBAS_ID,
    removal_id, YearRemoved, removal_type, removal_type_reduced,
    removal_HYRIV_ID, removal_MAIN_RIV, removal_line_distance_m,
    pre_community_id = community_id, pre_year = sample_year, pre_HYRIV_ID = HYRIV_ID,
    pre_MAIN_RIV = MAIN_RIV, pre_line_distance_m = community_line_distance_m,
    pre_total_richness = total_richness, pre_native_richness = native_richness,
    pre_prop_non_native_total = prop_non_native_total, pre_prop_non_native_known = prop_non_native_known,
    pre_total_known_taxa = total_known_taxa, pre_n_non_native_taxa = n_non_native_taxa,
    pre_n_native_taxa = n_native_taxa,
    years_pre_to_removal = YearRemoved - sample_year
  )

log_step(pre_baselines, "06a_pre_baselines", "last pre-removal community sample per candidate removal")

post_samples <- pre_baselines %>%
  #### FIX: community_samples_network ALSO has country/method columns,
  #### which would otherwise collide with pre_baselines' own
  #### country/method (aliased from country_sample/method_sample
  #### earlier) and get silently renamed to country.x/country.y,
  #### method.x/method.y by dplyr — leaving no plain "country"/"method"
  #### column downstream. Drop the duplicates from the RIGHT side of
  #### the join; pre_baselines' values are already correct (same site,
  #### so same country/method regardless of which side they came from).
  left_join(community_samples_network %>% select(-country, -method),
            by = c("site_id", "HYBAS_ID"), relationship = "many-to-many") %>%
  filter(sample_year > YearRemoved) %>%
  mutate(
    years_since_removal = sample_year - YearRemoved,
    log_years_since_removal = log1p(years_since_removal),
    removal_event_id = paste(site_id, removal_id, YearRemoved, sep = "_")
  )

additional_removal_counts <- post_samples %>%
  select(removal_event_id, site_id, HYBAS_ID, focal_removal_id = removal_id, YearRemoved, sample_year) %>%
  inner_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, other_removal_id = removal_id,
                                       other_YearRemoved = YearRemoved),
             by = "HYBAS_ID", relationship = "many-to-many") %>%
  filter(other_removal_id != focal_removal_id, other_YearRemoved > YearRemoved, other_YearRemoved <= sample_year) %>%
  group_by(removal_event_id, sample_year) %>%
  summarise(additional_removals_since_event = n_distinct(other_removal_id), .groups = "drop")

baseline_trajectories <- post_samples %>%
  left_join(additional_removal_counts, by = c("removal_event_id", "sample_year")) %>%
  mutate(additional_removals_since_event = replace_na(additional_removals_since_event, 0L))

log_step(baseline_trajectories, "06b_baseline_trajectories", "all post-removal samples paired with candidate removals")


#### Network distance/direction — cached, computed once per unique segment pair
unique_pairs <- baseline_trajectories %>% distinct(pre_HYRIV_ID, removal_HYRIV_ID)
unique_hyriv_needed <- unique(c(unique_pairs$pre_HYRIV_ID, unique_pairs$removal_HYRIV_ID))
unique_hyriv_needed <- unique_hyriv_needed[!is.na(unique_hyriv_needed) & unique_hyriv_needed %in% names(next_down_vec)]

cat("\nBuilding downstream-path cache for", length(unique_hyriv_needed), "unique HYRIV_IDs...\n")
cat("(now memoized — later calls reuse shared downstream trunks, so this should\n")
cat("speed up noticeably as it progresses, not stay at a constant pace)\n")
downstream_path_cache <- vector("list", length(unique_hyriv_needed))
names(downstream_path_cache) <- unique_hyriv_needed
for (idx in seq_along(unique_hyriv_needed)) {
  downstream_path_cache[[idx]] <- get_downstream_path(unique_hyriv_needed[idx])
  if (idx %% 2000 == 0) {
    cat("  ...", idx, "of", length(unique_hyriv_needed), "processed (",
        round(100 * idx / length(unique_hyriv_needed), 1), "%), cache size:",
        length(ls(path_cache_env)), "nodes\n")
  }
}
cat("Finished. Total unique nodes cached:", length(ls(path_cache_env)), "\n")

network_stats_pairs <- unique_pairs %>%
  mutate(network_info = pmap(list(pre_HYRIV_ID, removal_HYRIV_ID),
                             ~ get_network_connection_segments_cached(..1, ..2, downstream_path_cache))) %>%
  tidyr::unnest(network_info)

baseline_trajectories <- baseline_trajectories %>%
  left_join(network_stats_pairs, by = c("pre_HYRIV_ID", "removal_HYRIV_ID")) %>%
  mutate(
    same_segment_distance_m = if_else(pre_HYRIV_ID == removal_HYRIV_ID,
                                      abs(removal_line_distance_m - pre_line_distance_m), NA_real_),
    network_distance_m = case_when(pre_HYRIV_ID == removal_HYRIV_ID ~ same_segment_distance_m,
                                   TRUE ~ segment_path_distance_m),
    network_distance_km = network_distance_m / 1000
  )

path_barrier_counts <- baseline_trajectories %>%
  distinct(path_hyriv_ids) %>%
  mutate(n_amber_barriers_on_path = map_int(path_hyriv_ids, count_amber_barriers_on_path, amber_df = amber_network))

baseline_trajectories <- baseline_trajectories %>%
  left_join(path_barrier_counts, by = "path_hyriv_ids") %>%
  mutate(network_direction = factor(network_direction,
                                    levels = c("same_segment", "upstream", "downstream", "other_branch_same_microbasin", "not_connected_in_hydrorivers")))

#### Restrict to HydroRIVERS-connected pairs (excludes not_connected_in_hydrorivers)
#### per the original design decision.
baseline_connected <- baseline_trajectories %>%
  filter(network_direction != "not_connected_in_hydrorivers", !is.na(network_distance_km)) %>%
  droplevels()

log_step(baseline_connected, "06c_baseline_connected", "restricted to HydroRIVERS-connected community-removal pairs")

#### De-duplicate to nearest connected removal per post-observation
#### (a post-removal sample can be a valid candidate for >1 removal;
#### keep only the nearest one to avoid pseudo-replication).
baseline_connected <- baseline_connected %>%
  mutate(post_observation_id = paste(site_id, sample_year, HYBAS_ID, sep = "_"))

baseline_nearest <- baseline_connected %>%
  group_by(post_observation_id) %>%
  slice_min(order_by = network_distance_km, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  droplevels()

log_step(baseline_nearest, "06d_baseline_nearest", "de-duplicated to nearest connected removal per post-observation")

cat("\nBaseline structure — network direction counts:\n")
print(baseline_nearest %>% count(network_direction, sort = TRUE))


#### 6e. Removal-history fixed effects ##############################
#### Two NEW covariates, both requested explicitly:
####
#### n_removals_before_first_measurement — BASIN-level, time-invariant.
#### How many DRE removals happened in this basin BEFORE monitoring
#### ever started there. Matters because even the earliest "baseline"
#### sample in a basin may already postdate historical removals your
#### before/after design can't see as a pair — this is a confound on
#### how "pristine" any given baseline really is.
####
#### n_prior_removals_at_focal_event — PER-ROW. How many OTHER removals
#### had already happened in this basin before the FOCAL removal being
#### analysed in this row. This is the baseline-structure equivalent of
#### "number of removals overall" — cumulative disturbance history up
#### to the point of the event actually being modelled, as distinct
#### from z_additional_removals_since_event (which only counts removals
#### AFTER the focal one, already in the dataset).
first_sample_year_per_basin <- community_samples_network %>%
  group_by(HYBAS_ID) %>%
  summarise(first_sample_year = min(sample_year), .groups = "drop")

n_removals_before_first_measurement_df <- first_sample_year_per_basin %>%
  left_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, removal_id, YearRemoved),
            by = "HYBAS_ID", relationship = "many-to-many") %>%
  group_by(HYBAS_ID) %>%
  summarise(n_removals_before_first_measurement = sum(YearRemoved < first(first_sample_year), na.rm = TRUE),
            .groups = "drop")

cat("\nHistorical removal exposure (before monitoring started) by basin:\n")
print(summary(n_removals_before_first_measurement_df$n_removals_before_first_measurement))
cat("Basins with >=1 pre-monitoring removal:",
    sum(n_removals_before_first_measurement_df$n_removals_before_first_measurement > 0), "of",
    nrow(n_removals_before_first_measurement_df), "\n")

baseline_nearest <- baseline_nearest %>%
  left_join(n_removals_before_first_measurement_df, by = "HYBAS_ID") %>%
  mutate(n_removals_before_first_measurement = replace_na(n_removals_before_first_measurement, 0L))

#### Cumulative prior removals at the focal event (per row)
n_prior_removals_at_focal_df <- baseline_nearest %>%
  distinct(HYBAS_ID, removal_id, YearRemoved) %>%
  left_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, other_removal_id = removal_id,
                                      other_YearRemoved = YearRemoved),
            by = "HYBAS_ID", relationship = "many-to-many") %>%
  filter(other_removal_id != removal_id, other_YearRemoved < YearRemoved) %>%
  group_by(HYBAS_ID, removal_id, YearRemoved) %>%
  summarise(n_prior_removals_at_focal_event = n_distinct(other_removal_id), .groups = "drop")

baseline_nearest <- baseline_nearest %>%
  left_join(n_prior_removals_at_focal_df, by = c("HYBAS_ID", "removal_id", "YearRemoved")) %>%
  mutate(n_prior_removals_at_focal_event = replace_na(n_prior_removals_at_focal_event, 0L))

cat("\nCumulative prior removals at the focal event (per row) summary:\n")
print(summary(baseline_nearest$n_prior_removals_at_focal_event))

log_step(baseline_nearest, "06e_baseline_removal_history", "removal-history fixed effects added")


#### ============================================================
#### 7. Join AMBER-style covariates (dist_up, dist_mouth, climate,
#### population density, river length) onto the baseline dataset,
#### matching the AMBER predictor set exactly.
#### ============================================================
post_samples <- pre_baselines %>%
  left_join(community_samples_network %>% select(-country, -method),
            by = c("site_id", "HYBAS_ID"), relationship = "many-to-many") %>%
  filter(sample_year > YearRemoved) %>%
  mutate(
    years_since_removal = sample_year - YearRemoved,
    log_years_since_removal = log1p(years_since_removal),
    removal_event_id = paste(site_id, removal_id, YearRemoved, sep = "_")
  )

additional_removal_counts <- post_samples %>%
  select(removal_event_id, site_id, HYBAS_ID, focal_removal_id = removal_id, YearRemoved, sample_year) %>%
  inner_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, other_removal_id = removal_id,
                                       other_YearRemoved = YearRemoved),
             by = "HYBAS_ID", relationship = "many-to-many") %>%
  filter(other_removal_id != focal_removal_id, other_YearRemoved > YearRemoved, other_YearRemoved <= sample_year) %>%
  group_by(removal_event_id, sample_year) %>%
  summarise(additional_removals_since_event = n_distinct(other_removal_id), .groups = "drop")

baseline_trajectories <- post_samples %>%
  left_join(additional_removal_counts, by = c("removal_event_id", "sample_year")) %>%
  mutate(additional_removals_since_event = replace_na(additional_removals_since_event, 0L))

#### unique_pairs is cheap to rebuild — do NOT rebuild downstream_path_cache,
#### it's unaffected by the country/method fix and you already have it.
unique_pairs <- baseline_trajectories %>% distinct(pre_HYRIV_ID, removal_HYRIV_ID)

network_stats_pairs <- unique_pairs %>%
  mutate(network_info = pmap(list(pre_HYRIV_ID, removal_HYRIV_ID),
                             ~ get_network_connection_segments_cached(..1, ..2, downstream_path_cache))) %>%
  tidyr::unnest(network_info)

baseline_trajectories <- baseline_trajectories %>%
  left_join(network_stats_pairs, by = c("pre_HYRIV_ID", "removal_HYRIV_ID")) %>%
  mutate(
    same_segment_distance_m = if_else(pre_HYRIV_ID == removal_HYRIV_ID,
                                      abs(removal_line_distance_m - pre_line_distance_m), NA_real_),
    network_distance_m = case_when(pre_HYRIV_ID == removal_HYRIV_ID ~ same_segment_distance_m,
                                   TRUE ~ segment_path_distance_m),
    network_distance_km = network_distance_m / 1000
  )

path_barrier_counts <- baseline_trajectories %>%
  distinct(path_hyriv_ids) %>%
  mutate(n_amber_barriers_on_path = map_int(path_hyriv_ids, count_amber_barriers_on_path, amber_df = amber_network))

baseline_trajectories <- baseline_trajectories %>%
  left_join(path_barrier_counts, by = "path_hyriv_ids") %>%
  mutate(network_direction = factor(network_direction,
                                    levels = c("same_segment", "upstream", "downstream", "other_branch_same_microbasin", "not_connected_in_hydrorivers")))

baseline_connected <- baseline_trajectories %>%
  filter(network_direction != "not_connected_in_hydrorivers", !is.na(network_distance_km)) %>%
  droplevels() %>%
  mutate(post_observation_id = paste(site_id, sample_year, HYBAS_ID, sep = "_"))

baseline_nearest <- baseline_connected %>%
  group_by(post_observation_id) %>%
  slice_min(order_by = network_distance_km, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  droplevels()

first_sample_year_per_basin <- community_samples_network %>%
  group_by(HYBAS_ID) %>%
  summarise(first_sample_year = min(sample_year), .groups = "drop")

n_removals_before_first_measurement_df <- first_sample_year_per_basin %>%
  left_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, removal_id, YearRemoved),
            by = "HYBAS_ID", relationship = "many-to-many") %>%
  group_by(HYBAS_ID) %>%
  summarise(n_removals_before_first_measurement = sum(YearRemoved < first(first_sample_year), na.rm = TRUE),
            .groups = "drop")

baseline_nearest <- baseline_nearest %>%
  left_join(n_removals_before_first_measurement_df, by = "HYBAS_ID") %>%
  mutate(n_removals_before_first_measurement = replace_na(n_removals_before_first_measurement, 0L))

n_prior_removals_at_focal_df <- baseline_nearest %>%
  distinct(HYBAS_ID, removal_id, YearRemoved) %>%
  left_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, other_removal_id = removal_id,
                                      other_YearRemoved = YearRemoved),
            by = "HYBAS_ID", relationship = "many-to-many") %>%
  filter(other_removal_id != removal_id, other_YearRemoved < YearRemoved) %>%
  group_by(HYBAS_ID, removal_id, YearRemoved) %>%
  summarise(n_prior_removals_at_focal_event = n_distinct(other_removal_id), .groups = "drop")

baseline_nearest <- baseline_nearest %>%
  left_join(n_prior_removals_at_focal_df, by = c("HYBAS_ID", "removal_id", "YearRemoved")) %>%
  mutate(n_prior_removals_at_focal_event = replace_na(n_prior_removals_at_focal_event, 0L))

baseline_nearest <- baseline_nearest %>%
  join_amber_covariates() %>%
  mutate(
    z_dist_up = z_safe(log_dist_up), z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual), z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    z_network_distance_km = z_safe(log1p(network_distance_km)),
    z_n_amber_barriers_on_path = z_safe(log1p(n_amber_barriers_on_path)),
    z_additional_removals_since_event = z_safe(log1p(additional_removals_since_event)),
    z_years_pre_to_removal = z_safe(years_pre_to_removal),
    z_n_removals_before_first_measurement = z_safe(log1p(n_removals_before_first_measurement)),
    z_n_prior_removals_at_focal_event = z_safe(log1p(n_prior_removals_at_focal_event)),
    site_id = factor(site_id), HYBAS_ID = factor(HYBAS_ID), country = factor(country),
    method = factor(method), sample_year = factor(sample_year)
  )

common_baseline_cols <- c("site_id", "HYBAS_ID", "country", "method", "sample_year",
                          "years_since_removal", "log_years_since_removal", "years_pre_to_removal",
                          "z_years_pre_to_removal", "network_direction", "z_network_distance_km",
                          "z_n_amber_barriers_on_path", "z_additional_removals_since_event",
                          "z_n_removals_before_first_measurement", "z_n_prior_removals_at_focal_event",
                          "z_dist_up", "z_dist_mouth", "z_log_population_density",
                          "z_mean_temp_annual", "z_precip_annual", "z_log_river_km",
                          "removal_event_id", "post_observation_id")

#### These objects (climate_wide, hyde_full, hyde_1961, river_len_df,
#### site_river_position) were built during AMBER prep. Load them if
#### not already in this session.
if (!exists("climate_wide")) climate_wide <- readRDS(file.path(amber_out_dir, "..", "climate_wide.rds"))
if (!exists("river_len_df")) river_len_df <- readRDS(file.path(amber_out_dir, "..", "river_len_df.rds"))
if (!exists("hyde_full")) hyde_full <- readRDS(file.path(amber_out_dir, "..", "hyde_full.rds"))
if (!exists("hyde_1961")) hyde_1961 <- readRDS(file.path(amber_out_dir, "..", "hyde_1961.rds"))
if (!exists("site_river_position")) site_river_position <- readRDS(file.path(amber_out_dir, "..", "site_river_position.rds"))

cat("\nNOTE: if the four readRDS() calls above fail with 'file not found',\n")
cat("add explicit saveRDS() calls for climate_wide, river_len_df, hyde_full,\n")
cat("hyde_1961, and site_river_position at the end of 01_prepare_amber_data_v2.R\n")
cat("(they're built there but weren't saved to disk in the current version) —\n")
cat("or re-paste those five objects into this session before running Section 7.\n")

join_amber_covariates <- function(df) {
  df %>%
    mutate(site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year)),
           HYBAS_ID = as.character(HYBAS_ID)) %>%
    left_join(site_river_position %>%
                mutate(site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year))),
              by = c("site_id", "sample_year")) %>%
    left_join(climate_wide %>% mutate(HYBAS_ID = as.character(HYBAS_ID), sample_year = as.integer(as.character(sample_year))),
              by = c("HYBAS_ID", "sample_year")) %>%
    left_join(river_len_df %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
    left_join(hyde_full %>% mutate(HYBAS_ID = as.character(HYBAS_ID), year = as.integer(year)) %>%
                select(HYBAS_ID, sample_year = year, log_population_density),
              by = c("HYBAS_ID", "sample_year")) %>%
    left_join(hyde_1961 %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
    mutate(log_population_density = if_else(is.na(log_population_density) & sample_year < 1961,
                                            log_population_density_1961, log_population_density)) %>%
    select(-log_population_density_1961)
}

exists("hyde_full")
names(hyde_full)
"log_population_density" %in% names(hyde_full)
head(hyde_full)

hyde_full <- hyde_pop_density %>%
  rename(mean_population_density = population_density) %>%
  mutate(log_population_density = log1p(mean_population_density))

hyde_1961 <- hyde_full %>%
  filter(year == 1961) %>%
  select(HYBAS_ID, log_population_density_1961 = log_population_density)

"log_population_density" %in% names(baseline_nearest)

baseline_nearest <- baseline_nearest %>%
  select(-any_of(c(
    "HYRIV_ID", "dist_up_km", "dist_down_km", "dist_mouth_km",
    "log_dist_up", "log_dist_down", "log_dist_mouth",
    "mean_temp_annual", "precip_annual",
    "river_km", "log_river_km",
    "log_population_density", "log_population_density_1961"
  )))

#### Now safe to re-run the join
baseline_nearest <- baseline_nearest %>%
  join_amber_covariates() %>%
  mutate(
    z_dist_up = z_safe(log_dist_up), z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual), z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    z_network_distance_km = z_safe(log1p(network_distance_km)),
    z_n_amber_barriers_on_path = z_safe(log1p(n_amber_barriers_on_path)),
    z_additional_removals_since_event = z_safe(log1p(additional_removals_since_event)),
    z_years_pre_to_removal = z_safe(years_pre_to_removal),
    z_n_removals_before_first_measurement = z_safe(log1p(n_removals_before_first_measurement)),
    z_n_prior_removals_at_focal_event = z_safe(log1p(n_prior_removals_at_focal_event)),
    site_id = factor(site_id), HYBAS_ID = factor(HYBAS_ID), country = factor(country),
    method = factor(method), sample_year = factor(sample_year)
  )

c("country", "method", "log_population_density", "z_dist_up", "z_dist_mouth") %in% names(baseline_nearest)

log_step(baseline_nearest, "07_baseline_with_covariates", "AMBER-style covariates joined")


#### ============================================================
#### 8. Build 4 BASELINE model-ready datasets (one per response)
#### ============================================================

common_baseline_cols <- c("site_id", "HYBAS_ID", "country", "method", "sample_year",
                          "years_since_removal", "log_years_since_removal", "years_pre_to_removal",
                          "z_years_pre_to_removal", "network_direction", "z_network_distance_km",
                          "z_n_amber_barriers_on_path", "z_additional_removals_since_event",
                          "z_n_removals_before_first_measurement", "z_n_prior_removals_at_focal_event",
                          "z_dist_up", "z_dist_mouth", "z_log_population_density",
                          "z_mean_temp_annual", "z_precip_annual", "z_log_river_km",
                          "removal_event_id", "post_observation_id")

DRE1i_TotalRichnessBaseline_data <- baseline_nearest %>%
  filter(!is.na(total_richness), !is.na(pre_total_richness)) %>%
  transmute(value = as.integer(total_richness), pre_value = as.numeric(pre_total_richness),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

DRE2i_NativeRichnessBaseline_data <- baseline_nearest %>%
  filter(!is.na(native_richness), !is.na(pre_native_richness)) %>%
  transmute(value = as.integer(native_richness), pre_value = as.numeric(pre_native_richness),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

DRE3i_PropTotalBaseline_data <- baseline_nearest %>%
  filter(!is.na(prop_non_native_total), !is.na(pre_prop_non_native_total), total_richness > 0) %>%
  transmute(non_native = as.integer(n_non_native_taxa), total_taxa = as.integer(total_richness),
            prop_non_native = as.numeric(prop_non_native_total),
            pre_value = as.numeric(pre_prop_non_native_total),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

DRE4i_PropNativeBaseline_data <- baseline_nearest %>%
  filter(!is.na(prop_non_native_known), !is.na(pre_prop_non_native_known), total_known_taxa > 0) %>%
  transmute(non_native = as.integer(n_non_native_taxa), total_taxa = as.integer(total_known_taxa),
            prop_non_native = as.numeric(prop_non_native_known),
            pre_value = as.numeric(pre_prop_non_native_known),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

log_step(DRE1i_TotalRichnessBaseline_data, "08_DRE1i_baseline", "final baseline total richness rows")
log_step(DRE2i_NativeRichnessBaseline_data, "08_DRE2i_baseline", "final baseline native richness rows")
log_step(DRE3i_PropTotalBaseline_data, "08_DRE3i_baseline", "final baseline proportion (total denom) rows")
log_step(DRE4i_PropNativeBaseline_data, "08_DRE4i_baseline", "final baseline proportion (known denom) rows")


#### ============================================================
#### 9. YEAR-0 / CONTROL-BASIN STRUCTURE
#### ============================================================
#### Every community sample becomes one point on a continuous time
#### series. removal_status distinguishes:
####   "removed"      — site in a basin with >=1 DRE removal; time
####                    axis = years relative to the FOCAL removal
####                    (negative = pre, positive = post)
####   "control"      — site in a basin with ZERO DRE removals across
####                    the whole time series; time axis = years
####                    relative to that site's FIRST sample
####
#### For "removed" sites with multiple removals, the FIRST removal in
#### the basin is used as the reference point (simplest, most
#### defensible choice — using the nearest-in-time removal to each
#### sample would double-count information already captured by
#### additional_removals_since_event).
#### ============================================================

#### 9a. Identify control (never-removed) basins ###################
removed_basin_ids <- unique(dre_microbasin$removal_HYBAS_ID)

control_basin_check <- community_samples_network %>%
  distinct(HYBAS_ID) %>%
  mutate(is_control_basin = !(HYBAS_ID %in% removed_basin_ids))

n_control_basins <- sum(control_basin_check$is_control_basin)
n_removed_basins <- sum(!control_basin_check$is_control_basin)

cat("\n============================================================\n")
cat("CONTROL BASIN AVAILABILITY CHECK (was outstanding in the handoff — now run)\n")
cat("============================================================\n")
cat("Total basins with community data:", nrow(control_basin_check), "\n")
cat("Basins with >=1 DRE removal (used for removal reference):", n_removed_basins, "\n")
cat("Basins with ZERO DRE removals (control candidates):", n_control_basins, "\n")

#### Restrict control basins to ones with usable data (non-zero
#### total_known_taxa, i.e. the proportion responses are computable) —
#### this was drafted but never run in the original handoff.
control_basin_usability <- community_samples_network %>%
  filter(HYBAS_ID %in% control_basin_check$HYBAS_ID[control_basin_check$is_control_basin]) %>%
  group_by(HYBAS_ID) %>%
  summarise(n_samples = n(), n_usable_prop_samples = sum(total_known_taxa > 0, na.rm = TRUE), .groups = "drop")

cat("\nControl basins with >=1 usable proportion-model sample:",
    sum(control_basin_usability$n_usable_prop_samples > 0), "of", nrow(control_basin_usability), "\n")
cat("Control basins with >=2 samples (needed for any within-basin comparison):",
    sum(control_basin_usability$n_samples >= 2), "\n")

write_csv(control_basin_usability, file.path(diag_dir, "QC_control_basin_usability.csv"))

if (n_control_basins < 10) {
  warning("Fewer than 10 control basins available. The removal_status x time interaction ",
          "in the Year-0 models will be poorly identified — treat those results cautiously ",
          "and consider the baseline structure as primary.")
}

#### 9b. Site sample-count check (also outstanding in the handoff) ####
site_sample_counts <- community_samples_network %>%
  mutate(has_removal_in_basin = HYBAS_ID %in% removed_basin_ids) %>%
  group_by(site_id, has_removal_in_basin) %>%
  summarise(n_samples = n_distinct(sample_year), .groups = "drop")

cat("\nSamples per site, split by removal vs control basin:\n")
print(site_sample_counts %>% group_by(has_removal_in_basin) %>%
        summarise(median_samples = median(n_samples), mean_samples = round(mean(n_samples), 2),
                  max_samples = max(n_samples), n_sites = n()))

write_csv(site_sample_counts, file.path(diag_dir, "QC_site_sample_counts_removal_vs_control.csv"))

cat("\nDECISION based on the above: if median samples per site is <= 3-4,\n")
cat("use a LINEAR time term (log_years) rather than a spline in the Year-0\n")
cat("models below — a spline needs more repeated visits per site than a\n")
cat("linear term to be identifiable. Check the printed medians before\n")
cat("choosing df in Section 10's formula.\n")


#### 9c. Build the continuous Year-0 time series #####################
#### Reference removal per basin = FIRST removal chronologically.
first_removal_per_basin <- dre_microbasin %>%
  group_by(removal_HYBAS_ID) %>%
  slice_min(order_by = YearRemoved, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(HYBAS_ID = removal_HYBAS_ID, reference_removal_id = removal_id, reference_removal_year = YearRemoved)

year0_data <- community_samples_network %>%
  left_join(first_removal_per_basin, by = "HYBAS_ID") %>%
  group_by(site_id) %>%
  mutate(first_sample_year = min(sample_year)) %>%
  ungroup() %>%
  mutate(
    removal_status = if_else(!is.na(reference_removal_id), "removed", "control"),
    reference_year = if_else(removal_status == "removed", reference_removal_year, first_sample_year),
    years_relative_to_reference = sample_year - reference_year,
    log_years_relative = sign(years_relative_to_reference) * log1p(abs(years_relative_to_reference))
  )

#### Additional removals since the reference point (removed sites only;
#### 0 for control sites by construction).
year0_additional_removals <- year0_data %>%
  filter(removal_status == "removed") %>%
  select(site_id, HYBAS_ID, sample_year, reference_removal_id, reference_removal_year) %>%
  inner_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, other_removal_id = removal_id,
                                       other_YearRemoved = YearRemoved),
             by = "HYBAS_ID", relationship = "many-to-many") %>%
  filter(other_removal_id != reference_removal_id, other_YearRemoved <= sample_year,
         other_YearRemoved > reference_removal_year) %>%
  group_by(site_id, sample_year) %>%
  summarise(additional_removals_since_reference = n_distinct(other_removal_id), .groups = "drop")

year0_data <- year0_data %>%
  left_join(year0_additional_removals, by = c("site_id", "sample_year")) %>%
  mutate(additional_removals_since_reference = replace_na(additional_removals_since_reference, 0L))

log_step(year0_data, "09c_year0_data", "continuous time series, removal + control basins combined")

cat("\nYear-0 dataset — removal_status counts:\n")
print(year0_data %>% count(removal_status))


#### 9c-ii. Cumulative removals-to-date and years-since-last-removal ####
#### Two NEW covariates, both explicitly requested: "how many removals
#### overall" and "years since the last removal" — the latter uses the
#### MOST RECENT prior removal, not just the first, so it captures
#### multiple-removal exposure properly rather than treating a basin's
#### disturbance history as a single event.
cumulative_removal_stats <- year0_data %>%
  filter(removal_status == "removed") %>%
  distinct(site_id, HYBAS_ID, sample_year) %>%
  left_join(dre_microbasin %>% select(HYBAS_ID = removal_HYBAS_ID, removal_id, YearRemoved),
            by = "HYBAS_ID", relationship = "many-to-many") %>%
  filter(YearRemoved <= sample_year) %>%
  group_by(site_id, sample_year) %>%
  summarise(n_removals_to_date = n_distinct(removal_id),
            years_since_last_removal = first(sample_year) - max(YearRemoved), .groups = "drop")


year0_data <- year0_data %>%
  left_join(cumulative_removal_stats, by = c("site_id", "sample_year")) %>%
  mutate(
    n_removals_to_date = replace_na(n_removals_to_date, 0L),
    #### years_since_last_removal is genuinely undefined for control
    #### sites and for removed-basin samples that predate ANY removal
    #### (n_removals_to_date == 0) — left as NA, not 0, since 0 would
    #### wrongly imply "a removal happened this same year".
    monitoring_phase = case_when(
      removal_status == "control" ~ "control",
      n_removals_to_date == 0 ~ "pre_removal",
      TRUE ~ "post_removal"
    ) %>% factor(levels = c("control", "pre_removal", "post_removal"))
  )

cat("\nYear-0 monitoring_phase counts:\n")
print(year0_data %>% count(monitoring_phase))
cat("\nn_removals_to_date summary (post_removal phase only):\n")
print(summary(year0_data$n_removals_to_date[year0_data$monitoring_phase == "post_removal"]))
cat("\nyears_since_last_removal summary (post_removal phase only):\n")
print(summary(year0_data$years_since_last_removal[year0_data$monitoring_phase == "post_removal"]))

log_step(year0_data, "09c-ii_year0_removal_history", "cumulative removal count + recency added")


#### 9d. Join AMBER-style covariates onto Year-0 data ################
#### Also joins n_removals_before_first_measurement (basin-level,
#### built in Section 6e) — relevant here too, since a "control" basin
#### could still have had removals predating the whole monitoring
#### window, same confound as for the baseline structure.
year0_data <- year0_data %>%
  join_amber_covariates() %>%
  left_join(n_removals_before_first_measurement_df, by = "HYBAS_ID") %>%
  mutate(
    n_removals_before_first_measurement = replace_na(n_removals_before_first_measurement, 0L),
    z_dist_up = z_safe(log_dist_up), z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual), z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    z_additional_removals_since_reference = z_safe(log1p(additional_removals_since_reference)),
    z_n_removals_to_date = z_safe(log1p(n_removals_to_date)),
    z_n_removals_before_first_measurement = z_safe(log1p(n_removals_before_first_measurement)),
    removal_status = factor(removal_status, levels = c("control", "removed")),
    site_id = factor(site_id), HYBAS_ID = factor(HYBAS_ID), country = factor(country),
    method = factor(method), sample_year = factor(sample_year)
  )

#### CRITICAL FIX: years_since_last_removal is NA for control and
#### pre_removal rows by construction (nothing to be "since"). If left
#### as a plain NA-containing column, brms would silently drop EVERY
#### control-basin row via listwise deletion when this term is in the
#### formula — defeating the entire purpose of including controls.
#### Fix: z-score ONLY within the post_removal subset (so the scaling
#### is statistically meaningful there), then set it to a neutral 0
#### for all other rows. Because those rows are constant at 0
#### regardless of beta, they contribute nothing to the linear
#### predictor and aren't excluded — the coefficient is estimated
#### purely from variation within post_removal rows, which is exactly
#### the intended interpretation ("recency of last removal, among
#### samples that have experienced one").
z_years_since_last_removal_vec <- rep(0, nrow(year0_data))
post_removal_mask <- year0_data$monitoring_phase == "post_removal"
z_years_since_last_removal_vec[post_removal_mask] <- z_safe(year0_data$years_since_last_removal[post_removal_mask])
year0_data$z_years_since_last_removal <- z_years_since_last_removal_vec

log_step(year0_data, "09d_year0_with_covariates", "AMBER-style covariates + removal-history covariates joined")


#### ============================================================
#### 10. Build 4 YEAR-0 model-ready datasets
#### ============================================================

common_year0_cols <- c("site_id", "HYBAS_ID", "country", "method", "sample_year",
                       "removal_status", "monitoring_phase", "years_relative_to_reference", "log_years_relative",
                       "z_additional_removals_since_reference", "z_n_removals_to_date",
                       "z_years_since_last_removal", "z_n_removals_before_first_measurement",
                       "z_dist_up", "z_dist_mouth",
                       "z_log_population_density", "z_mean_temp_annual", "z_precip_annual", "z_log_river_km")

DRE1ii_TotalRichnessYear0_data <- year0_data %>%
  filter(!is.na(total_richness)) %>%
  transmute(value = as.integer(total_richness), across(all_of(common_year0_cols))) %>%
  droplevels()

DRE2ii_NativeRichnessYear0_data <- year0_data %>%
  filter(!is.na(native_richness)) %>%
  transmute(value = as.integer(native_richness), across(all_of(common_year0_cols))) %>%
  droplevels()

DRE3ii_PropTotalYear0_data <- year0_data %>%
  filter(!is.na(prop_non_native_total), total_richness > 0) %>%
  transmute(non_native = as.integer(n_non_native_taxa), total_taxa = as.integer(total_richness),
            prop_non_native = as.numeric(prop_non_native_total), across(all_of(common_year0_cols))) %>%
  droplevels()

DRE4ii_PropNativeYear0_data <- year0_data %>%
  filter(!is.na(prop_non_native_known), total_known_taxa > 0) %>%
  transmute(non_native = as.integer(n_non_native_taxa), total_taxa = as.integer(total_known_taxa),
            prop_non_native = as.numeric(prop_non_native_known), across(all_of(common_year0_cols))) %>%
  droplevels()

log_step(DRE1ii_TotalRichnessYear0_data, "10_DRE1ii_year0", "final Year-0 total richness rows")
log_step(DRE2ii_NativeRichnessYear0_data, "10_DRE2ii_year0", "final Year-0 native richness rows")
log_step(DRE3ii_PropTotalYear0_data, "10_DRE3ii_year0", "final Year-0 proportion (total denom) rows")
log_step(DRE4ii_PropNativeYear0_data, "10_DRE4ii_year0", "final Year-0 proportion (known denom) rows")


#### ============================================================
#### 11. Save all 8 datasets
#### ============================================================

save_and_report <- function(dat, name) {
  saveRDS(dat, file.path(dre_out_dir, paste0(name, ".rds")))
  write_csv(dat, file.path(dre_out_dir, paste0(name, ".csv")))
  cat(" -", name, " (n =", nrow(dat), ")\n")
}

cat("\n============================================================\n")
cat("SAVING 8 FINAL DRE DATASETS\n")
cat("============================================================\n")
save_and_report(DRE1i_TotalRichnessBaseline_data, "DRE1i_TotalRichnessBaseline_data")
save_and_report(DRE2i_NativeRichnessBaseline_data, "DRE2i_NativeRichnessBaseline_data")
save_and_report(DRE3i_PropTotalBaseline_data, "DRE3i_PropTotalBaseline_data")
save_and_report(DRE4i_PropNativeBaseline_data, "DRE4i_PropNativeBaseline_data")
save_and_report(DRE1ii_TotalRichnessYear0_data, "DRE1ii_TotalRichnessYear0_data")
save_and_report(DRE2ii_NativeRichnessYear0_data, "DRE2ii_NativeRichnessYear0_data")
save_and_report(DRE3ii_PropTotalYear0_data, "DRE3ii_PropTotalYear0_data")
save_and_report(DRE4ii_PropNativeYear0_data, "DRE4ii_PropNativeYear0_data")

write_csv(record_log, file.path(dre_out_dir, "TableS1_DRE_data_cleaning_record_counts.csv"))


#### ============================================================
#### 12. DIAGNOSTICS: correlations, VIF, snap-threshold sensitivity,
#### removal-vs-control balance check
#### ============================================================

pdf(file.path(diag_dir, "DRE_assumption_checking_figures.pdf"), width = 9, height = 6)

#### 12a. Predictor correlation matrix + heatmap (baseline structure) ####
baseline_predictor_cols <- c("z_network_distance_km", "z_n_amber_barriers_on_path",
                             "z_additional_removals_since_event", "z_years_pre_to_removal",
                             "z_n_removals_before_first_measurement", "z_n_prior_removals_at_focal_event",
                             "z_dist_up", "z_dist_mouth", "z_log_population_density",
                             "z_mean_temp_annual", "z_precip_annual", "z_log_river_km")

baseline_cor <- cor(DRE1i_TotalRichnessBaseline_data[baseline_predictor_cols], use = "pairwise.complete.obs")
write.csv(baseline_cor, file.path(diag_dir, "QC_DRE_baseline_predictor_correlation.csv"))

baseline_cor_long <- as.data.frame(baseline_cor) %>% mutate(var1 = rownames(.)) %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation")

print(
  ggplot(baseline_cor_long, aes(x = var1, y = var2, fill = correlation)) +
    geom_tile() + geom_text(aes(label = round(correlation, 2)), size = 2.8) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "DRE baseline structure: predictor correlation matrix", x = NULL, y = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

#### 12b. VIF check (full baseline predictor set, via lm on total richness) ####
sapply(DRE1i_TotalRichnessBaseline_data[c(
  "z_network_distance_km", "z_n_amber_barriers_on_path",
  "z_additional_removals_since_event", "z_years_pre_to_removal",
  "z_n_removals_before_first_measurement", "z_n_prior_removals_at_focal_event",
  "z_dist_up", "z_dist_mouth", "z_log_population_density",
  "z_mean_temp_annual", "z_precip_annual", "z_log_river_km"
)], function(x) length(unique(x)))

table(baseline_nearest$n_prior_removals_at_focal_event)
table(baseline_nearest$n_removals_before_first_measurement)

summary(baseline_nearest$log_dist_up)
summary(baseline_nearest$log_dist_mouth)
summary(baseline_nearest$mean_temp_annual)
summary(baseline_nearest$precip_annual)

class(site_river_position$sample_year); head(site_river_position$sample_year)
class(climate_wide$HYBAS_ID); head(climate_wide$HYBAS_ID)
class(baseline_nearest$HYBAS_ID); head(as.character(baseline_nearest$HYBAS_ID))



class(baseline_nearest$sample_year)

#### Repair sample_year if it got silently corrupted by as.integer(factor)
baseline_nearest <- baseline_nearest %>%
  mutate(sample_year = as.integer(as.character(sample_year)))

#### Confirm it looks like real years now, not small index numbers
summary(baseline_nearest$sample_year)

#### Now safe to rebuild the covariate join
baseline_nearest <- baseline_nearest %>%
  select(-any_of(c(
    "HYRIV_ID", "dist_up_km", "dist_down_km", "dist_mouth_km",
    "log_dist_up", "log_dist_down", "log_dist_mouth",
    "mean_temp_annual", "precip_annual", "river_km", "log_river_km",
    "log_population_density", "log_population_density_1961"
  ))) %>%
  join_amber_covariates() %>%
  mutate(
    z_dist_up = z_safe(log_dist_up), z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual), z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    z_network_distance_km = z_safe(log1p(network_distance_km)),
    z_n_amber_barriers_on_path = z_safe(log1p(n_amber_barriers_on_path)),
    z_additional_removals_since_event = z_safe(log1p(additional_removals_since_event)),
    z_years_pre_to_removal = z_safe(years_pre_to_removal),
    z_n_removals_before_first_measurement = z_safe(log1p(n_removals_before_first_measurement)),
    z_n_prior_removals_at_focal_event = z_safe(log1p(n_prior_removals_at_focal_event)),
    site_id = factor(site_id), HYBAS_ID = factor(HYBAS_ID), country = factor(country),
    method = factor(method), sample_year = factor(sample_year)
  )

#### Verify the fix actually worked before re-running VIF
summary(baseline_nearest$log_dist_up)
summary(baseline_nearest$mean_temp_annual)


baseline_nearest <- baseline_nearest %>%
  mutate(sample_year = YearRemoved + years_since_removal)

#### Should now show a real year range, not 1-26
summary(baseline_nearest$sample_year)

baseline_nearest <- baseline_nearest %>%
  select(-any_of(c(
    "HYRIV_ID", "dist_up_km", "dist_down_km", "dist_mouth_km",
    "log_dist_up", "log_dist_down", "log_dist_mouth",
    "mean_temp_annual", "precip_annual", "river_km", "log_river_km",
    "log_population_density", "log_population_density_1961"
  ))) %>%
  join_amber_covariates() %>%
  mutate(
    z_dist_up = z_safe(log_dist_up), z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual), z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    z_network_distance_km = z_safe(log1p(network_distance_km)),
    z_n_amber_barriers_on_path = z_safe(log1p(n_amber_barriers_on_path)),
    z_additional_removals_since_event = z_safe(log1p(additional_removals_since_event)),
    z_years_pre_to_removal = z_safe(years_pre_to_removal),
    z_n_removals_before_first_measurement = z_safe(log1p(n_removals_before_first_measurement)),
    z_n_prior_removals_at_focal_event = z_safe(log1p(n_prior_removals_at_focal_event)),
    site_id = factor(site_id), HYBAS_ID = factor(HYBAS_ID), country = factor(country),
    method = factor(method), sample_year = factor(sample_year)
  )

summary(baseline_nearest$log_dist_up)
summary(baseline_nearest$mean_temp_annual)



DRE1i_TotalRichnessBaseline_data <- baseline_nearest %>%
  filter(!is.na(total_richness), !is.na(pre_total_richness)) %>%
  transmute(value = as.integer(total_richness), pre_value = as.numeric(pre_total_richness),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

DRE2i_NativeRichnessBaseline_data <- baseline_nearest %>%
  filter(!is.na(native_richness), !is.na(pre_native_richness)) %>%
  transmute(value = as.integer(native_richness), pre_value = as.numeric(pre_native_richness),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

DRE3i_PropTotalBaseline_data <- baseline_nearest %>%
  filter(!is.na(prop_non_native_total), !is.na(pre_prop_non_native_total), total_richness > 0) %>%
  transmute(non_native = as.integer(n_non_native_taxa), total_taxa = as.integer(total_richness),
            prop_non_native = as.numeric(prop_non_native_total),
            pre_value = as.numeric(pre_prop_non_native_total),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

DRE4i_PropNativeBaseline_data <- baseline_nearest %>%
  filter(!is.na(prop_non_native_known), !is.na(pre_prop_non_native_known), total_known_taxa > 0) %>%
  transmute(non_native = as.integer(n_non_native_taxa), total_taxa = as.integer(total_known_taxa),
            prop_non_native = as.numeric(prop_non_native_known),
            pre_value = as.numeric(pre_prop_non_native_known),
            across(all_of(common_baseline_cols))) %>%
  droplevels()

summary(DRE1i_TotalRichnessBaseline_data$z_dist_up)
summary(DRE1i_TotalRichnessBaseline_data$z_mean_temp_annual)



vif_check_dre <- lm(value ~ z_network_distance_km + z_n_amber_barriers_on_path +
                      z_additional_removals_since_event + z_years_pre_to_removal +
                      z_n_removals_before_first_measurement + z_n_prior_removals_at_focal_event +
                      z_dist_up + z_dist_mouth + z_log_population_density +
                      z_mean_temp_annual + z_precip_annual + z_log_river_km,
                    data = DRE1i_TotalRichnessBaseline_data)
vif_results_dre <- car::vif(vif_check_dre)
vif_df_dre <- tibble(term = names(vif_results_dre), VIF = as.numeric(vif_results_dre))
write_csv(vif_df_dre, file.path(diag_dir, "QC_DRE_predictor_VIF.csv"))
cat("\nDRE baseline predictor VIF check:\n")
print(vif_df_dre)
cat("\nWATCH: z_n_prior_removals_at_focal_event and z_additional_removals_since_event\n")
cat("both derive from the same removal-count process (before vs. after the focal\n")
cat("event) and could be collinear in basins with many removals clustered in time —\n")
cat("check the VIF for these two specifically before trusting both coefficients.\n")
cat("\nNOTE: the earlier glmmTMB exploratory work found z_n_segment_steps\n")
cat("collinear with z_network_distance_km (VIF ~9-10) and dropped it — that\n")
cat("term is not included in this predictor set at all, consistent with\n")
cat("that finding. If z_network_distance_km itself now shows high VIF against\n")
cat("z_dist_up/z_dist_mouth, consider dropping one of the three.\n")

#### 12c. Snapping / distance-to-mouth / distance-upstream sensitivity ####
#### Re-uses AMBER's site_river_position at multiple hypothetical snap
#### thresholds is not directly re-derivable here without the raw
#### points (that lives in the earlier DRE spatial-join sections, not
#### shown). Instead, this checks SENSITIVITY of results to a
#### downstream analysis choice available at this stage: whether
#### excluding "other_branch_same_microbasin" connections (same basin,
#### different tributary — the least direct/most uncertain distance
#### measure) changes model-relevant summary statistics meaningfully.
snap_sensitivity_dre <- tibble(
  connection_set = c("all_connected", "direct_only_no_other_branch"),
  n_rows = c(nrow(DRE1i_TotalRichnessBaseline_data),
             sum(DRE1i_TotalRichnessBaseline_data$network_direction != "other_branch_same_microbasin")),
  median_network_distance_km = c(
    median(exp(DRE1i_TotalRichnessBaseline_data$z_network_distance_km), na.rm = TRUE),  # approx, on z-scale — see note
    NA_real_)
)
cat("\nDRE connection-type sensitivity (direct upstream/downstream/same-segment\n")
cat("vs. all HydroRIVERS-connected including other-branch):\n")
print(DRE1i_TotalRichnessBaseline_data %>% count(network_direction, sort = TRUE) %>%
        mutate(pct = round(100 * n / sum(n), 1)))
write_csv(DRE1i_TotalRichnessBaseline_data %>% count(network_direction, sort = TRUE),
          file.path(diag_dir, "QC_DRE_network_direction_sensitivity.csv"))
cat("\nFor a true snap-THRESHOLD sensitivity check (100m/250m/500m/1km/2km,\n")
cat("matching AMBER's Section 5b), re-run the underlying community/DRE point\n")
cat("snapping step (upstream of this script, in your DRE Sections 1-25) at\n")
cat("each threshold and re-run this whole script on each — that is the only\n")
cat("way to test it properly, since the snap itself happens before this script.\n")

#### 12d. Removal-vs-control basin baseline balance check ############
#### Are control basins systematically different from removal basins
#### BEFORE any removal happened? If yes, the removal_status effect in
#### the Year-0 models is confounded with pre-existing differences,
#### not a clean causal test.
balance_check <- community_samples_network %>%
  mutate(is_removed_basin = HYBAS_ID %in% removed_basin_ids) %>%
  group_by(site_id) %>%
  slice_min(order_by = sample_year, n = 1, with_ties = FALSE) %>%  # first sample per site only
  ungroup() %>%
  group_by(is_removed_basin) %>%
  summarise(
    n_sites = n(),
    mean_total_richness = mean(total_richness, na.rm = TRUE),
    mean_native_richness = mean(native_richness, na.rm = TRUE),
    mean_prop_non_native_total = mean(prop_non_native_total, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(balance_check, file.path(diag_dir, "QC_removal_vs_control_baseline_balance.csv"))
cat("\nRemoval vs. control basin balance check (using each site's FIRST\n")
cat("sample only, before any removal could have had an effect):\n")
print(balance_check)
cat("\nLarge differences here mean removal and control basins were not\n")
cat("comparable to begin with — the Year-0 removal_status effect would then\n")
cat("reflect pre-existing basin differences as much as the removal itself.\n")

#### 12e. Response distributions + dispersion/zero-inflation checks ####
#### Matches AMBER Section 16's standard: histograms for every response,
#### negbinomial dispersion checks for the two richness responses, and a
#### binomial deviance/dispersion check for the two proportion responses
#### (now that binomial, not beta-binomial, is the primary family).

print(ggplot(DRE1i_TotalRichnessBaseline_data, aes(x = value)) + geom_histogram(binwidth = 1) +
        labs(title = "DRE baseline: total richness — raw distribution", x = "Total richness", y = "Count") + theme_classic())
print(ggplot(DRE2i_NativeRichnessBaseline_data, aes(x = value)) + geom_histogram(binwidth = 1) +
        labs(title = "DRE baseline: native richness — raw distribution", x = "Native richness", y = "Count") + theme_classic())
print(ggplot(DRE3i_PropTotalBaseline_data, aes(x = prop_non_native)) + geom_histogram(bins = 30) +
        labs(title = "DRE baseline: proportion non-native (total denom.) — raw distribution",
             x = "Proportion non-native", y = "Count") + theme_classic())
print(ggplot(DRE4i_PropNativeBaseline_data, aes(x = prop_non_native)) + geom_histogram(bins = 30) +
        labs(title = "DRE baseline: proportion non-native (known denom.) — raw distribution",
             x = "Proportion non-native", y = "Count") + theme_classic())

#### Overdispersion (variance/mean) — negbinomial justification, same
#### logic as AMBER Section 16b.
dre_overdispersion_check <- tibble(
  response = c("DRE1i_total_richness", "DRE2i_native_richness"),
  mean_value = c(mean(DRE1i_TotalRichnessBaseline_data$value), mean(DRE2i_NativeRichnessBaseline_data$value)),
  var_value = c(var(DRE1i_TotalRichnessBaseline_data$value), var(DRE2i_NativeRichnessBaseline_data$value))
) %>% mutate(var_mean_ratio = var_value / mean_value, note = "ratio >> 1 supports negative binomial over Poisson")
write_csv(dre_overdispersion_check, file.path(diag_dir, "QC_DRE_overdispersion_variance_mean_ratio.csv"))
cat("\nDRE overdispersion check (variance/mean ratio):\n")
print(dre_overdispersion_check)

#### Zero-inflation + dispersion via DHARMa on quick glm.nb fits
run_dre_zi_check <- function(dat, response_name) {
  m <- MASS::glm.nb(value ~ z_network_distance_km + z_n_amber_barriers_on_path +
                      z_additional_removals_since_event + z_n_removals_before_first_measurement +
                      z_dist_up + z_dist_mouth + z_log_population_density +
                      z_mean_temp_annual + z_precip_annual + z_log_river_km, data = dat)
  sim <- DHARMa::simulateResiduals(m, n = 250)
  zi_test <- DHARMa::testZeroInflation(sim, plot = FALSE)
  disp_test <- DHARMa::testDispersion(sim, plot = FALSE)
  plot(sim, main = paste(response_name, "— DHARMa residual diagnostics"))
  tibble(response = response_name, n_zero_obs = sum(dat$value == 0), pct_zero_obs = round(100 * mean(dat$value == 0), 2),
         zi_ratio_obs_sim = as.numeric(zi_test$statistic), zi_p_value = zi_test$p.value,
         dispersion_ratio = as.numeric(disp_test$statistic), dispersion_p_value = disp_test$p.value)
}

dre_zi_results <- bind_rows(
  run_dre_zi_check(DRE1i_TotalRichnessBaseline_data, "DRE1i_total_richness"),
  run_dre_zi_check(DRE2i_NativeRichnessBaseline_data, "DRE2i_native_richness")
)
write_csv(dre_zi_results, file.path(diag_dir, "QC_DRE_zero_inflation_dispersion_tests.csv"))
cat("\nDRE zero-inflation/dispersion test results:\n")
print(dre_zi_results)

#### Binomial deviance/dispersion check for the proportion responses —
#### since binomial (not beta-binomial) is now the primary family,
#### check whether it's actually adequate (non-significant dispersion
#### test) rather than assuming it from the earlier glmmTMB finding alone.
run_dre_binom_dispersion_check <- function(dat, response_name) {
  m <- glm(cbind(non_native, total_taxa - non_native) ~ z_network_distance_km + z_n_amber_barriers_on_path +
             z_additional_removals_since_event + z_n_removals_before_first_measurement +
             z_dist_up + z_dist_mouth + z_log_population_density +
             z_mean_temp_annual + z_precip_annual + z_log_river_km,
           data = dat, family = binomial())
  disp_ratio <- sum(residuals(m, type = "pearson")^2) / m$df.residual
  tibble(response = response_name, pearson_dispersion_ratio = disp_ratio,
         note = "ratio near 1 supports binomial; ratio >> 1 means overdispersion the REs aren't absorbing — reconsider beta-binomial")
}

dre_binom_dispersion <- bind_rows(
  run_dre_binom_dispersion_check(DRE3i_PropTotalBaseline_data, "DRE3i_prop_total_denom"),
  run_dre_binom_dispersion_check(DRE4i_PropNativeBaseline_data, "DRE4i_prop_known_denom")
)
write_csv(dre_binom_dispersion, file.path(diag_dir, "QC_DRE_binomial_dispersion_check.csv"))
cat("\nDRE binomial dispersion check (Pearson residual ratio, ~1 = adequate):\n")
print(dre_binom_dispersion)
cat("\nIf either ratio is well above 1, that's evidence the binomial choice needs\n")
cat("revisiting for that specific response — don't assume it transfers automatically\n")
cat("from the glmmTMB finding on the OLD data structure (nearest-removal only, no\n")
cat("control basins, no removal-history covariates) to this expanded one.\n")

dev.off()

cat("\n============================================================\n")
cat("DRE v2 PREP COMPLETE\n")
cat("============================================================\n")
cat("8 model-ready datasets saved to:", dre_out_dir, "\n")
cat("Diagnostics saved to:", diag_dir, "\n")



#### FINAL VALIDATION — run before uploading anything to Kelvin ####
all_dre_datasets <- list(
  DRE1i_TotalRichnessBaseline = DRE1i_TotalRichnessBaseline_data,
  DRE2i_NativeRichnessBaseline = DRE2i_NativeRichnessBaseline_data,
  DRE3i_PropTotalBaseline = DRE3i_PropTotalBaseline_data,
  DRE4i_PropNativeBaseline = DRE4i_PropNativeBaseline_data,
  DRE1ii_TotalRichnessYear0 = DRE1ii_TotalRichnessYear0_data,
  DRE2ii_NativeRichnessYear0 = DRE2ii_NativeRichnessYear0_data,
  DRE3ii_PropTotalYear0 = DRE3ii_PropTotalYear0_data,
  DRE4ii_PropNativeYear0 = DRE4ii_PropNativeYear0_data
)

validation_summary <- purrr::imap_dfr(all_dre_datasets, function(dat, name) {
  numeric_cols <- names(dat)[sapply(dat, is.numeric)]
  n_constant <- sum(sapply(dat[numeric_cols], function(x) length(unique(na.omit(x))) <= 1))
  constant_names <- paste(numeric_cols[sapply(dat[numeric_cols], function(x) length(unique(na.omit(x))) <= 1)], collapse = ", ")
  
  sample_year_range <- if ("sample_year" %in% names(dat)) {
    yrs <- suppressWarnings(as.integer(as.character(dat$sample_year)))
    paste(range(yrs, na.rm = TRUE), collapse = "-")
  } else "no sample_year col"
  
  tibble(
    dataset = name,
    n_rows = nrow(dat),
    n_sites = n_distinct(dat$site_id),
    n_constant_numeric_cols = n_constant,
    constant_col_names = constant_names,
    sample_year_range = sample_year_range,
    n_country_levels = n_distinct(dat$country),
    any_na_in_key_cols = anyNA(dat[c("site_id", "HYBAS_ID", "sample_year")])
  )
})

print(validation_summary, width = Inf)

site_sample_counts %>% group_by(has_removal_in_basin) %>%
  summarise(median_samples = median(n_samples), mean_samples = round(mean(n_samples), 2), n_sites = n())

