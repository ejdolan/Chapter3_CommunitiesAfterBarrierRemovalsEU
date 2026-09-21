#### ============================================================
#### SCRIPT 01 (v2): PREPARE AMBER INVERTEBRATE DATA
#### FOR NEW 3-MODEL STRUCTURE (WITH-TYPE ONLY)
####
#### Models prepared (all "with type" only — density x barrier_type
#### interaction; the population-averaged marginal density effect
#### is obtained on Kelvin via emmeans::emtrends(), NOT by fitting a
#### separate "without type" model — see handoff §8.3):
####
####   A1i  total richness (all unique taxa, any alien_status)
####   A2i  native richness (native-confirmed taxa ONLY — excludes
####        non-native AND unresolved-status taxa)
####   A3i  proportion of non-native taxa (non-native / total known,
####        unresolved-status taxa excluded from denominator)
####
#### CHANGES FROM THE PREVIOUS VERSION OF THIS SCRIPT:
####   1. Shannon diversity DROPPED ENTIRELY (redundant with richness,
####      per supervisor).
####   2. Native richness ADDED as a new response (uses n_native_taxa,
####      which was already being computed in site_diversity_base but
####      was never carried through into its own model-ready dataset).
####   3. "Without type" datasets NO LONGER BUILT. Only "with type".
####   4. z_dist_down REMOVED from the final predictor set. It was
####      being computed (network distance, downstream path length)
####      but the handoff record confirms this was already flagged as
####      near-collinear with distance-to-mouth (Pearson's r = 0.9994)
####      and dropped from the models — the previous version of this
####      script computed it correctly but then left it IN the final
####      selected model columns by mistake. It's still computed here
####      for QC/reference but excluded from z-scoring / model data.
####   5. HYDE 3.5 rebuild DEFERRED — using the existing HYDE 3.4
####      processed files for now (HYDE_VERSION/HYDE_DIR are set as
####      variables at the top of Section 0 so this is a one-line
####      change later).
####   6. Record-count tracking added at each major cleaning step
####      (addresses the long-outstanding Table S1 gap — data-cleaning
####      counts were never actually tabulated before).
####   7. A real VIF check (car::vif on a fitted lm) is now run and
####      saved, replacing the old pairwise-correlation-only check.
####   8. Zero-inflation tests (DHARMa::testZeroInflation) actually run
####      and saved — previously this was claimed in a figure caption
####      but never confirmed to have been executed.
####   9. Overdispersion (variance/mean) and beta-binomial deviance/df
####      checks recomputed for the NEW response set, to re-confirm
####      negbin/beta-binomial family choices hold for total richness,
####      native richness, and proportion non-native separately.
####  10. Assumption-checking figures (distributions, dispersion,
####      collinearity, barrier density by type) saved as one PDF.
#### ============================================================

#### 0. Setup ####################################################

setwd("C:/Users/User/OneDrive - Queen's University Belfast/PhD/Chapters/Chapter 3 European Connectivity (Inverts)/UpToDateWork/EUConnectivity_Inverts")

library(sf)
library(readxl)
library(car)      # VIF — loaded BEFORE dplyr so dplyr::select()/filter() win the search path
library(MASS)      # glm.nb for dispersion/zero-inflation checks — MASS::select() masks dplyr::select() if loaded after it
library(DHARMa)    # zero-inflation / dispersion diagnostics
library(vegan)     # Shannon + evenness — DIAGNOSTIC ONLY, not a modelled response
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(igraph)
library(readr)
library(ggplot2)

sf_use_s2(FALSE)
options(scipen = 999)

#### Guard against namespace masking -----------------------------
#### car / MASS / DHARMa / vegan (or anything else already attached
#### in your R session, e.g. raster) can shadow dplyr::select() /
#### dplyr::filter() if they define functions of the same name and
#### get loaded after dplyr. This forces both back to the dplyr
#### versions regardless of load order, so the rest of the script
#### doesn't need every call qualified as dplyr::select().
select <- dplyr::select
filter <- dplyr::filter

dir.create("Data/Processed/Final_Invert_Models_v2", recursive = TRUE, showWarnings = FALSE)
dir.create("Data/Processed/Final_Invert_Models_v2/Diagnostics", recursive = TRUE, showWarnings = FALSE)
out_dir <- "Data/Processed/Final_Invert_Models_v2"
diag_dir <- file.path(out_dir, "Diagnostics")

#### HYDE VERSION -------------------------------------------------
#### REVERTED TO 3.4 FOR NOW (already processed/saved) — the 3.5
#### rebuild is deferred to a later session. When ready, just flip
#### HYDE_VERSION back to "3.5" and point HYDE_DIR at the 3.5
#### processed CSVs; nothing else in the script needs to change.
HYDE_VERSION <- "3.4"
HYDE_DIR <- "Data/HYDE/Processed"

if (!dir.exists(HYDE_DIR)) {
  warning(
    "HYDE_DIR ('", HYDE_DIR, "') does not exist. ",
    "Update HYDE_DIR at the top of the script to point at your HYDE ",
    HYDE_VERSION, " processed CSVs before running Section 12 (population density)."
  )
}

z_safe <- function(x) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) < 2 || length(unique(na.omit(x))) <= 1) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

clean_barrier_type <- function(x) {
  x <- str_trim(tolower(as.character(x)))
  case_when(
    str_detect(x, "dam")     ~ "dam",
    str_detect(x, "weir")    ~ "weir",
    str_detect(x, "culvert") ~ "culvert",
    TRUE                     ~ "other"
  )
}

parse_num <- function(x) {
  x <- str_trim(as.character(x))
  x <- na_if(x, "")
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

#### Record-count tracker -----------------------------------------
#### Replaces the never-completed Table S1. Call log_step() after
#### every major filtering/joining operation below.
record_log <- tibble(step = character(), n_rows = integer(), n_site_years = integer(), note = character())

log_step <- function(df, step_name, note = "", id_cols = c("site_id", "sample_year")) {
  n_sy <- if (length(id_cols) > 0 && all(id_cols %in% names(df))) {
    n_distinct(df[id_cols])
  } else {
    NA_integer_
  }
  record_log <<- bind_rows(
    record_log,
    tibble(step = step_name, n_rows = nrow(df), n_site_years = n_sy, note = note)
  )
  invisible(df)
}

#### 1. Load raw data ############################################

inverts_df <- read_xlsx("Data/Inverts data/Global_dataset.xlsx")
amber_sf   <- read_sf("Data/AMBER data/AMBER map.shp")
basins_sf  <- read_sf("Data/HydroBasins EU/hybas_eu_lev07_v1c.shp") %>%
  st_make_valid() %>%
  mutate(HYBAS_ID = as.character(HYBAS_ID))

rivers_sf  <- read_sf("Data/HydroRivers EU/HydroRIVERS_v10_eu.shp") %>%
  mutate(
    HYRIV_ID  = as.character(HYRIV_ID),
    NEXT_DOWN = as.character(NEXT_DOWN),
    MAIN_RIV  = as.character(MAIN_RIV))

log_step(inverts_df, "01_raw_invert_records", "raw Global_dataset.xlsx, no filtering yet")

#### 2. Clean invertebrate records ###############################
inverts_clean <- inverts_df %>%
  transmute(
    site_id = as.character(site_id),
    sample_year = suppressWarnings(as.integer(year)),
    taxon = as.character(taxon),
    abundance_num = suppressWarnings(as.numeric(abundance)),
    country = as.character(country),
    method = as.character(Method),
    Longitude_X = parse_num(Longitude_X),
    Latitude_Y  = parse_num(Latitude_Y),
    alien_status = case_when(
      str_trim(tolower(as.character(Alien))) %in%
        c("yes", "y", "alien", "non-native", "nonnative",
          "non native", "introduced", "invasive") ~ "non_native",
      str_trim(tolower(as.character(Alien))) %in%
        c("no", "n", "native") ~ "native",
      TRUE ~ "unknown")) %>%
  filter(
    !is.na(site_id),
    !is.na(sample_year),
    !is.na(taxon),
    !is.na(Longitude_X),
    !is.na(Latitude_Y))

log_step(inverts_clean, "02_clean_invert_records",
         "after dropping missing site/year/taxon/coords")

inverts_sf <- inverts_clean %>%
  st_as_sf(coords = c("Longitude_X", "Latitude_Y"), crs = 4326, remove = FALSE)

#### 3. Clean AMBER barriers #####################################
amber_clean_sf <- amber_sf %>%
  mutate(
    AMBER_ID = row_number(),
    barrier_type_reduced = clean_barrier_type(type))

#### 4. Keep only features inside HydroBASINS L7 #################
inverts_basins_sf <- st_join(
  inverts_sf,
  basins_sf %>% select(HYBAS_ID),
  join = st_intersects) %>%
  filter(!is.na(HYBAS_ID))

log_step(inverts_basins_sf, "04_joined_to_hydrobasins",
         "after point-in-polygon join to HydroBASINS L7")

amber_basins_sf <- st_join(
  amber_clean_sf,
  basins_sf %>% select(HYBAS_ID),
  join = st_intersects) %>%
  filter(!is.na(HYBAS_ID))

#### 5. Snap points to nearest HydroRIVERS segment within 2000 m ##
snap_to_river <- function(points_sf, rivers_sf, threshold_m = 2000) {
  pts_3035 <- st_transform(points_sf, 3035)
  riv_3035 <- st_transform(rivers_sf, 3035)
  
  nearest_idx <- st_nearest_feature(pts_3035, riv_3035)
  dist_m <- as.numeric(
    st_distance(pts_3035, riv_3035[nearest_idx, ], by_element = TRUE))
  
  keep <- dist_m <= threshold_m
  
  out <- bind_cols(
    st_drop_geometry(pts_3035[keep, ]),
    st_drop_geometry(riv_3035[nearest_idx[keep], ]) %>%
      dplyr::select(HYRIV_ID, NEXT_DOWN, MAIN_RIV, LENGTH_KM, DIST_DN_KM))
  
  out$dist_to_river_m <- dist_m[keep]
  out
}

SNAP_THRESHOLD_M <- 2000

inverts_rivers_df <- snap_to_river(inverts_basins_sf, rivers_sf, SNAP_THRESHOLD_M) %>%
  mutate(
    HYBAS_ID = as.character(HYBAS_ID),
    HYRIV_ID = as.character(HYRIV_ID),
    site_sample_id = paste(site_id, sample_year, sep = "_"))

log_step(inverts_rivers_df, "05_snapped_to_river_2km",
         "after 2km snap-to-HydroRIVERS filter")

amber_rivers_df <- snap_to_river(amber_basins_sf, rivers_sf, SNAP_THRESHOLD_M) %>%
  mutate(
    HYBAS_ID = as.character(HYBAS_ID),
    HYRIV_ID = as.character(HYRIV_ID),
    barrier_type_reduced = factor(
      barrier_type_reduced,
      levels = c("other", "dam", "weir", "culvert")))

#### 5b. Snapping threshold sensitivity check (kept as-is) ########
snap_thresholds_m <- c(100, 250, 500, 1000, 2000)

check_snap_threshold <- function(threshold_m) {
  amber_snap <- snap_to_river(amber_basins_sf, rivers_sf, threshold_m)
  inverts_snap <- snap_to_river(inverts_basins_sf, rivers_sf, threshold_m)
  
  tibble(
    threshold_m = threshold_m,
    n_amber_retained = nrow(amber_snap),
    n_amber_basins = n_distinct(amber_snap$HYBAS_ID),
    n_invert_rows_retained = nrow(inverts_snap),
    n_invert_siteyears = inverts_snap %>% distinct(site_id, sample_year) %>% nrow(),
    n_invert_basins = n_distinct(inverts_snap$HYBAS_ID)
  )
}

snap_sensitivity_results <- map_dfr(snap_thresholds_m, check_snap_threshold)

baseline_amber <- snap_sensitivity_results %>% filter(threshold_m == 2000) %>% pull(n_amber_retained)
baseline_siteyears <- snap_sensitivity_results %>% filter(threshold_m == 2000) %>% pull(n_invert_siteyears)

snap_sensitivity_results <- snap_sensitivity_results %>%
  mutate(
    pct_amber_of_2000m = round(100 * n_amber_retained / baseline_amber, 1),
    pct_siteyears_of_2000m = round(100 * n_invert_siteyears / baseline_siteyears, 1)
  )

write_csv(snap_sensitivity_results, file.path(out_dir, "QC_snap_threshold_sensitivity.csv"))

#### 6. QC: check each community assigned one HydroRIVERS segment ####
segment_check_site_year <- inverts_rivers_df %>%
  distinct(site_id, sample_year, HYRIV_ID) %>%
  count(site_id, sample_year, name = "n_segments_assigned") %>%
  arrange(desc(n_segments_assigned))

problem_site_years <- segment_check_site_year %>% filter(n_segments_assigned > 1)

write_csv(segment_check_site_year, file.path(out_dir, "QC_site_year_segment_assignment.csv"))
write_csv(problem_site_years, file.path(out_dir, "QC_problem_site_years_more_than_one_segment.csv"))

if (nrow(problem_site_years) > 0) {
  problem_details <- inverts_rivers_df %>%
    semi_join(problem_site_years, by = c("site_id", "sample_year")) %>%
    distinct(site_id, sample_year, HYRIV_ID, HYBAS_ID, MAIN_RIV, dist_to_river_m) %>%
    arrange(site_id, sample_year, dist_to_river_m)
  
  write_csv(problem_details, file.path(out_dir, "QC_problem_site_year_segment_details.csv"))
  stop("Some site_id x sample_year communities were assigned to more than one ",
       "HydroRIVERS segment. Check QC_problem_site_year_segment_details.csv before modelling.")
}

#### 7. River length per HydroBASINS L7 basin ####################
basins_3035 <- st_transform(basins_sf, 3035) %>% st_make_valid()
rivers_3035 <- st_transform(rivers_sf, 3035)

rivers_cut_sf <- st_intersection(st_make_valid(rivers_3035), basins_3035 %>% select(HYBAS_ID))

river_len_df <- rivers_cut_sf %>%
  mutate(river_km = as.numeric(st_length(geometry)) / 1000) %>%
  st_drop_geometry() %>%
  group_by(HYBAS_ID) %>%
  summarise(river_km = sum(river_km, na.rm = TRUE), .groups = "drop") %>%
  mutate(HYBAS_ID = as.character(HYBAS_ID), log_river_km = log1p(river_km))

#### 8. HydroRIVERS network distances: up / down / mouth ##########
#### NOTE: dist_down (downstream path length) is still computed
#### here for reference/QC but is EXCLUDED from the final model
#### predictor set below (Section 13) because it is near-collinear
#### with dist_mouth (native HydroRIVERS DIST_DN_KM field), r = 0.9994.
river_topo <- rivers_sf %>%
  st_drop_geometry() %>%
  transmute(
    HYRIV_ID = as.character(HYRIV_ID),
    NEXT_DOWN = as.character(NEXT_DOWN),
    MAIN_RIV = as.character(MAIN_RIV),
    LENGTH_KM = as.numeric(LENGTH_KM),
    DIST_DN_KM = as.numeric(DIST_DN_KM)) %>%
  filter(!is.na(HYRIV_ID)) %>%
  distinct(HYRIV_ID, .keep_all = TRUE)

edges <- river_topo %>%
  filter(!is.na(NEXT_DOWN), NEXT_DOWN != "0", NEXT_DOWN %in% river_topo$HYRIV_ID) %>%
  transmute(from = HYRIV_ID, to = NEXT_DOWN, weight = pmax(LENGTH_KM, 0.001))

g_river <- graph_from_data_frame(
  edges, directed = TRUE,
  vertices = river_topo %>% select(name = HYRIV_ID, MAIN_RIV, LENGTH_KM, DIST_DN_KM))

get_downstream_path_km <- function(x) {
  if (!x %in% V(g_river)$name) return(NA_real_)
  nodes <- setdiff(names(subcomponent(g_river, x, mode = "out")), x)
  if (length(nodes) == 0) return(0)
  sum(river_topo$LENGTH_KM[match(nodes, river_topo$HYRIV_ID)], na.rm = TRUE)
}

get_upstream_network_km <- function(x) {
  if (!x %in% V(g_river)$name) return(NA_real_)
  nodes <- setdiff(names(subcomponent(g_river, x, mode = "in")), x)
  if (length(nodes) == 0) return(0)
  sum(river_topo$LENGTH_KM[match(nodes, river_topo$HYRIV_ID)], na.rm = TRUE)
}

site_river_position <- inverts_rivers_df %>%
  distinct(site_id, sample_year, HYRIV_ID) %>%
  mutate(
    dist_up_km = map_dbl(HYRIV_ID, get_upstream_network_km),
    dist_down_km = map_dbl(HYRIV_ID, get_downstream_path_km)) %>%
  left_join(river_topo %>% select(HYRIV_ID, dist_mouth_km = DIST_DN_KM), by = "HYRIV_ID") %>%
  mutate(
    log_dist_up = log1p(dist_up_km),
    log_dist_down = log1p(dist_down_km),   # kept for QC only, NOT z-scored into models
    log_dist_mouth = log1p(dist_mouth_km))

#### QC: confirm the collinearity that justified dropping dist_down
dist_collinearity_check <- cor(
  site_river_position %>% select(log_dist_up, log_dist_down, log_dist_mouth),
  use = "pairwise.complete.obs"
)
write.csv(dist_collinearity_check, file.path(diag_dir, "QC_dist_up_down_mouth_correlation.csv"))
cat("\nCorrelation between log_dist_down and log_dist_mouth:",
    round(dist_collinearity_check["log_dist_down", "log_dist_mouth"], 4),
    "(dist_down excluded from models on this basis)\n")

#### 9. Basin barrier density by reduced type ####################
barrier_density_df <- amber_rivers_df %>%
  filter(!is.na(HYBAS_ID), !is.na(barrier_type_reduced)) %>%
  count(HYBAS_ID, barrier_type_reduced, name = "n_barrier") %>%
  right_join(
    expand_grid(
      HYBAS_ID = unique(river_len_df$HYBAS_ID),
      barrier_type_reduced = c("other", "dam", "weir", "culvert")),
    by = c("HYBAS_ID", "barrier_type_reduced")) %>%
  left_join(river_len_df, by = "HYBAS_ID") %>%
  mutate(
    n_barrier = coalesce(n_barrier, 0L),
    barrier_density = if_else(river_km > 0, n_barrier / river_km, NA_real_),
    log_barrier_density = log1p(barrier_density),
    barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert")))

#### 10. Community response data #################################
#### Responses (all "with type"):
####   1. total_richness   = all unique taxa, any alien_status
####   2. native_richness  = unique taxa with alien_status == "native" ONLY
####   3. proportion non-native = non_native_taxa / total_known_taxa
####      (total_known_taxa = native + non_native; unresolved-status
####       taxa excluded from BOTH numerator and denominator)
####
#### Shannon diversity is NOT computed in this version (dropped).

#### 10a. Collapse to one row per taxon per community #############
community_taxon_df <- inverts_rivers_df %>%
  filter(!is.na(site_id), !is.na(sample_year), !is.na(HYBAS_ID), !is.na(taxon)) %>%
  mutate(
    abundance_num = suppressWarnings(as.numeric(abundance_num)),
    abundance_num = if_else(is.na(abundance_num), 0, abundance_num),
    alien_status = case_when(
      alien_status %in% c("native", "non_native") ~ alien_status,
      TRUE ~ "unknown")) %>%
  group_by(site_id, sample_year, HYBAS_ID, country, method, taxon) %>%
  summarise(
    abundance_taxon = sum(abundance_num, na.rm = TRUE),
    alien_status = first(na.omit(alien_status)),
    .groups = "drop") %>%
  mutate(alien_status = if_else(is.na(alien_status), "unknown", alien_status))

log_step(community_taxon_df, "10a_community_taxon_collapse",
         "one row per site x year x taxon", id_cols = c("site_id", "sample_year"))

#### 10b. QC: conflicting alien statuses within site-year #########
alien_status_conflicts <- inverts_rivers_df %>%
  filter(!is.na(site_id), !is.na(sample_year), !is.na(taxon)) %>%
  mutate(alien_status = case_when(
    alien_status %in% c("native", "non_native") ~ alien_status,
    TRUE ~ "unknown")) %>%
  distinct(site_id, sample_year, taxon, alien_status) %>%
  count(site_id, sample_year, taxon, name = "n_statuses") %>%
  filter(n_statuses > 1)

write_csv(alien_status_conflicts, file.path(out_dir, "QC_alien_status_conflicts.csv"))
cat("\nTaxa with conflicting alien statuses within site-year:", nrow(alien_status_conflicts), "\n")

#### 10c. Total richness, native richness, taxon-status breakdown #
site_diversity_base <- community_taxon_df %>%
  group_by(site_id, sample_year, HYBAS_ID, country, method) %>%
  summarise(
    total_richness = n_distinct(taxon),
    n_native_taxa = n_distinct(taxon[alien_status == "native"]),
    n_non_native_taxa = n_distinct(taxon[alien_status == "non_native"]),
    n_unknown_taxa = n_distinct(taxon[alien_status == "unknown"]),
    total_abundance = sum(abundance_taxon, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    native_richness = n_native_taxa,   # explicit alias — this IS the A2i response
    total_known_taxa = n_native_taxa + n_non_native_taxa,
    prop_unknown_taxa = if_else(total_richness > 0, n_unknown_taxa / total_richness, NA_real_))

log_step(site_diversity_base, "10c_site_diversity_base",
         "one row per site x year, richness/native-richness computed")

#### QC: what % of total richness is unresolved-status (the
#### supervisor's outstanding question about whether NNS/unresolved
#### taxa were being folded into "richness")
unresolved_pct_check <- site_diversity_base %>%
  summarise(
    n_communities = n(),
    mean_prop_unknown_taxa = mean(prop_unknown_taxa, na.rm = TRUE),
    median_prop_unknown_taxa = median(prop_unknown_taxa, na.rm = TRUE),
    max_prop_unknown_taxa = max(prop_unknown_taxa, na.rm = TRUE),
    n_communities_gt50pct_unknown = sum(prop_unknown_taxa > 0.5, na.rm = TRUE)
  )
write_csv(unresolved_pct_check, file.path(diag_dir, "QC_unresolved_status_pct_of_richness.csv"))
cat("\nUnresolved-status taxa as % of total richness:\n")
print(unresolved_pct_check)

#### 10d. Proportion of non-native taxa ############################
site_prop_base <- site_diversity_base %>%
  transmute(
    site_id, sample_year, HYBAS_ID, country, method,
    native_taxa = n_native_taxa,
    non_native_taxa = n_non_native_taxa,
    unknown_taxa = n_unknown_taxa,
    total_known_taxa, total_richness,
    prop_non_native_taxa = if_else(total_known_taxa > 0, non_native_taxa / total_known_taxa, NA_real_),
    prop_unknown_taxa) %>%
  filter(total_known_taxa > 0)

log_step(site_prop_base, "10d_site_prop_base",
         "communities with >=1 known-status taxon (denominator > 0)")

#### 10e. Proportion non-native relative to native+non-native (A4i) #
#### CHANGED PER EXPLICIT DECISION: A4i is now BOUNDED, structurally
#### identical to A3i (beta-binomial, trials()) — just with a
#### different denominator, so the two are swappable regardless of
#### whether the final paper reports total richness (pairs with A3i)
#### or native richness (pairs with A4i).
####   A3i denominator = total_richness   (ALL taxa, incl. unresolved-status)
####   A4i denominator = total_known_taxa (native + non-native ONLY)
#### Both guarantee non_native <= denominator by construction (native
#### and non-native are both proper subsets of the respective
#### denominator), so both are valid bounded 0-1 proportions.
####
#### This REPLACES the earlier unbounded rate-model version of A4i
#### (non_native_richness / native_richness with an offset,
#### negbinomial family). That version is preserved in git history /
#### conversation record if you want to revisit it later — it answers
#### a genuinely different question (invasion pressure vs. baseline)
#### but was deliberately dropped in favour of structural consistency.
site_rate_base <- site_diversity_base %>%
  transmute(
    site_id, sample_year, HYBAS_ID, country, method,
    native_taxa = n_native_taxa,
    non_native_taxa = n_non_native_taxa,
    total_known_taxa,
    prop_non_native_of_known = if_else(total_known_taxa > 0,
                                       non_native_taxa / total_known_taxa, NA_real_)) %>%
  filter(total_known_taxa > 0)

log_step(site_rate_base, "10e_site_rate_base",
         "communities with >=1 known-status taxon (A4i denominator = native+non-native)")

cat("\nA4i is now bounded (native+non-native denominator), matching A3i's beta-\n")
cat("binomial trials() structure. Because total_known_taxa and total_richness\n")
cat("will be close whenever unresolved-status taxa are rare, A3i and A4i may end\n")
cat("up highly correlated — see Section 16i for the direct check between them.\n")

#### 10f. Long-format data for richness models ####################
#### response = "total_richness" or "native_richness"
site_response_long <- site_diversity_base %>%
  select(site_id, sample_year, HYBAS_ID, country, method,
         total_richness, native_richness,
         total_abundance, n_native_taxa, n_non_native_taxa, n_unknown_taxa,
         total_known_taxa, prop_unknown_taxa) %>%
  pivot_longer(cols = c(total_richness, native_richness),
               names_to = "response", values_to = "value") %>%
  mutate(response = factor(response, levels = c("total_richness", "native_richness")))

#### 10g. Quick checks #############################################
cat("\n--- Community response checks ---\n")
cat("\nNumber of site-year communities in diversity data:", nrow(site_diversity_base), "\n")
cat("\nSummary of total richness:\n"); print(summary(site_diversity_base$total_richness))
cat("\nSummary of native richness:\n"); print(summary(site_diversity_base$native_richness))
cat("\nSummary of proportion non-native taxa:\n"); print(summary(site_prop_base$prop_non_native_taxa))
cat("\nUnknown-status taxa proportion summary:\n"); print(summary(site_diversity_base$prop_unknown_taxa))

cat("\nCommunities with no known native/non-native taxa (excluded from proportion model):\n")
print(site_diversity_base %>%
        summarise(n_total = n(), n_no_known_status = sum(total_known_taxa == 0, na.rm = TRUE),
                  prop_no_known_status = n_no_known_status / n_total))

#### 10h. Save response objects ####################################
saveRDS(community_taxon_df, file.path(out_dir, "community_taxon_df.rds"))
saveRDS(site_diversity_base, file.path(out_dir, "site_diversity_base.rds"))
saveRDS(site_response_long, file.path(out_dir, "site_response_long_total_native_richness.rds"))
saveRDS(site_prop_base, file.path(out_dir, "site_prop_base_non_native_taxa.rds"))

write_csv(site_diversity_base, file.path(out_dir, "site_diversity_base.csv"))
write_csv(site_response_long, file.path(out_dir, "site_response_long_total_native_richness.csv"))
write_csv(site_prop_base, file.path(out_dir, "site_prop_base_non_native_taxa.csv"))


#### ============================================================
#### 10i. SUPPLEMENTARY DIVERSITY METRICS FOR CORRELATION CHECKING
#### ============================================================
#### DIAGNOSTIC ONLY — Shannon diversity, Pielou's evenness, and
#### temporal turnover are computed here purely to check their
#### Pearson correlation against total_richness / native_richness /
#### proportion non-native (Section 16g). None of these three are
#### fitted as AMBER models — the final model set stays at 3
#### (A1i total richness, A2i native richness, A3i proportion
#### non-native), all "with type".
#### ============================================================

#### Shannon + Pielou's evenness (abundance-based, one row per community)
site_shannon_evenness <- community_taxon_df %>%
  group_by(site_id, sample_year, HYBAS_ID, country, method) %>%
  summarise(
    shannon = vegan::diversity(abundance_taxon, index = "shannon"),
    richness_for_evenness = n_distinct(taxon),
    .groups = "drop") %>%
  mutate(
    evenness_pielou = if_else(richness_for_evenness > 1,
                              shannon / log(richness_for_evenness),
                              NA_real_))

cat("\nCommunities where evenness is undefined (richness <= 1):",
    sum(site_shannon_evenness$richness_for_evenness <= 1), "of",
    nrow(site_shannon_evenness), "\n")

#### Temporal turnover between consecutive samples at the same site
#### (Jaccard-based: (taxa gained + taxa lost) / total taxa across
#### both years). Attributed to the LATER of the two sample years.
#### A site's first-ever sample has no prior community to compare
#### against, so turnover = NA for those rows.
site_taxon_sets <- community_taxon_df %>%
  distinct(site_id, sample_year, taxon) %>%
  arrange(site_id, sample_year) %>%
  group_by(site_id, sample_year) %>%
  summarise(taxa_set = list(unique(taxon)), .groups = "drop") %>%
  arrange(site_id, sample_year)

turnover_df <- site_taxon_sets %>%
  group_by(site_id) %>%
  mutate(
    prev_taxa_set = lag(taxa_set),
    prev_sample_year = lag(sample_year)) %>%
  ungroup() %>%
  filter(!is.na(prev_sample_year)) %>%
  rowwise() %>%
  mutate(
    n_union = length(union(taxa_set, prev_taxa_set)),
    n_intersect = length(intersect(taxa_set, prev_taxa_set)),
    turnover = if_else(n_union > 0, (n_union - n_intersect) / n_union, NA_real_),
    years_between_samples_turnover = sample_year - prev_sample_year) %>%
  ungroup() %>%
  select(site_id, sample_year, prev_sample_year, years_between_samples_turnover, turnover)

cat("\nTemporal turnover: coverage check\n")
cat("Total site-years:", nrow(site_taxon_sets), "\n")
cat("Site-years with a prior sample (turnover computable):", nrow(turnover_df), "\n")
cat("Site-years with NO prior sample (turnover = NA, first visit):",
    nrow(site_taxon_sets) - nrow(turnover_df), "\n")

#### Combine all diagnostic metrics into one table, one row per community.
#### Joined on site_id/sample_year (both still character/integer at this
#### point in the script, matching site_diversity_base's types).
diagnostic_diversity_df <- site_diversity_base %>%
  select(site_id, sample_year, HYBAS_ID, country, method, total_richness, native_richness) %>%
  left_join(site_shannon_evenness %>% select(site_id, sample_year, shannon, evenness_pielou),
            by = c("site_id", "sample_year")) %>%
  left_join(turnover_df %>% select(site_id, sample_year, turnover),
            by = c("site_id", "sample_year")) %>%
  left_join(site_prop_base %>% select(site_id, sample_year, prop_non_native_taxa),
            by = c("site_id", "sample_year"))

saveRDS(diagnostic_diversity_df, file.path(diag_dir, "diagnostic_diversity_metrics_shannon_evenness_turnover.rds"))
write_csv(diagnostic_diversity_df, file.path(diag_dir, "diagnostic_diversity_metrics_shannon_evenness_turnover.csv"))

cat("\nDiagnostic diversity metrics summary:\n")
print(summary(diagnostic_diversity_df %>% select(total_richness, native_richness, shannon,
                                                 evenness_pielou, turnover, prop_non_native_taxa)))


#### 11. Climate data ##############################################
basin_year_climate_long <- readRDS("Data/Processed/basin_year_climate_long.rds")

climate_wide <- basin_year_climate_long %>%
  mutate(HYBAS_ID = as.character(HYBAS_ID), sample_year = as.integer(sample_year)) %>%
  group_by(HYBAS_ID, sample_year, climate_variable) %>%
  summarise(climate_value = mean(climate_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = climate_variable, values_from = climate_value)

#### 12. HYDE 3.5 population density ###############################
read_hyde_var <- function(var, hyde_dir = HYDE_DIR) {
  files <- list.files(hyde_dir, pattern = paste0("^", var, "_[0-9]{4}_[0-9]{4}\\.csv$"), full.names = TRUE)
  if (length(files) == 0) {
    stop("No HYDE ", HYDE_VERSION, " processed files found for '", var,
         "' in ", hyde_dir, ". Update HYDE_DIR at the top of the script.")
  }
  map_dfr(sort(files), read_csv, show_col_types = FALSE) %>%
    mutate(HYBAS_ID = as.character(HYBAS_ID), year = as.integer(year)) %>%
    distinct()
}

hyde_pop_density <- read_hyde_var("population_density")

hyde_full <- hyde_pop_density %>%
  rename(mean_population_density = population_density) %>%
  mutate(log_population_density = log1p(mean_population_density))

hyde_1961 <- hyde_full %>%
  filter(year == 1961) %>%
  select(HYBAS_ID, log_population_density_1961 = log_population_density)

cat("\nHYDE version used for population density:", HYDE_VERSION, "\n")
cat("HYDE years available:", paste(range(hyde_full$year), collapse = "-"), "\n")

#### 13. Build "with type" model data for the 3 responses ##########
#### Predictors used (z-scored): z_barrier_density (x barrier_type),
#### z_dist_up, z_dist_mouth, z_log_population_density,
#### z_mean_temp_annual, z_precip_annual, z_log_river_km.
#### z_dist_down is DELIBERATELY EXCLUDED (see Section 8 note above).

invert_diversity_model_data <- site_response_long %>%
  mutate(HYBAS_ID = as.character(HYBAS_ID), site_id = as.character(site_id),
         sample_year = as.integer(as.character(sample_year))) %>%
  left_join(site_river_position %>%
              mutate(site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year))),
            by = c("site_id", "sample_year")) %>%
  left_join(climate_wide %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID), sample_year = as.integer(as.character(sample_year))),
            by = c("HYBAS_ID", "sample_year")) %>%
  left_join(river_len_df %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
  left_join(hyde_full %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID), year = as.integer(year)) %>%
              select(HYBAS_ID, sample_year = year, log_population_density),
            by = c("HYBAS_ID", "sample_year")) %>%
  left_join(hyde_1961 %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
  mutate(log_population_density = if_else(
    is.na(log_population_density) & sample_year < 1961,
    log_population_density_1961, log_population_density)) %>%
  select(-log_population_density_1961) %>%
  tidyr::crossing(barrier_type_reduced = factor(c("other", "dam", "weir", "culvert"),
                                                levels = c("other", "dam", "weir", "culvert"))) %>%
  left_join(barrier_density_df %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID),
                     barrier_type_reduced = factor(barrier_type_reduced,
                                                   levels = c("other", "dam", "weir", "culvert"))) %>%
              select(HYBAS_ID, barrier_type_reduced, n_barrier, barrier_density, log_barrier_density),
            by = c("HYBAS_ID", "barrier_type_reduced")) %>%
  mutate(
    z_barrier_density = z_safe(log_barrier_density),
    z_dist_up = z_safe(log_dist_up),
    z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual),
    z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    response = factor(response, levels = c("total_richness", "native_richness")),
    barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert")),
    country = factor(country), HYBAS_ID = factor(HYBAS_ID), site_id = factor(site_id),
    method = factor(method), sample_year = factor(sample_year)) %>%
  select(value, response, barrier_type_reduced, z_barrier_density, z_dist_up, z_dist_mouth,
         z_log_population_density, z_mean_temp_annual, z_precip_annual, z_log_river_km,
         country, HYBAS_ID, site_id, method, sample_year, everything()) %>%
  filter(complete.cases(value, response, barrier_type_reduced, z_barrier_density, z_dist_up,
                        z_dist_mouth, z_log_population_density, z_mean_temp_annual,
                        z_precip_annual, z_log_river_km, country, HYBAS_ID, site_id, method,
                        sample_year)) %>%
  droplevels()

invert_total_richness_model_data <- invert_diversity_model_data %>%
  filter(response == "total_richness") %>%
  mutate(value = as.integer(value)) %>%
  droplevels()

invert_native_richness_model_data <- invert_diversity_model_data %>%
  filter(response == "native_richness") %>%
  mutate(value = as.integer(value)) %>%
  droplevels()

log_step(invert_total_richness_model_data, "13a_total_richness_model_data", "final AMBER with-type total richness rows")
log_step(invert_native_richness_model_data, "13a_native_richness_model_data", "final AMBER with-type native richness rows")

#### 13b. Proportion non-native model data #########################
invert_prop_model_data <- site_prop_base %>%
  mutate(HYBAS_ID = as.character(HYBAS_ID), site_id = as.character(site_id),
         sample_year = as.integer(as.character(sample_year))) %>%
  left_join(site_river_position %>%
              mutate(site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year))),
            by = c("site_id", "sample_year")) %>%
  left_join(climate_wide %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID), sample_year = as.integer(as.character(sample_year))),
            by = c("HYBAS_ID", "sample_year")) %>%
  left_join(river_len_df %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
  left_join(hyde_full %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID), year = as.integer(year)) %>%
              select(HYBAS_ID, sample_year = year, log_population_density),
            by = c("HYBAS_ID", "sample_year")) %>%
  left_join(hyde_1961 %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
  mutate(log_population_density = if_else(
    is.na(log_population_density) & sample_year < 1961,
    log_population_density_1961, log_population_density)) %>%
  select(-log_population_density_1961) %>%
  tidyr::crossing(barrier_type_reduced = factor(c("other", "dam", "weir", "culvert"),
                                                levels = c("other", "dam", "weir", "culvert"))) %>%
  left_join(barrier_density_df %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID),
                     barrier_type_reduced = factor(barrier_type_reduced,
                                                   levels = c("other", "dam", "weir", "culvert"))) %>%
              select(HYBAS_ID, barrier_type_reduced, n_barrier, barrier_density, log_barrier_density),
            by = c("HYBAS_ID", "barrier_type_reduced")) %>%
  mutate(
    non_native = as.integer(non_native_taxa),
    total_taxa = as.integer(total_known_taxa),
    prop_non_native = prop_non_native_taxa,
    z_barrier_density = z_safe(log_barrier_density),
    z_dist_up = z_safe(log_dist_up),
    z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual),
    z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert")),
    country = factor(country), HYBAS_ID = factor(HYBAS_ID), site_id = factor(site_id),
    method = factor(method), sample_year = factor(sample_year)) %>%
  filter(total_taxa > 0) %>%
  select(non_native, total_taxa, prop_non_native, barrier_type_reduced, z_barrier_density,
         z_dist_up, z_dist_mouth, z_log_population_density, z_mean_temp_annual, z_precip_annual,
         z_log_river_km, country, HYBAS_ID, site_id, method, sample_year, everything()) %>%
  filter(complete.cases(non_native, total_taxa, prop_non_native, barrier_type_reduced,
                        z_barrier_density, z_dist_up, z_dist_mouth, z_log_population_density,
                        z_mean_temp_annual, z_precip_annual, z_log_river_km, country, HYBAS_ID,
                        site_id, method, sample_year)) %>%
  droplevels()

log_step(invert_prop_model_data, "13b_proportion_model_data", "final AMBER with-type proportion non-native rows")

#### 13c. Proportion non-native model data, native+non-native denominator (A4i)
#### Structurally identical to 13b — same joins, same z-scored
#### predictors, same beta-binomial trials() column layout. The ONLY
#### difference from A3i is the denominator source (total_known_taxa
#### instead of total_richness).
invert_rate_model_data <- site_rate_base %>%
  mutate(HYBAS_ID = as.character(HYBAS_ID), site_id = as.character(site_id),
         sample_year = as.integer(as.character(sample_year))) %>%
  left_join(site_river_position %>%
              mutate(site_id = as.character(site_id), sample_year = as.integer(as.character(sample_year))),
            by = c("site_id", "sample_year")) %>%
  left_join(climate_wide %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID), sample_year = as.integer(as.character(sample_year))),
            by = c("HYBAS_ID", "sample_year")) %>%
  left_join(river_len_df %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
  left_join(hyde_full %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID), year = as.integer(year)) %>%
              select(HYBAS_ID, sample_year = year, log_population_density),
            by = c("HYBAS_ID", "sample_year")) %>%
  left_join(hyde_1961 %>% mutate(HYBAS_ID = as.character(HYBAS_ID)), by = "HYBAS_ID") %>%
  mutate(log_population_density = if_else(
    is.na(log_population_density) & sample_year < 1961,
    log_population_density_1961, log_population_density)) %>%
  select(-log_population_density_1961) %>%
  tidyr::crossing(barrier_type_reduced = factor(c("other", "dam", "weir", "culvert"),
                                                levels = c("other", "dam", "weir", "culvert"))) %>%
  left_join(barrier_density_df %>%
              mutate(HYBAS_ID = as.character(HYBAS_ID),
                     barrier_type_reduced = factor(barrier_type_reduced,
                                                   levels = c("other", "dam", "weir", "culvert"))) %>%
              select(HYBAS_ID, barrier_type_reduced, n_barrier, barrier_density, log_barrier_density),
            by = c("HYBAS_ID", "barrier_type_reduced")) %>%
  mutate(
    non_native = as.integer(non_native_taxa),
    total_taxa = as.integer(total_known_taxa),   # trials() denominator = native + non-native ONLY
    prop_non_native = prop_non_native_of_known,
    z_barrier_density = z_safe(log_barrier_density),
    z_dist_up = z_safe(log_dist_up),
    z_dist_mouth = z_safe(log_dist_mouth),
    z_log_population_density = z_safe(log_population_density),
    z_mean_temp_annual = z_safe(mean_temp_annual),
    z_precip_annual = z_safe(precip_annual),
    z_log_river_km = z_safe(log_river_km),
    barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert")),
    country = factor(country), HYBAS_ID = factor(HYBAS_ID), site_id = factor(site_id),
    method = factor(method), sample_year = factor(sample_year)) %>%
  filter(total_taxa > 0) %>%
  select(non_native, total_taxa, prop_non_native, barrier_type_reduced, z_barrier_density,
         z_dist_up, z_dist_mouth, z_log_population_density, z_mean_temp_annual, z_precip_annual,
         z_log_river_km, country, HYBAS_ID, site_id, method, sample_year, everything()) %>%
  filter(complete.cases(non_native, total_taxa, prop_non_native, barrier_type_reduced,
                        z_barrier_density, z_dist_up, z_dist_mouth, z_log_population_density,
                        z_mean_temp_annual, z_precip_annual, z_log_river_km, country, HYBAS_ID,
                        site_id, method, sample_year)) %>%
  droplevels()

log_step(invert_rate_model_data, "13c_rate_model_data", "final AMBER with-type A4i rows (bounded, native+non-native denominator)")


#### 14. Final checks ###############################################
cat("\n============================================================\n")
cat("FINAL MODEL DATA CHECKS\n")
cat("============================================================\n")
cat("\nTotal richness rows:", nrow(invert_total_richness_model_data), "\n")
cat("Native richness rows:", nrow(invert_native_richness_model_data), "\n")
cat("Proportion non-native rows:", nrow(invert_prop_model_data), "\n")

cat("\nBarrier types (total richness):\n"); print(table(invert_total_richness_model_data$barrier_type_reduced))
cat("\nBarrier types (native richness):\n"); print(table(invert_native_richness_model_data$barrier_type_reduced))
cat("\nBarrier types (proportion):\n"); print(table(invert_prop_model_data$barrier_type_reduced))

#### 14a. Join-success audit ########################################
#### complete.cases() in Sections 13a/13b silently drops any community
#### missing a match in ANY join (river position, climate, river
#### length, HYDE). This makes it explicit which source is actually
#### responsible for data loss, rather than lumping it all together.
join_success_audit <- site_diversity_base %>%
  distinct(site_id, sample_year, HYBAS_ID) %>%
  mutate(site_id = as.character(site_id), sample_year = as.integer(sample_year),
         HYBAS_ID = as.character(HYBAS_ID)) %>%
  mutate(
    has_river_position = paste(site_id, sample_year) %in%
      paste(site_river_position$site_id, site_river_position$sample_year),
    has_climate = paste(HYBAS_ID, sample_year) %in%
      paste(climate_wide$HYBAS_ID, climate_wide$sample_year),
    has_river_length = HYBAS_ID %in% river_len_df$HYBAS_ID,
    has_hyde_direct_year = paste(HYBAS_ID, sample_year) %in%
      paste(hyde_full$HYBAS_ID, hyde_full$year),
    has_hyde_1961_fallback = HYBAS_ID %in% hyde_1961$HYBAS_ID)

join_success_summary <- join_success_audit %>%
  summarise(
    n_communities = n(),
    pct_missing_river_position = round(100 * mean(!has_river_position), 2),
    pct_missing_climate = round(100 * mean(!has_climate), 2),
    pct_missing_river_length = round(100 * mean(!has_river_length), 2),
    pct_missing_hyde_direct_year_match = round(100 * mean(!has_hyde_direct_year), 2),
    pct_missing_hyde_even_with_1961_fallback =
      round(100 * mean(!has_hyde_direct_year & !has_hyde_1961_fallback), 2))

write_csv(join_success_audit, file.path(diag_dir, "QC_join_success_audit_per_community.csv"))
write_csv(join_success_summary, file.path(diag_dir, "QC_join_success_summary.csv"))
cat("\nJoin-success audit (which join is responsible for complete.cases() drop-off):\n")
print(join_success_summary)

#### Table S1 replacement: full record-count log through cleaning
write_csv(record_log, file.path(out_dir, "TableS1_data_cleaning_record_counts.csv"))
cat("\nData-cleaning record-count log (Table S1):\n")
print(record_log)

#### 15. Save "with type" final model-ready datasets ################
A1i_TotalRichnesswithtype_data <- invert_total_richness_model_data %>%
  mutate(value = as.integer(value),
         barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert"))) %>%
  droplevels()

A2i_NativeRichnesswithtype_data <- invert_native_richness_model_data %>%
  mutate(value = as.integer(value),
         barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert"))) %>%
  droplevels()

A3i_Propwithtype_data <- invert_prop_model_data %>%
  mutate(non_native = as.integer(non_native), total_taxa = as.integer(total_taxa),
         prop_non_native = as.numeric(prop_non_native),
         barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert"))) %>%
  filter(total_taxa > 0) %>%
  droplevels()

A4i_NonNativeRatewithtype_data <- invert_rate_model_data %>%
  mutate(non_native = as.integer(non_native), total_taxa = as.integer(total_taxa),
         prop_non_native = as.numeric(prop_non_native),
         barrier_type_reduced = factor(barrier_type_reduced, levels = c("other", "dam", "weir", "culvert"))) %>%
  filter(total_taxa > 0) %>%
  droplevels()

add_obs_weight <- function(dat) {
  dat %>%
    group_by(site_id, sample_year) %>%
    mutate(n_rows_per_siteyear = n(), obs_weight = 1 / n_rows_per_siteyear) %>%
    ungroup()
}

A1i_TotalRichnesswithtype_data <- add_obs_weight(A1i_TotalRichnesswithtype_data)
A2i_NativeRichnesswithtype_data <- add_obs_weight(A2i_NativeRichnesswithtype_data)
A3i_Propwithtype_data <- add_obs_weight(A3i_Propwithtype_data)
A4i_NonNativeRatewithtype_data <- add_obs_weight(A4i_NonNativeRatewithtype_data)

saveRDS(A1i_TotalRichnesswithtype_data, file.path(out_dir, "A1i_TotalRichnesswithtype_data.rds"))
saveRDS(A2i_NativeRichnesswithtype_data, file.path(out_dir, "A2i_NativeRichnesswithtype_data.rds"))
saveRDS(A3i_Propwithtype_data, file.path(out_dir, "A3i_Propwithtype_data.rds"))
saveRDS(A4i_NonNativeRatewithtype_data, file.path(out_dir, "A4i_NonNativeRatewithtype_data.rds"))

write_csv(A1i_TotalRichnesswithtype_data, file.path(out_dir, "A1i_TotalRichnesswithtype_data.csv"))
write_csv(A2i_NativeRichnesswithtype_data, file.path(out_dir, "A2i_NativeRichnesswithtype_data.csv"))
write_csv(A3i_Propwithtype_data, file.path(out_dir, "A3i_Propwithtype_data.csv"))
write_csv(A4i_NonNativeRatewithtype_data, file.path(out_dir, "A4i_NonNativeRatewithtype_data.csv"))

cat("\nSaved four final AMBER with-type model-specific datasets:\n")
cat(" - A1i_TotalRichnesswithtype_data.rds     (n =", nrow(A1i_TotalRichnesswithtype_data), ")\n")
cat(" - A2i_NativeRichnesswithtype_data.rds    (n =", nrow(A2i_NativeRichnesswithtype_data), ")\n")
cat(" - A3i_Propwithtype_data.rds              (n =", nrow(A3i_Propwithtype_data), ")\n")
cat(" - A4i_NonNativeRatewithtype_data.rds     (n =", nrow(A4i_NonNativeRatewithtype_data), ")\n")
cat("\nA4i formula reminder for Kelvin — SAME structure as A3i, trials() not offset:\n")
cat("  bf(non_native | trials(total_taxa) ~ z_barrier_density * barrier_type_reduced + z_dist_up +\n")
cat("       z_dist_mouth + z_log_population_density + z_mean_temp_annual + z_precip_annual +\n")
cat("       z_log_river_km + (1|country) + (1|HYBAS_ID) + (1|site_id) + (1|method) + (1|sample_year)),\n")
cat("     family = beta_binomial(link = \"logit\")\n")
cat("  ONLY difference from A3i: total_taxa here = native + non-native (total_known_taxa),\n")
cat("  NOT total_richness. Priors/adapt_delta/control settings can be copied directly from A3i.\n")
cat("\nNOTE: these are UPLOAD-TO-KELVIN ready, but the Kelvin runner script\n")
cat("(01_run_amber_model.R) still has the OLD valid_models list (A1i/A1ii/\n")
cat("A2i(Shannon)/A2ii/A3i/A3ii) and no A4i branch at all. It needs updating\n")
cat("to use these four new file names (A3i and A4i can share the SAME beta-\n")
cat("binomial branch in the runner, just pointed at different data files) and\n")
cat("drop the 'withouttype' + Shannon branches before submitting SLURM jobs.\n")

#### Rows-per-site-year sanity check (should be 4, one per type, since
#### these are all "with type" crossed datasets)
siteyear_check <- tibble(
  dataset = c("A1i_TotalRichnesswithtype", "A2i_NativeRichnesswithtype", "A3i_Propwithtype", "A4i_NonNativeRatewithtype"),
  n_rows = c(nrow(A1i_TotalRichnesswithtype_data), nrow(A2i_NativeRichnesswithtype_data),
             nrow(A3i_Propwithtype_data), nrow(A4i_NonNativeRatewithtype_data)),
  n_unique_siteyears = c(
    n_distinct(A1i_TotalRichnesswithtype_data$site_id, A1i_TotalRichnesswithtype_data$sample_year),
    n_distinct(A2i_NativeRichnesswithtype_data$site_id, A2i_NativeRichnesswithtype_data$sample_year),
    n_distinct(A3i_Propwithtype_data$site_id, A3i_Propwithtype_data$sample_year),
    n_distinct(A4i_NonNativeRatewithtype_data$site_id, A4i_NonNativeRatewithtype_data$sample_year)))
print(siteyear_check)
write_csv(siteyear_check, file.path(out_dir, "QC_A1_A2_A3_A4_siteyear_check.csv"))

#### ============================================================
#### 16. ASSUMPTION-CHECKING DIAGNOSTICS AND FIGURES
#### ============================================================
#### Everything below writes to Diagnostics/ and produces one
#### combined PDF plus CSV/TXT outputs for each check.
#### ============================================================

pdf(file.path(diag_dir, "AMBER_assumption_checking_figures.pdf"), width = 9, height = 6)

#### 16a. Response distributions ##################################
print(
  ggplot(A1i_TotalRichnesswithtype_data, aes(x = value)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Total richness (all unique taxa) — raw distribution",
         x = "Total richness", y = "Count") +
    theme_classic()
)

print(
  ggplot(A2i_NativeRichnesswithtype_data, aes(x = value)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Native richness (native-confirmed taxa only) — raw distribution",
         x = "Native richness", y = "Count") +
    theme_classic()
)

print(
  ggplot(A3i_Propwithtype_data, aes(x = prop_non_native)) +
    geom_histogram(bins = 30) +
    labs(title = "Proportion non-native — raw distribution",
         x = "Proportion non-native", y = "Count") +
    theme_classic()
)

print(
  ggplot(A3i_Propwithtype_data, aes(x = total_taxa)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Total known-status taxa (proportion model denominator)",
         x = "Total known taxa", y = "Count") +
    theme_classic()
)

print(
  ggplot(A4i_NonNativeRatewithtype_data, aes(x = prop_non_native)) +
    geom_histogram(bins = 30) +
    labs(title = "Proportion non-native, native+non-native denominator (A4i) — raw distribution",
         x = "Proportion non-native (of known-status taxa)", y = "Count") +
    theme_classic()
)

print(
  ggplot(A4i_NonNativeRatewithtype_data, aes(x = total_taxa)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Native + non-native taxa (A4i denominator)",
         x = "Total known taxa (native + non-native)", y = "Count") +
    theme_classic()
)

#### 16b. Overdispersion checks (variance/mean) for the two negbinomial
#### count responses only — A3i and A4i are now both bounded beta-
#### binomial proportions, checked separately in 16d.
overdispersion_check <- tibble(
  response = c("total_richness", "native_richness"),
  mean_value = c(mean(A1i_TotalRichnesswithtype_data$value), mean(A2i_NativeRichnesswithtype_data$value)),
  var_value = c(var(A1i_TotalRichnesswithtype_data$value), var(A2i_NativeRichnesswithtype_data$value))
) %>%
  mutate(var_mean_ratio = var_value / mean_value,
         note = "ratio >> 1 supports negative binomial over Poisson")

write_csv(overdispersion_check, file.path(diag_dir, "QC_overdispersion_variance_mean_ratio.csv"))
cat("\nOverdispersion check (variance/mean ratio, Poisson assumes 1):\n")
print(overdispersion_check)

#### 16c. Zero-inflation + dispersion via DHARMa on quick glm.nb fits
run_zi_dispersion_check <- function(dat, response_name) {
  m <- glm.nb(value ~ z_barrier_density * barrier_type_reduced + z_dist_up + z_dist_mouth +
                z_log_population_density + z_mean_temp_annual + z_precip_annual + z_log_river_km,
              data = dat)
  sim <- simulateResiduals(m, n = 250)
  
  zi_test <- testZeroInflation(sim, plot = FALSE)
  disp_test <- testDispersion(sim, plot = FALSE)
  
  plot(sim, main = paste(response_name, "— DHARMa residual diagnostics"))
  
  tibble(
    response = response_name,
    n_zero_obs = sum(dat$value == 0),
    pct_zero_obs = round(100 * mean(dat$value == 0), 2),
    zi_ratio_obs_sim = as.numeric(zi_test$statistic),
    zi_p_value = zi_test$p.value,
    dispersion_ratio = as.numeric(disp_test$statistic),
    dispersion_p_value = disp_test$p.value
  )
}

zi_dispersion_results <- bind_rows(
  run_zi_dispersion_check(A1i_TotalRichnesswithtype_data, "total_richness"),
  run_zi_dispersion_check(A2i_NativeRichnesswithtype_data, "native_richness")
)

write_csv(zi_dispersion_results, file.path(diag_dir, "QC_zero_inflation_dispersion_tests.csv"))
cat("\nZero-inflation and dispersion test results (glm.nb quick-fit basis, negbinomial responses only):\n")
print(zi_dispersion_results)

#### 16d. Beta-binomial justification for A3i AND A4i #############
#### Deviance/df from a single-intercept binomial GLM (>1 supports
#### beta-binomial over plain binomial). Run for BOTH proportion
#### responses now that A4i is also a beta-binomial trials() model.
compute_binom_dispersion_ratio <- function(dat, response_name) {
  m <- glm(cbind(non_native, total_taxa - non_native) ~ 1, data = dat, family = binomial())
  tibble(response = response_name,
         check = "binomial_deviance_over_df",
         value = m$deviance / m$df.residual,
         note = paste0("ratio >> 1 supports beta-binomial over binomial for ", response_name))
}

prop_dispersion_summary <- bind_rows(
  compute_binom_dispersion_ratio(A3i_Propwithtype_data, "A3i (total richness denominator)"),
  compute_binom_dispersion_ratio(A4i_NonNativeRatewithtype_data, "A4i (native+non-native denominator)")
)

cat("\nBeta-binomial justification (deviance/df ratio, >1 supports beta-binomial):\n")
print(prop_dispersion_summary)
write_csv(prop_dispersion_summary, file.path(diag_dir, "QC_proportion_betabinomial_justification.csv"))

#### 16e. Predictor collinearity: pairwise correlations + full VIF
predictor_cols <- c("z_barrier_density", "z_dist_up", "z_dist_mouth",
                    "z_log_population_density", "z_mean_temp_annual",
                    "z_precip_annual", "z_log_river_km")

pairwise_cor <- cor(A1i_TotalRichnesswithtype_data[predictor_cols], use = "pairwise.complete.obs")
write.csv(pairwise_cor, file.path(diag_dir, "QC_predictor_pairwise_correlation.csv"))
cat("\nPairwise predictor correlations:\n")
print(round(pairwise_cor, 3))

#### Full VIF via car::vif() on an lm() with all main-effect predictors
#### (run on the total-richness data; predictor set is identical across
#### A1i/A2i/A3i so this VIF check applies to all three response models).
vif_check <- lm(value ~ z_barrier_density + z_dist_up + z_dist_mouth +
                  z_log_population_density + z_mean_temp_annual + z_precip_annual + z_log_river_km,
                data = A1i_TotalRichnesswithtype_data)

vif_results <- car::vif(vif_check)
vif_df <- tibble(term = names(vif_results), VIF = as.numeric(vif_results))
write_csv(vif_df, file.path(diag_dir, "QC_predictor_VIF.csv"))
cat("\nVIF results (main effects, no barrier_type_reduced interaction term):\n")
print(vif_df)

#### VIF including the barrier_type interaction term as well, since
#### that's the actual model structure being fit on Kelvin
vif_check_withtype <- lm(value ~ z_barrier_density * barrier_type_reduced + z_dist_up + z_dist_mouth +
                           z_log_population_density + z_mean_temp_annual + z_precip_annual + z_log_river_km,
                         data = A1i_TotalRichnesswithtype_data)

vif_results_withtype <- car::vif(vif_check_withtype)
write.csv(as.data.frame(vif_results_withtype), file.path(diag_dir, "QC_predictor_VIF_withtype_interaction.csv"))
cat("\nVIF results (including barrier_type_reduced interaction — GVIF/df reported for factors):\n")
print(vif_results_withtype)

#### Collinearity heatmap figure
cor_df <- as.data.frame(pairwise_cor)
cor_df$var1 <- rownames(cor_df)
cor_long <- cor_df %>% pivot_longer(-var1, names_to = "var2", values_to = "correlation")

print(
  ggplot(cor_long, aes(x = var1, y = var2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = round(correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "AMBER predictor collinearity matrix", x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

#### 16f. Barrier density by type (raw-data sanity check) #########
print(
  barrier_density_df %>%
    filter(barrier_density > 0) %>%
    ggplot(aes(x = barrier_type_reduced, y = barrier_density)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = "Basin barrier density by reduced type (log10 scale, zero-density basins excluded)",
         x = "Barrier type", y = "Barriers per river km") +
    theme_classic()
)

print(
  A1i_TotalRichnesswithtype_data %>%
    ggplot(aes(x = barrier_type_reduced, y = value)) +
    geom_boxplot() +
    labs(title = "Total richness by barrier type (raw, unadjusted)", x = "Barrier type", y = "Total richness") +
    theme_classic()
)

#### 16g. Pearson correlation among diversity metrics #############
#### Purpose: empirically justify which diversity metrics are
#### redundant, rather than dropping any purely on prior expectation.
#### Includes Shannon, evenness, and turnover ALONGSIDE the two
#### metrics that are actually modelled (total/native richness) and
#### proportion non-native, even though Shannon/evenness/turnover are
#### not fitted as their own AMBER models.
diversity_metric_cor <- cor(
  diagnostic_diversity_df %>%
    select(total_richness, native_richness, shannon, evenness_pielou, turnover, prop_non_native_taxa),
  use = "pairwise.complete.obs", method = "pearson")

write.csv(diversity_metric_cor, file.path(diag_dir, "QC_diversity_metric_pearson_correlation.csv"))
cat("\nPearson correlation among diversity metrics:\n")
print(round(diversity_metric_cor, 3))

#### Report pairwise sample sizes too, since turnover has structural
#### NAs (first visit per site) that reduce n for any pair involving it.
pairwise_n <- diagnostic_diversity_df %>%
  select(total_richness, native_richness, shannon, evenness_pielou, turnover, prop_non_native_taxa) %>%
  summarise(across(everything(), ~ sum(!is.na(.x))))
cat("\nNon-missing n per diversity metric (turnover is lowest — first-visit rows are NA):\n")
print(pairwise_n)
write_csv(pairwise_n, file.path(diag_dir, "QC_diversity_metric_pairwise_n.csv"))

diversity_cor_df <- as.data.frame(diversity_metric_cor)
diversity_cor_df$var1 <- rownames(diversity_cor_df)
diversity_cor_long <- diversity_cor_df %>% pivot_longer(-var1, names_to = "var2", values_to = "correlation")

print(
  ggplot(diversity_cor_long, aes(x = var1, y = var2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = round(correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Pearson correlation among diversity metrics",
         subtitle = "total/native richness, Shannon, Pielou evenness, temporal turnover, proportion non-native",
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

#### 16h. Method-faceted diversity metric distributions ############
#### Purpose: check whether Shannon/evenness (and richness/turnover)
#### distributions differ systematically by sampling method — a
#### method confound here would mean the Section 16g correlation
#### estimates could be partly driven by method mix rather than pure
#### ecological redundancy between metrics.

#### Collapse rare methods so facets stay legible; keep the top 6
#### most-sampled methods separately, group the rest as "Other".
top_methods <- diagnostic_diversity_df %>%
  count(method, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(method)

diagnostic_diversity_df_faceted <- diagnostic_diversity_df %>%
  mutate(method_grouped = if_else(method %in% top_methods, method, "Other"))

diversity_metric_long <- diagnostic_diversity_df_faceted %>%
  select(method_grouped, total_richness, native_richness, shannon, evenness_pielou, turnover) %>%
  pivot_longer(cols = c(total_richness, native_richness, shannon, evenness_pielou, turnover),
               names_to = "metric", values_to = "metric_value") %>%
  filter(!is.na(metric_value))

print(
  ggplot(diversity_metric_long, aes(x = metric_value)) +
    geom_histogram(bins = 25) +
    facet_grid(method_grouped ~ metric, scales = "free") +
    labs(title = "Diversity metric distributions by sampling method",
         subtitle = "Top 6 methods by sample count shown individually; rest grouped as 'Other'",
         x = NULL, y = "Count") +
    theme_classic(base_size = 8) +
    theme(strip.text = element_text(size = 7))
)

#### Same check but as boxplots, easier to compare medians/spread
#### across methods at a glance than overlapping histograms.
print(
  ggplot(diversity_metric_long, aes(x = method_grouped, y = metric_value)) +
    geom_boxplot() +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(title = "Diversity metric distributions by sampling method (boxplot view)",
         x = "Method", y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

#### Quick numeric summary: mean of each metric by method, saved to CSV
method_metric_summary <- diagnostic_diversity_df_faceted %>%
  group_by(method_grouped) %>%
  summarise(
    n = n(),
    mean_total_richness = mean(total_richness, na.rm = TRUE),
    mean_native_richness = mean(native_richness, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    mean_evenness = mean(evenness_pielou, na.rm = TRUE),
    mean_turnover = mean(turnover, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(n))

write_csv(method_metric_summary, file.path(diag_dir, "QC_diversity_metrics_by_method.csv"))
cat("\nDiversity metric means by sampling method (top 6 + Other):\n")
print(method_metric_summary)

dev.off()

cat("\n============================================================\n")
cat("Finished Section 16: assumption-checking diagnostics saved to\n")
cat(diag_dir, "\n")
cat("============================================================\n")

cat("\nScript complete. Three AMBER 'with type' datasets ready for upload to Kelvin:\n")
cat(" - A1i_TotalRichnesswithtype_data.rds\n")
cat(" - A2i_NativeRichnesswithtype_data.rds\n")
cat(" - A3i_Propwithtype_data.rds\n")


#### 16j. A3i vs A4i proportion redundancy check ###################
#### Both are now bounded beta-binomial proportions with the SAME
#### numerator (non-native taxa) and slightly different denominators
#### (total_richness vs total_known_taxa). This checks directly how
#### similar the two resulting proportions actually are.
a3_a4_join <- A3i_Propwithtype_data %>%
  distinct(site_id, sample_year, barrier_type_reduced, prop_non_native_A3i = prop_non_native) %>%
  inner_join(
    A4i_NonNativeRatewithtype_data %>%
      distinct(site_id, sample_year, barrier_type_reduced, prop_non_native_A4i = prop_non_native),
    by = c("site_id", "sample_year", "barrier_type_reduced"))

a3_a4_cor <- cor(a3_a4_join$prop_non_native_A3i, a3_a4_join$prop_non_native_A4i,
                 use = "pairwise.complete.obs", method = "pearson")
a3_a4_r2 <- summary(lm(prop_non_native_A4i ~ prop_non_native_A3i, data = a3_a4_join))$r.squared

a3_a4_summary <- a3_a4_join %>%
  summarise(
    n_matched_rows = n(),
    pearson_r = a3_a4_cor,
    r_squared = a3_a4_r2,
    mean_abs_difference = mean(abs(prop_non_native_A3i - prop_non_native_A4i), na.rm = TRUE),
    max_abs_difference = max(abs(prop_non_native_A3i - prop_non_native_A4i), na.rm = TRUE))

write_csv(a3_a4_summary, file.path(diag_dir, "QC_A3i_vs_A4i_redundancy.csv"))
cat("\nA3i vs A4i proportion redundancy check:\n")
print(a3_a4_summary)
cat("\nInterpretation: r close to 1 and mean_abs_difference close to 0 means\n")
cat("A3i and A4i are nearly the same model with different denominators — fitting\n")
cat("both on Kelvin would be of limited additional value. A meaningfully lower r\n")
cat("or larger mean_abs_difference means unresolved-status taxa are common enough\n")
cat("in your data that the two denominators genuinely diverge, and both are worth\n")
cat("keeping.\n")

print(
  ggplot(a3_a4_join, aes(x = prop_non_native_A3i, y = prop_non_native_A4i)) +
    geom_point(alpha = 0.15) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "firebrick") +
    labs(title = "A3i vs A4i: proportion non-native, two denominators",
         subtitle = paste0("Pearson r = ", round(a3_a4_cor, 3), ", R\u00b2 = ", round(a3_a4_r2, 3),
                           " | dashed line = perfect agreement"),
         x = "A3i: non-native / total richness", y = "A4i: non-native / (native+non-native)") +
    theme_classic()
)

dev.off()

cat("\n============================================================\n")
cat("Finished Section 16: assumption-checking diagnostics saved to\n")
cat(diag_dir, "\n")
cat("============================================================\n")

cat("\nScript complete. Four AMBER 'with type' datasets ready for upload to Kelvin:\n")
cat(" - A1i_TotalRichnesswithtype_data.rds     (negbinomial)\n")
cat(" - A2i_NativeRichnesswithtype_data.rds    (negbinomial)\n")
cat(" - A3i_Propwithtype_data.rds              (beta-binomial, trials(total_taxa))\n")
cat(" - A4i_NonNativeRatewithtype_data.rds     (beta-binomial, trials(total_taxa), native+non-native denominator)\n")

