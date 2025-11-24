library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

mds_df <- read_delim("/Users/jaimeortiz/Library/CloudStorage/Box-Box/0_PhD_NINA/WildGenome/SEA_CUCUMBER/GAL_MAINLAND_Analysis/lostruct_per_chrom/Chr01/Chr01_MDS_results.tsv", delim = "\t")
chr <- "Chr01"
chrom <- "Chr01"

# ============================================================================
# LOSTRUCT OUTLIER REGION DETECTION AND VISUALIZATION
# ============================================================================
# This script identifies outlier genomic regions using Local PCA (lostruct)
# and creates a multi-panel visualization showing MDS axes across chromosomes
# Following methodology from Huang et al. (2020) and related studies
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# ----------------------------------------------------------------------------
# PARAMETERS
# ----------------------------------------------------------------------------
z_thr <- 4              # Standard deviation threshold for outlier detection (±4 SD)
gap_threshold <- 20     # Maximum number of windows between clusters to merge them
min_windows <- 5        # Minimum number of windows required for a cluster to be retained
max_MDS <- 10          # Number of MDS axes to analyze (MDS1 through MDS10)

# ----------------------------------------------------------------------------
# VALIDATE INPUT DATA
# ----------------------------------------------------------------------------
# Ensure the input dataframe (mds_df) contains required columns:
# - chrom: chromosome/linkage group identifier
# - start, end, mid: genomic coordinates of windows
# - window: sequential window index
# - MDS1, MDS2, ..., MDS10: MDS values for each axis
stopifnot(all(c("chrom","start","end","mid","window") %in% names(mds_df)))

# ----------------------------------------------------------------------------
# FUNCTION: CLUSTER OUTLIERS FOR ONE MDS AXIS AND ONE TAIL
# ----------------------------------------------------------------------------
# This function identifies and clusters outlier windows for either the
# positive tail (above +4 SD) or negative tail (below -4 SD) of one MDS axis
#
# Arguments:
#   df: dataframe with MDS scores and genomic positions
#   mds_col: name of the MDS column to analyze (e.g., "MDS1")
#   tail: either "positive" (upper outliers) or "negative" (lower outliers)
#   gap_threshold: max windows between clusters to merge
#   min_windows: minimum cluster size to retain
#   mds_idx: numeric index of MDS axis (1-10)
#
# Returns:
#   dataframe of clustered outlier regions with genomic coordinates
# ----------------------------------------------------------------------------
cluster_one_tail <- function(df, mds_col, tail = c("positive","negative"),
                             gap_threshold = 20, min_windows = 5, mds_idx = NA_integer_) {
  # Validate tail argument
  tail <- match.arg(tail)
  
  # Extract MDS values for this axis
  v <- df[[mds_col]]
  
  # Calculate mean and standard deviation for threshold determination
  m <- mean(v, na.rm = TRUE)
  s <- sd(v, na.rm = TRUE)
  
  # Define upper and lower thresholds (mean ± z_thr * SD)
  upper <- m + z_thr * s
  lower <- m - z_thr * s
  
  # Identify outlier windows based on tail direction
  # Positive: windows >= upper threshold
  # Negative: windows <= lower threshold
  idx_rows <- if (tail == "positive") which(v >= upper) else which(v <= lower)
  
  # Return NULL if no outliers found
  if (length(idx_rows) == 0) return(NULL)
  
  # Extract outlier windows with their genomic coordinates
  out <- df[idx_rows, c("chrom","start","end","mid","window")]
  
  # Sort by window index to enable clustering by adjacency
  out <- out[order(out$window), ]
  
  # ---- INITIALIZE FIRST CLUSTER ----
  # Start with the first outlier window
  res <- out[1,]
  names(res) <- c("LG","start","stop","mid","window_start")
  res$window_stop <- out$window[1]
  res$n_windows <- 1L
  res$mds <- paste0(tail, mds_idx)  # e.g., "positive1" or "negative3"
  res$mds_num <- if (tail == "positive") mds_idx else -mds_idx  # positive: 1-10, negative: -1 to -10
  res$cluster <- 1L
  cidx <- 1L  # Current cluster index
  
  # ---- ITERATE THROUGH REMAINING OUTLIERS TO BUILD CLUSTERS ----
  if (nrow(out) > 1) {
    for (i in 2:nrow(out)) {
      # Check if current window is within gap_threshold of previous window
      # If yes: extend current cluster
      # If no: start new cluster
      if (out$window[i] <= out$window[i-1] + 1 + gap_threshold) {
        # EXTEND CURRENT CLUSTER
        # Update the stop position and window count
        res$window_stop[cidx] <- out$window[i]
        res$stop[cidx] <- out$end[i]
        res$n_windows[cidx] <- res$n_windows[cidx] + 1L
      } else {
        # START NEW CLUSTER
        cidx <- cidx + 1L
        add <- out[i,]
        names(add) <- c("LG","start","stop","mid","window_start")
        add$window_stop <- out$window[i]
        add$n_windows <- 1L
        add$mds <- paste0(tail, mds_idx)
        add$mds_num <- if (tail == "positive") mds_idx else -mds_idx
        add$cluster <- cidx
        res <- rbind(res, add)
      }
    }
  }
  
  # ---- FILTER BY MINIMUM CLUSTER SIZE ----
  # Retain only clusters with at least min_windows windows
  res <- res[res$n_windows >= min_windows, , drop = FALSE]
  
  # Return NULL if all clusters were too small
  if (!nrow(res)) return(NULL)
  
  return(res)
}

# ----------------------------------------------------------------------------
# BUILD CLUSTER DATABASE ACROSS ALL MDS AXES
# ----------------------------------------------------------------------------
# Initialize empty dataframe to store all clusters
MDS_CLUSTER_ALL <- NULL

# Create list of MDS column names present in the data
mds_cols_present <- paste0("MDS", 1:max_MDS)
mds_cols_present <- mds_cols_present[mds_cols_present %in% names(mds_df)]

# Loop through each MDS axis
for (j in seq_along(mds_cols_present)) {
  mds_col <- mds_cols_present[j]
  
  # Detect POSITIVE outliers (above +4 SD)
  pos <- cluster_one_tail(mds_df, mds_col, "positive",
                          gap_threshold, min_windows, mds_idx = j)
  
  # Detect NEGATIVE outliers (below -4 SD)
  neg <- cluster_one_tail(mds_df, mds_col, "negative",
                          gap_threshold, min_windows, mds_idx = j)
  
  # Add clusters to master dataframe
  if (!is.null(pos)) MDS_CLUSTER_ALL <- rbind(MDS_CLUSTER_ALL, pos)
  if (!is.null(neg)) MDS_CLUSTER_ALL <- rbind(MDS_CLUSTER_ALL, neg)
}

# Handle case where no clusters were found
if (is.null(MDS_CLUSTER_ALL)) {
  MDS_CLUSTER_ALL <- data.frame(
    LG=character(), start=numeric(), stop=numeric(), mid=numeric(),
    window_start=integer(), window_stop=integer(), n_windows=integer(),
    mds=character(), mds_num=integer(), cluster=integer(), stringsAsFactors = FALSE
  )
}

# Calculate physical size of each cluster in base pairs
MDS_CLUSTER_ALL$cluster_size <- MDS_CLUSTER_ALL$stop - MDS_CLUSTER_ALL$start

# Write individual cluster file (one row per cluster per MDS axis)
write_delim(MDS_CLUSTER_ALL,
            paste0("cluster_outlier_allMDS_max", length(mds_cols_present),
                   "_merge", gap_threshold, "_filter", min_windows,
                   "_sdLim", z_thr, ".txt"),
            delim = "\t")

cat("\n=== Individual clusters identified ===\n")
cat("Total clusters:", nrow(MDS_CLUSTER_ALL), "\n")
print(table(MDS_CLUSTER_ALL$mds))

# ----------------------------------------------------------------------------
# FUNCTION: MERGE OVERLAPPING REGIONS ACROSS MDS AXES
# ----------------------------------------------------------------------------
# For each chromosome, merge regions that overlap in genomic position
# across different MDS axes. This ensures vertical lines span regions
# detected by ANY MDS axis at that genomic location.
#
# Arguments:
#   regions_df: dataframe of regions for one chromosome
#
# Returns:
#   dataframe of merged regions with extended boundaries
# ----------------------------------------------------------------------------
merge_overlapping_regions <- function(regions_df) {
  if (nrow(regions_df) == 0) return(NULL)
  
  # Sort regions by start position
  regions_df <- regions_df[order(regions_df$start), ]
  
  # Initialize merged regions with first region
  merged <- regions_df[1, ]
  merged$merged_region_id <- 1
  current_idx <- 1
  
  # Iterate through remaining regions
  if (nrow(regions_df) > 1) {
    for (i in 2:nrow(regions_df)) {
      # Check if current region overlaps with the last merged region
      if (regions_df$start[i] <= merged$stop[current_idx]) {
        # OVERLAPPING: extend the merged region to encompass both
        # Use the maximum stop position to get largest boundaries
        merged$stop[current_idx] <- max(merged$stop[current_idx], regions_df$stop[i])
        merged$n_windows[current_idx] <- merged$n_windows[current_idx] + regions_df$n_windows[i]
        merged$cluster_size[current_idx] <- merged$stop[current_idx] - merged$start[current_idx]
      } else {
        # NOT OVERLAPPING: start a new merged region
        current_idx <- current_idx + 1
        new_region <- regions_df[i, ]
        new_region$merged_region_id <- current_idx
        merged <- rbind(merged, new_region)
      }
    }
  }
  
  return(merged)
}

# ----------------------------------------------------------------------------
# MERGE OVERLAPPING REGIONS PER CHROMOSOME
# ----------------------------------------------------------------------------
# Apply merging function to each chromosome separately
if (nrow(MDS_CLUSTER_ALL) > 0) {
  merged_regions <- MDS_CLUSTER_ALL %>%
    group_by(LG) %>%
    group_modify(~ merge_overlapping_regions(.x)) %>%
    ungroup() %>%
    mutate(global_region_id = paste(LG, merged_region_id, sep = "_"))
  
  # Create boundary lines for plotting
  # Each merged region contributes TWO vertical lines (start and stop)
  boundary_lines_merged <- merged_regions %>%
    select(LG, start, stop, global_region_id) %>%
    distinct() %>%
    pivot_longer(c(start, stop), names_to = "boundary", values_to = "bp") %>%
    mutate(position_mb = bp / 1e6) %>%
    rename(chrom = LG)
  
} else {
  merged_regions <- data.frame()
  boundary_lines_merged <- data.frame()
}

cat("\n=== Merged regions across MDS axes ===\n")
if (nrow(merged_regions) > 0) {
  print(merged_regions %>% 
          select(LG, start, stop, cluster_size, merged_region_id, global_region_id) %>%
          arrange(LG, start))
}

# Write merged regions file
write_delim(merged_regions,
            paste0("merged_regions_allMDS_z", z_thr, ".txt"),
            delim = "\t")

# ----------------------------------------------------------------------------
# PREPARE DATA FOR PLOTTING
# ----------------------------------------------------------------------------
# Convert wide format (MDS1, MDS2, ... as columns) to long format
# (one row per window per MDS axis) for faceted plotting
mds_long <- mds_df %>%
  select(chrom, mid, all_of(mds_cols_present)) %>%
  pivot_longer(cols = all_of(mds_cols_present),
               names_to = "MDS_axis", values_to = "MDS_value")

# ---- FIX ORDERING ISSUE ----
# Convert MDS_axis to factor with proper numeric ordering
# This ensures MDS1, MDS2, ..., MDS9, MDS10 (not MDS1, MDS10, MDS2, ...)
mds_long <- mds_long %>%
  mutate(MDS_axis = factor(MDS_axis, 
                           levels = paste0("MDS", 1:max_MDS),
                           ordered = TRUE))

# Also fix ordering in boundary_lines_merged if regions exist
if (nrow(boundary_lines_merged) > 0) {
  # Extract numeric part from global_region_id for proper sorting
  boundary_lines_merged <- boundary_lines_merged %>%
    mutate(
      # Extract chromosome and numeric region ID
      region_num = as.numeric(sub(".*_", "", global_region_id)),
      # Create factor for proper ordering
      global_region_id = factor(global_region_id,
                                levels = unique(global_region_id[order(chrom, region_num)]),
                                ordered = TRUE)
    )
}

# ----------------------------------------------------------------------------
# CREATE FACETED GRID PLOT
# ----------------------------------------------------------------------------
# Grid layout: MDS axes (rows) × Chromosomes (columns)
# Vertical dotted lines mark merged outlier region boundaries
# Lines span all MDS rows within each chromosome column
p_grid <- ggplot(mds_long, aes(x = mid / 1e6, y = MDS_value)) +
  
  # ---- VERTICAL BOUNDARY LINES ----
# Draw dotted vertical lines at start/stop of each merged region
# These lines span all MDS rows due to facet_grid behavior
# Color distinguishes different regions
geom_vline(data = boundary_lines_merged,
           aes(xintercept = position_mb, color = global_region_id),
           linetype = "dotted", 
           linewidth = 0.8, 
           alpha = 0.7) +
  
  # ---- DATA POINTS ----
# Plot MDS values as black points
geom_point(size = 0.5, alpha = 0.55, color = "black") +
  
  # ---- FACETING ----
# Create grid with MDS axes as rows and chromosomes as columns
# scales = "free" allows each panel to have independent axis ranges
facet_grid(MDS_axis ~ chrom, 
           scales = "free",
           labeller = labeller(
             # Rename chromosomes from "Chr01" to "LG1"
             chrom = function(x) gsub("^Chr0?", "LG", x),
             # Rename MDS axes from "MDS1" to "mds1" (lowercase)
             MDS_axis = function(x) gsub("^MDS", "mds", x)
           )) +
  
  # ---- COLOR SCALE ----
# Discrete colors for different regions
scale_color_discrete(name = "Region") +
  
  # ---- THEME CUSTOMIZATION ----
theme_bw() +
  theme(
    # Remove background grid
    panel.grid = element_blank(),
    
    # Minimal spacing between panels
    panel.spacing = unit(0.06, "lines"),
    
    # Black borders around panels
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    # Font sizes for axes
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 10),
    
    # Facet strip formatting
    strip.text.x = element_text(size = 9, face = "bold"),  # Chromosome labels
    strip.text.y = element_text(size = 8, angle = 0),      # MDS axis labels
    strip.background = element_rect(fill = "white", color = "black"),
    
    # Legend formatting
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  ) +
  
  # ---- AXIS LABELS ----
labs(x = "position (Mb)", 
     y = "Value on each mds axis")

# Display plot
print(p_grid)

# Save plot as high-resolution PNG
ggsave("local_pca_grid_merged_z4.png", 
       plot = p_grid, 
       width = 14, 
       height = 18, 
       dpi = 300)

# ----------------------------------------------------------------------------
# SUMMARY OUTPUT
# ----------------------------------------------------------------------------
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Individual clusters detected:", nrow(MDS_CLUSTER_ALL), "\n")
cat("Merged regions across MDS axes:", 
    nrow(merged_regions %>% select(LG, merged_region_id) %>% distinct()), "\n")
cat("\nFiles created:\n")
cat("  1. cluster_outlier_allMDS_max", length(mds_cols_present), 
    "_merge", gap_threshold, "_filter", min_windows, 
    "_sdLim", z_thr, ".txt\n", sep="")
cat("  2. merged_regions_allMDS_z", z_thr, ".txt\n", sep="")
cat("  3. local_pca_grid_merged_z4.png\n")
