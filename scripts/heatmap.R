#############################################
# Sinkhole Metabolite-Taxonomy Correlation Heatmap
# Author: MDBoer
# Date: 2026-08-25
# Description: Plot correlation heatmaps with
#              metabolic and taxonomic annotations
#############################################


# ----------------------------------------------------------------------------
# 1. Auto-install & Load Dependencies
# ----------------------------------------------------------------------------
pkgs <- c("dplyr", "gridExtra", "grid", "glue", 
          "ComplexHeatmap", "circlize", "tidyr")
for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing", p, "..."))
    if (p %in% c("ComplexHeatmap", "circlize")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else {
      install.packages(p)
    }
  }
}
invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))


# ----------------------------------------------------------------------------
# 2. Configuration 
# ----------------------------------------------------------------------------

RANK <- 'phylum'  # Change taxonomy level as needed
SAMPLE_GROUPS <- c("Mesopelagic", "Acid Lake", "Deep Sinkhole", 
                   "Shallow Sinkhole", "Platform", "Surface")

# Color palette (consistent with previous analysis)
CLUSTER_COLORS <- c(
  'Surface'          = "#009E73",
  'Platform'         = "#ddcc77",
  'Shallow Sinkhole' = "#a2e7f2",
  'Deep Sinkhole'    = "#0072B2",
  'Acid Lake'        = "#E69F00",
  'Mesopelagic'      = "#CC79A7"
)

# Output directories
OUTPUT_DIR <- "../output"
FIGURE_DIR <- "../figures"
DATA_DIR   <- "../data"

# Directories for file paths
PATHS <- list(
  tax_abun = glue("{DATA_DIR}/mean_{RANK}_log.csv"),
  met_abun = "{DATA_DIR}/met_grouped_combined.csv",
  correlations = glue("{DATA_DIR}/{RANK}_correlations.csv"),
  components = "{DATA_DIR}/components.csv",
  group_data = glue("{DATA_DIR}/{RANK}_max_group_all.csv"),
)


# Number of replicate columns to add
N_REPLICATES <- 10


# ----------------------------------------------------------------------------
# 3. Helper Functions
# ----------------------------------------------------------------------------

rescale01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

create_color_gradient <- function(color, n = 50) {
  colorRampPalette(c("#FFFFFF", color))(n)
}

# ----------------------------------------------------------------------------
# 4. Load and Preprocess Data
# ----------------------------------------------------------------------------
cat("Loading data...\n")

tax_abun_data <- read.csv(PATHS$tax_abun, row.names = 1, check.names = FALSE)
met_abun_data <- read.csv(PATHS$met_abun, row.names = 1, check.names = FALSE)
hm_data       <- read.csv(PATHS$correlations, row.names = 1, check.names = FALSE)
comp_data     <- read.csv(PATHS$components, sep = '\t', row.names = 4)
group_data    <- read.csv(PATHS$group_data, row.names = 1)

cat(sprintf("Loaded: %d metabolites, %d taxa, %d correlation pairs\n",
            nrow(hm_data), nrow(tax_abun_data), ncol(hm_data)))

# Standardize row names
rownames(comp_data) <- paste('met_', rownames(comp_data), sep = '')
comp_data['id'] <- rownames(comp_data)
rownames(met_abun_data) <- paste('met_', rownames(met_abun_data), sep = '')
group_data$max_group <- group_data$`max_group`

# Merge metadata
met_full <- data.frame(id = rownames(met_abun_data), 
                       met_abun_data, 
                       max_group = group_data$`max_group`,
                       row.names = NULL)

# Apply scaling
met_scaled <- t(apply(met_abun_data, 1, rescale01))
tax_scaled <- t(apply(tax_abun_data, 1, rescale01))
colnames(met_scaled) <- colnames(met_abun_data)
colnames(tax_scaled) <- colnames(tax_abun_data)

# ----------------------------------------------------------------------------
# 5. Create Annotation Colors
# ----------------------------------------------------------------------------

cat("Creating color annotations...\n")

annotation_colors <- list()

# Gradient colors for each sample group
for (grp in SAMPLE_GROUPS) {
  grp_clean <- gsub("[^A-Za-z0-9\\.]", ".", grp)
  max_val <- max(met_abun_data[, grp], na.rm = TRUE)
  annotation_colors[[grp]] <- create_color_gradient(CLUSTER_COLORS[grp])(50)[1:as.integer(max_val * 50)]
}

# Degree gradient
annotation_colors$Degree <- create_color_gradient("#000000")(50)

# Discrete colors for max_group
annotation_colors$max_group <- CLUSTER_COLORS

# Match column names in annotation data
sample_cols <- intersect(colnames(met_abun_data), names(annotation_colors))

# --- Section 7: Prepare Replicated Data (Optional) ---
if (N_REPLICATES > 0) {
  cat("Adding replicate columns...\n")
  reps_to_add <- c('met_Halogenated_fatty_acids', 'met_Organoiodides', 
                   'met_Halobenzenes', 'met_Alkyl_iodides')
  
  for (rep_col in reps_to_add) {
    if (rep_col %in% colnames(hm_data)) {
      hm_data <- cbind(hm_data, hm_data[, rep(rep_col, N_REPLICATES)])
    }
  }
}
# ----------------------------------------------------------------------------
# 6. Cluster and Plot Heatmaps
# ----------------------------------------------------------------------------

cat("Computing hierarchical clustering...\n")

# Replace NA values for clustering
mat_for_clust <- hm_data
mat_for_clust[is.na(mat_for_clust)] <- 0

# Row and column clustering
row_dist <- dist(mat_for_clust)
row_clust <- hclust(row_dist, method = "average")

col_dist <- dist(t(mat_for_clust))
col_clust <- hclust(col_dist, method = "average")


cat("Generating ComplexHeatmap...\n")

top_ha <- HeatmapAnnotation(
  df = met_abun_data[, sample_cols, drop = FALSE],
  col = annotation_colors[sample_cols],
  annotation_name_gp = gpar(fontsize = 10)
)

left_ha <- rowAnnotation(
  df = tax_scaled,
  col = annotation_colors[sample_cols],
  annotation_name_gp = gpar(fontsize = 10)
)

ht <- Heatmap(
  hm_data,
  cluster_rows = row_clust,
  cluster_columns = col_clust,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = top_ha,
  left_annotation = left_ha,
  name = "Correlation"
)

# Optional: Add decorative rectangle
decorate_heatmap_body(ht, {
  grid.rect(x = unit((3 + 6) / 2, "native"),
            y = unit(0.5, "npc"),
            width = unit(4, "native"),
            height = unit(1, "npc"),
            gp = gpar(col = "red", lwd = 2, lty = 2))
})

pdf(file = paste0(FIGURE_DIR, "/correlation_heatmap_complex.png"), width = 14, height = 12)
draw(ht)

# ----------------------------------------------------------------------------
# 7. Diagnostic Plot
# ----------------------------------------------------------------------------
diagnostic_plot <- ggplot(data.frame(x = 1, y = 1, fill = "Deep Sinkhole"), 
                          aes(x = x, y = y, fill = fill)) +
  geom_point(size = 5) +
  scale_fill_gradient(low = "#FFFFFF", high = "#F00") +
  theme_void()

ggsave(filename = paste0(OUTPUT_DIR, "/legend.svg"), 
       plot = diagnostic_plot, width = 2, height = 2)