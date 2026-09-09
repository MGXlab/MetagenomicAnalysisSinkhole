# ============================================================================
# Sinkhole Saba Metagenomics Analysis Script
# Purpose: Ordination analysis (PCA, RDA, CCA) of microbial community data
# Author: MDBoer
# ============================================================================

library(dplyr)
library(ggplot2)
library(vegan)
library(compositions)

# ----------------------------------------------------------------------------
# 1. PATH CONFIGURATION + GLOBAL VARIABLES
# ----------------------------------------------------------------------------
# Directories
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

# Color palette (consistent with previous analysis)
CLUSTER_COLORS <- c(
  'Surface'          = "#009E73",
  'Platform'         = "#ddcc77",
  'Shallow Sinkhole' = "#a2e7f2",
  'Deep Sinkhole'    = "#0072B2",
  'Acid Lake'        = "#E69F00",
  'Mesopelagic'      = "#CC79A7"
)

ENV_VARIABLES <- c("pH (na)","Oxygen (μmol/kg)","Phosphate (μmol/kg)",
                   "Nitrate (μmol/kg)","Nitrite (μmol/kg)","Ammonium (μmol/kg)",
                   "Chlorophyll (μg/kg)","Silicate (μmol/kg)","DIC (μmol/kg)",
                   "Salinity (ppt)","Transmissivity (%)","Temperature (°C)")


# ----------------------------------------------------------------------------
# 2. DATA LOADING
# ----------------------------------------------------------------------------
cat("Loading taxonomic abundance data...\n")
tad_f <- read.csv(tax_loc, row.names = 1, check.names = FALSE, sep = ",")

cat("Loading metadata...\n")
met_f <- read.csv(meta_loc, row.names = 1, sep = ",", check.names = FALSE)

# ----------------------------------------------------------------------------
# 3. TAXONOMIC DATA TRANSFORMATION
# ----------------------------------------------------------------------------
cat("Transforming taxonomic data...\n")

df_tad <- as.data.frame(tad_f, check.names = FALSE)

# Filter low-abundance samples (sum >= 10)
df_tad_filt <- t(df_tad[rowSums(df_tad) >= 10, ])

# Centered log-ratio transformation (compositional data)
df_tad_scaled <- clr(df_tad_filt)
# Scale metadata variables
df_met <- as.data.frame(scale(met_f[,ENV_VARIABLES]))

# ----------------------------------------------------------------------------
# 4. PRINCIPAL COMPONENT ANALYSIS (PCA)
# ----------------------------------------------------------------------------
cat("Running PCA...\n")

pca_result <- prcomp(df_tad_scaled, center = TRUE, scale. = FALSE)
pca_data <- as.data.frame(pca_result$x)
pca_data$SampleID <- rownames(pca_data)

# Calculate variance explained by each PC
variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
x_label_pca <- paste0("PC1 (", round(variance_explained[1], 1), "%)")
y_label_pca <- paste0("PC2 (", round(variance_explained[2], 1), "%)")

# Merge with metadata
metadata_with_id <- met_f
metadata_with_id$SampleID <- rownames(metadata_with_id)
pca_plot_data <- merge(pca_data, metadata_with_id, by = "SampleID")
# PCA Plot
pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group, label = SampleID)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = sample_colors) +
  labs(
    title = "PCA of Microbial Communities",
    subtitle = paste0("Explains ", round(sum(variance_explained[1:2]), 1), "% variance (PC1+PC2)"),
    x = x_label_pca,
    y = y_label_pca
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

pca_plot
# ----------------------------------------------------------------------------
# 5. REDUNDANCY ANALYSIS (RDA)
# ----------------------------------------------------------------------------
cat("Running RDA...\n")

# Formula: all variables in df_met except 'group' as predictors
response_vars <- setdiff(names(df_met), "group")
if (length(response_vars) > 0) {
  rda_formula <- as.formula(paste("df_tad_scaled ~", paste(response_vars, collapse = " + ")))
  
  rda_result <- rda(rda_formula, data = df_met, na.action = na.omit)
  
  # Extract scores for plotting
  site_scores <- scores(rda_result, display = "sites", scaling = 1)
  biplot_scores <- scores(rda_result, display = "bp", scaling = 1)
  
  site_scores_df <- as.data.frame(site_scores)
  biplot_scores_df <- as.data.frame(biplot_scores)
  
  # Add group information to site scores
  site_scores_df$group <- df_met[rownames(site_scores_df), "group"]
  
  # Calculate variance explained by RDA axes
  if (!is.null(rda_result$CCA)) {
    eigenvalues <- rda_result$CCA$eig
    total_variance <- sum(rda_result$CCA$eig, rda_result$CA$eig)
    rda_var_explained <- (eigenvalues / total_variance) * 100
    x_label_rda <- paste0("RDA1 (", round(rda_var_explained[1], 2), "%)")
    y_label_rda <- paste0("RDA2 (", round(rda_var_explained[2], 2), "%)")
    
    # RDA Biplot
    rda_plot <- ggplot() +
      geom_segment(
        data = biplot_scores_df,
        aes(x = 0, y = 0, xend = RDA1 * 4, yend = RDA2 * 4),
        arrow = arrow(length = unit(0.15, "cm")),
        color = "red", size = 1
      ) +
      geom_text(
        data = biplot_scores_df,
        aes(x = RDA1 * 4, y = RDA2 * 4, label = rownames(biplot_scores_df)),
        color = "black", vjust = -0.5, hjust = 0.5
      ) +
      geom_point(
        data = site_scores_df,
        aes(x = RDA1, y = RDA2, color = group),
        size = 3
      ) +
      scale_color_manual(values = sample_colors) +
      labs(
        title = "RDA Biplot: Environmental Drivers of Community Structure",
        subtitle = paste0("Variance explained: RDA1=", round(rda_var_explained[1], 1), "%, RDA2=", round(rda_var_explained[2], 1), "%"),
        x = x_label_rda,
        y = y_label_rda
      ) +
      theme_minimal() +
      theme(legend.title = element_blank())
    
    print(rda_plot)
  }
  
# ----------------------------------------------------------------------------
# 6. STATISTICAL SIGNIFICANCE TESTS
# ----------------------------------------------------------------------------
  cat("\n--- RDA Significance Tests ---\n")
  
  cat("\nOverall model test:\n")
  print(anova.cca(rda_result, permutations = 999))
  
  cat("\nTest by term:\n")
  print(anova.cca(rda_result, by = "term", permutations = 999))
  
  cat("\nTest by axis:\n")
  print(anova.cca(rda_result, by = "axis", permutations = 999))
} else {
  warning("No response variables available for RDA.")
}

# ----------------------------------------------------------------------------
# 7. CANONICAL CORRESPONDENCE ANALYSIS (CCA)
# ----------------------------------------------------------------------------
cat("\nRunning CCA...\n")

cca_formula <- as.formula(paste("df_tad_scaled ~", paste(response_vars, collapse = " + ")))
cca_model <- cca(cca_formula, data = df_met, na.action = na.omit)

# Calculate variance explained by CCA axes
if (!is.null(cca_model$CCA)) {
  eig_vals <- cca_model$CCA$eig
  var_explained <- (eig_vals / sum(eig_vals)) * 100
  axis1_pct <- round(var_explained[1], 1)
  axis2_pct <- round(var_explained[2], 1)
  
  # Standard CCA plot
  cat("\nPlotting CCA results...\n")
  plot(cca_model, display = c("sites", "species"), scaling = 2,
       main = "Canonical Correspondence Analysis (CCA)", xlab = "", ylab = "")
  title(xlab = paste0("CCA1 (", axis1_pct, "%)"),
        ylab = paste0("CCA2 (", axis2_pct, "%)"))
  
  # Add environmental vectors
  envfit_result <- envfit(cca_model, df_met[, response_vars, drop = FALSE], permutations = 999)
  plot(envfit_result, p.max = 0.05, col = "blue")
}

# ----------------------------------------------------------------------------
# 8. MULTIVARIATE ANALYSIS OF VARIANCE (MANOVA)
# ----------------------------------------------------------------------------
cat("\nRunning MANOVA...\n")

if (length(response_vars) > 0 && "group" %in% colnames(df_met)) {
  manova_formula <- as.formula(
    paste("cbind(", paste(response_vars, collapse = ","), ") ~ group")
  )
  manova_model <- manova(manova_formula, data = df_met)
  
  cat("\nMANOVA Summary:\n")
  print(summary(manova_model))
  
  cat("\nPer-Response ANOVA:\n")
  print(summary.aov(manova_model))
}


# ----------------------------------------------------------------------------
# END OF SCRIPT
# ----------------------------------------------------------------------------
cat("\nAnalysis complete!\n")

