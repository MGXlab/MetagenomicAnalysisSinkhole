#############################################
# Sinkhole Environmental Data Analysis Script
# Author: MD Boer
# Date: 2026-08-25
# Description: Visualization and analysis of 
#              sinkhole environmental data
#############################################
# ----------------------------------------------------------------------------
# 1. Install and load required libraries
# ----------------------------------------------------------------------------
pkgs <- c("ggplot2","readxl","dplyr","tibble","tidyr","ggbreak","ggbeeswarm","ellipse",
          "RColorBrewer","GGally","phyloseq","fmsb")

# Install missing (phyloseq via Bioconductor)
to_install <- setdiff(pkgs, installed.packages()[,"Package"])
if (length(to_install)) {
  if ("phyloseq" %in% to_install && !require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  lapply(to_install, function(p) if(p=="phyloseq") BiocManager::install(p,ask=FALSE) else install.packages(p))
}
lapply(pkgs, library, character.only = TRUE, quietly = TRUE)

library(phyloseq)      # Microbiome data analysis
library(ggplot2)       # Graphics
library(readxl)        # Excel import
library(dplyr)         # Data manipulation
library(tibble)        # Tibble utilities
library(tidyr)         # Data tidying
library(ggbreak)       # Break axes
library(ggbeeswarm)    # Beeswarm plots
library(ellipse)       # Ellipse functions
library(RColorBrewer)  # Color palettes
library(GGally)        # Pairwise plots

# ----------------------------------------------------------------------------
# 2. Configuration + Global Variables
# ----------------------------------------------------------------------------

# Define color scheme for clusters (reusable constant)
CLUSTER_COLORS <- c(
  'Surface'          = "#009E73",
  'Platform'         = "#F0E442",
  'Shallow Sinkhole' = "#a2e7f2",
  'Deep Sinkhole'    = "#0072B2",
  'Acid Lake'        = "#E69F00",
  'Mesopelagic'      = "#CC79A7"
)

CLUSTER_LEVELS <- c('Surface','Platform','Shallow Sinkhole',
                    'Deep Sinkhole','Acid Lake','Mesopelagic')

# Define environmental variables in display order
ENV_VARIABLES <- c("pH (na)","Oxygen (μmol/kg)","Phosphate (μmol/kg)",
                   "Nitrate (μmol/kg)","Nitrite (μmol/kg)","Ammonium (μmol/kg)",
                   "Chlorophyll (μg/kg)","Silicate (μmol/kg)","DIC (μmol/kg)",
                   "Salinity (ppt)","Transmissivity (%)","Temperature (°C)")


# Theme function for consistent styling
theme_seaborn_style <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "#F7F7F7", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "black"),
      panel.border = element_rect(color = "black")
    )
}

# Rescale function for normalization
rescale01 <- function(x) (x - min(x)) / (max(x) - min(x))

# Radar chart helper function
plot_radar <- function(data, title, colors = CLUSTER_COLORS) {
  require(fmsb)  # Ensure fmsb package is installed
  
  # Add max and min rows for scaling
  radar_data <- rbind(
    max = rep(1, ncol(data)),
    min = rep(0, ncol(data)),
    data
  )
  
  radarchart(
    radar_data,
    axistype = 1,
    pcol = "black", pfcol = rgb(0,0,1,0.3), plwd = 2,
    cglcol = "black", cglty = 1, cglwd = 0.8,
    axislabcol = "black", vlcex = 0.8,
    title = title
  )
}
# ----------------------------------------------------------------------------
# 3. Load and Preprocess Data
# ----------------------------------------------------------------------------
df_path <- "../data/abiotic.csv"  # Replace with actual path
df      <- read.csv(df_path, row.names = 1, check.names = FALSE, sep=',')

# Factor ordering for consistent visualization
df$Area   <- factor(df$Area, levels = c('FS','E','S','W','N','O'))
df$Cluster <- factor(df$Cluster, levels = CLUSTER_LEVELS)


# ----------------------------------------------------------------------------
# 4. Plotting and Analysis
# ----------------------------------------------------------------------------

## 4.1 Depth Distribution by Area (Beeswarm)
p_depth <- ggplot(df, aes(x=Sinkhole, y=`Depth (m)`, color=Cluster)) +
  geom_beeswarm(cex=3, size=2) +
  geom_hline(yintercept=80, linetype="dashed", color="grey40") +
  scale_y_continuous(trans = "reverse", limits = c(700, 0), expand = c(0,0)) +
  scale_y_break(c(10, 50), scales = 2) +
  scale_y_break(c(300, 450), scales = 0.5) +
  scale_color_manual(values = CLUSTER_COLORS) +
  facet_grid(~Area, scales="free_x", space="free_x") +
  labs(y="Depth (m)", color="", caption = "Grey dashed line indicates 80m reference depth") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )
ggsave("../figures/depth_distribution.png", p_depth, width = 12, height = 6, dpi = 300)

## 4.2 Environmental Variables Comparison (Boxplot + Jitter)
df_long <- df %>%
  select(all_of(ENV_VARIABLES), Cluster) %>%
  pivot_longer(cols = all_of(ENV_VARIABLES), 
               names_to = "Variable", values_to = "Value")

p_env_vars <- ggplot(df_long, aes(x=Cluster, y=Value, color=Cluster, fill=Cluster)) +
  geom_boxplot(alpha = 0.2, outlier.alpha = 0) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~Variable, scales="free_y", ncol = 4) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  scale_color_manual(values = CLUSTER_COLORS) +
  scale_y_log10() +
  theme_seaborn_style() +
  labs(title = "Environmental Parameters by Cluster",
       y = "Value (log scale)", x = "")

ggsave("../figures/environmental_variables.png", p_env_vars, width = 14, height = 10, dpi = 300)

## 4.3 Correlation Heatmap (if correlation matrix exists)
if(exists("cord") && !is.null(cord)) {
  my_colors <- colorRampPalette(brewer.pal(5, "Spectral"))(100)
  plotcorr(cord, col=my_colors[abs(cord)*50+50], mar=c(1,1,1,1))
  ggsave("../figures/correlation_heatmap.png", width = 10, height = 10, dpi = 300)
} else {
  message("Correlation matrix 'cord' not found - skipping correlation plot")
}

## 4.4 Pairwise Relationships
if("MC" %in% names(df)) {
  p_pairs <- ggpairs(df, mapping = ggplot2::aes(colour=MC))
  ggsave("../figures/pairwise_relationships.png", p_pairs, width = 12, height = 12, dpi = 300)
}

## 4.5 Group Means (Normalized)
df_grouped <- df %>%
  select(where(is.numeric), Cluster) %>%
  group_by(Cluster) %>%
  summarise(across(ENV_VARIABLES, mean, na.rm = TRUE), .groups = 'drop')

df_scaled <- df_grouped %>%
  select(-Cluster) %>%
  mutate(across(ENV_VARIABLES, rescale01)) %>%
  mutate(group = df_grouped$Cluster)

## 4.6 Radar Charts by Cluster
cairo_pdf("../figures/radar_charts.pdf", width = 8, height = 10)
par(mfrow = c(2, 3))
for (g in CLUSTER_LEVELS) {
  group_data <- df_scaled %>%
    filter(group == g) %>%
    select(-group)
  
  plot_radar(group_data, g)
}
# ----------------------------------------------------------------------------
# 3. Save Output
# ----------------------------------------------------------------------------
write.csv(df_grouped, "../data/luster_summary_means.csv", row.names = FALSE)
write.csv(df_scaled, "../data/scaled_cluster_data.csv", row.names = FALSE)