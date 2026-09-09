# ============================================================================
# Sinkhole Saba Metagenomics Analysis Script
# Purpose: Plotting trees with heatmap data
# Author: MDBoer
# ============================================================================

library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(pheatmap)

# ----------------------------------------------------------------------------
# 1. PATH CONFIGURATION
# ----------------------------------------------------------------------------
ptt <- "../data/labeled_tree.tre"
hm_data_path <- '../data/hm_data.csv'
bin_info_path <- '../data/bin_info.csv'

# ----------------------------------------------------------------------------
# 2. DATA LOADING
# ----------------------------------------------------------------------------
hm_data <- read.csv(hm_data_path, row.names=1,stringsAsFactor=F)
bin_info <- read.csv(bin_info_path, row.names=1)
t <- read.tree(ptt)
p <- ggtree(rt) + geom_tiplab()

# ----------------------------------------------------------------------------
# 2. TREE PLOTTING
# ----------------------------------------------------------------------------

node_cyano <- 152
node_chloro <- 154
node_plancto <- 125
node_proteo <- 175
node_thems <- 145
node_ver <- 113
node_actino <- 120
node_marini <- 104
node_nip <- 141

p2 <- p + geom_tiplab(data=tips_to_label, mapping=aes(label=label))

comb <- p2 +
  geom_cladelabel(node=node_cyano, label="Cyanobacteria", color="black")+
  geom_cladelabel(node=node_chloro, label="Chloroflexota", color="black") +
  geom_cladelabel(node=node_plancto, label="Planctomycetota", color="black") +
  geom_cladelabel(node=node_proteo, label="Proteobacteria", color="black") +
  geom_cladelabel(node=164, label="Proteobacteria", color="black") +
  geom_cladelabel(node=node_thems, label="Thermoplasmatota", color="black")+
  geom_cladelabel(node=node_ver, label="Veruccomicrobiota", color="black")+
  geom_cladelabel(node=node_actino, label="Actinobacteriota", color="black")+
  geom_cladelabel(node=node_marini, label="Marinisomata", color="black")+
  scale_fill_gradientn(colors=c("lightgray", "blue", "darkblue"))

ggsave('../figures/Comb_plot.svg', comb, scale=1.8)

c_Surface = colorRampPalette(c("#FFFFFF", "#009E73"))(50)
c_Platform = colorRampPalette(c("#FFFFFF", "#ddcc77"))(50)
c_Shallow_Sinkhole = colorRampPalette(c("#FFFFFF", "#a2e7f2"))(50)
c_Deep_Sinkhole = colorRampPalette(c("#FFFFFF", "#0072B2"))(50)
c_Acid_Lake = colorRampPalette(c("#FFFFFF", "#E69F00"))(50)
c_Mesopelagic = colorRampPalette(c("#FFFFFF", "#CC79A7"))(50)
c_Black = colorRampPalette(c("#FFFFFF", "#000"))(50)

row_colors = list(Surface=c_Surface, Shallow.Sinkhole=c_Shallow_Sinkhole, Mesopelagic = c_Mesopelagic, Deep.Sinkhole = c_Deep_Sinkhole, Acid.Lake=c_Acid_Lake, Platform=c_Platform, contamination = c_Black, completeness=c_Black)

# ----------------------------------------------------------------------------
# 2. HEATMAP PLOTTING
# ----------------------------------------------------------------------------

pheatmap(select(data, c('completeness','contamination')), annotation_colors = row_colors, show_rownames = FALSE, annotation_row=select(data, c('Shallow.Sinkhole')))
data <- bin_info

t$tip.label
intersect(rownames(bin_info),t$tip.label)
data_colored <- data %>%
  mutate(
    A_col = scales::col_numeric(c("white", "#009E73"), domain = NULL)(rescale(Surface)),
    B_col = scales::col_numeric(c("white", "#ddcc77"), domain = NULL)(rescale(Platform)),
    C_col = scales::col_numeric(c("white", "#a2e7f2"), domain = NULL)(rescale(Shallow.Sinkhole)),
    D_col = scales::col_numeric(c("white", "#0072B2"), domain = NULL)(rescale(Deep.Sinkhole)),
    E_col = scales::col_numeric(c("white", "#E69F00"), domain = NULL)(rescale(Acid.Lake)),
    F_col = scales::col_numeric(c("white", "#CC79A7"), domain = NULL)(rescale(Mesopelagic)),
    G_col = scales::col_numeric(c("white", "#000000"), domain = NULL)(rescale(completeness)),
    H_col = scales::col_numeric(c("white", "#F00"), domain = NULL)(rescale(contamination)),
  )


data_for_gheatmap <- data_colored %>%
  select(A_col, B_col, C_col, D_col,E_col, F_col, G_col, H_col)
dfg <- data_for_gheatmap[t$tip.label, , drop = FALSE]

data_for_gheatmap <- data %>%
  select(Surface, Platform, Shallow.Sinkhole, Deep.Sinkhole, Acid.Lake, Mesopelagic, completeness, contamination)

plot <- gheatmap(p = p, data = data_for_gheatmap, offset = 0.5, width = 0.3)
hm_data <- read.csv(hm_data_path, row.names=1,stringsAsFactor=F)
my_colors <- brewer.pal(10, "Blues")
plot <- gheatmap(p = p, data = dfg, offset = 0.5, width = 0.3) + scale_fill_identity()
plot <- gheatmap(p = p, data = hm_data, offset = 0.5, width = 0.3) + scale_fill_gradientn(colors = my_colors, breaks = c(0,1,2,3,4,5,6,7,8,9,10))


