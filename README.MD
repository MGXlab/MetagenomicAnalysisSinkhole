# Code and Data for Sinkhole Paper

This repository contains code and data related to the paper "Specialized microbiome enriched in nitrification and halogen cycling potential resolved from multi-omics investigation of a Caribbean suboxic and acidic marine sinkhole" by Meine D. Boer *et al*.   

Scripts are located in ```scripts/```.   
- ```abiotic.R``` Visualization and analysis of sinkhole environmental data
- ```clust_data.py``` Script to cluster environmental data
- ```FA_general.ipynb``` Notebook for analysis and plotting of KEGG modules
- ```FA_plotting.ipynb``` Notebook for analysis and plotting of KEGG KOs
- ```Heatmap.R``` Script to calculate and plot microbe-metabolite heatmaps
- ```RDA.R``` Script for taxonomic analysis
- ```tree.R``` Script for generating trees with heatmaps

Files  are located in ```data/```:
- ```abiotic.csv``` Environmental data on samples collected
- ```bin_info.csv``` Information on dereplicated MAGs
- ```contig_abundance_example.tsv``` Example data for normalized contig abundance
- ```contig_annotations_example.tsv``` Example data for contig annotation by emapper
- ```contig_annotations_example.tsv``` Example data for contig annotation by emapper
- ```genus_frac.csv``` Relative abundances of taxa
- ```genus_nreads.csv``` Read counts of taxa
- ```KEGG_module_database.csv``` Translation table of KOs to modules
- Sub-folder ```drep_MAGs``` FASTA files of dereplicated MAGs
