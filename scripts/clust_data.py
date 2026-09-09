import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage

clust_data = pd.read_csv('../data/abiotic.csv', index_col=0, delimiter=',')

clust_data.replace(0, np.nan, inplace=True)
for i in clust_data.columns:
    clust_data[i].fillna(clust_data[i].min(), inplace=True)
for i in clust_data.columns[-7:]:
    clust_data[i] = (clust_data[i]).map(np.log10)

for i in clust_data.columns:
    clust_data[i].fillna(clust_data[i].min(), inplace=True)

clust_data2 = clust_data.apply(zscore, axis=0, nan_policy='omit')

link_rows = linkage(clust_data2, method='ward', metric='euclidean')
link_cols = linkage(clust_data2.T, method='ward', metric='euclidean')

clust_data2.to_csv('../data/abiotic_clust_data.csv')