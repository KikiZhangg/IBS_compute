#!/usr/bin/env python3

import os
import glob
import itertools
import pandas as pd
import numpy as np
from ibs_pairwise import compute_ibs_score

# For clustering and plotting
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

def get_short_name(file_path):
    base = os.path.basename(file_path)
    filename_no_ext, _ = os.path.splitext(base)
    parts = filename_no_ext.split("DNA23andMe", 1)
    short_name = parts[0] if len(parts) > 1 else filename_no_ext
    return short_name

def compute_all_pairwise_ibses(folder_path):
    genotype_files = sorted(glob.glob(os.path.join(folder_path, "*.txt")))
    results = []
    

    for file1, file2 in itertools.combinations(genotype_files, 2):
        ibs_val = compute_ibs_score(file1, file2)
        file1_short = get_short_name(file1)
        file2_short = get_short_name(file2)
        
        results.append({
            'file1': file1_short,
            'file2': file2_short,
            'ibs': ibs_val
        })
    
    df = pd.DataFrame(results)
    return df

if __name__ == "__main__":
    folder = '/Users/kikizhang/Downloads/genotype_files/'
    df_ibses = compute_all_pairwise_ibses(folder)

    print("\nPairwise IBS:")
    print(df_ibses)


    # 1) Create a square matrix of IBS scores from the pairwise data
    individuals = sorted(
        set(df_ibses['file1']).union(set(df_ibses['file2']))
    )

    # Initialize a DataFrame (N x N) with NaN
    ibs_matrix = pd.DataFrame(
        data=np.nan,
        index=individuals,
        columns=individuals
    )

    for ind in individuals:
        ibs_matrix.loc[ind, ind] = 1.0

    for _, row in df_ibses.iterrows():
        f1, f2, score = row['file1'], row['file2'], row['ibs']
        ibs_matrix.loc[f1, f2] = score
        ibs_matrix.loc[f2, f1] = score


    # 2) Convert IBS matrix to distance matrix = 1 - IBS
    distance_matrix = 1.0 - ibs_matrix
    # 3) Use scipy to perform hierarchical clustering
    condensed_dist = squareform(distance_matrix.values, checks=False)
    # Perform hierarchical clustering with Ward's method
    Z = linkage(condensed_dist, method='ward')

    plt.figure(figsize=(8, 6))
    dendrogram(Z, labels=individuals, leaf_rotation=90,leaf_font_size=10)
    plt.title("Hierarchical Clustering Dendrogram (Ward linkage)")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.show()
