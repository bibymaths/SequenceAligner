#!/usr/bin/env python3
"""
analysis_pipeline.py — Load a huge binary matrix, perform downsampling, PCA, t-SNE, clustering,
and generate related plots without exhausting memory.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import IncrementalPCA, PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans


def load_matrix(bin_path, stats_path):
    """
    Load binary matrix using numpy.memmap with automatic dimension inference.
    """
    # Read stats from JSON
    stats = json.load(open(stats_path))
    nx_json = stats.get('query_length')
    ny_json = stats.get('target_length')

    # Compute total floats from file size
    fsize = os.path.getsize(bin_path)
    n_floats = fsize // np.dtype('float64').itemsize

    # Infer dimensions
    if nx_json and ny_json and nx_json * ny_json == n_floats:
        n_cols, n_rows = int(nx_json), int(ny_json)
    elif nx_json and n_floats % int(nx_json) == 0:
        n_cols = int(nx_json)
        n_rows = int(n_floats // n_cols)
    elif ny_json and n_floats % int(ny_json) == 0:
        n_rows = int(ny_json)
        n_cols = int(n_floats // n_rows)
    else:
        # Fallback to square
        side = int(np.sqrt(n_floats))
        n_rows, n_cols = side, side

    print(f"Loading matrix with shape: ({n_rows}, {n_cols})")
    mat = np.memmap(bin_path, dtype='float64', mode='r', shape=(n_rows, n_cols))
    return mat, n_rows, n_cols


def downsample(mat, target_rows=1000, target_cols=1000):
    """
    Block-average downsampling to given resolution.
    """
    rows, cols = mat.shape
    row_block = max(1, rows // target_rows)
    col_block = max(1, cols // target_cols)
    new_rows = rows // row_block
    new_cols = cols // col_block

    trimmed = mat[:new_rows*row_block, :new_cols*col_block]
    ds = trimmed.reshape(new_rows, row_block, new_cols, col_block).mean(axis=(1,3))
    print(f"Downsampled to: ({new_rows}, {new_cols})")
    return ds


def plot_heatmap(mat, outpath, cmap='viridis'):
    plt.figure(figsize=(8, 6))
    plt.imshow(mat, aspect='auto', cmap=cmap)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"Saved heatmap: {outpath}")


def main():
    bin_path = 'results/local_dp_matrix.bin'
    stats_path = 'results/local_stats.json'

    # Load and inspect matrix
    mat, n_rows, n_cols = load_matrix(bin_path, stats_path)

    # 1) Heatmap of downsampled matrix
    ds = downsample(mat, target_rows=800, target_cols=800)
    plot_heatmap(ds, 'downsampled_heatmap.png')

    # 2) Incremental PCA (50 components)
    n_comp = 50
    batch_size = 1000
    ipca = IncrementalPCA(n_components=n_comp)
    for start in range(0, n_rows, batch_size):
        ipca.partial_fit(mat[start:start+batch_size])
    # Scree plot
    cumvar = np.cumsum(ipca.explained_variance_ratio_)
    plt.figure(figsize=(6,4))
    plt.plot(cumvar, marker='o')
    plt.xlabel('Number of components')
    plt.ylabel('Cumulative explained variance')
    plt.tight_layout()
    plt.savefig('scree_plot.png', dpi=150)
    plt.close()
    print("Saved scree_plot.png")

    # 3) 2D PCA for visualization
    pca2 = PCA(n_components=2)
    X_pca2 = pca2.fit_transform(ipca.transform(mat))
    plt.figure(figsize=(6,6))
    plt.scatter(X_pca2[:,0], X_pca2[:,1], s=1)
    plt.xlabel('PC1'); plt.ylabel('PC2')
    plt.tight_layout()
    plt.savefig('pca2d_scatter.png', dpi=150)
    plt.close()
    print("Saved pca2d_scatter.png")

    # 4) t-SNE on a random subset
    subset = np.random.choice(n_rows, min(5000, n_rows), replace=False)
    tsne = TSNE(n_components=2, init='random', random_state=42)
    X_tsne = tsne.fit_transform(ipca.transform(mat[subset]))
    plt.figure(figsize=(6,6))
    plt.scatter(X_tsne[:,0], X_tsne[:,1], s=5)
    plt.xlabel('t-SNE1'); plt.ylabel('t-SNE2')
    plt.tight_layout()
    plt.savefig('tsne_scatter.png', dpi=150)
    plt.close()
    print("Saved tsne_scatter.png")

    # 5) KMeans clustering on PCA50 reduced data
    kmeans = KMeans(n_clusters=10, random_state=42)
    clusters = kmeans.fit_predict(ipca.transform(mat))
    plt.figure(figsize=(6,6))
    plt.scatter(X_pca2[:,0], X_pca2[:,1], c=clusters, s=1, cmap='tab10')
    plt.xlabel('PC1'); plt.ylabel('PC2')
    plt.tight_layout()
    plt.savefig('cluster_pca2d.png', dpi=150)
    plt.close()
    print("Saved cluster_pca2d.png")


if __name__ == '__main__':
    main()