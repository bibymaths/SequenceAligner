#!/usr/bin/env python3
"""
analysis_pipeline.py — Pipeline for large DP matrices written with int32 header:
- Parses header of two int32s for rows & cols
- Memory-maps the DP matrix (int32 → float64)
- Downsamples heatmap
- Incremental PCA with auto component selection
- 2D projection via UMAP or PCA
- Clustering with MiniBatchKMeans
- Configurable via CLI, logs progress
"""
import os
import argparse
import logging
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
from sklearn.decomposition import IncrementalPCA
from sklearn.cluster import MiniBatchKMeans
matplotlib.use('Agg')

# Optional UMAP
try:
    import umap
    HAS_UMAP = True
except ImportError:
    from sklearn.manifold import TSNE
    HAS_UMAP = False


def parse_args():
    p = argparse.ArgumentParser(description="Analysis pipeline for DP matrix with int32 header")
    p.add_argument("--bin", required=True,
                   help="Path to binary DP matrix file with int32 header")
    p.add_argument("--down_rows", type=int, default=10000,
                   help="Heatmap downsample rows")
    p.add_argument("--down_cols", type=int, default=10000,
                   help="Heatmap downsample cols")
    p.add_argument("--batch", type=int, default=1000,
                   help="Batch size for IncrementalPCA/clustering")
    p.add_argument("--max_components", type=int, default=100,
                   help="Max PCA components to fit incrementally")
    p.add_argument("--var_threshold", type=float, default=0.90,
                   help="Explained variance threshold for component selection")
    p.add_argument("--clusters", type=int, default=10,
                   help="Number of clusters for KMeans")
    return p.parse_args()


def load_matrix(bin_path):
    # Read header: first 8 bytes = two little-endian int32s
    with open(bin_path, 'rb') as f:
        header = f.read(8)
    rows, cols = np.frombuffer(header, dtype='<i4')
    total_bytes = os.path.getsize(bin_path)
    expected = 8 + rows * cols * 4
    if total_bytes != expected:
        logging.warning(
            f"File size {total_bytes} != header size {expected}; extra padding ignored"
        )
    logging.info(f"Loaded header dims: {rows}×{cols}")
    # Memory-map int32 body starting at offset 8
    mat_int = np.memmap(bin_path, dtype='<i4', mode='r', offset=8, shape=(rows, cols))
    return mat_int, rows, cols


def downsample(mat, target_rows, target_cols):
    rows, cols = mat.shape
    rstep = max(1, rows // target_rows)
    cstep = max(1, cols // target_cols)
    trimmed = mat[:(rows // rstep) * rstep, :(cols // cstep) * cstep]
    ds = trimmed.reshape(-1, rstep, trimmed.shape[1] // cstep, cstep).mean(axis=(1,3))
    logging.info(f"Downsampled to: {ds.shape}")
    return ds


def plot_heatmap(mat, path, cmap='viridis'):
    plt.figure(figsize=(8,8))
    # plt.imshow(mat,
    #            aspect='auto',
    #            cmap='viridis',
    #            vmin=mat.min(),  # optional: fixes your color scale
    #            vmax=mat.max(),  # optional
    #            interpolation='nearest'  # no smoothing & minimal overhead
    #            )
    plt.matshow(mat)
    plt.axis('off')
    plt.savefig(path, dpi=150, bbox_inches='tight', pad_inches=0)
    plt.close()
    logging.info(f"Heatmap saved to {path}")


def incremental_pca(mat, batch, max_comp, var_thresh):
    rows = mat.shape[0]
    ipca = IncrementalPCA(n_components=max_comp)
    for i in tqdm(range(0, rows, batch), desc="PCA fitting"):
        ipca.partial_fit(mat[i:i+batch])
    cumvar = np.cumsum(ipca.explained_variance_ratio_)
    n_sel = np.searchsorted(cumvar, var_thresh) + 1
    logging.info(f"Selected {n_sel}/{max_comp} PCs for {var_thresh*100:.1f}% variance")
    return ipca, n_sel, cumvar


def project_2d(ipca, mat, batch, n_sel):
    rows = mat.shape[0]
    if HAS_UMAP:
        reducer = umap.UMAP(n_components=2, random_state=42)
        idx = np.random.choice(rows, min(5000, rows), replace=False)
        data = ipca.transform(mat[idx])[:,:n_sel]
        emb = reducer.fit_transform(data)
    else:
        ipca2 = IncrementalPCA(n_components=2)
        for i in tqdm(range(0, rows, batch), desc="PCA2 fitting"):
            block = ipca.transform(mat[i:i+batch])[:,:n_sel]
            ipca2.partial_fit(block)
        emb_parts = []
        for i in range(0, rows, batch):
            emb_parts.append(ipca2.transform(ipca.transform(mat[i:i+batch])[:,:n_sel]))
        emb = np.vstack(emb_parts)
    return emb


def cluster_data(reduced, batch, n_clusters):
    rows = reduced.shape[0]
    km = MiniBatchKMeans(n_clusters=n_clusters, batch_size=batch, random_state=42)
    for i in tqdm(range(0, rows, batch), desc="Clustering"):
        km.partial_fit(reduced[i:i+batch])
    labels = km.predict(reduced)
    return labels


def scatter_plot(points, labels, path, title, cmap='tab10'):
    plt.figure(figsize=(6,6))
    plt.scatter(points[:,0], points[:,1], c=labels, s=2, cmap=cmap)
    plt.title(title)
    plt.savefig(path, dpi=150)
    plt.close()
    logging.info(f"Scatter saved to {path}")


def scree_plot(cumvar, path):
    plt.figure(figsize=(6,4))
    plt.plot(cumvar, marker='o')
    plt.xlabel('PC #'); plt.ylabel('Cum. explained var')
    plt.savefig(path, dpi=150)
    plt.close()
    logging.info(f"Scree plot saved to {path}")


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
    args = parse_args()
    mat, rows, cols = load_matrix(args.bin)

    ds = downsample(mat, args.down_rows, args.down_cols)
    plot_heatmap(ds, 'downsampled_heatmap_1.png')

    # ipca, n_sel, cumvar = incremental_pca(mat, args.batch, args.max_components, args.var_threshold)
    # scree_plot(cumvar, 'scree_plot.png')
    #
    # emb = project_2d(ipca, mat, args.batch, n_sel)
    # scatter_plot(emb, np.zeros(len(emb)), 'projection2d.png', '2D Projection')
    #
    # reduced_full = ipca.transform(mat)[:,:n_sel]
    # labels = cluster_data(reduced_full, args.batch, args.clusters)
    # scatter_plot(emb, labels, 'cluster2d.png', 'Clustering on 2D')

if __name__ == '__main__':
    main()
