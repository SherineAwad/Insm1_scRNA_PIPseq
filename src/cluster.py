#!/usr/bin/env python3
import scanpy as sc
import argparse
import os

# ------------------------
# Argument parser
# ------------------------
parser = argparse.ArgumentParser(description="Clustering for scRNA PIP-seq")
parser.add_argument('--input', required=True, help='Input h5ad file (filtered & normalized)')
parser.add_argument('--output', required=True, help='Output h5ad file with clustering results')
parser.add_argument('--resolution', type=float, default=0.5, help='Leiden clustering resolution')
args = parser.parse_args()

# ------------------------
# Base name for saving plots
# ------------------------
base_name = os.path.splitext(os.path.basename(args.output))[0]

# ------------------------
# Load AnnData
# ------------------------
adata = sc.read(args.input)
print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ------------------------
# Check PCA
# ------------------------
if 'X_pca' not in adata.obsm.keys():
    sc.pp.pca(adata, svd_solver='arpack')
    print("PCA computed.")
else:
    print("PCA already exists, skipping PCA computation.")

# ------------------------
# Check neighbors graph
# ------------------------
if 'neighbors' not in adata.uns.keys():
    sc.pp.neighbors(adata)
    print("Neighborhood graph computed.")
else:
    print("Neighbors graph already exists, skipping computation.")

# ------------------------
# Leiden clustering
# ------------------------
sc.tl.leiden(adata, resolution=args.resolution)
print(f"Leiden clustering done. Number of clusters: {adata.obs['leiden'].nunique()}")

# ------------------------
# UMAP embedding
# ------------------------
if 'X_umap' not in adata.obsm.keys():
    sc.tl.umap(adata)
    print("UMAP embedding computed.")
else:
    print("UMAP already exists, skipping computation.")

# ------------------------
# Plot UMAP colored by clusters
# ------------------------
sc.pl.umap(adata, color='leiden', size=2, save=f"_{base_name}_clusters.png", show=False)
print(f"UMAP plot of clusters saved as {base_name}_clusters.png")

# ------------------------
# Save clustered AnnData
# ------------------------
adata.write(args.output)
print(f"Clustered AnnData saved to {args.output}")

