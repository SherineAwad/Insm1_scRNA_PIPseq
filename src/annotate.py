#!/usr/bin/env python3
import scanpy as sc
import argparse
import os
import pandas as pd

# ------------------------
# Arguments
# ------------------------
parser = argparse.ArgumentParser(description="Annotate clusters and plot UMAP")
parser.add_argument("--input", required=True, help="Input .h5ad file")
parser.add_argument("--output", required=True, help="Output .h5ad file")
parser.add_argument("--annotations", required=True, help="annotations.txt (cluster<Comma>label)")
parser.add_argument("--cluster_column", default="leiden", help="Column in adata.obs to use as cluster labels (default: leiden)")
args = parser.parse_args()

# ------------------------
# Load AnnData
# ------------------------
adata = sc.read(args.input)
print(f"Loaded AnnData: {adata.n_obs} cells x {adata.n_vars} genes")

# ------------------------
# Load annotations
# ------------------------
ann = pd.read_csv(
    args.annotations,
    sep=",",
    header=None,
    names=["cluster", "label"],
    dtype=str
)

# Map cluster -> label
cluster_to_label = dict(zip(ann["cluster"], ann["label"]))
adata.obs["annotation"] = adata.obs[args.cluster_column].map(cluster_to_label)

print("Annotations applied")

# ------------------------
# UMAP with legend
# ------------------------
sc.pl.umap(
    adata,
    color="annotation",
    size=2,
    show=False,
    save=f"_{os.path.splitext(os.path.basename(args.output))[0]}_annotations.png"
)

# ------------------------
# UMAP with labels on clusters
# ------------------------
sc.pl.umap(
    adata,
    color="annotation",
    size=2,
    legend_loc="on data",
    show=False,
    save=f"_{os.path.splitext(os.path.basename(args.output))[0]}_annotationsON.png"
)

# ------------------------
# Save
# ------------------------
adata.write(args.output)
print(f"Saved annotated AnnData to {args.output}")

