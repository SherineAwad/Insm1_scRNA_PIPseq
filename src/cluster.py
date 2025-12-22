#!/usr/bin/env python3
import scanpy as sc
import argparse
import os

# ------------------------
# Arguments
# ------------------------
parser = argparse.ArgumentParser(description="Leiden clustering only")
parser.add_argument("--input", required=True, help="Analysed .h5ad file")
parser.add_argument("--output", required=True, help="Clustered .h5ad file")
parser.add_argument("--resolution", type=float, default=1.0, help="Leiden resolution")
args = parser.parse_args()

# ------------------------
# Load AnnData
# ------------------------
adata = sc.read(args.input)
print(f"Loaded AnnData: {adata.n_obs} cells x {adata.n_vars} genes")

# ------------------------
# Leiden clustering ONLY
# ------------------------
sc.tl.leiden(adata, resolution=args.resolution)
print(f"Clusters found: {adata.obs['leiden'].nunique()}")

# ------------------------
# UMAP (reuse existing embedding)
# ------------------------
sc.pl.umap(
    adata,
    color="leiden",
    size=2,
    show=False,
    save=f"_{os.path.splitext(os.path.basename(args.output))[0]}_clusters.png"
)


sc.pl.umap(
    adata,
    color="leiden",
    size=2,
    legend_loc="on data",
    show=False,
    save=f"_{os.path.splitext(os.path.basename(args.output))[0]}_clustersON.png"
)


# ------------------------
# Save
# ------------------------
adata.write(args.output)
print(f"Saved clustered AnnData to {args.output}")

