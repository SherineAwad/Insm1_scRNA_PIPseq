#!/usr/bin/env python3
import scanpy as sc
import argparse
import os

# ------------------------
# Arguments
# ------------------------
parser = argparse.ArgumentParser(
    description="Plot marker expression from clustered AnnData"
)
parser.add_argument("--input", required=True, help="Clustered .h5ad file")
parser.add_argument("--markers", required=True,
                    help="Text file with one marker gene per line (symbols)")
args = parser.parse_args()

# ------------------------
# Load data
# ------------------------
adata = sc.read(args.input)
print(f"Loaded AnnData: {adata.n_obs} cells x {adata.n_vars} genes")

# ------------------------
# Make temporary adata for plotting with gene_symbols as var_names
# ------------------------
adata_plot = adata.raw.to_adata()
if 'gene_symbols' in adata_plot.var.columns:
    adata_plot.var_names = adata_plot.var['gene_symbols'].astype(str)
    print("Using gene_symbols for plotting")
else:
    print("Warning: 'gene_symbols' not found in raw.var. Markers must match var_names.")

# ------------------------
# Read markers
# ------------------------
with open(args.markers) as f:
    markers = [line.strip() for line in f if line.strip()]

# Keep only markers present in adata
markers_exist = [m for m in markers if m in adata_plot.var_names]
missing_markers = [m for m in markers if m not in adata_plot.var_names]

if missing_markers:
    print(f"Warning: {len(missing_markers)} markers not found and will be skipped: {', '.join(missing_markers)}")
if not markers_exist:
    raise ValueError("No markers found in the data. Aborting plotting.")

print(f"Markers to plot ({len(markers_exist)}): {', '.join(markers_exist)}")

# ------------------------
# Create figures directory
# ------------------------
fig_dir = "figures"
os.makedirs(fig_dir, exist_ok=True)

# Base name for plot files
base_name = os.path.splitext(os.path.basename(args.input))[0]

# ------------------------
# Feature plots: one per gene
# ------------------------
for gene in markers_exist:
    print(f"Plotting {gene}...")
    sc.pl.umap(
        adata_plot,
        color=gene,
        size=2,
        show=False,
        save=f"_{base_name}_{gene}.png"
    )

print(f"Feature plots saved to {fig_dir}/ (one per gene).")

