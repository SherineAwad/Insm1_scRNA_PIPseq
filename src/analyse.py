import scanpy as sc
import argparse
import matplotlib.pyplot as plt
import os

# ------------------------
# Parse command-line arguments
# ------------------------
parser = argparse.ArgumentParser(description='Preprocess filtered scRNA-seq AnnData object and plot UMAP.')
parser.add_argument('--input', required=True, help='Input .h5ad file (already filtered)')
parser.add_argument('--output', required=True, help='Output .h5ad file')
args = parser.parse_args()

# ------------------------
# Base name for saving plots
# ------------------------
base_name = os.path.splitext(os.path.basename(args.output))[0]

# ------------------------
# Read filtered AnnData
# ------------------------
adata = sc.read(args.input)
print(f"Loaded filtered AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ------------------------
# Normalize and log-transform
# ------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata 

# ------------------------
# Highly variable genes (HVG)
# ------------------------
# Use seurat_v3 if available, otherwise fall back to seurat
try:
    import skmisc
    flavor_used = 'seurat_v3'
except ImportError:
    flavor_used = 'seurat'
    print("Warning: skmisc not installed, falling back to 'seurat' flavor for HVG selection")

sc.pp.highly_variable_genes(adata, flavor=flavor_used, n_top_genes=2000)
print(f"Using HVG flavor: {flavor_used}")
print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")

# ------------------------
# Scale the data
# ------------------------
sc.pp.scale(adata, max_value=10)

# ------------------------
# PCA
# ------------------------
sc.tl.pca(adata, svd_solver='arpack')

# ------------------------
# Compute neighborhood graph
# ------------------------
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
print("Neighborhood graph computed.")

# ------------------------
# UMAP computation
# ------------------------
sc.tl.umap(adata)

# Plot UMAP colored by sample
if 'sample' in adata.obs.columns:
    sc.pl.umap(adata, color='sample', size=2, save=f"_{base_name}.png")
    print(f"UMAP plot colored by sample saved as {base_name}_UMAP.png")
else:
    print("Warning: 'sample' column not found in adata.obs. UMAP coloring skipped.")

# ------------------------
# Save processed AnnData
# ------------------------
adata.write(args.output)
print(f"Processed AnnData saved to {args.output}")

