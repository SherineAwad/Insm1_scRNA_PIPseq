import scanpy as sc
import argparse
import matplotlib.pyplot as plt
import os

# ------------------------
# Parse command-line arguments
# ------------------------
parser = argparse.ArgumentParser(description='Preprocess scRNA-seq AnnData object and plot UMAP by sample.')
parser.add_argument('--input', required=True, help='Input .h5ad file')
parser.add_argument('--output', required=True, help='Output .h5ad file')
args = parser.parse_args()

# ------------------------
# Base name for saving plots
# ------------------------
base_name = os.path.splitext(os.path.basename(args.output))[0]

# ------------------------
# Read AnnData
# ------------------------
adata = sc.read(args.input)
print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ------------------------
# QC metrics: mitochondrial genes
# ------------------------
if 'gene_symbols' in adata.var.columns:
    gene_symbols_str = adata.var['gene_symbols'].astype(str)
    mt_genes = gene_symbols_str.str.upper().str.startswith('MT-')
    n_mt_genes = mt_genes.sum()
    print(f"Detected {n_mt_genes} mitochondrial genes")
    if n_mt_genes > 0:
        print("MT genes examples:", gene_symbols_str.loc[mt_genes].tolist()[:10])
else:
    mt_genes = None
    print("Warning: 'gene_symbols' not found in adata.var. Percent mito will not be calculated.")

# Total counts and gene counts per cell
if hasattr(adata.X, "A"):
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
else:
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)

# Percent mitochondrial counts
if mt_genes is not None and mt_genes.sum() > 0:
    if hasattr(adata.X, "A"):
        adata.obs['percent_mito'] = adata[:, mt_genes].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1 * 100
    else:
        adata.obs['percent_mito'] = adata[:, mt_genes].X.sum(axis=1) / adata.X.sum(axis=1) * 100
else:
    adata.obs['percent_mito'] = 0.0

# ------------------------
# Normalize and log-transform
# ------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # store raw counts

# ------------------------
# Highly variable genes (HVG)
# ------------------------
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000)
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
print(f"Preprocessed AnnData saved to {args.output}")

