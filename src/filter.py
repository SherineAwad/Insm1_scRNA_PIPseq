import scanpy as sc
import argparse
import matplotlib.pyplot as plt

# ------------------------
# Parse command-line arguments
# ------------------------
parser = argparse.ArgumentParser(description='Filter AnnData based on QC metrics and plot metrics.')
parser.add_argument('--input', required=True, help='Input .h5ad file')
parser.add_argument('--output', required=True, help='Output .h5ad file')
args = parser.parse_args()

# ------------------------
# Read AnnData
# ------------------------
adata = sc.read(args.input)
print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ------------------------
# Calculate QC metrics
# ------------------------
# Convert gene_symbols to string to handle Categorical type
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

# Total counts per cell
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
# Plot QC metrics before filtering
# ------------------------
sc.pl.violin(
    adata,
    ['n_genes', 'n_counts', 'percent_mito'],
    jitter=0.4,
    multi_panel=True,
    save='_before_filtering.png'
)

# ------------------------
# Filter cells
# ------------------------
adata_filtered = adata[
    (adata.obs['n_genes'] >= 800) & (adata.obs['n_genes'] <= 8000) &
    (adata.obs['n_counts'] >= 1200) & (adata.obs['n_counts'] <= 30000) &
    (adata.obs['percent_mito'] < 25)
].copy()

print(f"Filtered AnnData: {adata_filtered.shape[0]} cells x {adata_filtered.shape[1]} genes")

# ------------------------
# Plot QC metrics after filtering
# ------------------------
sc.pl.violin(
    adata_filtered,
    ['n_genes', 'n_counts', 'percent_mito'],
    jitter=0.4,
    multi_panel=True,
    save='_after_filtering.png'
)

# ------------------------
# Save filtered AnnData
# ------------------------
adata_filtered.write(args.output)
print(f"Filtered AnnData saved to {args.output}")
print("QC violin plots saved as *_before_filtering.png and *_after_filtering.png")

