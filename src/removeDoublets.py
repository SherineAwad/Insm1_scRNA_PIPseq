import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import os
import re

sys.modules['importlib.metadata'] = importlib_metadata

# ---------------- Argument Parsing ----------------
parser = argparse.ArgumentParser(description="Remove doublets from AnnData object")
parser.add_argument("--input", required=True, help="Input .h5ad file")
parser.add_argument("--output", required=True, help="Output .h5ad file")
parser.add_argument("--markers", required=True,
                    help="Text file with one marker gene per line (symbols)")
parser.add_argument('--threshold', type=float, default=0.5, help='Doublet score threshold for filtering')
args = parser.parse_args()

markers_file = args.markers
threshold = args.threshold
input_file = args.input

print(f"Markers file: {markers_file}")
print(f"Doublet score threshold: {threshold}")

# ---------------- Extract original doublet threshold from filename ----------------
match = re.search(r'doubletDetected([0-9.]+)', os.path.basename(input_file))
input_doublet_score = match.group(1) if match else "NA"

# ---------------- Read in AnnData object ----------------
adata = sc.read(input_file)

# Filter cells by predicted_doublet score
adata = adata[adata.obs['predicted_doublet'] <= threshold].copy()
print(f"✅ Kept cells with predicted_doublet <= {threshold}")

# ---------------- Save UMAP plot of clustering ----------------
sc.pl.umap(
    adata,
    color=["leiden"],
    save=f"doubletsRemoved_{input_doublet_score}_{threshold}_clusters.png",
    legend_loc="on data",
    show=False
)

# ---------------- Prepare for plotting marker genes ----------------
# Use gene_symbols from raw if available
if adata.raw is not None and 'gene_symbols' in adata.raw.var.columns:
    adata_plot = adata.raw.to_adata()
    adata_plot.var_names = adata_plot.var['gene_symbols'].astype(str)
    print("Using gene_symbols from raw for plotting markers")
else:
    adata_plot = adata.copy()
    print("Warning: 'gene_symbols' not found in raw.var. Markers must match var_names exactly.")

# ---------------- Read marker genes ----------------
marker_genes = [line.strip() for line in open(markers_file) if line.strip()]

# Plot each marker gene separately
for gene in marker_genes:
    if gene in adata_plot.var_names:
        sc.pl.scatter(
            adata_plot,
            color=gene,
            title=gene,
            basis='umap',
            save=f"doubletRemoved_{input_doublet_score}_{threshold}_{gene}.png",
            show=False
        )
    else:
        print(f"⚠️ Gene not found in data: {gene}")

# ---------------- Finalize ----------------
adata.obs_names_make_unique()
adata.write(args.output)
print(f"✅ Saved filtered AnnData to {args.output}")
print(f"✅ All marker plots saved in 'figures/' (Scanpy default)")

