import scanpy as sc
import pandas as pd
import argparse

# ------------------------
# Parse command-line args
# ------------------------
parser = argparse.ArgumentParser(description='Read multiple 10X matrices and save combined AnnData.')
parser.add_argument('--samples', required=True, help='Path to samples.tsv')
parser.add_argument('--output', required=True, help='Output filename for combined AnnData (.h5ad)')
args = parser.parse_args()

# ------------------------
# Read sample names
# ------------------------
samples = pd.read_csv(args.samples, header=None)[0].tolist()

adatas = []

for sample in samples:
    matrix_file = f"{sample}_matrix.mtx.gz"
    barcodes_file = f"{sample}_barcodes.tsv.gz"
    features_file = f"{sample}_features.tsv.gz"

    # Read barcodes and features
    barcodes = pd.read_csv(barcodes_file, header=None)[0].tolist()
    features = pd.read_csv(features_file, sep='\t', header=None)
    gene_ids = features[0].tolist()       # ENSMUSG IDs, safe for merging
    gene_symbols = features[1].tolist()   # optional, can store in adata.var later

    # Read matrix and transpose (10X is genes x cells)
    adata = sc.read_mtx(matrix_file).T     # <- important transpose
    adata.var_names = gene_ids
    adata.obs_names = barcodes
    adata.obs['sample'] = sample

    # Optional: store gene symbols
    adata.var['gene_symbols'] = gene_symbols

    adatas.append(adata)

# ------------------------
# Combine all samples
# ------------------------
adata_combined = adatas[0].concatenate(
    adatas[1:],
    batch_key='sample',
    batch_categories=samples
)

# ------------------------
# Save combined AnnData
# ------------------------
adata_combined.write(args.output)
print(f"Combined AnnData saved to {args.output}")

