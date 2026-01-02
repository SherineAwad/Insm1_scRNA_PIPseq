import scanpy as sc
import doubletdetection
import numpy as np
from scipy import sparse
import argparse
import os
import matplotlib.pyplot as plt

# Parse input file
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--threshold', required=True)

args = parser.parse_args()

myObject = args.input
threshold = float(args.threshold)
newObject = args.output

base_name = os.path.splitext(os.path.basename(newObject))[0]

# Read AnnData object
adata = sc.read_h5ad(myObject)

# Use raw counts for doubletdetection if available
if adata.raw is not None:
    adata_use = adata.raw.X.copy()
else:
    adata_use = adata.X.copy()

# Clean NaNs in adata_use
if sparse.issparse(adata_use):
    print("Cleaning sparse matrix (adata_use)...")
    adata_use.data[np.isnan(adata_use.data)] = 0
    check_matrix = adata_use.toarray()
else:
    print("Cleaning dense matrix (adata_use)...")
    adata_use = np.nan_to_num(adata_use)
    check_matrix = adata_use

# Sanity check
if np.isnan(check_matrix).any():
    raise ValueError("NaNs still present in adata_use after cleanup.")

# Remove zero-count cells (required for doubletdetection)
if sparse.issparse(adata_use):
    cell_sums = np.array(adata_use.sum(axis=1)).ravel()
else:
    cell_sums = adata_use.sum(axis=1)

adata_use = adata_use[cell_sums > 0]

# Run doublet detection
print("Running doublet detection...")
clf = doubletdetection.BoostClassifier(n_iters=5, standard_scaling=True)
doublets = clf.fit(adata_use).predict(
    p_thresh=1e-16,
    voter_thresh=threshold
)

# Threshold plot
figname1 = "figures/" + args.threshold + "threshold_test.png"
print(figname1)
fig1 = doubletdetection.plot.threshold(clf, show=False, p_step=6)
fig1.savefig(figname1, dpi=300)
plt.close(fig1)

# Convergence plot
figname2 = "figures/" + args.threshold + "conversion_test.png"
print(figname2)
fig2 = doubletdetection.plot.convergence(
    clf,
    show=False,
    p_thresh=1e-16,
    voter_thresh=threshold
)
fig2.savefig(figname2, dpi=300)
plt.close(fig2)

# Save results to original AnnData
adata.obs['doublet_score'] = clf.doublet_score()
adata.obs['predicted_doublet'] = doublets

# Save updated AnnData object
adata.write(newObject, compression="gzip")
print("Doublet detection completed and results saved.")

