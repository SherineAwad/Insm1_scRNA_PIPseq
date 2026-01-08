#!/usr/bin/env python3
import scanpy as sc
import argparse
import os
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Filter, recluster, QC")
    parser.add_argument("--input", required=True, help="Input h5ad")
    parser.add_argument("--output", required=True, help="Output h5ad")
    args = parser.parse_args()

    out_base = os.path.splitext(os.path.basename(args.output))[0]
    input_base = os.path.splitext(os.path.basename(args.input))[0]

    # ----------------------------
    # Load
    # ----------------------------
    print("Loading h5ad...")
    adata = sc.read_h5ad(args.input)

    # ----------------------------
    # Filter cell types
    # ----------------------------
    remove_celltypes = {"Microglia", "Endothelial", "Astrocytes", "RPE"}

    print("Before filtering:")
    print(adata.obs["annotation"].value_counts())

    adata = adata[~adata.obs["annotation"].isin(remove_celltypes)].copy()

    print("After filtering:")
    print(adata.obs["annotation"].value_counts())

    # ----------------------------
    # UMAP after removing clusters (before reclustering)
    # Use the original UMAP coordinates
    # ----------------------------
    removed_umap_file = f"_{input_base}_removedClusters_UMAP.png"
    sc.pl.umap(
        adata,
        color="annotation",
        legend_loc="on data",
        title="After Removing Selected Cell Types",
        show=False,
        save=removed_umap_file
    )
    # shrink cluster numbers (on data)
    for txt in plt.gca().texts:
        txt.set_fontsize(6)
    plt.draw()
    plt.close()

    # ----------------------------
    # Reclustering (use existing log counts)
    # ----------------------------
    print("Reclustering...")
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(
        adata,
        resolution=1.0,
        key_added="leiden_recluster",
        flavor="igraph",
        n_iterations=2,
        directed=False
    )

    # ----------------------------
    # Reclustered UMAP
    # ----------------------------
    sc.pl.umap(
        adata,
        color="leiden_recluster",
        legend_loc="on data",
        title="Reclustered Leiden",
        show=False,
        save=f"_{out_base}_UMAP_recluster.png"
    )
    for txt in plt.gca().texts:
        txt.set_fontsize(6)
    plt.draw()
    plt.close()

    # ----------------------------
    # QC violin
    # ----------------------------
    qc_vars = ["n_genes", "n_counts", "percent_mito"]

    sc.pl.violin(
        adata,
        keys=qc_vars,
        groupby="leiden_recluster",
        rotation=90,
        stripplot=False,
        multi_panel=True,
        show=False,
        save=f"_{out_base}_QC_violin_by_cluster.png"
    )
    for ax in plt.gcf().axes:
        ax.tick_params(axis='x', labelsize=8)
    plt.draw()
    plt.close()

    # ----------------------------
    # Save h5ad
    # ----------------------------
    print(f"Saving reclustered h5ad: {args.output}")
    adata.write_h5ad(args.output)

    print("DONE.")
    print(f"Figures:")
    print(f" {removed_umap_file}")
    print(f" _{out_base}_UMAP_recluster.png")
    print(f" _{out_base}_QC_violin_by_cluster.png")

if __name__ == "__main__":
    main()

