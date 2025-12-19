# Analysing INSM1 Overexpression: scRNA-seq (PIP-seq) Data


# Benchmark studies comparing Seurat and Scanpy
- Rich, J. M., Moses, L., Einarsson, P. H., Jackson, K., Luebbert, L., Booeshaghi, A. S., et al. (2024).  
  *The impact of package selection and versioning on single-cell RNA-seq analysis.* **bioRxiv**.

- Billato, I., Pages, H., Carey, V., Waldron, L., Sales, G., Romualdi, C., & Risso, D. (2025).  
  *Benchmarking large-scale single-cell RNA-seq analysis.* **bioRxiv**.

- Zappia, L., Richter, S., Ramírez-Suástegui, C., Kfuri-Rubens, R., Vornholz, L., Wang, W., et al. (2025).  
  *Feature selection methods affect the performance of scRNA-seq data integration and querying.* **Nature Methods**, 1–11.

# 💡💡💡
# Thanh chose to proceed with this analysis using **Scanpy**.
# 💡💡💡


## Filtering: Same Neurog2 filtering as below

```python
adata_filtered = adata[
    (adata.obs['n_genes'] >= 800) & (adata.obs['n_genes'] <= 8000) &
    (adata.obs['n_counts'] >= 1200) & (adata.obs['n_counts'] <= 30000) &
    (adata.obs['percent_mito'] < 25)
].copy()

```
Basically: 

- Retain cells with **800–8,000 detected genes** to remove low-quality cells and potential doublets.
- Keep cells with **1,200–30,000 total UMI counts** to exclude low-depth captures and overly complex libraries.
- Exclude cells with **>25% mitochondrial reads**, as this indicates poor cell quality or stress.


### Before filtering
![Violin plot before filtering](figures/violin_before_filtering.png?v=2)

### After filtering
![Violin plot after filtering](figures/violin_after_filtering.png?v=2)


## Umap 

We took the filtered AnnData object, normalized it with log-transformation, and stored the log-normalized counts in `adata.raw`. We then selected highly variable genes using the **Seurat v3 flavor**, which closely mimics Seurat’s method of mean-variance trend fitting with Loess regression, making it the closest approach in Scanpy to Seurat’s HVG selection. After that, we scaled the data, performed PCA, computed the neighborhood graph, and finally calculated UMAP embeddings for visualization.


![UMAP](figures/umap_Insm1_analysed.png?v=3)



