# Analysing INSM1 Overexpression: scRNA-seq (PIP-seq) Data


# 🔬🔬 Benchmark studies of Seurat vs Scanpy: Shows Outputs Can Differ on the Same Input and Similar Parameters


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
![Violin plot before filtering](ofigures/violin_before_filtering.png?v=2)

### After filtering
![Violin plot after filtering](ofigures/violin_after_filtering.png?v=2)


## ## UMAP (Generated with the Parameters Below — Tweaking May Change Results)
 
- **Normalization:** `sc.pp.normalize_total(adata, target_sum=1e4)`  
- **Log-transformation:** `sc.pp.log1p(adata)`  
- **Raw counts stored:** `adata.raw = adata`  

**Highly Variable Genes (HVG) Selection:**  
- Flavor: `seurat_v3` (falls back to `seurat` if `skmisc` not installed)  
- Number of top genes: `2000`  

**Scaling:**  
- `sc.pp.scale(adata, max_value=10)`  

**PCA:**  
- `sc.tl.pca(adata, svd_solver='arpack')`  

**Neighborhood Graph:**  
- `sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)`  

**UMAP Computation:**  
- `sc.tl.umap(adata)`  

*Note: Tweaking any of these parameters may produce slightly different UMAP embeddings — this is expected and does not indicate errors.*


# Rerunning the analysis using the modified genome 

### Before filtering
![Violin plot before filtering](figures/violin_before_filtering.png?v=2)

### After filtering
![Violin plot after filtering](figures/violin_after_filtering.png?v=2)


### UMAP 

![UMAP](figures/umap_mInsm1_analysed.png?v=3)

## Clustering 

![CLUSTERS](figures/umap_mInsm1_clustered_clusters.png?v=2)


## Marker genes feature plots 

<img src="figures/umap_mInsm1_clustered_Prdx6.png?v=2" alt="Prdx6" width="33%"><img src="figures/umap_mInsm1_clustered_Sox9.png?v=2" alt="Sox9" width="33%"><img src="figures/umap_mInsm1_clustered_Rho.png?v=2" alt="Rho" width="33%">
<img src="figures/umap_mInsm1_clustered_Pax6.png?v=2" alt="Pax6" width="33%"><img src="figures/umap_mInsm1_clustered_Slc17a7.png?v=2" alt="Slc17a7" width="33%"><img src="figures/umap_mInsm1_clustered_Slc1a3.png?v=2" alt="Slc1a3" width="33%">
<img src="figures/umap_mInsm1_clustered_Apoe.png?v=2" alt="Apoe" width="33%"><img src="figures/umap_mInsm1_clustered_clusters.png?v=2" alt="clusters" width="33%"><img src="figures/umap_mInsm1_clustered_Lhx2.png?v=2" alt="Lhx2" width="33%">
<img src="figures/umap_mInsm1_clustered_Aqp4.png?v=2" alt="Aqp4" width="33%"><img src="figures/umap_mInsm1_clustered_Glul.png?v=2" alt="Glul" width="33%"><img src="figures/umap_mInsm1_clustered_Malat1.png?v=2" alt="Malat1" width="33%">
<img src="figures/umap_mInsm1_clustered_Rlbp1.png?v=2" alt="Rlbp1" width="33%"><img src="figures/umap_mInsm1_clustered_mt-Atp6.png?v=2" alt="mt-Atp6" width="33%"><img src="figures/umap_mInsm1_clustered_Otx2.png?v=2" alt="Otx2" width="33%">
<img src="figures/umap_mInsm1_clustered_Abca8a.png?v=2" alt="Abca8a" width="33%"><img src="figures/umap_mInsm1_clustered_Nrl.png?v=2" alt="Nrl" width="33%"><img src="figures/umap_mInsm1_clustered_Vim.png?v=2" alt="Vim" width="33%">
<img src="figures/umap_mInsm1_clustered_Notch1.png?v=2" alt="Notch1" width="33%"><img src="figures/umap_mInsm1_clustered_Slc6a9.png?v=2" alt="Slc6a9" width="33%"><img src="figures/umap_mInsm1_clustered_Arr3.png?v=2" alt="Arr3" width="33%">
<img src="figures/umap_mInsm1_clustered_Hes1.png?v=2" alt="Hes1" width="33%"><img src="figures/umap_mInsm1_clustered_Lhx4.png?v=2" alt="Lhx4" width="33%"><img src="figures/umap_mInsm1_clustered_Neurog2.png?v=2" alt="Neurog2" width="33%">
<img src="figures/umap_mInsm1_clustered_Hes5.png?v=2" alt="Hes5" width="33%"><img src="figures/umap_mInsm1_clustered_Cabp5.png?v=2" alt="Cabp5" width="33%"><img src="figures/umap_mInsm1_clustered_Bsn.png?v=2" alt="Bsn" width="33%">
<img src="figures/umap_mInsm1_clustered_Rbfox3.png?v=2" alt="Rbfox3" width="33%"><img src="figures/umap_mInsm1_clustered_Elavl3.png?v=2" alt="Elavl3" width="33%"><img src="figures/umap_mInsm1_clustered_Elavl4.png?v=2" alt="Elavl4" width="33%">
<img src="figures/umap_mInsm1_clustered_Isl1.png?v=2" alt="Isl1" width="33%"><img src="figures/umap_mInsm1_clustered_Prdm1.png?v=2" alt="Prdm1" width="33%"><img src="figures/umap_mInsm1_clustered_Ascl1.png?v=2" alt="Ascl1" width="33%">
<img src="figures/umap_mInsm1_clustered_Gfap.png?v=2" alt="Gfap" width="33%"><img src="figures/umap_mInsm1_clustered_Insm1.png?v=2" alt="Insm1" width="33%"><img src="figures/umap_mInsm1_clustered_Atoh7.png?v=2" alt="Atoh7" width="33%">
<img src="figures/umap_mInsm1_clustered_Gad1.png?v=2" alt="Gad1" width="33%"><img src="figures/umap_mInsm1_clustered_Chat.png?v=2" alt="Chat" width="33%"><img src="figures/umap_mInsm1_clustered_Calb2.png?v=2" alt="Calb2" width="33%">
<img src="figures/umap_mInsm1_clustered_Slc18a3.png?v=2" alt="Slc18a3" width="33%"><img src="figures/umap_mInsm1_clustered_Sox11.png?v=2" alt="Sox11" width="33%"><img src="figures/umap_mInsm1_clustered_Sebox.png?v=2" alt="Sebox" width="33%">
<img src="figures/umap_mInsm1_clustered_Olig2.png?v=2" alt="Olig2" width="33%"><img src="figures/umap_mInsm1_clustered_Ccr2.png?v=2" alt="Ccr2" width="33%"><img src="figures/umap_mInsm1_clustered_Calb1.png?v=2" alt="Calb1" width="33%">
<img src="figures/umap_mInsm1_clustered_Emx1.png?v=2" alt="Emx1" width="33%"><img src="figures/umap_mInsm1_clustered_Csf1r.png?v=2" alt="Csf1r" width="33%"><img src="figures/umap_mInsm1_clustered_Pax2.png?v=2" alt="Pax2" width="33%">
<img src="figures/umap_mInsm1_clustered_Tfap2a.png?v=2" alt="Tfap2a" width="33%"><img src="figures/umap_mInsm1_clustered_Kcnj8.png?v=2" alt="Kcnj8" width="33%"><img src="figures/umap_mInsm1_clustered_Pou4f2.png?v=2" alt="Pou4f2" width="33%">
<img src="figures/umap_mInsm1_clustered_Rpe65.png?v=2" alt="Rpe65" width="33%"><img src="figures/umap_mInsm1_clustered_Lhx1.png?v=2" alt="Lhx1" width="33%"><img src="figures/umap_mInsm1_clustered_Tie1.png?v=2" alt="Tie1" width="33%">
<img src="figures/umap_mInsm1_clustered_Foxn4.png?v=2" alt="Foxn4" width="33%"><img src="figures/umap_mInsm1_clustered_Acta2.png?v=2" alt="Acta2" width="33%">




