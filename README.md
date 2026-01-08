# Analysing INSM1 Overexpression: scRNA-seq (PIP-seq) Data


```
 ___                     _ 
|_ _|_ __  ___ _ __ ___ / |
 | || '_ \/ __| '_ ` _ \| |
 | || | | \__ \ | | | | | |
|___|_| |_|___/_| |_| |_|_|
                           
  ___                                                  _             
 / _ \__   _____ _ __ _____  ___ __  _ __ ___  ___ ___(_) ___  _ __  
| | | \ \ / / _ \ '__/ _ \ \/ / '_ \| '__/ _ \/ __/ __| |/ _ \| '_ \ 
| |_| |\ V /  __/ | |  __/>  <| |_) | | |  __/\__ \__ \ | (_) | | | |
 \___/  \_/ \___|_|  \___/_/\_\ .__/|_|  \___||___/___/_|\___/|_| |_|
                              |_|                                   

```

# 🔬🔬 Benchmark studies of Seurat vs Scanpy: Shows Outputs Can Differ on the Same Input and Similar Parameters


- Rich, J. M., Moses, L., Einarsson, P. H., Jackson, K., Luebbert, L., Booeshaghi, A. S., et al. (2024).  
  *The impact of package selection and versioning on single-cell RNA-seq analysis.* **bioRxiv**.

- Billato, I., Pages, H., Carey, V., Waldron, L., Sales, G., Romualdi, C., & Risso, D. (2025).  
  *Benchmarking large-scale single-cell RNA-seq analysis.* **bioRxiv**.

- Zappia, L., Richter, S., Ramírez-Suástegui, C., Kfuri-Rubens, R., Vornholz, L., Wang, W., et al. (2025).  
  *Feature selection methods affect the performance of scRNA-seq data integration and querying.* **Nature Methods**, 1–11.

# 💡💡💡
# Thanh chose to proceed with this analysis using **Scanpy**.
# Rerunning the analysis using the modified genome: old results can still be restored from github 
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


### Before filtering
![Violin plot before filtering](figures/violin_before_filtering.png?v=2)

### After filtering
![Violin plot after filtering](figures/violin_after_filtering.png?v=2)


### UMAP 

![UMAP](figures/umap_mInsm1_analysed.png?v=3)

## Clustering 

![CLUSTERS](figures/umap_mInsm1_clustered_clusters.png?v=2)

![CLUSTERS ON](figures/umap_mInsm1_clustered_clustersON.png?v=3)

## Marker genes feature plots 

<img src="figures/umap_mInsm1_clustered_Gli1.png?v=3" alt="Gli1" width="33%"><img src="figures/umap_mInsm1_clustered_Malat1.png?v=3" alt="Malat1" width="33%"><img src="figures/umap_mInsm1_clustered_Arr3.png?v=3" alt="Arr3" width="33%">
<img src="figures/umap_mInsm1_clustered_Isl1.png?v=3" alt="Isl1" width="33%"><img src="figures/umap_mInsm1_clustered_Chat.png?v=3" alt="Chat" width="33%"><img src="figures/umap_mInsm1_clustered_Csf1r.png?v=3" alt="Csf1r" width="33%">
<img src="figures/umap_mInsm1_clustered_Ptprc.png?v=3" alt="Ptprc" width="33%"><img src="figures/umap_mInsm1_clustered_Prdx6.png?v=3" alt="Prdx6" width="33%"><img src="figures/umap_mInsm1_clustered_Rlbp1.png?v=3" alt="Rlbp1" width="33%">
<img src="figures/umap_mInsm1_clustered_Trpm1.png?v=3" alt="Trpm1" width="33%"><img src="figures/umap_mInsm1_clustered_Prdm1.png?v=3" alt="Prdm1" width="33%"><img src="figures/umap_mInsm1_clustered_Calb2.png?v=3" alt="Calb2" width="33%">
<img src="figures/umap_mInsm1_clustered_Onecut1.png?v=3" alt="Onecut1" width="33%"><img src="figures/umap_mInsm1_clustered_Rpe65.png?v=3" alt="Rpe65" width="33%"><img src="figures/umap_mInsm1_clustered_Sox9.png?v=3" alt="Sox9" width="33%">
<img src="figures/umap_mInsm1_clustered_Sox2.png?v=3" alt="Sox2" width="33%"><img src="figures/umap_mInsm1_clustered_Hes1.png?v=3" alt="Hes1" width="33%"><img src="figures/umap_mInsm1_clustered_Ascl1.png?v=3" alt="Ascl1" width="33%">
<img src="figures/umap_mInsm1_clustered_Slc18a3.png?v=3" alt="Slc18a3" width="33%"><img src="figures/umap_mInsm1_clustered_Cbln4.png?v=3" alt="Cbln4" width="33%"><img src="figures/umap_mInsm1_clustered_Isl2.png?v=3" alt="Isl2" width="33%">
<img src="figures/umap_mInsm1_clustered_Prox1.png?v=3" alt="Prox1" width="33%"><img src="figures/umap_mInsm1_clustered_mt-Atp6.png?v=3" alt="mt-Atp6" width="33%"><img src="figures/umap_mInsm1_clustered_Lhx4.png?v=3" alt="Lhx4" width="33%">
<img src="figures/umap_mInsm1_clustered_Gad2.png?v=3" alt="Gad2" width="33%"><img src="figures/umap_mInsm1_clustered_Sox11.png?v=3" alt="Sox11" width="33%"><img src="figures/umap_mInsm1_clustered_Pax2.png?v=3" alt="Pax2" width="33%">
<img src="figures/umap_mInsm1_clustered_Pou4f1.png?v=3" alt="Pou4f1" width="33%"><img src="figures/umap_mInsm1_clustered_Rho.png?v=3" alt="Rho" width="33%"><img src="figures/umap_mInsm1_clustered_Thy1.png?v=3" alt="Thy1" width="33%">
<img src="figures/umap_mInsm1_clustered_Neurog2.png?v=3" alt="Neurog2" width="33%"><img src="figures/umap_mInsm1_clustered_Cx3cr1.png?v=3" alt="Cx3cr1" width="33%"><img src="figures/umap_mInsm1_clustered_Onecut2.png?v=3" alt="Onecut2" width="33%">
<img src="figures/umap_mInsm1_clustered_Tfap2a.png?v=3" alt="Tfap2a" width="33%"><img src="figures/umap_mInsm1_clustered_Lhx1.png?v=3" alt="Lhx1" width="33%"><img src="figures/umap_mInsm1_clustered_Pax6.png?v=3" alt="Pax6" width="33%">
<img src="figures/umap_mInsm1_clustered_Otx2.png?v=3" alt="Otx2" width="33%"><img src="figures/umap_mInsm1_clustered_Hes5.png?v=3" alt="Hes5" width="33%"><img src="figures/umap_mInsm1_clustered_Grm6.png?v=3" alt="Grm6" width="33%">
<img src="figures/umap_mInsm1_clustered_Sebox.png?v=3" alt="Sebox" width="33%"><img src="figures/umap_mInsm1_clustered_Kcnj8.png?v=3" alt="Kcnj8" width="33%"><img src="figures/umap_mInsm1_clustered_Sall1.png?v=3" alt="Sall1" width="33%">
<img src="figures/umap_mInsm1_clustered_Slc17a7.png?v=3" alt="Slc17a7" width="33%"><img src="figures/umap_mInsm1_clustered_Abca8a.png?v=3" alt="Abca8a" width="33%"><img src="figures/umap_mInsm1_clustered_Prkca.png?v=3" alt="Prkca" width="33%">
<img src="figures/umap_mInsm1_clustered_Rbpms.png?v=3" alt="Rbpms" width="33%"><img src="figures/umap_mInsm1_clustered_Olig2.png?v=3" alt="Olig2" width="33%"><img src="figures/umap_mInsm1_clustered_Nefm.png?v=3" alt="Nefm" width="33%">
<img src="figures/umap_mInsm1_clustered_Tie1.png?v=3" alt="Tie1" width="33%"><img src="figures/umap_mInsm1_clustered_Slc1a3.png?v=3" alt="Slc1a3" width="33%"><img src="figures/umap_mInsm1_clustered_Nrl.png?v=3" alt="Nrl" width="33%">
<img src="figures/umap_mInsm1_clustered_Cabp5.png?v=3" alt="Cabp5" width="33%"><img src="figures/umap_mInsm1_clustered_Gfap.png?v=3" alt="Gfap" width="33%"><img src="figures/umap_mInsm1_clustered_Sncg.png?v=3" alt="Sncg" width="33%">
<img src="figures/umap_mInsm1_clustered_Ebf3.png?v=3" alt="Ebf3" width="33%"><img src="figures/umap_mInsm1_clustered_Foxn4.png?v=3" alt="Foxn4" width="33%"><img src="figures/umap_mInsm1_clustered_Apoe.png?v=3" alt="Apoe" width="33%">
<img src="figures/umap_mInsm1_clustered_Vim.png?v=3" alt="Vim" width="33%"><img src="figures/umap_mInsm1_clustered_Bsn.png?v=3" alt="Bsn" width="33%"><img src="figures/umap_mInsm1_clustered_Insm1.png?v=3" alt="Insm1" width="33%">
<img src="figures/umap_mInsm1_clustered_Nefl.png?v=3" alt="Nefl" width="33%"><img src="figures/umap_mInsm1_clustered_Pou4f2.png?v=3" alt="Pou4f2" width="33%"><img src="figures/umap_mInsm1_clustered_Sfrp2.png?v=3" alt="Sfrp2" width="33%">
<img src="figures/umap_mInsm1_clustered_Emx1.png?v=3" alt="Emx1" width="33%"><img src="figures/umap_mInsm1_clustered_Notch1.png?v=3" alt="Notch1" width="33%"><img src="figures/umap_mInsm1_clustered_Rbfox3.png?v=3" alt="Rbfox3" width="33%">
<img src="figures/umap_mInsm1_clustered_Pecam1.png?v=3" alt="Pecam1" width="33%"><img src="figures/umap_mInsm1_clustered_Pdgfra.png?v=3" alt="Pdgfra" width="33%"><img src="figures/umap_mInsm1_clustered_Pou4f3.png?v=3" alt="Pou4f3" width="33%">
<img src="figures/umap_mInsm1_clustered_Acta2.png?v=3" alt="Acta2" width="33%"><img src="figures/umap_mInsm1_clustered_Lhx2.png?v=3" alt="Lhx2" width="33%"><img src="figures/umap_mInsm1_clustered_Igf2.png?v=3" alt="Igf2" width="33%">
<img src="figures/umap_mInsm1_clustered_Elavl3.png?v=3" alt="Elavl3" width="33%"><img src="figures/umap_mInsm1_clustered_Bhlhe23.png?v=3" alt="Bhlhe23" width="33%"><img src="figures/umap_mInsm1_clustered_Ccr2.png?v=3" alt="Ccr2" width="33%">
<img src="figures/umap_mInsm1_clustered_Vsx1.png?v=3" alt="Vsx1" width="33%"><img src="figures/umap_mInsm1_clustered_Aqp4.png?v=3" alt="Aqp4" width="33%"><img src="figures/umap_mInsm1_clustered_Slc6a9.png?v=3" alt="Slc6a9" width="33%">
<img src="figures/umap_mInsm1_clustered_Elavl4.png?v=3" alt="Elavl4" width="33%"><img src="figures/umap_mInsm1_clustered_Atoh7.png?v=3" alt="Atoh7" width="33%"><img src="figures/umap_mInsm1_clustered_Calb1.png?v=3" alt="Calb1" width="33%">
<img src="figures/umap_mInsm1_clustered_Fgf15.png?v=3" alt="Fgf15" width="33%"><img src="figures/umap_mInsm1_clustered_Glul.png?v=3" alt="Glul" width="33%"><img src="figures/umap_mInsm1_clustered_Pcp4.png?v=3" alt="Pcp4" width="33%">
<img src="figures/umap_mInsm1_clustered_Gad1.png?v=3" alt="Gad1" width="33%"><img src="figures/umap_mInsm1_clustered_Top2a.png?v=3" alt="Top2a" width="33%"><img src="figures/umap_mInsm1_clustered_Pcna.png?v=3"  alt="Pcna" width="33%">
<img src="figures/umap_mInsm1_clustered_Pecam1.png?v=3" alt="Pecam1" width="33%">

## Annotations 

![Annotations](figures/umap_mInsm1_annotated_annotations.png?v=3)

![Annotations ON](figures/umap_mInsm1_annotated_annotationsON.png?v=3)


# Going back a step to remove doublets from original clustered object

### DoubletDetection threshold plot explanation

- **Y-axis (Voting threshold)**  
  Fraction of classifier runs that must label a cell as a doublet.  
  Higher values = stricter, fewer cells called doublets.  
  Example: **0.8** means a cell must be called a doublet in ≥80% of runs.

- **X-axis (log10 p-value cutoff)**  
  Statistical stringency within each run.  
  More negative values = more stringent (fewer doublets).  
  Less negative values = more permissive (more doublets).

- **Color scale**  
  Number of cells predicted as doublets for each parameter combination.

- **Interpretation of 0.8 cutoff**  
  Selecting **0.8** yields a conservative and stable set of doublets  
  (≈ ~2,000–2,500 cells), avoiding over-filtering seen at lower thresholds.


![doublet_threshold 0.8](figures/0.8threshold_test.png?v=1)



### Trying different doublets thresholds 

| Threshold | Total Classified Cells | Singlets | Singlets (%) | Doublets | Doublets (%) |
|-----------|----------------------|----------|--------------|----------|--------------|
| 0.7       | 57,607               | 53,918   | 93.6%        | 3,689    | 6.4%         |
| 0.8       | 57,607               | 54,114   | 93.9%        | 3,493    | 6.1%         |
| 0.9       | 57,611               | 54,697   | 94.9%        | 2,914    | 5.1%         |


<img src="figures/mInsm1_doubletDetected0.7_doublets_umap.png?v=1" alt="P7" width ="75%"> 
<img src="figures/mInsm1_doubletDetected0.8_doublets_umap.png?v=1" alt="P8" width ="75%">
<img src="figures/mInsm1_doubletDetected0.9_doublets_umap.png?v=1" alt="P9" width ="75%">


### UMAP for clusters using doublet threshold 0.8

![doublet_threshold 0.8 UMAP](figures/umapdoubletsRemoved_0.8._0.5_clusters.png?v=2)

### Marker genes UMAP using doublet threshold 0.8 

<img src="figures/umapdoubletRemoved_0.8._0.5_Rho.png?v=2" alt="Rho" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Thy1.png?v=2" alt="Thy1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Opn1mw.png?v=2" alt="Opn1mw" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Calb2.png?v=2" alt="Calb2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Pou4f3.png?v=2" alt="Pou4f3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_EGFP.png?v=2" alt="EGFP" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Nrl.png?v=2" alt="Nrl" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Rbfox3.png?v=2" alt="Rbfox3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Chat.png?v=2" alt="Chat" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Ebf3.png?v=2" alt="Ebf3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Prdx6.png?v=2" alt="Prdx6" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Vim.png?v=2" alt="Vim" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Elavl3.png?v=2" alt="Elavl3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sall3.png?v=2" alt="Sall3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Cbln4.png?v=2" alt="Cbln4" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Pax6.png?v=2" alt="Pax6" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Abca8a.png?v=2" alt="Abca8a" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Neurod2.png?v=2" alt="Neurod2" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Slc18a3.png?v=2" alt="Slc18a3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Kcnj8.png?v=2" alt="Kcnj8" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Prox1.png?v=2" alt="Prox1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Notch1.png?v=2" alt="Notch1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Elavl4.png?v=2" alt="Elavl4" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sebox.png?v=2" alt="Sebox" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Nefm.png?v=2" alt="Nefm" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sox9.png?v=2" alt="Sox9" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Otx2.png?v=2" alt="Otx2" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Isl1.png?v=2" alt="Isl1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Olig2.png?v=2" alt="Olig2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Gli1.png?v=2" alt="Gli1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Slc1a3.png?v=2" alt="Slc1a3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Igf2.png?v=2" alt="Igf2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Ascl1.png?v=2" alt="Ascl1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Pdgfra.png?v=2" alt="Pdgfra" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Vsx1.png?v=2" alt="Vsx1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Apoe.png?v=2" alt="Apoe" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Neurod4.png?v=2" alt="Neurod4" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Prdm1.png?v=2" alt="Prdm1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sncg.png?v=2" alt="Sncg" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Cdk1.png?v=2" alt="Cdk1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Pcna.png?v=2" alt="Pcna" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Hes1.png?v=2" alt="Hes1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Gfap.png?v=2" alt="Gfap" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Calb1.png?v=2" alt="Calb1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Ptprc.png?v=2" alt="Ptprc" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Malat1.png?v=2" alt="Malat1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Pcp4.png?v=2" alt="Pcp4" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Gad2.png?v=2" alt="Gad2" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Emx1.png?v=2" alt="Emx1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Isl2.png?v=2" alt="Isl2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Glul.png?v=2" alt="Glul" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Arr3.png?v=2" alt="Arr3" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Rbpms.png?v=2" alt="Rbpms" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Nefl.png?v=2" alt="Nefl" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Top2a.png?v=2" alt="Top2a" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_mt-Atp6.png?v=2" alt="mt-Atp6" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Slc6a9.png?v=2" alt="Slc6a9" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Grm6.png?v=2" alt="Grm6" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Tfap2a.png?v=2" alt="Tfap2a" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Pou4f1.png?v=2" alt="Pou4f1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Rlbp1.png?v=2" alt="Rlbp1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Trpm1.png?v=2" alt="Trpm1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Cx3cr1.png?v=2" alt="Cx3cr1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Pax2.png?v=2" alt="Pax2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Mki67.png?v=2" alt="Mki67" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Hes6.png?v=2" alt="Hes6" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Gnat2.png?v=2" alt="Gnat2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Insm1.png?v=2" alt="Insm1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Csf1r.png?v=2" alt="Csf1r" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Foxn4.png?v=2" alt="Foxn4" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Lhx2.png?v=2" alt="Lhx2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Hes5.png?v=2" alt="Hes5" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Pecam1.png?v=2" alt="Pecam1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Pou4f2.png?v=2" alt="Pou4f2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sall1.png?v=2" alt="Sall1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Aqp4.png?v=2" alt="Aqp4" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Lhx4.png?v=2" alt="Lhx4" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Bhlhe23.png?v=2" alt="Bhlhe23" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Onecut1.png?v=2" alt="Onecut1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Tie1.png?v=2" alt="Tie1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Slc17a7.png?v=2" alt="Slc17a7" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Prkca.png?v=2" alt="Prkca" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Atoh7.png?v=2" alt="Atoh7" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Rpe65.png?v=2" alt="Rpe65" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Acta2.png?v=2" alt="Acta2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Neurod1.png?v=2" alt="Neurod1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Neurog2.png?v=2" alt="Neurog2" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Gad1.png?v=2" alt="Gad1" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Ptf1a.png?v=2" alt="Ptf1a" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Lhx1.png?v=2" alt="Lhx1" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Opn1sw.png?v=2" alt="Opn1sw" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Cabp5.png?v=2" alt="Cabp5" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Onecut2.png?v=2" alt="Onecut2" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Fgf15.png?v=2" alt="Fgf15" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sfrp2.png?v=2" alt="Sfrp2" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sox2.png?v=2" alt="Sox2" width="33%">
<img src="figures/umapdoubletRemoved_0.8._0.5_Bsn.png?v=2" alt="Bsn" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Sox11.png?v=2" alt="Sox11" width="33%"><img src="figures/umapdoubletRemoved_0.8._0.5_Ccr2.png?v=2" alt="Ccr2" width="33%">



### Annotaions 

![Annotations doubelt removed 0.8](figures/umap_mInsm1_annotated_annotationsON.png?v=1)
![Annotations doubelt removed 0.8](figures/umap_mInsm1_annotated_annotations.png?v=1)

### Remove non retinal neurons and reCluster
We removed:
- RPE
- Astrocytes
- Microglia
- Endothelia


Then we reClustered


