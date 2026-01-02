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


## Remove non retinal neurons and reCluster 


### We removed: 
  
  - RPE 
  - Astrocytes
  - Microglia 
  - Endothelia 

Then we reClustered 

### reClustering UMAP 

![reClustered](figures/umap_mInsm1_reClustered_UMAP_recluster.png?v=3)

### reClustered UMAP split by sample 

![reClustered by sample](figures/umap_mInsm1_reClustered_by_sample.png?v=3)


### reClustered UMAP QC violin plot 

![QC](figures/violin_mInsm1_reClustered_QC_violin_by_cluster.png?v=3)



### Marker genes' feature plots 

<img src="figures/umap_mInsm1_reClustered_Abca8a.png?v=3" alt="Abca8a" width="33%"><img src="figures/umap_mInsm1_reClustered_Csf1r.png?v=3" alt="Csf1r" width="33%"><img src="figures/umap_mInsm1_reClustered_Hes5.png?v=3" alt="Hes5" width="33%">
<img src="figures/umap_mInsm1_reClustered_Neurod4.png?v=3" alt="Neurod4" width="33%"><img src="figures/umap_mInsm1_reClustered_Pou4f2.png?v=3" alt="Pou4f2" width="33%"><img src="figures/umap_mInsm1_reClustered_Slc17a7.png?v=3" alt="Slc17a7" width="33%">
<img src="figures/umap_mInsm1_reClustered_Acta2.png?v=3" alt="Acta2" width="33%"><img src="figures/umap_mInsm1_reClustered_Cx3cr1.png?v=3" alt="Cx3cr1" width="33%"><img src="figures/umap_mInsm1_reClustered_Hes6.png?v=3" alt="Hes6" width="33%">
<img src="figures/umap_mInsm1_reClustered_Neurog2.png?v=3" alt="Neurog2" width="33%"><img src="figures/umap_mInsm1_reClustered_Pou4f3.png?v=3" alt="Pou4f3" width="33%"><img src="figures/umap_mInsm1_reClustered_Slc18a3.png?v=3" alt="Slc18a3" width="33%">
<img src="figures/umap_mInsm1_reClustered_Apoe.png?v=3" alt="Apoe" width="33%"><img src="figures/umap_mInsm1_reClustered_Ebf3.png?v=3" alt="Ebf3" width="33%"><img src="figures/umap_mInsm1_reClustered_Igf2.png?v=3" alt="Igf2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Notch1.png?v=3" alt="Notch1" width="33%"><img src="figures/umap_mInsm1_reClustered_Prdm1.png?v=3" alt="Prdm1" width="33%"><img src="figures/umap_mInsm1_reClustered_Slc1a3.png?v=3" alt="Slc1a3" width="33%">
<img src="figures/umap_mInsm1_reClustered_Aqp4.png?v=3" alt="Aqp4" width="33%"><img src="figures/umap_mInsm1_reClustered_EGFP.png?v=3" alt="EGFP" width="33%"><img src="figures/umap_mInsm1_reClustered_Insm1.png?v=3" alt="Insm1" width="33%">
<img src="figures/umap_mInsm1_reClustered_Nrl.png?v=3" alt="Nrl" width="33%"><img src="figures/umap_mInsm1_reClustered_Prdx6.png?v=3" alt="Prdx6" width="33%"><img src="figures/umap_mInsm1_reClustered_Slc6a9.png?v=3" alt="Slc6a9" width="33%">
<img src="figures/umap_mInsm1_reClustered_Arr3.png?v=3" alt="Arr3" width="33%"><img src="figures/umap_mInsm1_reClustered_Elavl3.png?v=3" alt="Elavl3" width="33%"><img src="figures/umap_mInsm1_reClustered_Isl1.png?v=3" alt="Isl1" width="33%">
<img src="figures/umap_mInsm1_reClustered_Olig2.png?v=3" alt="Olig2" width="33%"><img src="figures/umap_mInsm1_reClustered_Prkca.png?v=3" alt="Prkca" width="33%"><img src="figures/umap_mInsm1_reClustered_Sncg.png?v=3" alt="Sncg" width="33%">
<img src="figures/umap_mInsm1_reClustered_Ascl1.png?v=3" alt="Ascl1" width="33%"><img src="figures/umap_mInsm1_reClustered_Elavl4.png?v=3" alt="Elavl4" width="33%"><img src="figures/umap_mInsm1_reClustered_Isl2.png?v=3" alt="Isl2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Onecut1.png?v=3" alt="Onecut1" width="33%"><img src="figures/umap_mInsm1_reClustered_Prox1.png?v=3" alt="Prox1" width="33%"><img src="figures/umap_mInsm1_reClustered_Sox11.png?v=3" alt="Sox11" width="33%">
<img src="figures/umap_mInsm1_reClustered_Atoh7.png?v=3" alt="Atoh7" width="33%"><img src="figures/umap_mInsm1_reClustered_Emx1.png?v=3" alt="Emx1" width="33%"><img src="figures/umap_mInsm1_reClustered_Kcnj8.png?v=3" alt="Kcnj8" width="33%">
<img src="figures/umap_mInsm1_reClustered_Onecut2.png?v=3" alt="Onecut2" width="33%"><img src="figures/umap_mInsm1_reClustered_Ptf1a.png?v=3" alt="Ptf1a" width="33%"><img src="figures/umap_mInsm1_reClustered_Sox2.png?v=3" alt="Sox2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Bhlhe23.png?v=3" alt="Bhlhe23" width="33%"><img src="figures/umap_mInsm1_reClustered_Fgf15.png?v=3" alt="Fgf15" width="33%"><img src="figures/umap_mInsm1_reClustered_Lhx1.png?v=3" alt="Lhx1" width="33%">
<img src="figures/umap_mInsm1_reClustered_Opn1mw.png?v=3" alt="Opn1mw" width="33%"><img src="figures/umap_mInsm1_reClustered_Ptprc.png?v=3" alt="Ptprc" width="33%"><img src="figures/umap_mInsm1_reClustered_Sox9.png?v=3" alt="Sox9" width="33%">
<img src="figures/umap_mInsm1_reClustered_Bsn.png?v=3" alt="Bsn" width="33%"><img src="figures/umap_mInsm1_reClustered_Foxn4.png?v=3" alt="Foxn4" width="33%"><img src="figures/umap_mInsm1_reClustered_Lhx2.png?v=3" alt="Lhx2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Opn1sw.png?v=3" alt="Opn1sw" width="33%"><img src="figures/umap_mInsm1_reClustered_Rbfox3.png?v=3" alt="Rbfox3" width="33%"><img src="figures/umap_mInsm1_reClustered_Tfap2a.png?v=3" alt="Tfap2a" width="33%">
<img src="figures/umap_mInsm1_reClustered_Gad1.png?v=3" alt="Gad1" width="33%"><img src="figures/umap_mInsm1_reClustered_Lhx4.png?v=3" alt="Lhx4" width="33%"><img src="figures/umap_mInsm1_reClustered_Otx2.png?v=3" alt="Otx2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Rbpms.png?v=3" alt="Rbpms" width="33%"><img src="figures/umap_mInsm1_reClustered_Thy1.png?v=3" alt="Thy1" width="33%"><img src="figures/umap_mInsm1_reClustered_Cabp5.png?v=3" alt="Cabp5" width="33%">
<img src="figures/umap_mInsm1_reClustered_Gad2.png?v=3" alt="Gad2" width="33%"><img src="figures/umap_mInsm1_reClustered_Malat1.png?v=3" alt="Malat1" width="33%"><img src="figures/umap_mInsm1_reClustered_Pax2.png?v=3" alt="Pax2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Rho.png?v=3" alt="Rho" width="33%"><img src="figures/umap_mInsm1_reClustered_Tie1.png?v=3" alt="Tie1" width="33%"><img src="figures/umap_mInsm1_reClustered_Calb1.png?v=3" alt="Calb1" width="33%">
<img src="figures/umap_mInsm1_reClustered_Gfap.png?v=3" alt="Gfap" width="33%"><img src="figures/umap_mInsm1_reClustered_Mki67.png?v=3" alt="Mki67" width="33%"><img src="figures/umap_mInsm1_reClustered_Pax6.png?v=3" alt="Pax6" width="33%">
<img src="figures/umap_mInsm1_reClustered_Rlbp1.png?v=3" alt="Rlbp1" width="33%"><img src="figures/umap_mInsm1_reClustered_Top2a.png?v=3" alt="Top2a" width="33%"><img src="figures/umap_mInsm1_reClustered_Calb2.png?v=3" alt="Calb2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Gli1.png?v=3" alt="Gli1" width="33%"><img src="figures/umap_mInsm1_reClustered_mt-Atp6.png?v=3" alt="mt-Atp6" width="33%"><img src="figures/umap_mInsm1_reClustered_Pcna.png?v=3" alt="Pcna" width="33%">
<img src="figures/umap_mInsm1_reClustered_Rpe65.png?v=3" alt="Rpe65" width="33%"><img src="figures/umap_mInsm1_reClustered_Trpm1.png?v=3" alt="Trpm1" width="33%"><img src="figures/umap_mInsm1_reClustered_Cbln4.png?v=3" alt="Cbln4" width="33%">
<img src="figures/umap_mInsm1_reClustered_Glul.png?v=3" alt="Glul" width="33%"><img src="figures/umap_mInsm1_reClustered_Nefl.png?v=3" alt="Nefl" width="33%"><img src="figures/umap_mInsm1_reClustered_Pcp4.png?v=3" alt="Pcp4" width="33%">
<img src="figures/umap_mInsm1_reClustered_Sall1.png?v=3" alt="Sall1" width="33%"><img src="figures/umap_mInsm1_reClustered_Ccr2.png?v=3" alt="Ccr2" width="33%"><img src="figures/umap_mInsm1_reClustered_Gnat2.png?v=3" alt="Gnat2" width="33%">
<img src="figures/umap_mInsm1_reClustered_Nefm.png?v=3" alt="Nefm" width="33%"><img src="figures/umap_mInsm1_reClustered_Pdgfra.png?v=3" alt="Pdgfra" width="33%"><img src="figures/umap_mInsm1_reClustered_Sall3.png?v=3" alt="Sall3" width="33%">
<img src="figures/umap_mInsm1_reClustered_Vim.png?v=3" alt="Vim" width="33%"><img src="figures/umap_mInsm1_reClustered_Cdk1.png?v=3" alt="Cdk1" width="33%"><img src="figures/umap_mInsm1_reClustered_Grm6.png?v=3" alt="Grm6" width="33%">
<img src="figures/umap_mInsm1_reClustered_Neurod1.png?v=3" alt="Neurod1" width="33%"><img src="figures/umap_mInsm1_reClustered_Pecam1.png?v=3" alt="Pecam1" width="33%"><img src="figures/umap_mInsm1_reClustered_Sebox.png?v=3" alt="Sebox" width="33%">
<img src="figures/umap_mInsm1_reClustered_Vsx1.png?v=3" alt="Vsx1" width="33%"><img src="figures/umap_mInsm1_reClustered_Chat.png?v=3" alt="Chat" width="33%"><img src="figures/umap_mInsm1_reClustered_Hes1.png?v=3" alt="Hes1" width="33%">
<img src="figures/umap_mInsm1_reClustered_Neurod2.png?v=3" alt="Neurod2" width="33%"><img src="figures/umap_mInsm1_reClustered_Pou4f1.png?v=3" alt="Pou4f1" width="33%"><img src="figures/umap_mInsm1_reClustered_Sfrp2.png?v=3" alt="Sfrp2" width="33%">
i


# Going back a step to remove doublets from original clustered object


### Using threshold 0.8 

![doublet_threshold 0.8](figures/0.8threshold_test.png?v=1)



### Using threshold 0.9



![doublet_threshold 0.9](figures/0.9threshold_test.png?v=1) 
