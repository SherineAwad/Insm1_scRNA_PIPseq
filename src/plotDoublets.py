import scanpy as sc
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

def main():
    # Parse input file
    parser = argparse.ArgumentParser(description='Create UMAP plot highlighting doublets')
    parser.add_argument('--input', required=True, help='Path to h5ad file with doublet detection results')
    
    args = parser.parse_args()
    
    # Read AnnData object
    print(f"Reading AnnData object from: {args.input}")
    adata = sc.read_h5ad(args.input)
    
    # Get base name for output prefix
    base_name = os.path.splitext(os.path.basename(args.input))[0]
    
    # Check if doublet predictions exist
    if 'predicted_doublet' not in adata.obs.columns:
        raise ValueError("No 'predicted_doublet' column found in the AnnData object")
    
    # Create a copy with only classified cells (remove NaN)
    is_classified = ~adata.obs['predicted_doublet'].isna()
    classified_adata = adata[is_classified].copy()
    
    print(f"Total cells: {adata.shape[0]}")
    print(f"Classified cells: {classified_adata.shape[0]}")
    print(f"Unclassified cells (ignored): {adata.shape[0] - classified_adata.shape[0]}")
    
    # Check if UMAP exists, compute if not
    if 'X_umap' not in classified_adata.obsm:
        print("UMAP not found in object. Computing UMAP...")
        sc.pp.neighbors(classified_adata, n_neighbors=15, n_pcs=50)
        sc.tl.umap(classified_adata)
    
    # Create figures directory if it doesn't exist
    os.makedirs('figures', exist_ok=True)
    
    # FIX: Convert float to int for proper mapping
    # Create doublet_category column with proper mapping for FLOAT values
    classified_adata.obs['doublet_category'] = classified_adata.obs['predicted_doublet'].apply(
        lambda x: 'Doublet' if float(x) == 1.0 else 'Singlet'
    )
    
    # Verify the mapping worked
    doublet_count = (classified_adata.obs['doublet_category'] == 'Doublet').sum()
    singlet_count = (classified_adata.obs['doublet_category'] == 'Singlet').sum()
    
    print(f"\nVerification:")
    print(f"  Doublet cells: {doublet_count}")
    print(f"  Singlet cells: {singlet_count}")
    
    # Set color palette
    color_dict = {
        'Singlet': '#cccccc',  # light grey
        'Doublet': '#1f77b4'   # blue
    }
    
    # Create figure
    plt.figure(figsize=(8, 6))
    
    # Plot UMAP
    ax = sc.pl.umap(
        classified_adata,
        color='doublet_category',
        palette=color_dict,
        size=20,
        title=f'{base_name}\nDoublets: {doublet_count} ({doublet_count/len(classified_adata)*100:.1f}%)',
        frameon=False,
        show=False,
        legend_loc='right margin',
        return_fig=False
    )
    
    # Save figure
    output_file = f"figures/{base_name}_doublets_umap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nUMAP plot saved to: {output_file}")
    plt.close()
    
    # Print statistics
    print(f"\nDoublet Statistics:")
    print(f"  Total classified cells: {len(classified_adata)}")
    print(f"  Singlets: {singlet_count} ({singlet_count/len(classified_adata)*100:.1f}%)")
    print(f"  Doublets: {doublet_count} ({doublet_count/len(classified_adata)*100:.1f}%)")
    print(f"Done.")

if __name__ == '__main__':
    main()
