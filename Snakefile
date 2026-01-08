configfile: "config.yaml"

rule all:
    input:
        f"{config['output_prefix']}_clustered.h5ad",
        "figures/umap_mInsm1_clustered_by_sample.png",
        f"{config['output_prefix']}_doubletDetected{config['DOUBLET_THRESHOLD']}.h5ad",
        f"{config['output_prefix']}_doubletRemoved{config['DOUBLET_THRESHOLD']}_{config['DOUBLET_SCORE']}.h5ad",
        f"{config['output_prefix']}_annotated.h5ad",

rule read_matrices:
    input:
        samples=config["samples"]
    output:
        h5ad=f"{config['output_prefix']}.h5ad"
    conda:
        "sknb.yaml"
    shell:
        """
        python src/read_matrices.py \
            --samples {input.samples} \
            --output {output.h5ad}
        """

rule filter:
    input:
        h5ad=f"{config['output_prefix']}.h5ad"
    output:
        h5ad=f"{config['output_prefix']}_filtered.h5ad"
    conda:
        "sknb.yaml"
    shell:
        """
        python src/filter.py \
            --input {input.h5ad} \
            --output {output.h5ad}
        """

rule analyse:
    input:
        h5ad=f"{config['output_prefix']}_filtered.h5ad"
    output:
        h5ad=f"{config['output_prefix']}_analysed.h5ad"
    conda:
        "sknb.yaml"
    shell:
        """
        python src/analyse.py \
            --input {input.h5ad} \
            --output {output.h5ad}
        """

rule cluster:
    input: 
       h5ad=f"{config['output_prefix']}_analysed.h5ad"
    output: 
       h5ad=f"{config['output_prefix']}_clustered.h5ad"
    conda: 
      "sknb.yaml"
    params: 
        resolution = config['resolution'] 
    shell: 
        """
        python src/cluster.py --input {input} --output {output} --resolution {params} 
        """ 

rule plotMarkers: 
      input: 
         h5ad=f"{config['output_prefix']}_clustered.h5ad"
      params: 
        markers = config['MARKERS'] 
      output: 
       "figures/umap_mInsm1_clustered_by_sample.png"
      shell: 
          """
           python src/plotMarkers.py --input {input} --markers {params} 
          """ 


rule detectDoublets:
    input:
        h5ad=f"{config['output_prefix']}_clustered.h5ad"
    params:
        threshold=config['DOUBLET_THRESHOLD']  # this goes to the shell command
    output:
        f"{config['output_prefix']}_doubletDetected{config['DOUBLET_THRESHOLD']}.h5ad"  # output must be static
    shell:
        """
        python src/detectDoublets.py --input {input.h5ad} --output {output} --threshold {params.threshold}
        """

rule removeDoublets: 
     input: 
         f"{config['output_prefix']}_doubletDetected{config['DOUBLET_THRESHOLD']}.h5ad"
     output: 
        f"{config['output_prefix']}_doubletRemoved{config['DOUBLET_THRESHOLD']}_{config['DOUBLET_SCORE']}.h5ad" 
     params: 
        score=config['DOUBLET_SCORE']
     shell: 
        """
          python src/removeDoublets.py --input {input} --output {output} --markers markers.txt --threshold {params.score} 
        """


rule annotate:
     input:
        f"{config['output_prefix']}_doubletRemoved{config['DOUBLET_THRESHOLD']}_{config['DOUBLET_SCORE']}.h5ad"
     params:
       annotations = config['ANNOTATIONS']
     output:
       h5ad=f"{config['output_prefix']}_annotated.h5ad"
     shell:
       """
       python src/annotate.py --input {input} --output {output} --annotations {params}
       """


rule reCluster:
    input:
      f"{config['output_prefix']}_annotated.h5ad"
    output:
      f"{config['output_prefix']}_reClustered.h5ad"
    shell:
       """
       python src/reCluster.py --input {input} --output {output}
       """

rule rePlot:
    input:
       f"{config['output_prefix']}_reClustered.h5ad"
    output:
       "figures/umap_mInsm1_reClustered_by_sample.png"
    params:
       markers = config['MARKERS']
    shell:
       """
       python src/plotMarkers.py --input {input} --markers {params}
       """
