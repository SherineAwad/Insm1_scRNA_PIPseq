configfile: "config.yaml"

rule all:
    input:
        f"{config['output_prefix']}_analysed.h5ad"

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

