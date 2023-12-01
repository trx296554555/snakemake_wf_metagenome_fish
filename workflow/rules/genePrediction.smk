rule predict_coding_genes:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    output:
        genes=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gene.fa",
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.faa",
    message:
        "10: Predict protein coding genes using prodigal -------------------------"
    threads:
        1
    params:
        gff=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gff",
        stats=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.stats",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "checkm.yaml"
    log:
        config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.log"
    shell:
        """
        prodigal -i {input.contigs} -d {output.genes} -a {output.proteins} \
        -f gff -o {params.gff} -s {params.stats} -p meta > {log} 2>&1
        """