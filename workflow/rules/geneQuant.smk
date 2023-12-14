rule index_gene:
    input:
        config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gene.prune.fa"
    output:
        directory(config["root"] + "/" + config["folder"]["gene_quant"] + "/{sample}/{sample}.salmon.idx"),
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        12
    log:
        config["root"] + "/" + config["folder"]["gene_quant"] + "/{sample}.index.log"
    params:
        tmpdir=temp(directory(config["root"] + "/" + config["folder"]["gene_quant"] + "/{sample}/tmp"))
    message:
        " 9. Index gene expression using salmon"
    shell:
        """
        salmon index -t {input} -i {output} -p {threads} --tmpdir {params.tmpdir} 2>&1 > {log}
        """

rule quant_gene:
    input:
        idx=directory(config["root"] + "/" + config["folder"]["gene_quant"] + "/{sample}/{sample}.salmon.idx"),
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        dir=temp(directory(config["root"] + "/" + config["folder"]["gene_quant"] + "/{sample}")),
        sf=config["root"] + "/" + config["folder"]["gene_quant"] + "/{sample}/{sample}.sf",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        12
    log:
        config["root"] + "/" + config["folder"]["gene_quant"] + "/{sample}.quant.log"
    message:
        " 10. Quantify gene expression using salmon"
    shell:
        """
        salmon quant -i {input.idx} -l A --meta -1 {input.fq1} -2 {input.fq2} \
         -o {output.dir} -p {threads} --validateMappings --seqBias --quiet 2>&1 > {log}
         mv {output.dir}/quant.sf {output.sf}
        """
