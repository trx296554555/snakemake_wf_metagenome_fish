rule predict_coding_genes:
    input:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa",
    output:
        genes=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gene.fa",
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.faa",
    message:
        "10: Predict protein coding genes using prodigal -------------------------"
    threads:
        12
    params:
        tmp=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/tmp",
        gff=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gff",
        stats=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.stats",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "checkm.yaml"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["gene_prediction"] + "/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.log"
    shell:
        """
        mkdir -p {params.tmp}
        seqkit split -p {threads} {input.contigs} -O {params.tmp} > /dev/null 2>&1
        for i in {params.tmp}/*.fa; do
            prodigal -i $i -d $i.genes -a $i.proteins -p meta -m -q > /dev/null 2>&1 &
        done
        wait
        cat {params.tmp}/*.genes > {output.genes}
        cat {params.tmp}/*.proteins > {output.proteins}
        rm -rf {params.tmp}
        seqkit stats -T {output.proteins} > {log}
        """