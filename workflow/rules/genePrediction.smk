import os

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
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prediction.yaml"
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


rule pruning_redundancy:
    input:
        genes=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gene.fa",
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.faa",
    output:
        genes=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gene.prune.fa",
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
    message:
        "11: Pruning redundancy genes -------------------------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prediction.yaml"
    threads:
        6
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["gene_prediction"] + "/{sample}.log"
    log:
        config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}_prune.log"
    shell:
        """
        cd-hit-est -T {threads} -M 0 -i {input.genes} -o {output.genes} -c 0.95 -aS 0.9 -g 1 -sc 1 -sf 1 -d 0 > /dev/null 2>&1
        seqkit seq -n {output.genes} |seqkit grep -n -f - {input.proteins} > {output.proteins}
        sed -i 's/*//g' {output.proteins}
        seqkit stats -T {output.proteins} > {log}
        """


rule report_predict_coding_genes:
    input:
        megahit_logs=expand(
            config["root"] + "/" + config["folder"][
                "gene_prediction"] + "/{sample}/{sample}.log",sample=get_run_sample()),
        co_megahit_logs=expand(
            config["root"] + "/" + config["folder"][
                "gene_prediction"] + "/{item}/{item}.log",item=get_co_item()) if config["co_assemble"]["enable"] else []
    output:
        config["root"] + "/" + config["folder"]["reports"] + "/05_gene_prediction.report"
    run:
        log_list = input.megahit_logs + input.co_megahit_logs
        res_dict = {}
        header = "file\tinfos"
        for l_file in log_list:
            sample = os.path.basename(os.path.dirname(l_file))
            with open(l_file, "r") as f:
                lines = f.readlines()
                header = lines[0].strip()
                info = lines[-1].strip().split("\t")[1:]
                res_dict[sample] = '\t'.join(info)
        with open(output[0], "w") as f:
            if header:
                f.write(header + "\n")
            for sample, info in res_dict.items():
                f.write(sample + "\t" + info + "\n")