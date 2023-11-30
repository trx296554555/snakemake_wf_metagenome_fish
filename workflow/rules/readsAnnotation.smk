rule cat_reads:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
        # host_fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_1.fq",
        # host_fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_2.fq",
    output:
        all_fq=config["root"] + "/" + config["folder"]["reads_annotation_humann3"] + "/{sample}/{sample}_all.fq"
    message:
        "08 : Cat all reads together for humann3 ------------------------------------------"
    threads:
        12
    shell:
        """
        mkdir -p params.opt
        cat {input.fq1} {input.fq2} > {output.all_fq}
        """

rule generate_humann3_report:
    input:
        all_fq=config["root"] + "/" + config["folder"]["reads_annotation_humann3"] + "/{sample}/{sample}_all.fq"
    output:
        res_dir=directory(config["root"] + "/" + config["folder"]["reads_annotation_humann3"] + "/{sample}/{sample}_res")
    message:
        "09 : Run humann3 to generate reads functional annotation -------------------------"
    threads:
        24
    log:
        config["root"] + "/" + config["folder"]["reads_annotation_humann3"] + "/{sample}/{sample}.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "humann3.yaml"
    shell:
        """
        humann --input {input.all_fq} --output {output.res_dir} --threads {threads} \
         --verbose --remove-temp-output  > {log} 2>&1
        """
