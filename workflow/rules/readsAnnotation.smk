rule generate_kraken2_report:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        kraken2_report=config["root"] + "/" + config["folder"]["kraken2"] + "/{sample}.kraken",
        kraken2_out=temp(config["root"] + "/" + config["folder"]["kraken2"] + "/{sample}.out")
    message:
        "05: Execute Kraken2 for assigning taxonomic labels to DNA sequences ------------------------------------------"
    threads:
        12
    params:
        db=config["db_root"] + "/" + config["db"]["kraken2"]
    conda:
        "kraken2"
    log:
        config["root"] + "/" + config["folder"]["kraken2"] + "/{sample}.log"
    shell:
        """
        kraken2 --paired -db {params.db} --threads {threads} --memory-mapping \
        --report {output.kraken2_report} --out {output.kraken2_out} {input.fq1} {input.fq2} > {log} 2>&1
        """