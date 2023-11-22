rule generate_kraken2_report:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        kraken2_report=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.kraken",
        kraken2_out=temp(config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.out")
    message:
        "05: Execute Kraken2 for assigning taxonomic labels to DNA sequences ------------------------------------------"
    threads:
        24
    params:
        db=config["db_root"] + "/" + config["db"]["kraken2"]
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kraken2.yaml"
    log:
        config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.log"
    shell:
        """
        kraken2 --paired -db {params.db} --threads {threads} \
        --report {output.kraken2_report} --out {output.kraken2_out} {input.fq1} {input.fq2} > {log} 2>&1
        """


rule generate_bracken_report:
    input:
        kraken2_report=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.kraken"
    output:
        bracken_report=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.braken",
        bracken_out=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.braken.out"
    message:
        "06: Reestimate Kraken2 results at the species level using Bracken ------------------------------------------"
    threads:
        1
    params:
        db=config["db_root"] + "/" + config["db"]["kraken2"],
        reads_len=config["reads_length"],
        tmp_log=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.bracken.log",
        final_log=config["root"] + "/" + config["folder"]["reads_classify"] + "/{sample}.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kraken2.yaml"
    shell:
        """
        sed -i '1i ------- kraken2 log --------' {params.final_log}
        bracken -d {params.db} -i {input.kraken2_report} -o {output.bracken_out} \
        -w {output.bracken_report} -r {params.reads_len} > {params.tmp_log} 2>&1
        echo '------- bracken log --------' | cat - {params.tmp_log} >> {params.final_log} && rm {params.tmp_log}
        """