rule assemble_contigs:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq"
    output:
        contigs=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}.contigs.fa"
    message:
        "08: Assemble contigs using megahit -------------------------"
    threads:
        24
    params:
        opt=config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}_megahit"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "megahit.yaml"
    log:
        config["root"] + "/" + config["folder"]["assemble_contigs"] + "/{sample}/{sample}_run.log"
    shell:
        """
        rm -rf {params.opt}
        megahit -1 {input.fq1} -2 {input.fq2} -m 0.95 --min-contig-len 300 \
        --k-list 79,99,119,139 -t {threads} \
        --out-dir {params.opt} --out-prefix {wildcards.sample} > {log} 2>&1
        mv {params.opt}/{wildcards.sample}.contigs.fa {output.contigs}
        """