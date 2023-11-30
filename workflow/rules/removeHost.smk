def get_species(wildcards):
    species = all_sample_df[all_sample_df['Sample_ID'] == wildcards.sample]['Ref'].values[0]
    species = str(species).replace(' ','_')
    res = config["root"] + "/" + config["folder"]["index"] + "/" + f"{species}_bowtie2_index/{species}"
    return res


rule build_index:
    input:
        config["root"] + "/" + config["folder"]["index"] + "/{species}.fna.gz"
    output:
        config["root"] + "/" + config["folder"]["index"] + "/{species}_bowtie2_index"
    threads:
        48
    message:
        "03: Building the index for the reference genome ------------------------------------------"
    log:
        config["root"] + "/" + config["folder"]["index"] + "/{species}_build_index.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kneaddata.yaml"
    shell:
        """
        bowtie2-build --threads {threads} -f {input} {output} > {log} 2>&1
        """

## TODO
# remove_host shell中的缩进有点问题，项目完成后修改
rule remove_host:
    input:
        fq1=config["root"] + "/" + config["folder"]["data"] + "/{sample}/{sample}_1.clean.fq.gz",
        fq2=config["root"] + "/" + config["folder"]["data"] + "/{sample}/{sample}_2.clean.fq.gz",
    output:
        rename_fq1=temp(config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_1.rename.clean.fq"),
        rename_fq2=temp(config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_2.rename.clean.fq"),
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
        host_fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_1.fq",
        host_fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_2.fq",
    message:
        "04: Execute kneaddata to remove host reads ------------------------------------------"
    threads:
        12
    params:
        path = config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}",
        index=get_species
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kneaddata.yaml"
    log:
        config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_run.log"
    shell:
        """
            zcat {input.fq1} | sed 's/ 1.*/\/1/g' > {output.rename_fq1}
            zcat {input.fq2} | sed 's/ 2.*/\/2/g' > {output.rename_fq2}
            kneaddata -i1 {output.rename_fq1} -i2 {output.rename_fq2} --reference-db {params.index} \
            --threads {threads} --max-memory 50g \
            --bypass-trf --reorder --remove-intermediate-output --decontaminate-pairs lenient \
            -o {params.path} --output-prefix {wildcards.sample} --log {params.path}/kneaddata.log \
            > {log} 2>&1

            mv {params.path}/{wildcards.sample}_paired_1.fastq {output.fq1}
            mv {params.path}/{wildcards.sample}_paired_2.fastq {output.fq2}
            mv {params.path}/{wildcards.sample}_*paired*contam_1.fastq {output.host_fq1}
            mv {params.path}/{wildcards.sample}_*paired*contam_2.fastq {output.host_fq2}
            rm -f {params.path}/*unmatched* 
        """


rule report_rm_host:
    input:
        expand(config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/kneaddata.log", sample=get_run_sample())
    output:
        config["root"] + "/" + config["folder"]["reports"] + "/01_rm_host.report"
    params:
        script=config["root"] + "/" + config["scripts"] + "/statistics_kneaddata_res.py"
    shell:
        """
        python {params.script} {input} {output}
        """
