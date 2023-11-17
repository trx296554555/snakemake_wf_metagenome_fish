import os
import pandas as pd

def get_species(wildcards):
    species = all_sample_df[all_sample_df['Sample_ID'] == wildcards.sample]['Species'].values[0]
    res = config["root"] + "/" + config["folder"]["index"] + "/" + str(species).replace(' ','_') + "_bowtie2_index"
    print(res)
    return res

1111
# echo `which python`
# echo `which bowtie2-build`
rule build_index:
    input:
        config["root"] + "/" + config["folder"]["index"] + "/" + "{species}.fna.gz"
    output:
        config["root"] + "/" + config["folder"]["index"] + "/" + "{species}_bowtie2_index"
    threads:
        20
    message:
        "01 make index ------------------------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kneaddata.yaml"
    log:
        config["root"] + "/" + config["folder"]["rm_host"] + "{species}" + "/" + "make_index.log"
    shell:
        "bowtie2-build --threads {threads} {input} {output} > {log} 2>&1"


rule remove_host:
    input:
        fq1=config["root"] + "/" + config["folder"]["data"] + "/{sample}/{sample}_1.clean.fq.gz",
        fq2=config["root"] + "/" + config["folder"]["data"] + "/{sample}/{sample}_2.clean.fq.gz",
        index=get_species
    output:
        config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}_paired_2.fq",
    message:
        "03 run kneaddata to remove host reads ------------------------------------------"
    threads:
        40
    params:
        outdir=config["root"] + "/" + config["folder"]["rm_host"] + "/" + "{sample}",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "kneaddata.yaml"
    log:
        "logs/kneaddata/{sample}.log"
    shell:
        """
            mkdir -p {params.outdir} && cd {params.outdir}
            kneaddata -i1 {input.fq1} -i2 {input.fq2} --reference-db {input.index} \

            --reorder --remove-intermediate-output --threads {threads} \
            -o {params.outdir} --output-prefix {sample} --log {output}/kneaddata.log \
            > {log} 2>&1

            mv {input.fq1} {output.fq1}
            mv {input.fq2} {output.fq2}
            rm -rf *unmatched*
        """
