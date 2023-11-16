import os
import pandas as pd

all_sample_df = pd.read_csv(f'{config["root"]}/{config["meta"]["sampleList"]}')
# def get_species_ref(sample):
#     return all_sample_df[all_sample_df['Sample_ID'] == sample][['Species', 'Download_link']].values[0]
#
# download = [
#     'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/422/415/GCA_027422415.1_SKLP_Slan_1.0/GCA_027422415.1_SKLP_Slan_1.0_genomic.fna.gz',
#     'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/029/165/GCA_023029165.1_ASM2302916v1/GCA_023029165.1_ASM2302916v1_genomic.fna.gz']


# rule download_ref:
#     output:
#         config["root"] + "/" + config["folder"]["index"] + "/" + "{species}.fna.gz"
#     run:
#         wildcards.species = get_species_ref(wildcards.sample)[0]
#         shell("wget -c -O {output} {download}")
#
#
# rule build_index:
#     input:
#         config["root"] + "/" + config["folder"]["index"] + "/" + "{species}.fna.gz"
#     output:
#         config["root"] + "/" + config["folder"]["index"] + "/" + "{species}_bowtie2_index"
#     threads:
#         40
#     message:
#         "01 make index ------------------------------------------"
#     conda:
#         config["root"] + "/" + config["envs"] + "/" + "kneaddata.yaml"
#     log:
#         config["root"] + "/" + config["folder"]["rm_host"] + "{species}" + "/" + "make_index.log"
#     shell:
#         "bowtie2-build --threads {threads} {input} {output} > {log} 2>&1"


# def test(wildcards):
#     print(wildcards.sample)
#     return get_species_ref(wildcards.sample)[0]


rule remove_host:
    input:
        fq1=config["root"] + "/" + config["folder"]["data"] + "/{sample}_1.clean.fq.gz",
        # fq2=config["root"] + "/" +  config["folder"]["data"] + "/{sample}_2.clean.fq.gz",
        # index="123"
        # index=test
    output:
        fq3=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}_paired_1.fq",
        # fq4=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}_paired_2.fq",
        # host_fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_host_1.fq",
        # host_fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_host_2.fq",
    message:
        "03 run kneaddata to remove host reads ------------------------------------------"
    # threads:
    #     40
    # params:
    #     outdir=config["root"] + "/" + config["folder"]["rm_host"] + "/" + "{sample}",
    conda: "../config/envs/kneaddata.yaml"
    # log:
    #     "logs/kneaddata/{sample}.log"
    shell:
        """
        head {input.fq1} > {output.fq3}
        """
    # shell:
    #     """
    #         mkdir -p {params.outdir} && cd {params.outdir}
    #
    #         kneaddata -i1 {input.fq1} -i2 {input.fq2} -o {params.outdir} --output-prefix {sample} --reference-db {input.index}\
    #         --serial --cat-final-output --reorder --remove-intermediate-output\
    #          --threads {threads} --log {output}/kneaddata.log  \
    #         --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" --bowtie2-options "--very-sensitive" > {log} 2>&1
    #
    #         mv {params.outdir}/*_paired_1.fastq {output.fq1}
    #         mv {params.outdir}/*_paired_2.fastq {output.fq2}
    #         mv {params.outdir}/*paired_contam_1.fastq {output.host_fq1}
    #         mv {params.outdir}/*paired_contam_2.fastq {output.host_fq2}
    #         rm -rf *unmatched*
    #     """
