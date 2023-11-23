# If the reference genome contains MT sequence, can use only host reads
rule get_mitogenome:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
        host_fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_1.fq",
        host_fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_host_2.fq",
    output:
        all_fq1=temp(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}_all_1.fq"),
        all_fq2=temp(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}_all_2.fq"),
        annotation=temp(directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/mt_annotation")),
        assembly=temp(directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/mt_assembly")),
        res_dir=directory(config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}.result")
    message:
        "07: Use mitoZ to obtain the mitochondrial genome of the host from the host sequence ------------------------------------------"
    threads:
        24
    params:
        out_path=config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "mitozEnv.yaml"
    log:
        config["root"] + "/" + config["folder"]["mitogenome"] + "/{sample}/{sample}_mitoz.log"
    shell:
        """
        cat {input.fq1} {input.host_fq1} > {output.all_fq1}
        cat {input.fq2} {input.host_fq2} > {output.all_fq2}
        
        mitoz all --workdir {params.out_path} --outprefix {wildcards.sample} --fq1 {output.all_fq1} --fq2 {output.all_fq2} \
        --skip_filter --data_size_for_mt_assembly 0 --assembler megahit --clade Chordata --requiring_taxa Chordata \
        --insert_size 350 --kmers_megahit 139 --thread_number {threads} > {log} 2>&1 || true
        
        if [ ! -d {params.out_path}/{wildcards.sample}.result ]; then
            mkdir {params.out_path}/{wildcards.sample}.result
            mkdir {params.out_path}/mt_annotation
            touch {params.out_path}/task_failed.error
        else
            mv {output.res_dir}/{wildcards.sample}.megahit.result {params.out_path}
            mv {output.res_dir}/{wildcards.sample}.{wildcards.sample}.megahit.mitogenome.fa.result/* {output.res_dir}
            rm -rf {output.res_dir}/{wildcards.sample}.{wildcards.sample}.megahit.mitogenome.fa.result
        fi
        """