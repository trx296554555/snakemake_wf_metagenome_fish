import glob


rule generate_humann3_report:
    input:
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        all_fq=temp(config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_all.fq"),
        gf_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_genefamilies.tsv",
        path_abd_tsv=config["root"] + "/" + config["folder"][
            "reads_anno_humann3"] + "/{sample}/{sample}_pathabundance.tsv",
        path_cov_tsv=config["root"] + "/" + config["folder"][
            "reads_anno_humann3"] + "/{sample}/{sample}_pathcoverage.tsv",
    message:
        "12 : Run humann3 to generate reads functional annotation -------------------------"
    threads:
        24
    params:
        res_dir=directory(config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}"),
    log:
        config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "humann3.yaml"
    shell:
        """
        mkdir -p {params.res_dir}
        cat {input.fq1} {input.fq2} > {output.all_fq}
        humann --input {output.all_fq} --output {params.res_dir} --threads {threads} \
        --verbose --remove-temp-output > {log} 2>&1
        mv {params.res_dir}/{wildcards.sample}_all_genefamilies.tsv {output.gf_tsv}
        mv {params.res_dir}/{wildcards.sample}_all_pathabundance.tsv {output.path_abd_tsv}
        mv {params.res_dir}/{wildcards.sample}_all_pathcoverage.tsv {output.path_cov_tsv}
        """


# rule regroup_gene_family:
#     input:
#         gf_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_genefamilies.tsv",
#     output:
#         rxn_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_rxn.tsv",
#         ko_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_ko.tsv",
#         go_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_go.tsv",
#     message:
#         "12 : Regroup humann3 reports -------------------------"
#     shell:
#         """
#         """


rule merge_humann3_res:
    input:
        all_gf_tsv=expand(
            config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_genefamilies.tsv",
            sample=get_run_sample()),
    output:
        merger_gf_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/all_genefamilies.tsv",
    message:
        "12 : Merge humann3 results -------------------------"
    shell:
        """
        humann_join_tables -i humann_out/ -o ./humann3_genefamilies.tsv \
        --file_name genefamilies; humann_join_tables -i humann_out/ -o ./
        """



def gather_dbcan_report(wildcards):
    opt_dir = checkpoints.run_dbcan.get(**wildcards).output.opt
    opt_reports = glob.glob(opt_dir + "/*overview.txt")
    return opt_reports


checkpoint run_dbcan:
    input:
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
    output:
        opt=directory(config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}")
    message:
        "12 : Run dbcan to generate CAZy annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "annotation.yaml"
    threads:
        24
    params:
        db=config["db_root"] + "/" + config["db"]["dbcan"],
        type="protein",
        # tools: choose from 'hmmer', 'diamond', 'dbcansub', 'all';
        # use two or more tools, use ' ' to separate them, for example: tools="hmmer diamond"
        # dbcansub will take a lot of space uncontrolled, so we don't use it
        tools="hmmer diamond"
    log:
        config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.log"
    shell:
        """
        run_dbcan --out_dir {output} --db_dir {params.db}  {input.proteins} {params.type} \
        --hmm_cpu {threads} --dia_cpu {threads} --tf_cpu {threads} --stp_cpu {threads} -dt {threads} \
        --tools {params.tools} --out_dir {output.opt} --out_pre {wildcards.sample}_ 2>&1 > {log}
        """

rule dbcan_merge:
    input:
        gather_dbcan_report
    output:
        config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.done"
    message:
        "12 : Merge dbcan reports -------------------------"
    params:
        dir=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/results"
    shell:
        """
        mkdir -p {params.dir}
        mv {input} {params.dir}
        touch {output}
        """


rule run_rgi:
    input:
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
    output:
        opt=config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}_rgi.txt"
    message:
        "12 : Run rgi to generate ARGs annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "annotation.yaml"
    threads:
        24
    params:
        opt_prefix=config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}",
        # alignment tool: choose from 'diamond', 'blast'
        tool="diamond",
    log:
        config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}/{sample}.rgi.log"
    shell:
        """
        rgi main -i {input.proteins} -o {params.opt_prefix} \
        -t protein -a {params.tool} -n {threads} --clean --debug --include_loose 2>&1 > {log}
        rm -rf {params.opt_prefix}.json
        mv {params.opt_prefix}.txt {output.opt}
        """


rule run_vfdb:
    input:
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
    output:
        opt=config["root"] + "/" + config["folder"]["reads_anno_vfdb"] + "/{sample}_vfdb.txt"
    message:
        "12 : Run vfdb to generate virulence factors annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "annotation.yaml"
    threads:
        24
    params:
        db=config["db_root"] + "/" + config["db"]["vfdb"],
    log:
        config["root"] + "/" + config["folder"]["reads_anno_vfdb"] + "/{sample}/{sample}.vfdb.log"
    shell:
        """
        diamond blastp -d {params.db} -q {input.proteins} -o {output.opt}/{wildcards.sample}.diamond.out \
        -f 6 -p {threads} --id 80% --query-cover 70% --evalue 1e-3 2>&1 > {log}
        """
