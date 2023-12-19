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
        "12.1: Run humann3 to generate reads functional annotation -------------------------"
    threads:
        12
    params:
        res_dir=directory(config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}"),
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_anno_humann3"] + "/{sample}.humann3.log"
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


rule merge_humann3_report:
    input:
        all_gf_tsv=expand(
            config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/{sample}/{sample}_genefamilies.tsv",
            sample=get_run_sample()),
        all_path_abd_tsv=expand(config["root"] + "/" + config["folder"][
            "reads_anno_humann3"] + "/{sample}/{sample}_pathabundance.tsv",sample=get_run_sample()),
    output:
        merge_gf_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/all_genefamilies.tsv",
        merge_path_abd_tsv=config["root"] + "/" + config["folder"][
            "reads_anno_humann3"] + "/all_pathabundance.tsv",
        rxn_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/all_genefamilies_rxn.tsv",
        ko_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/all_genefamilies_ko.tsv",
        go_tsv=config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/all_genefamilies_go.tsv",
    message:
        "12.2: Merge humann3 results -------------------------"
    threads:
        1
    params:
        res_dir=directory(config["root"] + "/" + config["folder"]["reads_anno_humann3"])
    log:
        config["root"] + "/" + config["folder"]["reads_anno_humann3"] + "/merge.log"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "humann3.yaml"
    shell:
        """
        humann_join_tables -s -i {params.res_dir} -o {output.merge_gf_tsv} --file_name genefamilies > {log} 2>&1
        humann_join_tables -s -i {params.res_dir} -o {output.merge_path_abd_tsv} --file_name pathabundance >> {log} 2>&1
        sed -i '1s/_all_Abundance-RPKs//g' {output.merge_gf_tsv}
        sed -i '1s/_all_Abundance//g' {output.merge_path_abd_tsv}
        humann_regroup_table -i {output.merge_gf_tsv} -o {output.rxn_tsv} --groups uniref90_rxn 2>&1 >> {log} 2>&1
        humann_regroup_table -i {output.merge_gf_tsv} -o {output.ko_tsv} --groups uniref90_ko 2>&1 >> {log} 2>&1
        humann_regroup_table -i {output.merge_gf_tsv} -o {output.go_tsv} --groups uniref90_go 2>&1 >> {log} 2>&1
        """


rule salmon_index:
    input:
        config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.gene.prune.fa"
    output:
        directory(config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.salmon.idx"),
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        6
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_gene_quant"] + "/{sample}.index.log"
    log:
        config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.index.log"
    params:
        tmpdir=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/tmp"
    message:
        "12.3: Index prune gene from contig using salmon"
    shell:
        """
        salmon index -t {input} -i {output} -p {threads} --tmpdir {params.tmpdir} > {log} 2>&1
        """


rule quantify_gene_expression:
    input:
        idx=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.salmon.idx",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        expression=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.sf",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        12
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_gene_quant"] + "/{sample}.quant.log"
    log:
        config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.quant.log"
    message:
        "12.4: Quantify gene expression using salmon"
    params:
        dir=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}"
    shell:
        """
        salmon quant -i {input.idx} -l A --meta -1 {input.fq1} -2 {input.fq2} \
        -o {params.dir} -p {threads} --validateMappings --seqBias --quiet > {log} 2>&1
        mv {params.dir}/quant.sf {output.expression}
        """


rule generate_dbcan_report:
    input:
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
        expression=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.sf"
    output:
        opt_dir=directory(config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}"),
        dbcan_anno=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.anno",
        dbcan_tsv=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.tsv"
    message:
        "12.5: Run dbcan to generate CAZy annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "annotation.yaml"
    threads:
        12
    params:
        # tools: choose from 'hmmer', 'diamond', 'dbcansub', 'all';
        # use two or more tools, use ' ' to separate them, for example: tools="hmmer diamond"
        # dbcansub will take a lot of space uncontrolled, use with caution
        tools="all",
        db=config["db_root"] + "/" + config["db"]["dbcan"],
        type="protein",
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_anno_dbcan"] + "/{sample}.dbcan.log"
    log:
        config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.log"
    shell:
        """
        run_dbcan {input.proteins} {params.type} --out_dir {output.opt_dir} --out_pre {wildcards.sample}_ \
        --hmm_cpu {threads} --dia_cpu {threads} --tf_cpu {threads} --stp_cpu {threads} -dt {threads} \
        --db_dir {params.db} --tools {params.tools} > {log} 2>&1
        mv {output.opt_dir}/{wildcards.sample}_overview.txt {output.dbcan_anno}
        """


# rule dbcan_merge:
#     input:
#         gather_dbcan_report
#     output:
#         config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.done"
#     message:
#         "12 : Merge dbcan reports -------------------------"
#     params:
#         dir=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/results"
#     shell:
#         """
#         mkdir -p {params.dir}
#         mv {input} {params.dir}
#         touch {output}
#         """


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
