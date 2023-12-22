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
        "12.2: Merge humann3 reports -------------------------"
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


rule get_rgi_annotation:
    input:
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
    output:
        rgi_anno=config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}/{sample}.rgi.anno",
    message:
        "12.5: Run rgi to get ARGs annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "rgi.yaml"
    threads:
        12
    params:
        # alignment tool: choose from 'diamond', 'blast'
        tool="diamond",
        opt_prefix=config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}/{sample}",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_anno_rgi"] + "/{sample}.rgi.log"
    log:
        config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}/{sample}.rgi.log"
    shell:
        """
        rgi main -i {input.proteins} -o {params.opt_prefix} \
        -t protein -a {params.tool} -n {threads} --clean --debug --include_loose > {log} 2>&1
        rm -rf {params.opt_prefix}.json
        mv {params.opt_prefix}.txt {output.rgi_anno}
        """


rule generate_rgi_report:
    input:
        rgi_anno=config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}/{sample}.rgi.anno",
        expression=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.sf"
    output:
        rgi_tsv=config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/{sample}/{sample}.rgi.tsv"
    message:
        "12.6: Combine quantification and annotation results to generate rgi reports------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -q {input.expression} -a {input.rgi_anno} -t rgi -o {output.rgi_tsv}
        """


rule merge_rgi_report:
    input:
        all_rgi_tsv=expand(config["root"] + "/" + config["folder"][
            "reads_anno_rgi"] + "/{sample}/{sample}.rgi.tsv",sample=get_run_sample()),
    output:
        merge_rgi_tsv=config["root"] + "/" + config["folder"]["reads_anno_rgi"] + "/all.rgi.tsv",
    message:
        "12.7: Merge rgi reports -------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -a {input.all_rgi_tsv} -t rgi -o {output.merge_rgi_tsv}
        """


rule get_dbcan_annotation:
    input:
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
    output:
        dbcan_anno=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.anno",
    message:
        "12.8: Run dbcan to get CAZy annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "dbcan4.yaml"
    threads:
        12
    params:
        # tools: choose from 'hmmer', 'diamond', 'dbcansub', 'all';
        # use two or more tools, use ' ' to separate them, for example: tools="hmmer diamond"
        # dbcansub will take a lot of space uncontrolled, use with caution
        opt_dir=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}",
        tools="all",
        db=config["db_root"] + "/" + config["db"]["dbcan"],
        type="protein",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_anno_dbcan"] + "/{sample}.dbcan.log"
    log:
        config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.log"
    shell:
        """
        run_dbcan {input.proteins} {params.type} --out_dir {params.opt_dir} --out_pre {wildcards.sample}_ \
        --hmm_cpu {threads} --dia_cpu {threads} --tf_cpu {threads} --stp_cpu {threads} -dt {threads} \
        --db_dir {params.db} --tools {params.tools} > {log} 2>&1
        mv {params.opt_dir}/{wildcards.sample}_overview.txt {output.dbcan_anno}
        """


rule generate_dbcan_report:
    input:
        dbcan_anno=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.anno",
        expression=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.sf"
    output:
        dbcan_tsv=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.tsv"
    message:
        "12.9: Combine quantification and annotation results to generate dbcan reports------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -q {input.expression} -a {input.dbcan_anno} -t dbcan -o {output.dbcan_tsv}
        """


rule merge_dbcan_report:
    input:
        all_dbcan_tsv=expand(config["root"] + "/" + config["folder"][
            "reads_anno_dbcan"] + "/{sample}/{sample}.dbcan.tsv",sample=get_run_sample()),
    output:
        merge_dbcan_tsv=config["root"] + "/" + config["folder"]["reads_anno_dbcan"] + "/all.dbcan.tsv",
    message:
        "12.10: Merge dbcan reports -------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -a {input.all_dbcan_tsv} -t dbcan -o {output.merge_dbcan_tsv}
        """


rule get_vfdb_annotation:
    input:
        proteins=config["root"] + "/" + config["folder"]["gene_prediction"] + "/{sample}/{sample}.protein.prune.faa",
    output:
        vfdb_anno=config["root"] + "/" + config["folder"]["reads_anno_vfdb"] + "/{sample}/{sample}.vfdb.anno",
    message:
        "12.11: Run vfdb to generate virulence factors annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "dbcan4.yaml"
    threads:
        12
    params:
        db=config["db_root"] + "/" + config["db"]["vfdb"],
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["reads_anno_vfdb"] + "/{sample}.vfdb.log"
    log:
        config["root"] + "/" + config["folder"]["reads_anno_vfdb"] + "/{sample}/{sample}.vfdb.log"
    shell:
        """
        diamond blastp -d {params.db} -q {input.proteins} -o {output.vfdb_anno} \
        -f 6 -p {threads} --id 80% --query-cover 70% --evalue 1e-5 > {log} 2>&1
        """


rule generate_vfdb_report:
    input:
        vfdb_anno=config["root"] + "/" + config["folder"]["reads_anno_vfdb"] + "/{sample}/{sample}.vfdb.anno",
        expression=config["root"] + "/" + config["folder"]["reads_gene_quant"] + "/{sample}/{sample}.sf"
    output:
        vfdb_tsv=config["root"] + "/" + config["folder"]["reads_anno_vfdb"] + "/{sample}/{sample}.vfdb.tsv"
    message:
        "12.12: Combine quantification and annotation results to generate vfdb reports------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -q {input.expression} -a {input.vfdb_anno} -t vfdb -o {output.vfdb_tsv}
        """


rule merge_vfdb_report:
    input:
        all_vfdb_tsv=expand(config["root"] + "/" + config["folder"][
            "reads_anno_vfdb"] + "/{sample}/{sample}.vfdb.tsv",sample=get_run_sample()),
    output:
        merge_vfdb_tsv=config["root"] + "/" + config["folder"]["reads_anno_vfdb"] + "/all.vfdb.tsv",
    message:
        "12.13: Merge vfdb reports -------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -a {input.all_vfdb_tsv} -t vfdb -o {output.merge_vfdb_tsv}
        """
