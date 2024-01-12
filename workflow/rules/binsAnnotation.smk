## TODO 统计MAGs 大小，N50，编码蛋白数量
# new_MAGs_report = config["root"] + "/" + config["folder"]["reports"] + "/08_new_MAGs.report",
rule get_MAGs_prokka_anno:
    input:
        mag=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins/{bin}.fa",
        all_classify_tsv=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.all.summary.tsv"
    output:
        genes=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.gene.fa",
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.protein.faa",
        prokka_anno=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}/{bin}.prokka.txt",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prokka.yaml"
    message:
        "20: Run prokka to get bins annotation"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_anno_prokka"] + "/{bin}.prokka.log"
    threads:
        12
    params:
        mags_dir=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins",
        out_dir=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/{bin}",
    shell:
        """
        mkdir -p {params.out_dir}
        raw_id=$(basename {input.mag} .fa)
        new_id=$(printf "%04d" $(ls $(dirname {input.mag})| nl | awk '$2=="'$raw_id.fa'" {{print $1}}'))
        kingdom=$(cat {input.all_classify_tsv} | grep $raw_id$'\t' | cut -f3 | cut -d \; -f1 | sed 's/d__//')
        genus=$(cat {input.all_classify_tsv} | grep $raw_id$'\t' | cut -f3 | cut -d \; -f6 | sed 's/g__//')
        species=$(cat {input.all_classify_tsv} | grep $raw_id$'\t' | cut -f3 | cut -d \; -f7 | sed 's/s__//'|sed 's/ /_/g')
        if [ $genus ]; then genus=$genus;else genus="''"; fi
        if [ $species ]; then species=$species;else species="''"; fi
        ln -sf {input.mag} {params.out_dir}/$new_id.fa
        
        prokka --quiet --cpus {threads} --outdir {params.out_dir} --prefix {wildcards.bin} \
        --metagenome --force --gcode 11 --kingdom ${{kingdom}} --genus ${{genus}} --species ${{species}} \
        --locustag fish_bin_${{new_id}} --centre '' --compliant {input.mag} > /dev/null 2>&1
        seqkit replace --quiet -p ".*" -r ${{raw_id}}_{{nr}} -w 0 {params.out_dir}/{wildcards.bin}.ffn > {output.genes}
        seqkit replace --quiet -p ".*" -r ${{raw_id}}_{{nr}} -w 0 {params.out_dir}/{wildcards.bin}.faa > {output.proteins}
        ln -s {params.out_dir}/{wildcards.bin}.txt {output.prokka_anno}
        """


def collect_MAGs(wildcards):
    mag_bins_dir = checkpoints.get_MAGs.get(**wildcards).output.mags_dir
    mag_bins = glob.glob(mag_bins_dir + "/*.f*a")
    mag_names = [os.path.basename(x).replace('.fa','') for x in mag_bins]
    prokka_gene = [os.path.join(config["root"],config["folder"]["bins_anno_prokka"],bin,bin + ".gene.fa") for bin in
                   mag_names]
    prokka_prptein = [os.path.join(config["root"],config["folder"]["bins_anno_prokka"],bin,bin + ".protein.faa")
                      for bin in mag_names]
    return {"gene": prokka_gene, "protein": prokka_prptein}


rule get_MAGs_prokka_anno_res:
    input:
        unpack(collect_MAGs),
    output:
        all_mags_gene=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_gene.fa",
        all_mags_protein=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_protein.faa",
    shell:
        """
        cat {input.gene} > {output.all_mags_gene}
        cat {input.protein} > {output.all_mags_protein}
        """


rule pruning_MAGs_gene_redundancy:
    input:
        all_mags_gene=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_gene.fa",
        all_mags_protein=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_protein.faa",
    output:
        genes=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_gene.prune.fa",
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_protein.prune.faa",
        cluster_file=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_gene.prune.fa.clstr",
    message:
        "21: Pruning redundancy genes from all mag bins"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "prediction.yaml"
    threads:
        80
    params:
        clean_all_mags_gene=config["root"] + "/" + config["folder"][
            "bins_anno_prokka"] + "/all_mags/all_mags_gene.clean.fa",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_anno_prokka"] + "/pruning_MAGs_gene_redundancy.log"
    log:
        config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/pruning_MAGs_gene_redundancy.log"
    shell:
        """
        seqkit seq --threads {threads} -n {input.all_mags_protein} |seqkit grep --threads {threads} -n -f - {input.all_mags_gene} > {params.clean_all_mags_gene}
        cd-hit-est -T {threads} -M 0 -i {params.clean_all_mags_gene} -o {output.genes} -c 0.95 -aS 0.9 -g 1 -sc 1 -sf 1 -d 0 > /dev/null 2>&1
        seqkit seq --threads {threads} -n {output.genes} |seqkit grep --threads {threads} -n -f - {input.all_mags_protein} > {output.proteins}
        sed -i 's/*//g' {output.proteins}
        seqkit stats --threads {threads} -T {output.proteins} > {log}
        """


rule get_MAGs_gene_salmon_index:
    input:
        config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_gene.prune.fa"
    output:
        directory(config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags.salmon.idx"),
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        80
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_gene_quant"] + "/all_mags.salmon_index.log"
    log:
        config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags.salmon_index.log"
    params:
        tmpdir=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/tmp"
    message:
        "22: Index prune gene from all MAGs using salmon"
    shell:
        """
        salmon index -t {input} -i {output} -p {threads} --tmpdir {params.tmpdir} > {log} 2>&1
        """


rule quantify_MAGs_gene_expression:
    input:
        idx=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags.salmon.idx",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        expression=config["root"] + "/" + config["folder"]["bins_gene_quant"] + "/{sample}/{sample}.sf",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        12
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_gene_quant"] + "/{sample}.quant.log"
    log:
        config["root"] + "/" + config["folder"]["bins_gene_quant"] + "/{sample}/{sample}.quant.log"
    message:
        "23: Quantify gene expression using salmon"
    params:
        dir=config["root"] + "/" + config["folder"]["bins_gene_quant"] + "/{sample}"
    shell:
        """
        salmon quant -i {input.idx} -l A --meta -1 {input.fq1} -2 {input.fq2} \
        -o {params.dir} -p {threads} --validateMappings --seqBias --quiet > {log} 2>&1
        mv {params.dir}/quant.sf {output.expression}
        """


rule get_MAGs_rgi_annotation:
    input:
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_protein.prune.faa",
    output:
        rgi_anno=config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/annotation/all_mags.rgi.anno",
    message:
        "24.1: Run rgi to get ARGs annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "rgi.yaml"
    threads:
        80
    params:
        # alignment tool: choose from 'diamond', 'blast'
        tool="diamond",
        opt_prefix=config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/annotation/all_mags",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_anno_rgi"] + "/all_mags.rgi.annotation.log"
    log:
        config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/annotation/all_mags.rgi.annotation.log"
    shell:
        """
        rgi main -i {input.proteins} -o {params.opt_prefix} \
        -t protein -a {params.tool} -n {threads} --clean --debug --include_loose > {log} 2>&1
        rm -rf {params.opt_prefix}.json
        mv {params.opt_prefix}.txt {output.rgi_anno}
        """


rule generate_MAGs_rgi_report:
    input:
        rgi_anno=config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/annotation/all_mags.rgi.anno",
        expression=config["root"] + "/" + config["folder"]["bins_gene_quant"] + "/{sample}/{sample}.sf",
    output:
        rgi_tsv=config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/{sample}/{sample}.rgi.tsv"
    message:
        "24.2: Combine quantification and annotation results to generate rgi reports------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -q {input.expression} -a {input.rgi_anno} -t rgi -o {output.rgi_tsv}
        """


rule merge_MAGs_rgi_report:
    input:
        all_rgi_tsv=expand(config["root"] + "/" + config["folder"][
            "bins_anno_rgi"] + "/{sample}/{sample}.rgi.tsv",sample=get_run_sample()),
    output:
        merge_rgi_tsv=config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/all.rgi.tsv",
    message:
        "24.3: Merge rgi reports -------------------------"
    params:
        script=config["root"] + "/workflow/scripts/join_tables.py",
        tmp_dir=config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/tmp"
    shell:
        """
        mkdir -p {params.tmp_dir}
        for res in {input.all_rgi_tsv}; do ln -sf $res {params.tmp_dir}; done
        python {params.script} -i {params.tmp_dir} -o {output.merge_rgi_tsv}
        rm -rf {params.tmp_dir}
        """


rule get_MAGs_dbcan_annotation:
    input:
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_protein.prune.faa",
    output:
        dbcan_anno=config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/annotation/all_mags.dbcan.anno",
    message:
        "25.1: Run dbcan to get CAZy annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "dbcan4.yaml"
    threads:
        80
    params:
        # tools: choose from 'hmmer', 'diamond', 'dbcansub', 'all';
        # use two or more tools, use ' ' to separate them, for example: tools="hmmer diamond"
        # dbcansub will take a lot of space uncontrolled, use with caution
        opt_dir=config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/annotation",
        opt_prefix="all_mags",
        tools="all",
        db=config["db_root"] + "/" + config["db"]["dbcan"],
        type="protein",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_anno_dbcan"] + "/all_mags.dbcan.annotation.log"
    log:
        config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/annotation/all_mags.dbcan.annotation.log"
    shell:
        """
        run_dbcan {input.proteins} {params.type} --out_dir {params.opt_dir} --out_pre {params.opt_prefix}_ \
        --hmm_cpu {threads} --dia_cpu {threads} --tf_cpu {threads} --stp_cpu {threads} -dt {threads} \
        --db_dir {params.db} --tools {params.tools} > {log} 2>&1
        mv {params.opt_dir}/{params.opt_prefix}_overview.txt {output.dbcan_anno}
        """


rule generate_MAGs_dbcan_report:
    input:
        dbcan_anno=config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/annotation/all_mags.dbcan.anno",
        expression=config["root"] + "/" + config["folder"]["bins_gene_quant"] + "/{sample}/{sample}.sf",
    output:
        dbcan_tsv=config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/{sample}/{sample}.dbcan.tsv"
    message:
        "25.2: Combine quantification and annotation results to generate dbcan reports------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -q {input.expression} -a {input.dbcan_anno} -t dbcan -o {output.dbcan_tsv}
        """


rule merge_MAGs_dbcan_report:
    input:
        all_dbcan_tsv=expand(config["root"] + "/" + config["folder"][
            "bins_anno_dbcan"] + "/{sample}/{sample}.dbcan.tsv",sample=get_run_sample()),
    output:
        merge_dbcan_tsv=config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/all.dbcan.tsv",
    message:
        "25.3: Merge dbcan reports -------------------------"
    params:
        script=config["root"] + "/workflow/scripts/join_tables.py",
        tmp_dir=config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/tmp"
    shell:
        """
        mkdir -p {params.tmp_dir}
        for res in {input.all_dbcan_tsv}; do ln -sf $res {params.tmp_dir}; done
        python {params.script} -i {params.tmp_dir} -o {output.merge_dbcan_tsv}
        rm -rf {params.tmp_dir}
        """


rule get_MAGs_vfdb_annotation:
    input:
        proteins=config["root"] + "/" + config["folder"]["bins_anno_prokka"] + "/all_mags/all_mags_protein.prune.faa",
    output:
        vfdb_anno=config["root"] + "/" + config["folder"]["bins_anno_vfdb"] + "/annotation/all_mags.vfdb.anno",
    message:
        "26.1: Run vfdb to generate virulence factors annotation -------------------------"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "dbcan4.yaml"
    threads:
        80
    params:
        db=config["db_root"] + "/" + config["db"]["vfdb"],
        diamond="VF_full_db.dmnd",
        id="80%", cover="70%", evalue="1e-5",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_anno_vfdb"] + "/all_mags.vfdb.annotation.log"
    log:
        config["root"] + "/" + config["folder"]["bins_anno_vfdb"] + "/annotation/all_mags.vfdb.annotation.log"
    shell:
        """
        diamond blastp -d {params.db}/{params.diamond} -q {input.proteins} -o {output.vfdb_anno} \
        -f 6 -p {threads} --id {params.id} --query-cover {params.cover} --evalue {params.evalue} > {log} 2>&1
        """


rule generate_MAGs_vfdb_report:
    input:
        vfdb_anno=config["root"] + "/" + config["folder"]["bins_anno_vfdb"] + "/annotation/all_mags.vfdb.anno",
        expression=config["root"] + "/" + config["folder"]["bins_gene_quant"] + "/{sample}/{sample}.sf"
    output:
        vfdb_tsv=config["root"] + "/" + config["folder"]["bins_anno_vfdb"] + "/{sample}/{sample}.vfdb.tsv"
    message:
        "26.2: Combine quantification and annotation results to generate vfdb reports------------------------"
    params:
        script=config["root"] + "/workflow/scripts/salmon_anno_integration.py"
    shell:
        """
        python {params.script} -q {input.expression} -a {input.vfdb_anno} -t vfdb -o {output.vfdb_tsv}
        """


rule merge_MAGs_vfdb_report:
    input:
        all_vfdb_tsv=expand(config["root"] + "/" + config["folder"][
            "bins_anno_vfdb"] + "/{sample}/{sample}.vfdb.tsv",sample=get_run_sample()),
    output:
        merge_vfdb_tsv=config["root"] + "/" + config["folder"]["bins_anno_vfdb"] + "/all.vfdb.tsv",
    message:
        "26.3: Merge vfdb reports -------------------------"
    params:
        script=config["root"] + "/workflow/scripts/join_tables.py",
        tmp_dir=config["root"] + "/" + config["folder"]["bins_anno_vfdb"] + "/tmp"
    shell:
        """
        mkdir -p {params.tmp_dir}
        for res in {input.all_vfdb_tsv}; do ln -sf $res {params.tmp_dir}; done
        python {params.script} -i {params.tmp_dir} -o {output.merge_vfdb_tsv}
        rm -rf {params.tmp_dir}
        """


def get_MAGs_annotation_res():
    anno_res = {}
    if config["bins_anno"]["rgi_enable"]:
        anno_res['rgi'] = config["root"] + "/" + config["folder"]["bins_anno_rgi"] + "/all.rgi.tsv"
    if config["bins_anno"]["dbcan_enable"]:
        anno_res['dbcan'] = config["root"] + "/" + config["folder"]["bins_anno_dbcan"] + "/all.dbcan.tsv"
    if config["bins_anno"]["vfdb_enable"]:
        anno_res['vfdb'] = config["root"] + "/" + config["folder"]["bins_anno_vfdb"] + "/all.vfdb.tsv"
    return anno_res


rule report_MAGs_annotation:
    input:
        **get_MAGs_annotation_res()
    output:
        config["root"] + "/" + config["folder"]["reports"] + "/10_bin_annotation.report"
    run:
        anno_info_path = config["root"] + "/" + config["folder"]["reports"] + "/10_bin_annotation_info"
        info_dict = {
            'rgi': config["db_root"] + "/" + config["db"]["rgi"] + "/rgi_card_info.zip",
            'dbcan': config["db_root"] + "/" + config["db"]["dbcan"] + "/fam-substrate-mapping.tsv",
            'vfdb': config["db_root"] + "/" + config["db"]["vfdb"] + "/core_dataset/core_pro_info.txt"
        }


        def statistical_other_res(input_file, db_name):
            df = pd.read_csv(input_file,sep="\t",index_col=0)
            tpm = str(db_name + "_anno_total_TPM")
            num = str(db_name + "_anno_term_num")
            tmp_df = pd.DataFrame()
            tmp_df[tpm] = df.sum(axis=0).round(6)
            tmp_df[num] = df.apply(lambda x: len(x[x > 0]),axis=0)
            return tmp_df.T

        def statistical_humann_res(input_file, db_name):
            df = pd.read_csv(input_file,sep="\t",index_col=0)
            unclassified = str(db_name + "_unclassified_COUNT")
            count = str(db_name + "_anno_total_COUNT")
            num = str(db_name + "_anno_term_num")
            tmp_df = pd.DataFrame()
            tmp_df[unclassified] = pd.DataFrame(df.iloc[0:2].sum(axis=0).round(6))
            tmp_df[count] = pd.DataFrame(df.iloc[2:].sum(axis=0).round(6))
            tmp_df[num] = df.apply(lambda x: len(x[x > 0]) - 2,axis=0)
            return tmp_df.T


        if not os.path.exists(anno_info_path):
            os.makedirs(anno_info_path)
        res_df = pd.DataFrame()
        for db in input.keys():
            shutil.copy(info_dict[db], anno_info_path)
            if db == "humann3":
                for i in input[db]:
                    sub_db = os.path.basename(i).replace("_clean.tsv","").replace("all_","")
                    df_stat = statistical_humann_res(i,sub_db)
                    res_df = pd.concat([res_df, df_stat],axis=0)
            else:
                df_stat = statistical_other_res(input[db], db)
                res_df = pd.concat([res_df, df_stat],axis=0)
        res_df.to_csv(output[0],sep="\t",header=True,index=True)