rule concatenation_MAGs:
    input:
        mags_dir=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins",
    output:
        concatenation_MAGs_fa=config["root"] + "/" + config["folder"][
            "bins_dereplication"] + "/concatenation_MAGs/concatenation_MAGs.fa"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "tools.yaml"
    message:
        "18.1: Concatenation MAGs"
    params:
        concatenation_MAGs_dir=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/concatenation_MAGs",
    shell:
        """
        mkdir -p {params.concatenation_MAGs_dir}/tmp
        for f in {input.mags_dir}/*.f*a; do
            file_name=$(basename $f)
            seqkit replace -p '.*' -r ${{file_name/.fa/}}={{nr}} -w 0 $f > {params.concatenation_MAGs_dir}/tmp/${{file_name}}
        done
        cat {params.concatenation_MAGs_dir}/tmp/*.fa > {output.concatenation_MAGs_fa}
        rm -rf {params.concatenation_MAGs_dir}/tmp
        """


rule get_MAGs_salmon_index:
    input:
        concatenation_MAGs_fa=config["root"] + "/" + config["folder"][
            "bins_dereplication"] + "/concatenation_MAGs/concatenation_MAGs.fa"
    output:
        index=directory(
            config["root"] + "/" + config["folder"]["bins_dereplication"] + "/concatenation_MAGs/salmon_index"),
        done=touch(
            config["root"] + "/" + config["folder"]["bins_dereplication"] + "/concatenation_MAGs/salmon_index.done")
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    message:
        "18.2 Indexing of MAGs"
    threads:
        12
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_dereplication"] + "/mags_salmon_index.log"
    log:
        config["root"] + "/" + config["folder"]["bins_dereplication"] + "/concatenation_MAGs/mags_salmon_index.log"
    shell:
        """
        salmon index -t {input.concatenation_MAGs_fa} -i {output.index} -p {threads} > {log} 2>&1
        """


rule quantify_MAGs_expression:
    input:
        index=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/concatenation_MAGs/salmon_index",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        expression=config["root"] + "/" + config["folder"]["bins_classify"] + "/quant_MAGs/{sample}/{sample}.sf",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        12
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_classify"] + "/quant_MAGs/{sample}.quant.log"
    log:
        config["root"] + "/" + config["folder"]["bins_classify"] + "/quant_MAGs/{sample}/{sample}.quant.log"
    message:
        "18.3 Quantify MAGs expression using salmon"
    params:
        dir=config["root"] + "/" + config["folder"]["bins_classify"] + "/quant_MAGs/{sample}"
    shell:
        """
        salmon quant -i {input.index} -l A --meta -1 {input.fq1} -2 {input.fq2} \
        -o {params.dir} -p {threads} --validateMappings --seqBias --quiet > {log} 2>&1
        mv {params.dir}/quant.sf {output.expression}
        """


rule gtdbtk_classify_wf:
    input:
        mags_dir=ancient(config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins"),
    output:
        ar53_tsv=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf/res.ar53.summary.tsv",
        bac120_tsv=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.bac120.summary.tsv",
        ar53_msa=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.ar53.user_msa.fasta.gz",
        bac120_msa=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.bac120.user_msa.fasta.gz",
        classify_done=touch(config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify.done")
    conda:
        config["root"] + "/" + config["envs"] + "/" + "gtdbtk.yaml"
    message:
        "18.4: GTDB-Tk classify_wf"
    threads:
        80
    params:
        # -f/--full_tree : use the unsplit bacterial tree for the classify step;
        # this is the original GTDB-Tk approach (version < 2) and requires more than 403.74 GB of RAM to load the reference tree
        # force : continue processing if an error occurs on a single genome (default: False)
        out_dir=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf",
        prefix="gtdbtk"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf.log"
    log:
        config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf.log"
    shell:
        """
        gtdbtk classify_wf --cpus {threads} --genome_dir {input.mags_dir} --out_dir {params.out_dir} \
        --full_tree --extension fa --pplacer_cpus {threads} --force \
        --prefix {params.prefix} --skip_ani_screen > {log} 2>&1
        if [ -f {params.out_dir}/{params.prefix}.ar53.summary.tsv ]; then
            cp -f {params.out_dir}/{params.prefix}.ar53.summary.tsv {output.ar53_tsv}
        else
            touch {output.ar53_tsv} {output.ar53_tsv}.is_null
        fi
        if [ -f {params.out_dir}/{params.prefix}.bac120.summary.tsv ]; then
            cp -f {params.out_dir}/{params.prefix}.bac120.summary.tsv {output.bac120_tsv}
        else
            touch {output.bac120_tsv} {output.bac120_tsv}.is_null
        fi
        if [ -f {params.out_dir}/align/{params.prefix}.ar53.user_msa.fasta.gz ]; then
            cp -f {params.out_dir}/align/{params.prefix}.ar53.user_msa.fasta.gz {output.ar53_msa}
        else
            touch {output.ar53_msa} {output.ar53_msa}.is_null
        fi
        if [ -f {params.out_dir}/align/{params.prefix}.bac120.user_msa.fasta.gz ]; then
            cp -f {params.out_dir}/align/{params.prefix}.bac120.user_msa.fasta.gz {output.bac120_msa}
        else
            touch {output.bac120_msa} {output.bac120_msa}.is_null
        fi
        """


rule gtdbtk_infer_MAGs_tree:
    input:
        ar53_msa=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.ar53.user_msa.fasta.gz",
        bac120_msa=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.bac120.user_msa.fasta.gz",
    output:
        user_ar53_tree=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_mag_tree/user_ar53_unrooted.tree",
        user_bac120_tree=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_mag_tree/user_bac120_unrooted.tree",
        user_all_tree=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_mag_tree/user_all_unrooted.tree"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "gtdbtk.yaml"
    message:
        "19: GTDB-Tk infer user mag tree"
    threads:
        80
    params:
        user_tree_dir=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_mag_tree/user_tree",
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_classify"] + "/gtdbtk_infer_mag_tree.log"
    log:
        config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_infer_mag_tree.log"
    shell:
        """
        mkdir -p {params.user_tree_dir}
        if [ ! -f {input.ar53_msa}.is_null ]; then
            gunzip -kc {input.ar53_msa} > {params.user_tree_dir}/ar53_msa.fa
            gtdbtk infer --msa_file {params.user_tree_dir}/ar53_msa.fa --out_dir {params.user_tree_dir} --prefix ar53 --cpus {threads} > {log} 2>&1
            gtdbtk convert_to_itol --input_tree {params.user_tree_dir}/ar53.unrooted.tree --output_tree {params.user_tree_dir}/ar53.unrooted.itol.tree
            cp -f {params.user_tree_dir}/ar53.unrooted.tree {output.user_ar53_tree}
        else
            touch {output.user_ar53_tree} {output.user_ar53_tree}.is_null
        fi
        if [ ! -f {input.bac120_msa}.is_null ]; then
            gunzip -kc {input.bac120_msa} > {params.user_tree_dir}/bac120_msa.fa
            gtdbtk infer --msa_file {params.user_tree_dir}/bac120_msa.fa --out_dir {params.user_tree_dir} --prefix bac120 --cpus {threads} >> {log} 2>&1
            gtdbtk convert_to_itol --input_tree {params.user_tree_dir}/bac120.unrooted.tree --output_tree {params.user_tree_dir}/bac120.unrooted.itol.tree
            cp -f {params.user_tree_dir}/bac120.unrooted.tree {output.user_bac120_tree}
        else
            touch {output.user_bac120_tree} {output.user_bac120_tree}.is_null
        fi
        if [ ! -f {input.ar53_msa}.is_null ] -a [ ! -f {input.bac120_msa}.is_null ]; then
            cat {input.ar53_msa} {input.bac120_msa} > {params.user_tree_dir}/user_msa_all.fa.gz
            gunzip {params.user_tree_dir}/user_msa_all.fa.gz
            gtdbtk infer --msa_file {params.user_tree_dir}/user_msa_all.fa --out_dir {params.user_tree_dir} --prefix all --cpus {threads} >> {log} 2>&1
            gtdbtk convert_to_itol --input_tree {params.user_tree_dir}/all.unrooted.tree --output_tree {params.user_tree_dir}/all.unrooted.itol.tree
            cp -f {params.user_tree_dir}/all.unrooted.tree {output.user_all_tree}
        else
            touch {output.user_all_tree} {output.user_all_tree}.is_null
        fi
        """


rule get_MAGs_classify:
    input:
        ar53_tsv=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf/res.ar53.summary.tsv",
        bac120_tsv=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.bac120.summary.tsv",
        user_all_tree=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_mag_tree/user_all_unrooted.tree"
    output:
        all_classify_tsv=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.all.summary.tsv",
    params:
        species_bins=config["root"] + "/" + config["folder"][
            "bins_dereplication"] + "/species_bins/dereplicated_genomes",
    run:
        table = pd.DataFrame()
        if os.path.getsize(input.bac120_tsv) > 0:
            bac_table = pd.read_csv(input.bac120_tsv,sep="\t")
            bac_table = bac_table[["user_genome", "classification", "fastani_ani", "closest_placement_ani"]]
            bac_table["final_ani"] = bac_table[["fastani_ani", "closest_placement_ani"]].max(axis=1)
            table = pd.concat([table, bac_table])
        if os.path.getsize(input.ar53_tsv) > 0:
            ar53_table = pd.read_csv(input.ar53_tsv,sep="\t")
            ar53_table = ar53_table[["user_genome", "classification", "fastani_ani", "closest_placement_ani"]]
            ar53_table["final_ani"] = ar53_table[["fastani_ani", "closest_placement_ani"]].max(axis=1)
            table = pd.concat([table, ar53_table])

        species_list = os.listdir(params.species_bins)
        species_list = [x.replace('.fa','') for x in species_list]
        table["species_bins"] = table["user_genome"].apply(lambda x: x if x in species_list else "")
        table.fillna(0,inplace=True)
        table.rename(columns={'user_genome': 'strain_bins'},inplace=True)
        table = table[["strain_bins", "species_bins", "classification", "final_ani"]]
        table.to_csv(output.all_classify_tsv,sep="\t",index=False)


# This rule needs to be modified according to the user's needs
rule gtdbtk_de_novo_tree:
    input:
        bins_derep2=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/rep_bins_user",
    output:
        tree=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/denovo_tree/bac120.decorated.tree"
    conda:
        config["root"] + "/" + config["envs"] + "/" + "gtdbtk.yaml"
    message:
        " GTDB-Tk de_novo_wf"
    threads:
        80
    params:
        out_dir=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf/denovo_tree",
        prefix="fish",
        force=True,
        skip_gtdb_refs=False
    log:
        config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_de_novo_tree.log"
    shell:
        """
        gtdbtk de_novo_wf --cpus {threads} --genome_dir {input.bins_derep2} --out_dir {params.out_dir} \
        --bacteria --skip_ani_screen --force {params.force} --skip_gtdb_refs {params.skip_gtdb_refs} \
        --extension fa --pplacer_cpus  {threads} --prefix {params.prefix} \
        --outgroup_taxon p__Patescibacteria > {log} 2>&1
        mv {params.out_dir}/{params.prefix}.bac120.decorated.tree {output.tree}
        """


rule generate_MAGs_classify_report:
    input:
        all_expression=expand(config["root"] + "/" + config["folder"][
            "bins_classify"] + "/quant_MAGs/{sample}/{sample}.sf",sample=get_run_sample()),
        all_classify_tsv=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.all.summary.tsv",
    output:
        merge_classify_report=config["root"] + "/" + config["folder"]["bins_classify"] + "/MAGs_classify_report.tsv"
    run:
        def get_bins_counts_from_expression(sf, sample_name):
            raw_df = pd.read_csv(sf,sep="\t",header=0)
            raw_df['MAGs'] = raw_df['Name'].apply(lambda x: x.split("=")[0])
            # Group by MAGs name, and then sum each set of NumReads
            raw_df = raw_df.groupby(['MAGs']).agg({'NumReads': 'sum'}).reset_index()
            raw_df.rename(columns={'NumReads': sample_name},inplace=True)
            return raw_df


        table = pd.DataFrame(columns=["MAGs"])
        for raw_sf in input.all_expression:
            sample = os.path.basename(raw_sf).split(".")[0]
            # Merge by MAGs column
            table = pd.merge(table,get_bins_counts_from_expression(raw_sf,sample),on='MAGs',how='outer')

        classify_anno_table = pd.read_csv(input.all_classify_tsv,sep="\t",header=0)
        table = pd.merge(table,classify_anno_table,left_on='MAGs',right_on='strain_bins',how='left')
        table['Taxon'] = table['classification'].apply(lambda x: x.replace('d__','k__').replace(';','|'))
        table = table[["MAGs", "Taxon"] + get_run_sample()]
        table.rename(columns={'MAGs': 'TaxID'},inplace=True)
        table.sort_values(by=['Taxon'],inplace=True)
        table.to_csv(output.merge_classify_report,sep="\t",index=False)


rule get_other_MAGs_salmon_index:
    input:
        concatenation_MAGs_fa=config["db_root"] + "/" + config["db"]["mag_db"] + "/{other_mags}_concatenation.fa",
    output:
        index=directory(
            config["root"] + "/" + config["folder"]["bins_recall"] + "/all_salmon_index/{other_mags}_index"),
        done=touch(
            config["root"] + "/" + config["folder"]["bins_recall"] + "/all_salmon_index/{other_mags}_index.done")
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    message:
        "18.5 Indexing of other MAGs"
    threads:
        24
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_recall"] + "/{other_mags}_salmon_index.log"
    log:
        config["root"] + "/" + config["folder"]["bins_recall"] + "/all_salmon_index/{other_mags}_index.log"
    shell:
        """
        salmon index -t {input.concatenation_MAGs_fa} -i {output.index} -p {threads} > {log} 2>&1
        """


rule recall_to_other_MAGs:
    input:
        index=config["root"] + "/" + config["folder"]["bins_recall"] + "/all_salmon_index/{other_mags}_index",
        fq1=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_1.fq",
        fq2=config["root"] + "/" + config["folder"]["rm_host"] + "/{sample}/{sample}_paired_2.fq",
    output:
        expression=config["root"] + "/" + config["folder"][
            "bins_recall"] + "/quant_MAGs/{other_mags}/{sample}/{sample}.sf",
    conda:
        config["root"] + "/" + config["envs"] + "/" + "salmon.yaml"
    threads:
        12
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_recall"] + "/quant_MAGs/{other_mags}_{sample}.quant.log"
    log:
        config["root"] + "/" + config["folder"]["bins_recall"] + "/quant_MAGs/{other_mags}/{sample}/{sample}.quant.log"
    message:
        "18.6 Quantify other MAGs expression using salmon"
    params:
        dir=config["root"] + "/" + config["folder"]["bins_recall"] + "/quant_MAGs/{other_mags}/{sample}"
    shell:
        """
        salmon quant -i {input.index} -l A --meta -1 {input.fq1} -2 {input.fq2} \
        -o {params.dir} -p {threads} --validateMappings --seqBias --quiet > {log} 2>&1
        mv {params.dir}/quant.sf {output.expression}
        """


rule generate_other_MAGs_classify_report:
    input:
        all_expression=expand(config["root"] + "/" + config["folder"][
            "bins_recall"] + "/quant_MAGs/{other_mags}/{sample}/{sample}.sf",other_mags=config["bins_recall"][
            "other_mags_list"],sample=get_run_sample()),
    output:
        other_classify_report_list=expand(
            config["root"] + "/" + config["folder"][
                "bins_recall"] + "/{other_mags}_MAGs_classify_report.tsv",other_mags=
            config["bins_recall"]["other_mags_list"])
    params:
        classify_tsv_db=config["db_root"] + "/" + config["db"]["mag_db"]
    run:
        def get_bins_counts_from_expression(sf, sample_name):
            raw_df = pd.read_csv(sf,sep="\t",header=0)
            raw_df['MAGs'] = raw_df['Name'].apply(lambda x: x.split("_")[0])
            # Group by MAGs name, and then sum each set of NumReads
            raw_df = raw_df.groupby(['MAGs']).agg({'NumReads': 'sum'}).reset_index()
            raw_df.rename(columns={'NumReads': sample_name},inplace=True)
            return raw_df


        for out_rpt in output.other_classify_report_list:
            other_mag = os.path.basename(out_rpt).split("_MAGs_")[0]
            classify_tsv = params.classify_tsv_db + "/" + other_mag + "_metadata.tsv"

            table = pd.DataFrame(columns=["MAGs"])
            for raw_sf in input.all_expression:
                if other_mag in raw_sf:
                    sample = os.path.basename(raw_sf).split(".")[0]
                    # Merge by MAGs column
                    table = pd.merge(table,get_bins_counts_from_expression(raw_sf,sample),on='MAGs',how='outer')

            classify_anno_table = pd.read_csv(classify_tsv,sep="\t",header=0)
            table = pd.merge(table,classify_anno_table,left_on='MAGs',right_on='Genome',how='left')
            table['Taxon'] = table['Lineage'].apply(lambda x: x.replace('d__','k__').replace(';','|'))
            table = table[["MAGs", "Taxon"] + get_run_sample()]
            table.rename(columns={'MAGs': 'TaxID'},inplace=True)
            table.sort_values(by=['Taxon'],inplace=True)
            table.to_csv(out_rpt,sep="\t",index=False)


rule report_MAGs_classify:
    input:
        merge_classify_report=config["root"] + "/" + config["folder"]["bins_classify"] + "/MAGs_classify_report.tsv",
        other_classify_report_list=expand(config["root"] + "/" + config["folder"][
            "bins_recall"] + "/{other_mags}_MAGs_classify_report.tsv",other_mags=config["bins_recall"][
            "other_mags_list"]) if config["bins_recall"]["enable"] else []
    output:
        MAGs_classify_report=config["root"] + "/" + config["folder"]["reports"] + "/08_MAGs_classify.report"
    params:
        total_reads=config["root"] + "/" + config["folder"]["reports"] + "/01_rm_host.report"
    run:
        # Statistical recall rate
        table = pd.DataFrame(columns=["Sample"])
        reads_table = pd.read_csv(params.total_reads,sep="\t",header=0)
        table["Sample"] = get_run_sample()
        table = pd.merge(table,reads_table[['Unnamed: 0', 'final pair1']],left_on='Sample',right_on='Unnamed: 0')
        table.drop(columns=['Unnamed: 0'],inplace=True)
        table.rename(columns={'final pair1': 'Total'},inplace=True)

        all_report = [input.merge_classify_report] + input.other_classify_report_list
        for classify_report in all_report:
            prefix = os.path.basename(classify_report).split("_classify_")[0]
            counts_table = pd.read_csv(classify_report,sep="\t",header=0)
            counts_table['k__'] = counts_table['Taxon'].apply(lambda x: prefix + '_' + x.split("|")[0].replace('k__',''))
            counts_table = counts_table[["k__"] + get_run_sample()]
            counts_table = counts_table.groupby(['k__']).agg('sum').T.reset_index()
            microbe_list = counts_table.columns.tolist()[1:]
            counts_table.rename(columns=counts_table.iloc[0]).reset_index(drop=True)
            table = pd.merge(table,counts_table,left_on='Sample',right_on='index')
            table.drop(columns=['index'],inplace=True)

        table.to_csv(output.MAGs_classify_report,sep="\t",index=False)

        collect_result(input.merge_classify_report)
        if input.other_classify_report_list:
            collect_result(input.other_classify_report_list,'other_bins_classify')
