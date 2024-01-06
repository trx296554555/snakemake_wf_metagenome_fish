import os.path

import pandas as pd

rule gtdbtk_classify_wf:
    input:
        mags_dir=config["root"] + "/" + config["folder"]["bins_dereplication"] + "/mag_bins",
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
        "18: GTDB-Tk classify_wf"
    threads:
        80
    params:
        # -f/--full_tree : use the unsplit bacterial tree for the classify step;
        # this is the original GTDB-Tk approach (version < 2) and requires more than 403.74 GB of RAM to load the reference tree
        # force : continue processing if an error occurs on a single genome (default: False)
        out_dir=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf",
        prefix="gtdbtk"
    benchmark:
        config["root"] + "/benchmark/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf.benchmark.txt"
    log:
        config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify.log"
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


rule gtdbtk_infer_mag_tree:
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
        config["root"] + "/benchmark/" + config["folder"]["bins_classify"] + "/gtdbtk_infer_mag_tree.benchmark.txt"
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
        outdir=directory(config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf/denovo_tree"),
        prefix="fish",
        force=True,
        skip_gtdb_refs=False
    log:
        config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_de_novo_tree.log"
    shell:
        """
        gtdbtk de_novo_wf --cpus {threads} --genome_dir {input.bins_derep2} --out_dir {params.outdir} \
        --bacteria --skip_ani_screen --force {params.force} --skip_gtdb_refs {params.skip_gtdb_refs} \
        --extension fa --pplacer_cpus  {threads} --prefix {params.prefix} \
        --outgroup_taxon p__Patescibacteria > {log} 2>&1
        mv {params.outdir}/{params.prefix}.bac120.decorated.tree {output.tree}
        """


rule report_bins_classify:
    input:
        ar53_tsv=config["root"] + "/" + config["folder"]["bins_classify"] + "/gtdbtk_classify_wf/res.ar53.summary.tsv",
        bac120_tsv=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_classify_wf/res.bac120.summary.tsv",
        user_all_tree=config["root"] + "/" + config["folder"][
            "bins_classify"] + "/gtdbtk_mag_tree/user_all_unrooted.tree"
    output:
        bins_classify_report=config["root"] + "/" + config["folder"]["reports"] + "/08_bins_classify.report"
    run:
        if os.path.getsize(input.bac120_tsv) > 0:
            bac_table = pd.read_csv(input.bac120_tsv, sep="\t")
            print(bac_table.head)
            bac_table = bac_table[["user_genome", "classification", "fastani_ani", "closest_placement_ani"]]
            # final_ani 为 fastani_ani 和 closest_placement_ani 中的最大值
            bac_table["final_ani"] = bac_table[["fastani_ani", "closest_placement_ani"]].max(axis=1)